import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
const runtimeModules=process.env.RUNTIME_NODE_MODULES || '/Users/kathrynwalters/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.RUNTIME_NODE_MODULES=runtimeModules;
const {Presentation,PresentationFile,FileBlob}=await import(pathToFileURL(path.join(runtimeModules,'@oai/artifact-tool/dist/artifact_tool.mjs')));
const root=process.cwd(),tmp=path.join(root,'.codex-progression-build'),out=path.join(root,'output/MSI2_analysis_progression_final');
const skill='/Users/kathrynwalters/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations';
const {resolvePresentationFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')));
const family=resolvePresentationFont(), data=JSON.parse(await fs.readFile(path.join(out,'slide_content.json'),'utf8'));
const p=Presentation.create({slideSize:{width:1600,height:900}}),native=[];
const ink='#193647',muted='#51636D',teal='#087E82';
function text(s,t,x,y,w,h,size=25,bold=false,color=ink){let z=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});z.text=t;z.text.style={typeface:family,fontSize:size,bold,color,alignment:'left',verticalAlignment:'top',autoFit:'none'};return z;}
function lines(s,items,x,y,w,h,size=26){let step=h/Math.max(items.length,1);items.forEach((v,i)=>text(s,v,x,y+i*step,w,step-12,size));}
function table(s,values,x,y,w,h,size=24){const n=values[0].length;let widths=n===3?[w*.26,w*.29,w*.45]:Array(n).fill(w/n);if(values[0][0]==='Figure ID')widths=[w*.12,w*.29,w*.59];const t=s.tables.add({rows:values.length,columns:n,left:x,top:y,width:w,height:h,values,columnWidths:widths});t.borders.assign({style:'solid',fill:'#DCE3E6',width:1});for(let r=0;r<values.length;r++)for(let c=0;c<n;c++){let cell=t.getCell(r,c);cell.fill=r===0?'#E8EFF1':(r%2?'#FFFFFF':'#F5F7F8');cell.text.style={typeface:family,fontSize:size,bold:r===0,color:ink,alignment:'left',verticalAlignment:'middle',autoFit:'shrinkText'};}return t;}
async function pic(s,f,x,y,w,h){s.images.add({blob:await fs.readFile(path.join(root,f)),contentType:'image/png',alt:f,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
for(const d of data.slides){const s=p.slides.add();text(s,d.tag,55,24,1300,28,17,true,teal);text(s,d.title,55,66,1490,d.layout==='cover'?180:104,d.layout==='cover'?57:43,true);if(d.sign)text(s,'ΔΔ < 0: stronger relative 328 CLIP/SMI',840,25,700,27,19,true,'#A65428');if(d.why)text(s,d.why,58,172,1480,61,24,false,muted);
const y=d.why?244:194;
if(d.layout==='cover'){lines(s,d.findings,65,345,1370,230,35);}
else if(d.layout==='figure_full'){await pic(s,d.figures[0],50,235,1500,490);text(s,d.findings.join('  •  '),58,730,1480,64,22);}
else if(d.layout==='figure_large'){await pic(s,d.figures[0],48,228,1015,520);lines(s,d.findings,1100,270,435,400,25);}
else if(d.layout==='figure'){await pic(s,d.figures[0],45,y,1030,540);lines(s,d.findings,1110,y+30,435,490,25);}
else if(d.table){native.push(d.number);let full=['table_full','inventory'].includes(d.layout),has=d.findings.length>0;table(s,d.table,55,y,full||!has?1490:1000,full?540:510,d.layout==='inventory'?18:(full?20:24));if(has&&!full)lines(s,d.findings,1100,y+25,435,455,25);if(has&&full)text(s,d.findings.join('  '),60,755,1470,54,20);}
else{lines(s,d.findings,70,y+30,1420,Math.min(510,d.findings.length*97),d.layout==='progression'?28:31);}
if(d.changed)text(s,d.changed,58,794,1480,48,25,true,teal);
const src=(d.ids.length?d.ids.join(' · ')+' | ':'')+(d.figure_aliases[0]||d.scripts[0]||'Source details in speaker notes');text(s,src.length>180?src.slice(0,177)+'…':src,58,850,1450,21,13,false,muted);if(d.foot)text(s,d.foot,58,874,1460,22,13,false,muted);text(s,String(d.number),1535,855,45,25,17,false,muted);s.speakerNotes.textFrame.setText(d.note+'\n\nCaption: '+d.foot+'\n\nMain result: '+d.changed);
}
await fs.mkdir(path.join(tmp,'renders'),{recursive:true}); await fs.mkdir(path.join(root,'.codex-progression-checks'),{recursive:true}); const stamp=Date.now();
const candidate=path.join(tmp,'candidate.pptx'),final=path.join(tmp,`validated-${stamp}.pptx`);await(await PresentationFile.exportPptx(p)).save(candidate);
await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:final,pythonExecutable:'/Users/kathrynwalters/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','15240000,8572500','--validate-heading-fit',...native.flatMap(n=>['--require-native-table-slide',String(n)])],explicitTotalSlideCount:data.slides.length,requiredNativeTableOwnerSlides:native,fontPolicy:{basis:'design',families:[family]},verifyArtifactToolImport:true,receiptPath:path.join(root,'.codex-progression-checks',`validation-${stamp}.json`)});
const checked=await PresentationFile.importPptx(await FileBlob.load(final));for(let i=0;i<checked.slides.items.length;i++){const z=await checked.export({slide:checked.slides.items[i],format:'png',scale:1.5});await fs.writeFile(path.join(tmp,'renders',`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await z.arrayBuffer()));console.log('Rendered',i+1);}
await fs.copyFile(final,path.join(out,'MSI2_analysis_progression_final.pptx'));console.log('DONE',data.slides.length);
