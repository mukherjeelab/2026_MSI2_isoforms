from pathlib import Path
from reportlab.pdfgen import canvas
from pypdf import PdfReader
r=Path('.codex-progression-build/renders');out=Path('output/MSI2_analysis_progression_final/MSI2_analysis_progression_final.pdf')
files=sorted(r.glob('slide-*.png'));assert len(files)==59
c=canvas.Canvas(str(out),pagesize=(1152,648))
c.setTitle('MSI2 Isoform CLIP Analysis: From Initial Signal to Robustness Testing')
for f in files:
 c.drawImage(str(f),0,0,1152,648);c.showPage()
c.save();assert len(PdfReader(str(out)).pages)==59
print(out)
