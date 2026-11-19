#include <Rcpp.h>
#include <algorithm>
#include <random>
using namespace Rcpp;
// [[Rcpp::export]]
List dinuc_null(CharacterVector seqs, int iterations=200, int seed=260924) {
 std::mt19937 rng(seed);int n=seqs.size(),L=51;NumericMatrix positional(L-2,iterations);NumericMatrix per_sequence(n,2);NumericVector observed(L-2);
 for(int i=0;i<n;i++) {std::string s=as<std::string>(seqs[i]);for(int j=0;j<L-2;j++)if(s.substr(j,3)=="TAG"){observed[j]++;per_sequence(i,0)++;}
 for(int it=0;it<iterations;it++){std::vector<char> edges[256];for(int j=0;j<L-1;j++)edges[(unsigned char)s[j]].push_back(s[j+1]);for(int v=0;v<256;v++)std::shuffle(edges[v].begin(),edges[v].end(),rng);std::vector<char> stack(1,s[0]),trail;while(!stack.empty()){char v=stack.back();auto &e=edges[(unsigned char)v];if(e.empty()){trail.push_back(v);stack.pop_back();}else{char w=e.back();e.pop_back();stack.push_back(w);}}std::reverse(trail.begin(),trail.end());if(trail.size()!=51)stop("Invalid shuffled sequence");for(int j=0;j<L-2;j++)if(trail[j]=='T'&&trail[j+1]=='A'&&trail[j+2]=='G'){positional(j,it)++;per_sequence(i,1)+=1.0/iterations;}}
 }
 return List::create(_["observed"]=observed,_["null"]=positional,_["per_sequence"]=per_sequence);
}
