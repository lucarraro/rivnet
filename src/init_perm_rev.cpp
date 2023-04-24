 #include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List init_perm_rev_cpp(List downNode_rev, IntegerVector outlet){
int nNodes = downNode_rev.length();
IntegerVector nodesToExplore = outlet;
IntegerVector reverse_perm (nNodes);
IntegerVector upNodes(0);
IntegerVector tmp(10);
int node;
bool flag = false;
int k = 0;
while(nodesToExplore.length()>0){
node = nodesToExplore[0];
reverse_perm[k] = node;
nodesToExplore.erase(0,1);
flag = false;
if (downNode_rev[node-1] != R_NilValue) {
upNodes = downNode_rev[node-1];
flag = true;
}
k++;
while(flag){
node = upNodes[0];
reverse_perm[k] = node;
if (upNodes.size() > 1){
tmp = rep(0,upNodes.size() + nodesToExplore.size() - 1);
int i=0;
for ( ; i<(upNodes.size()-1); i++) tmp[i] = upNodes[i+1];
for (int j{0}; j<nodesToExplore.size(); i++, j++) tmp[i] = nodesToExplore[j];
nodesToExplore = tmp;
//for (int i=1; i<upNodes.size(); i++) nodesToExplore.push_front(upNodes[i]); // not faster
}
flag = false;
if (downNode_rev[node-1] != R_NilValue) {
upNodes = downNode_rev[node-1];
flag = true;
}
k++;
}
}
IntegerVector perm = rev(reverse_perm);
int noDAG = 0;
List Lout = List::create(_["perm"]=perm, _["noDAG"]=noDAG);
return(Lout);
}