#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List path_vel_cpp(S4 OCN, String str = "RN", bool includeDownstreamNode = false){

  Environment pkg = Environment::namespace_env("spam");
  Function asMatrixSpam = pkg["as.matrix.spam"];
  //
  List L = OCN.slot(str);
  int nNodes = L["nNodes"];
  NumericVector foo = NumericVector::create(1000,0.1*nNodes);
  foo = ceiling(foo);
  int siz = nNodes*max(foo);
  IntegerVector set_col ( siz );
  IntegerVector set_row ( siz );
  NumericVector set_values ( siz );
  List downstreamPath = L["downstreamPath"];
  S4 PathLength_S4 = L["downstreamPathLength"];
  NumericMatrix PathLength = asMatrixSpam(PathLength_S4);
  int outlet = L["outlet"];
  NumericVector leng = L["leng"];
  NumericVector velocity = L["velocity"];
  int k = 0;
  for (int i{0}; i<nNodes; ++i)
  {
    List dP_i = downstreamPath[i];
    for (int j{0}; j<nNodes; ++j)
    {
      if (!(dP_i[j] == R_NilValue ))
      {
        set_row(k) = i;
        set_col(k) = j;
        NumericVector path = dP_i[j];
        if ( !((i == outlet - 1) && (j == outlet - 1)))
        {
          NumericVector l = leng[path-1];
          NumericVector v = velocity[path-1];
          l = l/v;
          if (includeDownstreamNode)
          {
            double foo = sum(l);
            set_values(k) = PathLength(i, j) / foo;
          }
          else
          {
            if (path.length()>1)
            {
              // l.erase(l.end()); // why not working?
              l = l[Range(0,l.length()-2)];
              double foo = sum(l);
              set_values(k) = PathLength(i, j) / foo;
            }
            else
              set_values(k) = velocity[i];
          }
        }
        else
          set_values(k) = velocity[i];
        k++;
      }
    }
  }
  set_values = set_values[ Range(0, k-1) ];
  set_row = set_row[ Range(0, k-1) ] + 1;
  set_col = set_col[ Range(0, k-1) ] + 1;

  List Lout = List::create( _["i"] = set_row, _["j"] = set_col, _["values"] = set_values);
  return(Lout);
}
