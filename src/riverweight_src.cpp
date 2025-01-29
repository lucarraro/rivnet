#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double eval_weight_cpp(double dst, List weightList, double A){
  //
  String distFunction = weightList["func"];
  double weightDist = weightList["scale.length"];
  bool FA = weightList["FA"];
  double weight = 1; // default value when distFunction is not a standard name
  if ((distFunction=="exponential") | (distFunction=="gexponential") ){
    weight = exp(-dst/weightDist);
  } else if (distFunction=="cauchy"){
    weight = pow(weightDist,2)/(pow(weightDist,2) + pow(dst,2));
  } else if (distFunction=="power"){
    weight = 1/(pow(dst + 1,weightDist));			
  } else if (distFunction=="linear"){
    weight = 1-dst/weightDist;
    if (weight<0) weight = 0;
  } 
  if (FA==true) weight=weight*A;
  return(weight);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix eval_wu_exp_cpp(NumericMatrix val, S4 OCN, NumericVector wl, 
								bool FA = false, bool unweighted = false){
  //
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector A = FD["A"];
  IntegerVector downNode = FD["downNode"];
  IntegerVector perm = FD["perm"];
  NumericVector weights (nNodes);
  for (int i{ 0 }; i<nNodes; ++i)  weights[i] = 1;
  NumericMatrix val_wu = clone(val); 
  if (FA == true) {
    for (int i{ 0 }; i<nNodes; ++i) {
      for (int nf{0}; nf<nFields; ++nf) val_wu(i,nf) = val_wu(i,nf)*A[i];
      weights[i] = weights[i]*A[i];
    }
  }
  for (int i {0}; i<nNodes; ++i){
    int j = perm[i];
    int d = downNode[j-1];
	if (d!=0){
    for (int nf{0}; nf<nFields; ++nf) val_wu(d-1,nf) +=  val_wu(j-1,nf)*wl[j-1];
    if (unweighted==true){
		weights[d-1] += 1;
	} else {
	weights[d-1] +=  weights[j-1]*wl[j-1];
	}
	}
  }
  for (int i {0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) val_wu(i,nf) = val_wu(i,nf)/weights[i];
  }
  return(val_wu);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix eval_wu_euclidean_cpp(NumericMatrix val, S4 OCN, 
									List weightNum, List weightDen,
									IntegerVector hw){
  //
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector X = FD["X"];
  NumericVector Y = FD["Y"];
  IntegerVector downNode = FD["downNode"];
  NumericVector A = FD["A"];
  NumericMatrix num = clone(val);
  NumericVector den (nNodes);
  for (int i{ 0 }; i<nNodes; ++i) den[i] = 1;
  NumericMatrix ll (nNodes,nFields);
  double dst (1);
  double weight_num (1);
  double weight_den (1);
  IntegerVector alreadyProcessed (nNodes+1);
  IntegerVector path;
  int lhw = hw.length();
  // loop once on each headwater
  for (int i{ 0 }; i<lhw; ++i) {
    int j = hw[i];
    path = j;
    alreadyProcessed[j] = 1; // alreadyProcessed has nNodes + 1 components
    int k = downNode[j-1];
    bool flag = true;
    while (k!=0) { // move down
      int pathlength = path.length();
      // copy path elements to slots of k
      for (int ii{ 0 }; ii<pathlength; ++ii){
        int ip = path[ii];
        // here add to num and den
        dst = sqrt(pow(X[ip-1]-X[k-1],2) + pow(Y[ip-1]-Y[k-1],2));
		
		weight_num = eval_weight_cpp(dst,weightNum,A[ip-1]);
        weight_den = eval_weight_cpp(dst,weightDen,A[ip-1]);

        for (int nf{0}; nf<nFields; ++nf) num(k-1,nf) += val(ip-1,nf)*weight_num;
        den[k-1] += weight_den;
      }
      // check if node already processed;
      // if not, add to path
      // (if node already processed, only nodes upstream 
      // of the confluence with already explored path 
      // are copied)
      if (flag==true)  flag = (alreadyProcessed[k] == 0);
      if (flag==true){
        alreadyProcessed[k] = 1;
        path.push_back(k);
      }
      k = downNode[k-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix eval_wu_euclidean_stream_cpp(NumericMatrix val, S4 OCN, 
										List weightNum, List weightDen,
										IntegerVector hw, NumericVector distRiver, 
                                        IntegerVector rvr){
										
  //
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector X = FD["X"];
  NumericVector Y = FD["Y"];
  IntegerVector downNode = FD["downNode"];
  NumericVector A = FD["A"];
  NumericMatrix num = clone(val);
  NumericVector den (nNodes);
  for (int i{ 0 }; i<nNodes; ++i) den[i] = 1;
  LogicalVector isRiver (nNodes);
  for (int i{ 0 }; i<rvr.length(); ++i) {
    int j = rvr[i];
    isRiver[j-1] = true;
  }
  NumericMatrix ll (nNodes,nFields);
  double dst (1);
  double weight_num (1);
  double weight_den (1);
  
  IntegerVector alreadyProcessed (nNodes+1);
  IntegerVector path;
  int lhw = hw.length();
  // loop once on each headwater
  for (int i{ 0 }; i<lhw; ++i) {
    int j = hw[i];
    path = j;
    alreadyProcessed[j] = 1; // alreadyProcessed has nNodes + 1 components
    int k = downNode[j-1];
    bool flag = true;
    while (k!=0) { // move down
      int pathlength = path.length();
      // copy path elements to slots of k
      for (int ii{ 0 }; ii<pathlength; ++ii){
        int ip = path[ii];
        
        if (isRiver[ip-1] && isRiver[k-1]) dst = 0;
        if (~isRiver[ip-1] && ~isRiver[k-1]) dst = sqrt(pow(X[ip-1]-X[k-1],2) + pow(Y[ip-1]-Y[k-1],2));
        if (~isRiver[ip-1] && isRiver[k-1]) dst = distRiver[ip-1];
        
		weight_num = eval_weight_cpp(dst,weightNum,A[ip-1]);
        weight_den = eval_weight_cpp(dst,weightDen,A[ip-1]);
		
        for (int nf{0}; nf<nFields; ++nf) num(k-1,nf) += val(ip-1,nf)*weight_num;
        den[k-1] += weight_den;
      }
      // check if node already processed;
      // if not, add to path
      // (if node already processed, only nodes upstream 
      // of the confluence with already explored path 
      // are copied)
      if (flag==true)  flag = (alreadyProcessed[k] == 0);
      if (flag==true){
        alreadyProcessed[k] = 1;
        path.push_back(k);
      }
      k = downNode[k-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix eval_wu_euclidean_stream_cpp_equalND(NumericMatrix val, S4 OCN, 
										List weightNum, IntegerVector hw,
										NumericVector distRiver, 
                                        IntegerVector rvr,
										bool unweighted=false){
										
  //
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector X = FD["X"];
  NumericVector Y = FD["Y"];
  IntegerVector downNode = FD["downNode"];
  NumericVector A = FD["A"];
  NumericMatrix num = clone(val);
  NumericVector den (nNodes);
  for (int i{ 0 }; i<nNodes; ++i) den[i] = 1;
  LogicalVector isRiver (nNodes);
  for (int i{ 0 }; i<rvr.length(); ++i) {
    int j = rvr[i];
    isRiver[j-1] = true;
  }
  NumericMatrix ll (nNodes,nFields);
  double dst (1);
  double weight_num (1);
  
  IntegerVector alreadyProcessed (nNodes+1);
  IntegerVector path;
  int lhw = hw.length();
  // loop once on each headwater
  for (int i{ 0 }; i<lhw; ++i) {
    int j = hw[i];
    path = j;
    alreadyProcessed[j] = 1; // alreadyProcessed has nNodes + 1 components
    int k = downNode[j-1];
    bool flag = true;
    while (k!=0) { // move down
      int pathlength = path.length();
      // copy path elements to slots of k
      for (int ii{ 0 }; ii<pathlength; ++ii){
        int ip = path[ii];
        
        if (isRiver[ip-1] && isRiver[k-1]) dst = 0;
        if (~isRiver[ip-1] && ~isRiver[k-1]) dst = sqrt(pow(X[ip-1]-X[k-1],2) + pow(Y[ip-1]-Y[k-1],2));
        if (~isRiver[ip-1] && isRiver[k-1]) dst = distRiver[ip-1];
        
		weight_num = eval_weight_cpp(dst,weightNum,A[ip-1]);
		
        for (int nf{0}; nf<nFields; ++nf) num(k-1,nf) += val(ip-1,nf)*weight_num;
        if (unweighted==true){
			den[k-1] += 1;
		} else {
			den[k-1] += weight_num;
		}
      }
      // check if node already processed;
      // if not, add to path
      // (if node already processed, only nodes upstream 
      // of the confluence with already explored path 
      // are copied)
      if (flag==true)  flag = (alreadyProcessed[k] == 0);
      if (flag==true){
        alreadyProcessed[k] = 1;
        path.push_back(k);
      }
      k = downNode[k-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector eval_wu_generic_flow_cpp_equalND(NumericMatrix val, S4 OCN, 
									 List weightNum, IntegerVector hw,
									 bool unweighted = false){
  
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  IntegerVector downNode = FD["downNode"];
  IntegerVector toRN = FD["toRN"];
  NumericVector leng = FD["leng"];
  NumericVector A = FD["A"];
  bool stream = weightNum["stream"]; // but this doesn't differentiate num vs den
  NumericMatrix num (nNodes,nFields);
  NumericVector den (nNodes);	
  NumericMatrix ll (nNodes,nFields);
  LogicalVector alreadyProcessed (nNodes);
  LogicalVector isRN (nNodes);  
  for (int i{0}; i<nNodes; ++i){
    if (toRN[i] > 0) isRN[i] = true;
  }
  int source (1);
  NumericVector sourceval (nFields);
  double weight_num (1);
  for (int i{0}; i<hw.length(); ++i){ // 
    int j = hw[i];
    IntegerVector path (0);
    NumericVector dst (0);
    while (j != 0){
      if (alreadyProcessed[j-1] == false){
        path.push_back(j);
        dst.push_back(0);
      }
      for (int pp{0}; pp<path.length(); ++pp){
        if (((stream==true) & (isRN[j-1]==false)) | (stream==false)) dst[pp] +=  leng[j-1];
      }
      
      NumericVector innerProd (nFields);
      double sumweight = 0;
      for (int pp{0}; pp<path.length(); ++pp){
        source = path[pp];
		weight_num = eval_weight_cpp(dst[pp],weightNum,A[source-1]);
        sourceval = val(source-1,_);
        for (int nf{0}; nf<nFields; ++nf) innerProd[nf] += weight_num*sourceval[nf];
        if (unweighted==true){
			sumweight += 1;	
		} else {
		sumweight += weight_num;		
		}
		 
      }
      
      for (int nf{0}; nf<nFields; ++nf) num(j-1, nf) += innerProd[nf];
      den[j-1] += sumweight; 
      alreadyProcessed[j-1] = true;
      j = downNode[j-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector eval_wu_generic_flow_cpp(NumericMatrix val, S4 OCN, 
									 List weightNum, List weightDen,
                                     IntegerVector hw){
  
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  IntegerVector downNode = FD["downNode"];
  IntegerVector toRN = FD["toRN"];
  NumericVector leng = FD["leng"];
  NumericVector A = FD["A"];
  bool stream = weightNum["stream"]; // but this doesn't differentiate num vs den
  NumericMatrix num (nNodes,nFields);
  NumericVector den (nNodes);	
  NumericMatrix ll (nNodes,nFields);
  LogicalVector alreadyProcessed (nNodes);
  LogicalVector isRN (nNodes);  
  for (int i{0}; i<nNodes; ++i){
    if (toRN[i] > 0) isRN[i] = true;
  }
  int source (1);
  NumericVector sourceval (nFields);
  double weight_num (1);
  double weight_den (1);
  for (int i{0}; i<hw.length(); ++i){ // 
    int j = hw[i];
    IntegerVector path (0);
    NumericVector dst (0);
    while (j != 0){
      if (alreadyProcessed[j-1] == false){
        path.push_back(j);
        dst.push_back(0);
      }
      for (int pp{0}; pp<path.length(); ++pp){
        if (((stream==true) & (isRN[j-1]==false)) | (stream==false)) dst[pp] +=  leng[j-1];
      }
      
      NumericVector innerProd (nFields);
      double sumweight = 0;
      for (int pp{0}; pp<path.length(); ++pp){
        source = path[pp];
		weight_num = eval_weight_cpp(dst[pp],weightNum,A[source-1]);
		weight_den = eval_weight_cpp(dst[pp],weightDen,A[source-1]);
        sourceval = val(source-1,_);
        for (int nf{0}; nf<nFields; ++nf) innerProd[nf] += weight_num*sourceval[nf];
        sumweight += weight_den;	 
      }
      
      for (int nf{0}; nf<nFields; ++nf) num(j-1, nf) += innerProd[nf];
      den[j-1] += sumweight; 
      alreadyProcessed[j-1] = true;
      j = downNode[j-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector dist_to_river_cpp(S4 OCN, IntegerVector no_river, IntegerVector rvr){ 
  //
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector X = FD["X"];
  NumericVector Y = FD["Y"];
  NumericVector distRiver (nNodes);
  NumericVector eucledianDist (rvr.length());
  
  for (int i{ 0 }; i<no_river.length(); ++i) {
    int j = no_river[i];
    for (int k{ 0 }; k<rvr.length(); ++k) {
      int kk = rvr[k] - 1;
      eucledianDist[k] = sqrt(pow(X[j-1] - X[kk],2) + pow(Y[j-1] - Y[kk],2));
    }
    distRiver[j-1] = min(eucledianDist);
  }
  
  return(distRiver);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix eval_wu_euclidean_cpp_equalND(NumericMatrix val, S4 OCN, 
											List weightList, IntegerVector hw,
											bool unweighted = false){
  //
  int nFields = val.ncol();
  List FD = OCN.slot("FD");
  int nNodes = FD["nNodes"];
  NumericVector X = FD["X"];
  NumericVector Y = FD["Y"];
  IntegerVector downNode = FD["downNode"];
  NumericVector A = FD["A"];
  NumericMatrix num = clone(val);
  NumericVector den (nNodes);
  for (int i{ 0 }; i<nNodes; ++i) den[i] = 1;
  NumericMatrix ll (nNodes,nFields);
  double dst (1);
  double weight (1);
  
  IntegerVector alreadyProcessed (nNodes+1);
  IntegerVector path;
  int lhw = hw.length();
  // loop once on each headwater
  for (int i{ 0 }; i<lhw; ++i) {
    int j = hw[i];
    path = j;
    alreadyProcessed[j] = 1; // alreadyProcessed has nNodes + 1 components
    int k = downNode[j-1];
    bool flag = true;
    while (k!=0) { // move down
      int pathlength = path.length();
      // copy path elements to slots of k
      for (int ii{ 0 }; ii<pathlength; ++ii){
        int ip = path[ii];
        // here add to num and den
        dst = sqrt(pow(X[ip-1]-X[k-1],2) + pow(Y[ip-1]-Y[k-1],2));
		
		weight = eval_weight_cpp(dst,weightList,A[ip-1]);

		
        for (int nf{0}; nf<nFields; ++nf) num(k-1,nf) += val(ip-1,nf)*weight;
		if (unweighted==true){
			den[k-1] += 1;	
		} else {
		den[k-1] += weight;	
		}
        
      }
      // check if node already processed;
      // if not, add to path
      // (if node already processed, only nodes upstream 
      // of the confluence with already explored path 
      // are copied)
      if (flag==true)  flag = (alreadyProcessed[k] == 0);
      if (flag==true){
        alreadyProcessed[k] = 1;
        path.push_back(k);
      }
      k = downNode[k-1];
    }
  }
  for (int i{0}; i<nNodes; ++i){
    for (int nf{0}; nf<nFields; ++nf) ll(i,nf) = num(i,nf)/den[i];
  }
  
  return(ll);
}
