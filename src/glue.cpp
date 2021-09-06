#include<RcppCommon.h>
#include <Rcpp.h>
#include <string>
#include "LocatingArray.h"
#include "ConstraintGroup.h"
#include "FactorData.h"
#include "CSMatrix.h"
#include "VectorXf.h"
#include "Noise.h"
#include "Model.h"
using namespace std;
using namespace Rcpp;

//make the exposed class visable 
RCPP_EXPOSED_CLASS(LocatingArray);
RCPP_EXPOSED_CLASS(CSMatrix);

LocatingArray* makeLA(std::string file, std::string factorDataFile){
  LocatingArray* la = new LocatingArray(file, factorDataFile);
  return la;
}

CSMatrix* makeCSMatrix(LocatingArray *locatingArray){
  CSMatrix* cs = new CSMatrix(locatingArray);
  return cs;
}


RCPP_MODULE(LocatingArray_module){
  class_<LocatingArray>("LocatingArray")
  .factory<std::string,std::string>(makeLA, "makes a locating array")
  .method("getTests", &LocatingArray::getTests);
}

RCPP_MODULE(CSMatrix_module){
  class_<CSMatrix>("CSMatrix")
  .factory<LocatingArray*>(makeCSMatrix, "makes a cs matrix");
}

RCPP_MODULE(VectorXf_module){
  class_<VectorXf>("VectorXf")
  .constructor<int>()
  .method("writeToFile", &VectorXf::writeToFile)
  .method("loadResponse", &VectorXf::loadResponse);
}
