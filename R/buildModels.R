library('inline')
library('Rcpp')

#' Builds models using locating arrays
#' 
#' @param la_path Path to locating array TSV.
#' @param factor_data Path to factor data file in TSV format.
#' @param response_dir Path to the response directory in TSV format.
#' @param response_column Response column as a string
#' @param log_data TRUE, FALSE Apply logarithm to response data
#' @param max_terms Number of terms in each model
#' @param models_n Number of models
#' @param new_models_n Number of models to be created at each iteration
#' @param logit TRUE, FALSE Create logistic regression models rather than linear
#' @return A text output of all top models and iterations as well as summary statistics.
#' @export
#' @examples
#' locatingArray_path = system.file("extdata", "programRepair_logit/LA.tsv", package="LATools")
#' factorData_path = system.file("extdata", "programRepair_logit/factor_data.tsv", package="LATools")
#' response_path = system.file("extdata", "programRepair_logit/Response", package="LATools")
#' ##EXAMPLE 1: Use logisitc regression model on software repair data set.
#' buildModels(locatingArray_path, factorData_path, response_path, response_column = "better_than_original", logit = TRUE)
#' ##EXAMPLE 2: Use linear regression model on network optimization data set.
#' locatingArray_path = system.file("extdata", "network_linear/LA.tsv", package="LATools")
#' factorData_path = system.file("extdata", "network_linear/factor_data.tsv", package="LATools")
#' response_path = system.file("extdata", "network_linear/Response", package="LATools")
#' buildModels(locatingArray_path, factorData_path, response_path, response_column = "MOS", logit = FALSE)
buildModels <-function(la_path, factor_data_path, response_dir, response_column, log_data, max_terms, models_n, new_models_n, logit){
  
  if(missing(la_path)){ 
    stop("Please supply a locating array file.")
  }
  if(missing(factor_data_path)){
    stop("Please supply a factor data file.")
  }
  if(missing(response_dir)){
    stop("Please supply a response directory.")
  }
  if(missing(response_column)){
    stop("Please supply the name of your response.")
  }
  if(missing(log_data)){
    log_data = FALSE
  }
  if(missing(max_terms)){
    max_terms = 20
  }
  if(missing(models_n)){
    models_n = 10
  }
  if(missing(new_models_n)){
    new_models_n = 5
  }
  if(missing(logit)){
    logit = FALSE
  }
  
  ##load la module and make local la
  la_module <- Module("LocatingArray_module")
  LocatingArray <- la_module$LocatingArray
  (la_ptr <- new(LocatingArray, la_path, factor_data_path))
  ##get tests for response vector
  tests = la_ptr$getTests()
  
  la2 <- makeLA2(la_path, factor_data_path)
  
  ##load cs module and make local csMatrix
  cs_module <- Module("CSMatrix_module")
  CSMatrix <- cs_module$CSMatrix
  cs_local <- new(CSMatrix, la_ptr)
  
  cs2 <- makeCSMatrix2(la2)
  
  ##load vec module and make local vec
  vec_module <- Module("VectorXf_module")
  VectorXf <- vec_module$VectorXf
  vec_local <- new(VectorXf, tests)
  vec_local$loadResponse(response_dir, response_column, log_data)
  
  vec2 <- makeVectorXF2(tests, response_dir, response_column, log_data)
  
  #createModels(LocatingArray *locatingArray, VectorXf *response, CSMatrix *csMatrix,int maxTerms, int models_n, int newModels_n)
  models <- createModels_wrapper(la2, vec2, cs2, max_terms, models_n, new_models_n, logit)
  
}
