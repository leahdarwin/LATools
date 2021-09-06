#include <cmath>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <Rcpp.h>
#include<RcppCommon.h>
#include "CSMatrix.h"
#include "Model.h"
#include "LocatingArray.h"
#include "Noise.h"
#include "Occurrence.h"
#include "VectorXf.h"

using namespace std;
using namespace Rcpp;

// static member of the class
WorkSpace *Model::workSpace = NULL;

// loader section

VectorXf* loadResponseVector(VectorXf *response, string directory, string column, bool performLog) {
	struct dirent *dp;

	// open the directory
	DIR *dir = opendir(directory.c_str());
	if (dir == NULL) {
		cout << "Cannot open directory \"" << directory << "\"" << endl;
		exit(0);
	}

	int cols, col_i;
	string header;

	int rows = response->getLength();

	string tempString;
	float tempFloat;
	int tempInt;

	// the response count vector
	VectorXf *responseCount = new VectorXf(rows);

	// grab data vectors
	float *data = response->getData();
	float *responseData = responseCount->getData();

	// initialize vectors
	for (int row_i = 0; row_i < rows; row_i++) {
		data[row_i] = 0;
		responseData[row_i] = 0;
	}

	while ((dp = readdir(dir)) != NULL) {

		// go through not hidden files
		if (dp->d_name[0] != '.') {

			string file = directory + "/" + dp->d_name;
			cout << "File: \"" << file << "\"" << endl;

			// initialize the file input stream
			ifstream ifs(file.c_str(), ifstream::in);

			// read the headers
			ifs >> tempInt;

			// make sure the rows match
			if (tempInt != rows) {
				cout << "Row mismatch (LA vs RE): " << file << endl;
				cout << "Expected " << rows << " but received " << tempInt << endl;
				exit(0);
			}

			// find the relevant column
			col_i = -1;		// this is the column index we are lookin for
			string line = "";
			string value = "";
			getline(ifs, line); // terminating new line of rows
			getline(ifs, line); // headers
			stringstream colStream(line);
			for (int tcol_i = 0; colStream >> value; tcol_i++) {
				if (value == column) {
					col_i = tcol_i;
					break;
				}
			}

			// ensure we found the relevant column
			if (col_i == -1) {
				cout << "Could not find column \"" << column << "\" in \"" << file << "\"" << endl;
			} else {
				// read into the response vector struct
				data = response->getData();
				for (int row_i = 0; row_i < rows; row_i++) {
					getline(ifs, line);
					stringstream myStream(line);

					for (int tcol_i = 0; getline(myStream, value, '\t'); tcol_i++) {
						if (tcol_i == col_i) {
							if (value == "") {
								cout << "No value found" << endl;
							} else {
								tempFloat = atof(value.c_str());

								if (tempFloat == 0) {
									//cout << "Warning: Float value response was 0 in " << file << " in row " << row_i << endl;
								}

								// process the data
								responseData[row_i] += 1;
								if (performLog) {
									data[row_i] += log(tempFloat);
								} else {
									data[row_i] += tempFloat;
								}

							}

						}
					}

				}

			}

			// close the file
			ifs.close();

		}


	}

	// average the responses
	data = response->getData();
	float minResponse = data[0] / responseData[0];
	float maxResponse = data[0] / responseData[0];
	for (int row_i = 0; row_i < rows; row_i++) {
		//data[row_i] /= responseData[row_i];
		//cout << data[row_i] << " from " << responseData[row_i] << " responses " << endl;

		if (data[row_i] < minResponse) minResponse = data[row_i];
		if (data[row_i] > maxResponse) maxResponse = data[row_i];
	}

	// add noise as necessary
	float range = maxResponse - minResponse;
	//if (noise != NULL) {
	//	for (int row_i = 0; row_i < rows; row_i++) {
	//		data[row_i] = noise->addNoise(data[row_i], range);
	//	}
	//}

	// delete response count vector
	delete responseCount;

	// close the directory
	closedir(dir);

	// calculate SStot for r-squared calculations
	response->calculateSStot();

	cout << "Loaded responses" << endl;

	return response;
}


// comparison, not case sensitive.
bool compareOccurrence(const Occurrence *first, const Occurrence *second) {
	return (first->count < second->count);
}

void allocateOccurrences(Occurrence *occurrence, int t, int factors, list<Occurrence*> *occurrenceLists) {

	// no more factors or interactions to allocate for
	if (factors == 0 || t == 0) {
		occurrence->list = NULL;
		return;
	}

	// allocate memory for the next dimension of occurrences
	occurrence->list = new Occurrence[factors];

	for (int factor_i = 0; factor_i < factors; factor_i++) {
		// initialize occurrence
		occurrence->list[factor_i].factorList = new int[occurrence->factorList_n + 1];
		occurrence->list[factor_i].factorList_n = occurrence->factorList_n + 1;
                occurrence->list[factor_i].rSquaredContribution = 0;
		occurrence->list[factor_i].count = 0;
		occurrence->list[factor_i].magnitude = 0;

		// populate the factor list
		for (int factorList_i = 0; factorList_i < occurrence->factorList_n; factorList_i++) {
			occurrence->list[factor_i].factorList[factorList_i] = occurrence->factorList[factorList_i];
		}
		occurrence->list[factor_i].factorList[occurrence->factorList_n] = factor_i;

		// push new occurrence into list
		occurrenceLists[occurrence->factorList_n].push_back(&occurrence->list[factor_i]);

		// allocate the next dimension
		allocateOccurrences(&occurrence->list[factor_i], t - 1, factor_i, occurrenceLists);
	}
}

void deallocateOccurrences(Occurrence *occurrence, int factors) {

	if (occurrence != NULL && occurrence->list != NULL) {

		for (int factor_i = 0; factor_i < factors; factor_i++) {
			// deallocate list of factors for occurrences
			delete[] occurrence->list[factor_i].factorList;

			// deallocate the next dimension
			deallocateOccurrences(&occurrence->list[factor_i], factor_i);
		}
		delete[] occurrence->list;

	}

}

/*
	ANALYSIS Procedure:

	Priority queue begins with 1 model: the model with no terms but an intercept.
	At every iteration, the models are pulled out of the queue, one by one, and for each,
	the top (n) terms are taken, one at a time, based on the distance to residuals.
	For every term, it is added to the current model and r-squared is calculated.
	The new model is then placed in the next priority queue. This process stops
	when all models have as many terms as maxTerms. Each model taken out of the queue
	produces newModels_n new models that are added to the next queue. The queues are
	then priority queues (by model R^2) with a maximum number of models. The same model
	(with the same terms) can be generated in multiple ways, and in this case, duplicates
	will not be added and the function below prints "Duplicate Model!!!".

*/

int createModels(LocatingArray *locatingArray, VectorXf *response, CSMatrix *csMatrix,
					int maxTerms, int models_n, int newModels_n, bool logit) {
	cout << "Creating Models..." << endl;
	Model::setupWorkSpace(response->getLength(), maxTerms);

	// work variables
	Model *model;		// current model we are working on

	Model **topModels = new Model*[models_n];
	Model **nextTopModels = new Model*[models_n];

	// populate initial top models
	topModels[0] = new Model(response, maxTerms, csMatrix);
	for (int model_i = 1; model_i < models_n; model_i++)
		topModels[model_i] = NULL;

	// struct for linked list of top columns of CS matrix
	struct ColDetails {
		int termIndex;
		float dotProduct;
		bool used;

	} *colDetails;

	// allocate memory for top columns
	colDetails = new ColDetails[csMatrix->getCols()];

	while (topModels[0] != NULL && topModels[0]->getTerms() < maxTerms) {
		// LOOP HERE


		// make sure all next top models are NULL
		for (int model_i = 0; model_i < models_n; model_i++) {
			nextTopModels[model_i] = NULL;
		}

		// grab the models from the topModels priority queue
		for (int model_i = 0; model_i < models_n; model_i++) {

			// we are done finding the next top models if we hit a NULL model
			if (topModels[model_i] == NULL) break;


			// grab the model from top models queue
			model = topModels[model_i];

//			cout << "Grabbing model" << endl;
//			model->printModelFactors();
//			cout << endl;


			// grab the distances to columns in cs matrix
			for (int col_i = 0; col_i < csMatrix->getCols(); col_i++) {
				colDetails[col_i].dotProduct =
					csMatrix->getProductWithCol(col_i, model->getResiVec());
				colDetails[col_i].termIndex = col_i;
				colDetails[col_i].used = model->termExists(col_i);
			}


			// find the columns with the largest dot products (at most as many as we have models)
			for (int colsUsed = 0; colsUsed < newModels_n; colsUsed++) {

				float largestDotProduct = 0;
				int bestCol_i = -1;

				// find the unused column with smallest distance
				for (int col_i = 0; col_i < csMatrix->getCols(); col_i++) {

					// check if this column has larger dot product and should be marked as best
					if (!colDetails[col_i].used &&
						(bestCol_i == -1 || colDetails[col_i].dotProduct > largestDotProduct)) {
						bestCol_i = col_i;
						largestDotProduct = colDetails[col_i].dotProduct;
					}

				}

				// check if all columns have been used already
				if (bestCol_i == -1) {
					break;
				}

				// mark the term as used
				colDetails[bestCol_i].used = true;

        // create a new model and add the term to the model
				Model *newModel = new Model(model);
				newModel->addTerm(bestCol_i, logit);

				// add the term to the model temporarily
				//model->addTerm(bestCol_i);
				//model->leastSquares();

//				cout << "Ran least squares on:" << endl;
//				model->printModelFactors();
//				cout << endl;

				// check if the model is a duplicate
				bool isDuplicate = false;
				for (int model_i = 0; model_i < models_n; model_i++) {
					if (nextTopModels[model_i] == NULL) {
						break;
					} else if (nextTopModels[model_i]->isDuplicate(newModel, true)) {
						//cout << "Duplicate Model!!! Merged!" << endl;
						isDuplicate = true;
						break;
					}
				}

				if (!isDuplicate) {
					// find a possible next top model to replace
					if (nextTopModels[models_n - 1] == NULL ||
						nextTopModels[models_n - 1]->getRSquared() < newModel->getRSquared()) {

						// make sure we deallocate the older next top model
						if (nextTopModels[models_n - 1] != NULL) {
							delete nextTopModels[models_n - 1];
							nextTopModels[models_n - 1] = NULL;
						}

						// insert the new next top model
						nextTopModels[models_n - 1] = newModel;
                                                } else {
						delete newModel;
					}

					// perform swapping to maintain sorted list
					for (int model_i = models_n - 2; model_i >= 0; model_i--) {
						if (nextTopModels[model_i] == NULL ||
							nextTopModels[model_i]->getRSquared() < nextTopModels[model_i + 1]->getRSquared()) {

							Model *temp = nextTopModels[model_i];
							nextTopModels[model_i] = nextTopModels[model_i + 1];
							nextTopModels[model_i + 1] = temp;
						}
					}

                                    } else {
					delete newModel;
				    }
				}

				// remove the temporarily added term
				//model->removeTerm(bestCol_i);



			// delete the model from the top models queue since it has been processed
			delete model;

		}

		// copy next top models to top models
		for (int model_i = 0; model_i < models_n; model_i++) {
			topModels[model_i] = nextTopModels[model_i];
		}

		// find the top model
		if (topModels[0] != NULL) {
			cout << "Top Model (" << topModels[0]->getRSquared() << "):" << endl;
			topModels[0]->printModelFactors();
		} else {
			cout << "No Model" << endl;
		}

	}

	delete[] colDetails;

	// count occurrences
	Occurrence *occurrence = new Occurrence;
	occurrence->factorList = new int[0];
	occurrence->factorList_n = 0;
        occurrence->rSquaredContribution = 0;
	occurrence->count = 0;
	occurrence->magnitude = 0;
	list<Occurrence*> *occurrenceLists = new list<Occurrence*>[locatingArray->getT()];
	allocateOccurrences(occurrence, locatingArray->getT(), locatingArray->getFactors(), occurrenceLists);

	cout << endl;
	cout << "Final Models Ranking: " << endl;
	for (int model_i = 0; model_i < models_n; model_i++) {
		if (topModels[model_i] != NULL) {
			cout << "Model " << (model_i + 1) << " (" << topModels[model_i]->getRSquared() << "):" << endl;
			topModels[model_i]->printModelFactors();
			cout << endl;

			topModels[model_i]->countOccurrences(occurrence);
			delete topModels[model_i];
			topModels[model_i] = NULL;
		} else {
			break;
		}
	}

	// sort number of occurrences and print
	cout << "Occurrence Counts" << endl;
	for (int t = 0; t < locatingArray->getT(); t++) {

		occurrenceLists[t].sort(compareOccurrence);
		occurrenceLists[t].reverse();

                // print out the headers
		cout << right << setw(15) << "R^2 Contr." << " | " <<
		                 setw(15) << "Count" << " | " <<
				         setw(15) << "Magnitude" << " | " << setw(15) << "AvgMag." << " | " << "Factor Combination" << endl;
		cout << setw(15) << setfill('-') << "" << " | " <<
				setw(15) << setfill('-') << "" << " | " <<
				setw(15) << setfill('-') << "" << " | " <<
                                setw(15) << setfill('-') << "" << " | " <<
				setw(20) << setfill('-') << "" << setfill(' ') << endl;


		// print out each count
		for (std::list<Occurrence*>::iterator it = occurrenceLists[t].begin(); it != occurrenceLists[t].end(); it++) {

			// only print the count if it has a count
                        if ((*it)->count > 0) {
				cout << setw(15) << right << (*it)->rSquaredContribution << " | ";
				cout << setw(15) << right << (*it)->count << " | ";
				cout << setw(15) << right << (*it)->magnitude << " | ";
                                cout << setw(15) << right << (*it)->magnitude/(*it)->count << " | ";

				// print out the factor combination names
				for (int factorList_i = 0; factorList_i < (*it)->factorList_n; factorList_i++) {
					if (factorList_i != 0) cout << " & ";
					cout << locatingArray->getFactorData()->getFactorName((*it)->factorList[factorList_i]);
				}
				cout << endl;
			}

		}

		cout << endl;

	}

	deallocateOccurrences(occurrence, locatingArray->getFactors());
	delete[] occurrence->factorList;
	delete occurrence;

	delete[] occurrenceLists;
	//delete[] topModels;
	delete[] nextTopModels;

  return 0;

}

//[[Rcpp::export]]
RcppExport SEXP createModels_wrapper(SEXP la_, SEXP response_, SEXP csMatrix_, SEXP maxTerms_, SEXP models_n_, SEXP newModels_n_, SEXP logit_){
  
  //grab the objects as XPtrs to make la, response, and csMatrix
  Rcpp::XPtr<LocatingArray> la_ptr(la_);
  Rcpp::XPtr<VectorXf> response_ptr(response_);
  Rcpp::XPtr<CSMatrix> cs_ptr(csMatrix_);
  
  //convert the remaining parameters
  int maxTerms = as<int>(maxTerms_);
  int models_n = as<int>(models_n_);
  int newModels_n = as<int>(newModels_n_);
  bool logit = as<bool>(logit_);
  
  int models;
  
  //invoke the function
  models <- createModels(la_ptr, response_ptr, cs_ptr, maxTerms, models_n, newModels_n, logit);
  
  return wrap(models);
}

void reorderrows_la(LocatingArray *locatingArray, VectorXf *response, CSMatrix *csMatrix, int k, int c, string newla_path, string response_col, string response_dir){

	//set to zero as this can be done later in the model building step
	bool performLog = 0;

	//call reorder function
	csMatrix->reorderRows(k, c, &response);

	locatingArray->writeToFile(newla_path);
	response->writeToFile(response_col, response_dir);

}



int main(int argc, char **argv) {

	//LocatingArray *array = new LocatingArray("LA_SMALL.tsv");
	//FactorData *factorData = new FactorData("Factors_SMALL.tsv");
	//VectorXf *response = loadResponseVector("responses_SMALL", "Exposure", false);
	//LocatingArray *array = new LocatingArray("LA_LARGE.tsv");
	//FactorData *factorData = new FactorData("Factors_LARGE.tsv");
	//VectorXf *response = loadResponseVector("responses_LARGE", "Throughput", true);

	// ./Search LA_LARGE.tsv Factors_LARGE.tsv analysis responses_LARGE Throughput 1 13 50 50

	long long int seed = time(NULL);
	cout << "Seed:\t" << seed << endl;
	srand(seed);

	if (argc < 3) {
		cout << "Usage: " << argv[0] << " [LocatingArray.tsv] ([FactorData.tsv]) ..." << endl;
		return 0;
	}

	Noise *noise = NULL;
	LocatingArray *array = new LocatingArray(argv[1], argv[2]);

	CSMatrix *matrix = new CSMatrix(array);

	for (int arg_i = 3; arg_i < argc; arg_i++) {
		if (strcmp(argv[arg_i], "memchk") == 0) {
			int exit;
			cout << "Check memory and press ENTER" << endl;
			cin >> exit;
		} else if (strcmp(argv[arg_i], "analysis") == 0) {
			if (arg_i + 6 < argc) {
				bool performLog = atoi(argv[arg_i + 3]);
				int terms_n = atoi(argv[arg_i + 4]);
				int models_n = atoi(argv[arg_i + 5]);
				int newModels_n = atoi(argv[arg_i + 6]);

				VectorXf *response = new VectorXf(array->getTests());
				loadResponseVector(response, argv[arg_i + 1], argv[arg_i + 2], performLog);
				cout << "Response range: " << response->getData()[0] << " to " << response->getData()[response->getLength() - 1] << endl;

				createModels(array, response, matrix, terms_n, models_n, newModels_n, false);
				delete response;

				arg_i += 6;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [ResponsesDirectory] [response_column] [1/0 - perform log on responses] [nTerms] [nModels] [nNewModels]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "autofind") == 0) {
			if (arg_i + 3 < argc) {
				int k = atoi(argv[arg_i + 1]);
				int c = atoi(argv[arg_i + 2]);
				int startRows = atoi(argv[arg_i + 3]);

				matrix->autoFindRows(k, c, startRows);

				arg_i += 3;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [k Separation] [c Minimum Count] [Start Rows]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "fixla") == 0) {
			if (arg_i + 1 < argc) {
				matrix->exactFix();
				array->writeToFile(argv[arg_i + 1]);

				arg_i += 1;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [FixedOutputLA.tsv]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "model") == 0) {
			if (arg_i + 3 < argc) {
				string responseDir = string(argv[arg_i + 1]);
				string responseCol = string(argv[arg_i + 2]);
				int terms = atoi(argv[arg_i + 3]);

				float *coefficients = new float[terms];
				int *columns = new int[terms];

				arg_i += 3;
				for (int term_i = 0; term_i < terms; term_i++) {
					if (arg_i + 2 < argc) {
						coefficients[term_i] = atof(argv[arg_i + 1]);
						columns[term_i] = atoi(argv[arg_i + 2]);
					} else {
						cout << "Terms are missing" << endl;
						coefficients[term_i] = 0;
						columns[term_i] = 0;
					}
					arg_i += 2;
				}

				matrix->writeResponse(responseDir, responseCol, terms, coefficients, columns);

				delete [] coefficients;
				delete [] columns;

			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [ResponsesDirectory] [response_column] [Terms] [Term 0 coefficient] [Term 0 column] ..." << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "mtfixla") == 0) {
			if (arg_i + 4 < argc) {
				int k = atoi(argv[arg_i + 1]);
				int c = atoi(argv[arg_i + 2]);
				int totalRows = atoi(argv[arg_i + 3]);

				matrix->randomFix(k, c, totalRows);
				array->writeToFile(argv[arg_i + 4]);

				arg_i += 4;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [k Separation] [c Minimum Count] [Total Rows] [FixedOutputLA.tsv]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "sysfixla") == 0) {
			if (arg_i + 5 < argc) {
				int k = atoi(argv[arg_i + 1]);
				int c = atoi(argv[arg_i + 2]);
				int initialRows = atoi(argv[arg_i + 3]);
				int minChunk = atoi(argv[arg_i + 4]);

				matrix->systematicRandomFix(k, c, initialRows, minChunk);
				array->writeToFile(argv[arg_i + 5]);

				arg_i += 5;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [k Separation] [c Minimum Count] [Initial Rows] [Minimum Chunk] [FixedOutputLA.tsv]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "checkla") == 0) {
			if (arg_i + 2 < argc) {
				int k = atoi(argv[arg_i + 1]);
				int c = atoi(argv[arg_i + 2]);

				matrix->performCheck(k, c);

				arg_i += 2;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [k Separation] [c Minimum Count]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "noise") == 0) {
			if (arg_i + 1 < argc) {
				float ratio = atof(argv[arg_i + 1]);

				noise = new Noise(ratio);

				arg_i += 1;
			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [ratio]" << endl;
				arg_i = argc;
			}
		} else if (strcmp(argv[arg_i], "printcs") == 0) {
			cout << "CS Matrix:" << endl;
			matrix->print();
		} else if (strcmp(argv[arg_i], "reorderrowsla") == 0) {
			if (arg_i + 3 < argc && argc != 9) {

				cout << "arg_i: " << arg_i << endl;
				cout << "argc: " << argc << endl;

				int k = atoi(argv[arg_i + 1]);
				int c = atoi(argv[arg_i + 2]);

				matrix->reorderRows(k, c);

				//LEAH: check this function
				array->writeToFile(argv[arg_i + 3]);

				arg_i += 3;
			} else if (arg_i + 5 < argc){

				cout << "ENTERS CORRECT IF!" << endl;

				bool performLog = 0;

				VectorXf *response = new VectorXf(array->getTests());
				loadResponseVector(response, argv[arg_i + 1], argv[arg_i + 2], performLog);
				cout << "Response range: " << response->getData()[0] << " to " << response->getData()[response->getLength() - 1] << endl;


				int k = atoi(argv[arg_i + 3]);
				int c = atoi(argv[arg_i + 4]);

				matrix->reorderRows(k, c, &response);

				array->writeToFile(argv[arg_i + 5]);
				response->writeToFile(argv[arg_i + 1], argv[arg_i + 2]);

				arg_i += 5;

			} else {
				cout << "Usage: ... " << argv[arg_i];
				cout << " [k Separation] [c Minimum Count] [ReorderedOutputLA.tsv]" << endl;
				cout << " or " << endl;
				cout << "[ResponsesDirectory] [response_column] [k Separation] [c Minimum Count] [ReorderedOutputLA.tsv] " <<endl;
				arg_i = argc;
			}
		}
	}

	cout << endl;
	cout << "Other Stuff:" << endl;
	cout << "Columns in CSMatrix: " << matrix->getCols() << endl;

	delete matrix;
	delete array;

	return 0;

}
