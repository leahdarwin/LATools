#include <iostream>
#include <iomanip>

#include "LocatingArray.h"
#include "Model.h"
#include <math.h>

using namespace std;


Model::Model(VectorXf *response, int maxTerms, CSMatrix *csMatrix) {
	this->response = response;
	this->maxTerms = maxTerms;
	this->csMatrix = csMatrix;
        this->occurrences = 1;
	this->tests = response->getLength();

	// allocate memory for the coefficients vector
	coefVec = new float[maxTerms];

	// allocate memory for the residuals vector
	resiVec = new float[tests];

	// allocate memory for the model response
	modelResponse = new float[tests];

	// add the intercept
	terms = 1;
	hTermIndex = new TermIndex;
	hTermIndex->termIndex = 0;
	hTermIndex->next = NULL;

	// this is a one line way to get the same intercept above
	leastSquares();

        // this should always be 0
	hTermIndex->rSquared = this->rSquared;

}



Model::Model(Model *model) {
	this->response = model->response;
	this->maxTerms = model->maxTerms;
	this->csMatrix = model->csMatrix;
	this->rSquared = model->rSquared;
	this->terms = model->terms;
        this->occurrences = model->occurrences;
	this->tests = response->getLength();

	// allocate memory for the coefficients vector and copy
	coefVec = new float[maxTerms];
	for (int coef_i = 0; coef_i < maxTerms; coef_i++) {
		coefVec[coef_i] = model->coefVec[coef_i];
	}

	// allocate memory for the residuals vector
	resiVec = new float[tests];
	for (int resi_i = 0; resi_i < tests; resi_i++) {
		resiVec[resi_i] = model->resiVec[resi_i];
	}

	// allocate memory for the model response
	modelResponse = new float[tests];
	for (int resp_i = 0; resp_i < tests; resp_i++) {
		modelResponse[resp_i] = model->modelResponse[resp_i];
	}

	// copy term index list
	TermIndex **destTermIndex = &hTermIndex;
	for (TermIndex *pTermIndex = model->hTermIndex; pTermIndex != NULL;
			pTermIndex = pTermIndex->next) {
		(*destTermIndex) = new TermIndex;
		(*destTermIndex)->termIndex = pTermIndex->termIndex;
                (*destTermIndex)->rSquared = pTermIndex->rSquared;
		(*destTermIndex)->next = NULL;
		destTermIndex = &(*destTermIndex)->next;
	}
}

void Model::printModelFactors() {
	int factor1i, factor2i, level1, level2;

	// print out the number of terms
	cout << "Model with " << terms << " terms" << endl;

        // print out the headers
	cout << setw(15) << right << "R^2 Contr." << " | " << setw(15) << right << "Coefficient" << " | " << "Term" << endl;
	cout << setw(15) << setfill('-') << "" << " | " <<
			setw(15) << setfill('-') << "" << " | " <<
			setw(15) << setfill('-') << "" << setfill(' ') << endl;

	// print out each term
	int term_i = 0;
	for (TermIndex *pTermIndex = hTermIndex; pTermIndex != NULL; pTermIndex = pTermIndex->next) {
		// print out the coefficient
		int termIndex = pTermIndex->termIndex;
                cout << setw(15) << right << pTermIndex->rSquared << " | ";
		cout << setw(15) << right << coefVec[term_i++] << " | ";

		// print out the factor names and level names
		cout << csMatrix->getColName(csMatrix->getCol(termIndex)) << endl;

	}

	// calculate adjusted r-squared (terms - 2 to not include the intercept)
	float adjustedRSquared = 1 - (1 - rSquared) * (tests - 1) / (tests - terms);

	// print model r-squared
	if (rSquared == 1) cout << "Perfect Model!!!" << endl;
        cout << "Occurrences: " << occurrences << endl;
	cout << "R-Squared: " << rSquared << endl;
	cout << "Adjusted R-Squared: " << adjustedRSquared << endl;
}


void Model::logisticRegression(){

  int NUM_ITERATIONS = 1000;

  float weights[terms];
  float scores[tests];
  float gradients[terms];


  //init w zeros
  for(int i = 0; i<terms; i++){
    weights[i] = 0;
    gradients[i] = 0;
  }

  for(int i = 0; i<tests; i++){
    scores[i]=0;
  }


  // used when accessing CS Matrix
  CSCol *csCol;

  // used when looping through term indices
  TermIndex *pTermIndex;

  //start at the head of index list
  pTermIndex = hTermIndex;

  //get response vector
  float *responseData = response->getData();

  long totalCorrect = 0;

  for(int i = 0; i<NUM_ITERATIONS; i++){

    // (scores <- trainingSetLabels - sigmoid(trainingSetFeatures%*%weights))
    for (int row_i=0; row_i < tests; row_i++){

      //start at the head of index list
      pTermIndex = hTermIndex;

      float score = 0.0f;

      for (int col_i=0; col_i < terms; col_i++){

        csCol = csMatrix->getCol(pTermIndex->termIndex);

        score += csCol->dataVector[row_i] * weights[col_i];

        pTermIndex = pTermIndex->next;

      }

      scores[row_i] = responseData[row_i] - 1/(1+exp(-score));

    }

    //start at the head of index list
    pTermIndex = hTermIndex;

    // loop 2: (gradients <- t(trainingSetFeatures)%*%scores)
    for (int col_i=0; col_i < terms; col_i++){

      csCol = csMatrix->getCol(pTermIndex->termIndex);

      gradients[col_i] = 0.0f;

      for(int row_i=0; row_i < tests; row_i++){

        gradients[col_i] += csCol->dataVector[row_i] * scores[row_i];

      }

      //next column
      pTermIndex = pTermIndex->next;

    }

    //  loop 3: (weights <- weights*(1-LEARNING_RATE*REGULARIZATION) + LEARNING_RATE*gradients)
    for (int col_i=0; col_i < terms; col_i++){

      weights[col_i] += gradients[col_i]*0.1;

    }
  }

  //fill coeff vector
  for (int i = 1; i<terms; i++){
    coefVec[i]=weights[i];
  }


  for (int row_i = 0; row_i < tests; row_i++){

    double predictedValue, z=0;

    //start at the head of index list
    pTermIndex = hTermIndex;

    for (int col_i=0; col_i < terms; col_i++){

      //get col
      csCol = csMatrix->getCol(pTermIndex->termIndex);

      z += csCol->dataVector[row_i] * weights[col_i];

      //next col
      pTermIndex = pTermIndex->next;

    }

    //calculate predicted val
    predictedValue =  1/(1 + exp(-1 * z));

    if(predictedValue < 0.5) {
      predictedValue = 0;
    }
    else {
      predictedValue = 1;
    }

    if(predictedValue == responseData[row_i]) {
      totalCorrect++;
    }

  }

  rSquared =totalCorrect/double(tests);

}

void Model::leastSquares() {

	// used when accessing CS Matrix
	CSCol *csCol;

	// used when looping through term indices
	TermIndex *pTermIndex;

	// check dimensions
	if (terms > tests) {
		cout << "Cols (terms) cannot be more than rows (tests) to perform QR" << endl;
		return;
	}

	// do QR here (column by column)
	pTermIndex = hTermIndex;
	for (int col_i = 0; col_i < terms; col_i++) {

		csCol = csMatrix->getCol(pTermIndex->termIndex);
		// assign initial column A[col_i] to work vector
		for (int row_i = 0; row_i < tests; row_i++) {
			workSpace->workVec[row_i] = csCol->dataP[row_i];
		}

		// subtract appropriate other vectors
		float dotProd;
		for (int row_i = 0; row_i < col_i; row_i++) {
			// find the dot product of A[:][col_i] and Q[:][row_i]
			dotProd = 0;
			for (int dotrow_i = 0; dotrow_i < tests; dotrow_i++)
				dotProd += csCol->dataP[dotrow_i] * workSpace->dataQ[dotrow_i][row_i];

			// assign the dot product to the R matrix
			workSpace->dataR[row_i][col_i] = dotProd;

			// perform the necessary subtraction
			for (int subrow_i = 0; subrow_i < tests; subrow_i++)
				workSpace->workVec[subrow_i] -= dotProd * workSpace->dataQ[subrow_i][row_i];

		}

		// find the dot product of workVec and workVec
		dotProd = 0;
		for (int dotrow_i = 0; dotrow_i < tests; dotrow_i++)
			dotProd += workSpace->workVec[dotrow_i] * workSpace->workVec[dotrow_i];
		dotProd = sqrt(dotProd);

		// assign the dot product to the R matrix
		workSpace->dataR[col_i][col_i] = dotProd;

		// assign the normalized column to Q
		if (dotProd == 0) {
			for (int row_i = 0; row_i < tests; row_i++)
				workSpace->dataQ[row_i][col_i] = 0;
		} else {
			for (int row_i = 0; row_i < tests; row_i++)
				workSpace->dataQ[row_i][col_i] = workSpace->workVec[row_i] / dotProd;
		}

		pTermIndex = pTermIndex->next;
	}

	/*
	cout << "Q matrix" << endl;
	for (int row_i = 0; row_i < tests; row_i++) {
		for (int col_i = 0; col_i < terms; col_i++) {
			cout << workSpace->dataQ[row_i][col_i] << "\t";
		}
		cout << endl;
	}

	cout << "R matrix" << endl;
	for (int row_i = 0; row_i < terms; row_i++) {
		for (int col_i = 0; col_i < terms; col_i++) {
			cout << workSpace->dataR[row_i][col_i] << "\t";
		}
		cout << endl;
	}
	*/

	/* m >> n
	** A is (m by n)
	** Q is (m by n)
	** R is (n by n)
	** b is (m by 1)
	** Q-transpose is (n by m)
	**
	*/
	// We now have QRx = b

	// perform workVec = Q-transpose * b

	// first zero out the work vector up until n (or cols)
	for (int row_i = 0; row_i < terms; row_i++)
		workSpace->workVec[row_i] = 0;

	// we go column 1st because Q is transposed
	float *responseData = response->getData();
	for (int col_i = 0; col_i < tests; col_i++) {
		for (int row_i = 0; row_i < terms; row_i++) {
			// row and col are swapped because of the transpose
			workSpace->workVec[row_i] += workSpace->dataQ[col_i][row_i] * responseData[col_i];
		}
	}

	// We now have Rx = workVec but R is upper triangular (n by n)
	float rowSolution;
	for (int row_i = terms - 1; row_i >= 0; row_i--) {

		// initialize row solution
		rowSolution = workSpace->workVec[row_i];

		// subtract other parts of LHS of equation
		for (int col_i = terms - 1; col_i > row_i; col_i--) {
			rowSolution -= workSpace->dataR[row_i][col_i] * coefVec[col_i];
		}

		// divide for final row solution
		if (workSpace->dataR[row_i][row_i] == 0) {
			coefVec[row_i] = 0;
		} else {
			coefVec[row_i] = rowSolution / workSpace->dataR[row_i][row_i];
		}

	}

	// coefVec holds the least squares coefficients
	//cout << "coefficient vector" << endl;
	//for (int term_i = 0; term_i < terms; term_i++) {
	//	cout << coefVec[term_i] << " ";
	//}
	//cout << endl;

	// calculate residuals and r-squared
	float SSres = 0;
	for (int row_i = 0; row_i < tests; row_i++) {
		modelResponse[row_i] = 0;
	}

	// find model response
	pTermIndex = hTermIndex;
	for (int term_i = 0; term_i < terms; term_i++) {
		csCol = csMatrix->getCol(pTermIndex->termIndex);
		for (int row_i = 0; row_i < tests; row_i++) {
			modelResponse[row_i] += csCol->dataP[row_i] * coefVec[term_i];
		}
		pTermIndex = pTermIndex->next;
	}

	// find residuals and SSres
	for (int row_i = 0; row_i < tests; row_i++) {
		resiVec[row_i] = responseData[row_i] - modelResponse[row_i];

		// add the square of residual to SSres
		SSres += resiVec[row_i] * resiVec[row_i];
	}

	// find r-squared
	rSquared = 1 - (SSres / response->getSStot());

}

// get the residuals vector for this model
float *Model::getResiVec() {
	return resiVec;
}

float Model::getRSquared() {
	return rSquared;
}

// check if term is part of the model
bool Model::termExists(int col_i) {
	for (TermIndex *pTermIndex = hTermIndex; pTermIndex != NULL; pTermIndex = pTermIndex->next) {
		if (pTermIndex->termIndex == col_i) {
			// term was found!
			return true;
		}
	}
	// term never found
	return false;
}

// add a term (col_i of csMatrix) to the model. Returns true if added successfully
bool Model::addTerm(int col_i, bool logit) {
	TermIndex *termIndex = NULL;

	TermIndex **pTermIndex = &hTermIndex;
	for (; ; pTermIndex = &(*pTermIndex)->next) {

		if (*pTermIndex == NULL || (*pTermIndex)->termIndex > col_i) {
			// insert within (or at the end of) the list
			termIndex = new TermIndex;
			termIndex->termIndex = col_i;
			termIndex->next = (*pTermIndex);
			(*pTermIndex) = termIndex;
			terms++;

			// calculate r-squared contribution for this term
			float oldRSquared = this->rSquared;
			
			if(logit)
			  logisticRegression();
			else
			  leastSquares();

			termIndex->rSquared = this->occurrences * (this->rSquared - oldRSquared);

			return true;
		} else if ((*pTermIndex)->termIndex == col_i) {
			// already exists so cannot insert
			return false;
		}

	}
}

// remove a term from the model
bool Model::removeTerm(int col_i) {
	// loop through the term index list
	for (TermIndex **pTermIndex = &hTermIndex; *pTermIndex != NULL; pTermIndex = &(*pTermIndex)->next) {
		// check if it equals the term to be removed (col_i)
		if ((*pTermIndex)->termIndex == col_i) {
			// take the current term index out of the list
			TermIndex *removed = (*pTermIndex);
			(*pTermIndex) = removed->next;

			// decrease the number of terms
			terms--;

			// free its memory
			delete removed;

			return true;
		}
	}

	// term to remove was not found
	return false;
}

int Model::getTerms() {
	return terms;
}

// static
void Model::setupWorkSpace(int rows, int cols) {
	// allocate initial memory
	workSpace = new WorkSpace;

	// allocate memory for Q
	workSpace->dataQ = new float*[rows];
	for (int row_i = 0; row_i < rows; row_i++)
		workSpace->dataQ[row_i] = new float[cols];

	// allocate memory for R
	workSpace->dataR = new float*[cols];
	for (int row_i = 0; row_i < cols; row_i++)
		workSpace->dataR[row_i] = new float[cols];

	// allocate memory for work vector
	workSpace->workVec = new float[rows];
}

void Model::countOccurrences(Occurrence *occurrence) {
	// count occurrences for each term
	int term_i = 0;
	for (TermIndex *pTermIndex = hTermIndex; pTermIndex != NULL; pTermIndex = pTermIndex->next) {
		int termIndex = pTermIndex->termIndex;
                float rSquared = pTermIndex->rSquared;
		float magnitude = coefVec[term_i++];
		csMatrix->countOccurrences(csMatrix->getCol(termIndex), occurrence, 0, magnitude, rSquared);
	}
}

bool Model::isDuplicate(Model *model, bool merge) {
	// check if all terms match
	TermIndex *p1TermIndex = hTermIndex;
	TermIndex *p2TermIndex = model->hTermIndex;

	while (!(p1TermIndex == NULL && p2TermIndex == NULL)) {
		if (p1TermIndex->termIndex != p2TermIndex->termIndex) return false;

		p1TermIndex = p1TermIndex->next;
		p2TermIndex = p2TermIndex->next;
	}

        // now merge the duplicate models if requested
	if (merge) {
		p2TermIndex = model->hTermIndex;
		for (p1TermIndex = hTermIndex; p1TermIndex != NULL; p1TermIndex = p1TermIndex->next) {
			p1TermIndex->rSquared += p2TermIndex->rSquared;

			p2TermIndex = p2TermIndex->next;
		}

		this->occurrences += model->occurrences;
	}

	return true;
}

Model::~Model() {

	// delete term index list
	while (hTermIndex != NULL) {
		TermIndex *removed = hTermIndex;
		hTermIndex = hTermIndex->next;
		delete removed;
	}

	// delete coefficients vector
	delete[] coefVec;

	// delete residuals vector
	delete[] resiVec;

	// delete model response
	delete[] modelResponse;

}
