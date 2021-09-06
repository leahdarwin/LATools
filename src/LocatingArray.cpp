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
#include "Model.h"

using namespace Rcpp;


LocatingArray::LocatingArray(int factors, int *levelCounts) {
	tests = 0;

	this->factors = factors;

	t = 1;

	factorGrouping = new GroupingInfo*[factors];

	// load the level counts and grouping for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		factorGrouping[factor_i] = new GroupingInfo();
		factorGrouping[factor_i]->levels = levelCounts[factor_i];

		factorGrouping[factor_i]->grouped = false;
		factorGrouping[factor_i]->levelGroups = NULL;

		factorGrouping[factor_i]->conGroup = NULL;
		factorGrouping[factor_i]->conGroupIndex = -1;
	}

	nConGroups = 0;
	conGroups = new ConstraintGroup*[nConGroups];

	factorData = NULL;

}

LocatingArray::LocatingArray(string file, string factorDataFile) {
	int tempData;
	string tempString;

	ifstream ifs(file.c_str(), ifstream::in);

	// verify version
	ifs >> tempString;
	if (tempString != laVersion) {
		cout << "LA version must be " << laVersion << " for this software" << endl;
		cout << "But instead found: " << tempString << endl;
		cout << "Exiting" << endl;
		exit(0);
	}

	// read the first 2 lines of locating array
	ifs >> tests;
	ifs >> factors;

	t = 2;

	cout << tests << ", " << factors << endl;
	factorGrouping = new GroupingInfo*[factors];

	// load the level counts for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		factorGrouping[factor_i] = new GroupingInfo();

		ifs >> tempData;
		factorGrouping[factor_i]->levels = tempData;
	}

	// load the grouping for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		// load grouped bool
		ifs >> tempData;
		factorGrouping[factor_i]->grouped = tempData;

		// check grouping
		if (factorGrouping[factor_i]->grouped) {
			// allocate array for level grouping
			factorGrouping[factor_i]->levelGroups = new char[factorGrouping[factor_i]->levels];

			for (int level_i = 0; level_i < factorGrouping[factor_i]->levels; level_i++) {
				// load level group
				ifs >> tempData;
				factorGrouping[factor_i]->levelGroups[level_i] = tempData;
			}
		} else {
			factorGrouping[factor_i]->levelGroups = NULL;
		}

		factorGrouping[factor_i]->conGroup = NULL;
		factorGrouping[factor_i]->conGroupIndex = -1;
	}

	// load factor data (must be done before constraint groups)
	if (factorDataFile == "") {
		factorData = new FactorData(getGroupingInfo(), getFactors());
	} else {
		factorData = new FactorData(factorDataFile);
	}

	// load constraint groups
	ifs >> nConGroups;
	conGroups = new ConstraintGroup*[nConGroups];

	// load each constraint group
	for (int iConGroup = 0; iConGroup < nConGroups; iConGroup++) {
		conGroups[iConGroup] = new ConstraintGroup(this, ifs);
	}

	// load the tests now
	for (int test_i = 0; test_i < tests; test_i++) {
		char *levelRow = new char[factors];

		// now load each factor level for this specific test
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			ifs >> tempData;
			levelRow[factor_i] = tempData;
		}

		addLevelRow(levelRow);
		tests--;
	}

	ifs.close();



}





void LocatingArray::addLevelRow(char *levelRow) {
	levels.push_back(levelRow);
	tests++;
}

char *LocatingArray::remLevelRow() {
	char *levelRow = levels.back();
	levels.pop_back();
	tests--;
	return levelRow;
}

GroupingInfo **LocatingArray::getGroupingInfo() {
	return factorGrouping;
}

char **LocatingArray::getLevelMatrix() {
	return &levels[0];
}

int LocatingArray::getFactors() {
	return factors;
}

int LocatingArray::getTests() {
	return tests;
}

int LocatingArray::getT() {
	return t;
}

int LocatingArray::getNConGroups() {
	return nConGroups;
}

ConstraintGroup **LocatingArray::getConGroups() {
	return conGroups;
}

void LocatingArray::writeToFile(string file) {

	cout << "Writing to " << file << "...";

	ofstream ofs(file.c_str());

	// write version
	ofs << laVersion << endl;

	// initial tests and factors
	ofs << tests << "\t" << factors << endl;

	// write the level counts for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		ofs << (int)factorGrouping[factor_i]->levels << "\t";
	}
	ofs << endl;

	// write the grouping for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		// load grouped bool
		ofs << factorGrouping[factor_i]->grouped << "\t";

		// check grouping
		if (factorGrouping[factor_i]->grouped) {

			for (int level_i = 0; level_i < factorGrouping[factor_i]->levels; level_i++) {
				// write level group
				ofs << (int)factorGrouping[factor_i]->levelGroups[level_i] << "\t";
			}
		}
		ofs << endl;
	}

	// write constraint groups
	ofs << nConGroups << endl;
	for (int iConGroup = 0; iConGroup < nConGroups; iConGroup++) {
		conGroups[iConGroup]->writeToStream(ofs);
	}

	// write the tests now
	char **levelMatrix = getLevelMatrix();
	for (int test_i = 0; test_i < tests; test_i++) {
		// now write each factor level for this specific test
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			ofs << (int)levelMatrix[test_i][factor_i] << "\t";
		}
		ofs << endl;
	}

	ofs.close();

}

FactorData *LocatingArray::getFactorData() {
	return factorData;
}


LocatingArray::~LocatingArray() {
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		if (factorGrouping[factor_i]->grouped) {
			delete[] factorGrouping[factor_i]->levelGroups;
		}

		delete factorGrouping[factor_i];
	}
	delete[] factorGrouping;

	// deallocate the factor data if it exists
	if (factorData != NULL) {
		delete factorData;
	}

	for (int iConGroup = 0; iConGroup < nConGroups; iConGroup++) {
		delete conGroups[iConGroup];
	}
	delete[] conGroups;

	while (tests > 0) {
		char *levelRow = remLevelRow();

		delete[] levelRow;
	}
}

// [[Rcpp::export]]
SEXP makeLA2(std::string file, std::string factorDataFile){
  LocatingArray* la = new LocatingArray(file, factorDataFile);
  Rcpp::XPtr<LocatingArray> ptr(la);
  return ptr;
}

