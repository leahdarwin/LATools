#ifndef LOCATINGARRAY_H
#define LOCATINGARRAY_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Rcpp.h>
#include <RcppCommon.h>
#include "ConstraintGroup.h"
#include "FactorData.h"
#include "VectorXf.h"

using namespace std;

const string laVersion = "v2.0";

class ConstraintGroup;
class FactorData;

struct GroupingInfo {
	int levels;	// total levels for this factor
	
	bool grouped;
	char *levelGroups;
	
	ConstraintGroup *conGroup; // constraint group to which this factor belongs
	int conGroupIndex; // this factor's index in the constraint group above
};

class LocatingArray {
private:
	GroupingInfo **factorGrouping;
	
	vector <char*>levels;	// pointer to array of test levels (the main locating array)
	
	int tests;		// count of tests in locating array
	int factors;	// count of factors in locating array
	
	int t;			// covers t-way interactions
	
	int nConGroups;
	ConstraintGroup **conGroups;
	
	FactorData *factorData;
public:
	LocatingArray(int factors, int *levelCounts);
	LocatingArray(string file, string factorDataFile);
	
	LocatingArray makeLA(string file, string factorDataFile);
	
	void addLevelRow(char *levelRow);
	char *remLevelRow();
	
	GroupingInfo **getGroupingInfo();
	
	char **getLevelMatrix();
	
	int getFactors();
	int getTests();
	
	int getT();
	
	int getNConGroups();
	ConstraintGroup **getConGroups();
	
	FactorData *getFactorData();
	
	void writeToFile(string file);
	
	~LocatingArray();
};



#endif