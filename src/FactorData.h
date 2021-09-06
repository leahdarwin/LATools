#ifndef FACTORDATA_H
#define FACTORDATA_H

#include <string>

#include "LocatingArray.h"

using namespace std;

struct GroupingInfo;

class Factor {
public:

	// name of the factor
	string name;
	
	// level count
	int levels;
	
	// is numeric
	bool numeric;
	
	// names of the levels
	string *levelNames;
	
	// values of the levels (NULL if numeric is false)
	float *levelValues;
	
};

class FactorData {
private:
	// total factors (75-100ish for large dataset)
	int factorCount;
	
	// array of factor pointers
	Factor **factors;
	
public:
	FactorData(GroupingInfo **arrayInfo, int factorCount);
	FactorData(string file);
	
	Factor *getFactor(int factor_i);
	
	string getFactorName(int factor_i);
	
	string getFactorLevelName(int factor_i, int level_i);
	
	float getNumericFactorLevel(int factor_i, int level_i);
	
	~FactorData();
};

#endif
