
#include <fstream>
#include <sstream>

#include "FactorData.h"

using namespace std;

FactorData::FactorData(string file) {
	ifstream ifs(file.c_str(), ifstream::in);
	
	// grab the total number of factors
	ifs >> factorCount;
	
	// allocate necessary memory
	factors = new Factor*[factorCount];
	
	// read in line by line
	for (int factor_i = 0; factor_i < factorCount; factor_i++) {
		// allocate the factor
		Factor *factor = new Factor;
		
		// read the factor name
		ifs >> factor->name;
		
		// read the number of levels
		ifs >> factor->levels;
		
		// read whether the factor has numeric levels
		ifs >> factor->numeric;
		
		// allocate memory for each level
		factor->levelNames = new string[factor->levels];
		
		// read individual level names
		for (int level_i = 0; level_i < factor->levels; level_i++) {
			ifs >> factor->levelNames[level_i];
		}
		
		// read the level values if they are numeric
		if (factor->numeric) {
			// allocate memory for each level value
			factor->levelValues = new float[factor->levels];
			
			// read individual level values
			for (int level_i = 0; level_i < factor->levels; level_i++) {
				ifs >> factor->levelValues[level_i];
			}
		} else {
			factor->levelValues = NULL;
		}
		
		factors[factor_i] = factor;
	}
}

FactorData::FactorData(GroupingInfo **arrayInfo, int factorCount) {
	// grab the total number of factors
	this->factorCount = factorCount;
	
	// allocate necessary memory
	factors = new Factor*[factorCount];
	
	// read in line by line
	for (int factor_i = 0; factor_i < factorCount; factor_i++) {
		// allocate the factor
		Factor *factor = new Factor;
		
		// read the factor name
		ostringstream oss;
		oss << "F" << factor_i;
		factor->name = oss.str();
		
		// read the number of levels
		factor->levels = arrayInfo[factor_i]->levels;
		
		// read whether the factor has numeric levels
		factor->numeric = false;
		
		// allocate memory for each level
		factor->levelNames = new string[factor->levels];
		
		// read individual level names
		for (int level_i = 0; level_i < factor->levels; level_i++) {
			ostringstream oss;
			oss << "L" << level_i;
			factor->levelNames[level_i] = oss.str();
		}
		
		factor->levelValues = NULL;
		
		factors[factor_i] = factor;
	}
}

Factor *FactorData::getFactor(int factor_i) {
	return factors[factor_i];
}

string FactorData::getFactorName(int factor_i) {
	return factors[factor_i]->name;
}

string FactorData::getFactorLevelName(int factor_i, int level_i) {
	return factors[factor_i]->levelNames[level_i];
}

float FactorData::getNumericFactorLevel(int factor_i, int level_i) {
	if (factors[factor_i]->numeric) {
		return factors[factor_i]->levelValues[level_i];
	} else {
		return 0;
	}
}

FactorData::~FactorData() {
	// read in line by line
	for (int factor_i = 0; factor_i < factorCount; factor_i++) {
		Factor *factor = factors[factor_i];
		
		delete[] factor->levelNames;
		
		if (factor->numeric) {
			delete[] factor->levelValues;
		}
		
		delete factor;
	}
	
	delete[] factors;
}
