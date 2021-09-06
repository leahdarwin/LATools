#include "ConstraintGroup.h"

ConstraintGroup::ConstraintGroup(LocatingArray *array, ifstream &ifs) {
	
	// load constraint group factors
	ifs >> factors;
	factorIndeces = new int[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		ifs >> factorIndeces[factor_i];
	}
	
	// allocate the level counts for the factors
	int *levelCounts = new int[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		levelCounts[factor_i] = array->getGroupingInfo()[factorIndeces[factor_i]]->levels;
	}
	
	// create constraint group LA
	groupLA = new LocatingArray(factors, levelCounts);
	
	ifs >> constraints;
	
	// allocate memory for BoolResult objects
	boolConstraints = new BoolResult*[constraints];
	
	for (int constraint_i = 0; constraint_i < constraints; constraint_i++) {
		boolConstraints[constraint_i] = BoolResult::readBoolResult(array, ifs);
	}
	
	// link to factors in array
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		array->getGroupingInfo()[factorIndeces[factor_i]]->conGroup = this;
		array->getGroupingInfo()[factorIndeces[factor_i]]->conGroupIndex = factor_i;
	}
	
	// populate groupLA
	char *levelRow = new char[array->getFactors()];
	char *factorLevels = new char[factors];
	
	// loop through all locating array factors and set all to 0
	for (int factor_i = 0; factor_i < array->getFactors(); factor_i++) {
		levelRow[factor_i] = 0;
	}
	
	// add row of 0s to the array for use in recursive function
	array->addLevelRow(levelRow);
	
	// count of each level occurrence (used for weighting)
	int **settingCount = new int*[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		settingCount[factor_i] = new int[levelCounts[factor_i]];
		for (int level_i = 0; level_i < levelCounts[factor_i]; level_i++) {
			settingCount[factor_i][level_i] = 0;
		}
	}
	
	// count of full-factorial combinations (used for weighting)
	int fullFactorialCount = 0;
	
	// call recursive populator
	populateGroupLA(array, levelRow, factorLevels, 0, settingCount, fullFactorialCount);
	
	// assign row weights based on groupLA entries
	char **groupLevelMatrix = groupLA->getLevelMatrix();
	float *rowWeights = new float[groupLA->getTests()];
	weightMin = new int[groupLA->getTests()];
	weightMax = new int[groupLA->getTests()];
	float prevWeight = 0;
	
	for (int row_i = 0; row_i < groupLA->getTests(); row_i++) {
		rowWeights[row_i] = prevWeight;
		for (int col_i = 0; col_i < factors; col_i++) {
			rowWeights[row_i] += (float)fullFactorialCount / settingCount[col_i][groupLevelMatrix[row_i][col_i]];
//			cout << (int)groupLevelMatrix[row_i][col_i] << "\t";
		}
		
		weightMin[row_i] = (int)prevWeight;
		weightMax[row_i] = (int)rowWeights[row_i] - 1;
		
//		cout << rowWeights[row_i] << "\t" << rowWeights[row_i] - prevWeight << "\t" << weightMin[row_i] << "\t" << weightMax[row_i] << endl;
		
		prevWeight = rowWeights[row_i];
	}
	
	weightRandMax = (int)prevWeight;
	
	delete[] rowWeights;
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		delete[] settingCount[factor_i];
	}
	delete[] settingCount;
	delete[] levelCounts;
	levelRow = array->remLevelRow();
	delete[] levelRow;
	delete[] factorLevels;
}

void ConstraintGroup::populateGroupLA(LocatingArray *array, char *levelRow, char *factorLevels, int fixedFactors,
		int **settingCount, int &fullFactorialCount) {
	if (fixedFactors == factors) {
//		for (int factor_i = 0; factor_i < factors; factor_i++) {
//			cout << (int)factorLevels[factor_i] << "\t";
//		}
//		cout << endl;
		
		// increment full-factorial count
		fullFactorialCount++;
		
		// loop through all factors in constraint group
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			levelRow[factorIndeces[factor_i]] = factorLevels[factor_i];
		}
		
		// add factorLevels to groupLA if group constraints are satisfied
		if (getResult(0)) {
			char *groupLevelRow = new char[factors];
			for (int factor_i = 0; factor_i < factors; factor_i++) {
				groupLevelRow[factor_i] = factorLevels[factor_i];
				
				// increment the setting count
				settingCount[factor_i][factorLevels[factor_i]]++;
			}
			groupLA->addLevelRow(groupLevelRow);
		}
	} else {
		// call recursively for full-factorial design
		for (int level_i = 0; level_i < array->getGroupingInfo()[factorIndeces[fixedFactors]]->levels; level_i++) {
			factorLevels[fixedFactors] = level_i;
			populateGroupLA(array, levelRow, factorLevels, fixedFactors + 1, settingCount, fullFactorialCount);
		}
	}
}

bool ConstraintGroup::satisfiableInGroupLA(char *requireLevelRow, char *avoidLevelRow) {
	char **levelMatrix = groupLA->getLevelMatrix();
	
	for (int test_i = 0; test_i < groupLA->getTests(); test_i++) {
		bool match = true;
		
		// check requireLevelRow
		for (int factor_i = 0; factor_i < groupLA->getFactors() && match; factor_i++) {
			// check if the required level disagrees with the LA row
			if (requireLevelRow[factorIndeces[factor_i]] != -1 &&
				requireLevelRow[factorIndeces[factor_i]] != levelMatrix[test_i][factor_i]) {
				match = false;
			}
		}
		
		// check avoidLevelRow
		for (int factor_i = 0; factor_i < groupLA->getFactors() && match; factor_i++) {
			// check if the avoided level disagrees with the LA row
			if (avoidLevelRow[factorIndeces[factor_i]] != -1 &&
				avoidLevelRow[factorIndeces[factor_i]] == levelMatrix[test_i][factor_i]) {
				match = false;
			}
		}
		
		// check if the row satisfies
		if (match) return true;
	}
	
	return false;
}

bool ConstraintGroup::getResult(int test) {
	for (int constraint_i = 0; constraint_i < this->constraints; constraint_i++) {
		if (!boolConstraints[constraint_i]->getResult(test)) return false;
	}
	return true;
}

void ConstraintGroup::randPopulateLevelRow(char *levelRow) {
	GroupingInfo **groupingInfo = groupLA->getGroupingInfo();
	char **levelMatrix = groupLA->getLevelMatrix();
	
	int weightRand = rand() % weightRandMax;
	
	// use binary search to find weight window with weightRand
	int botRow = 0;
	int topRow = groupLA->getTests();
	
	while (topRow != botRow) {
		int row_i = (topRow + botRow) / 2;
		
		if (weightRand <= weightMax[row_i] && weightRand >= weightMin[row_i]) {
			topRow = row_i;
			botRow = row_i;
		} else if (weightRand < weightMin[row_i]) {
			topRow = row_i - 1;
		} else if (weightRand > weightMax[row_i]) {
			botRow = row_i + 1;
		} else {
			// should never reach here
			cout << "Mistake in binary search" << endl;
		}
	}
	
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		levelRow[factorIndeces[factor_i]] = levelMatrix[topRow][factor_i];
	}
}

void ConstraintGroup::writeToStream(ofstream &ofs) {
	// write constraint group factors
	ofs << factors;
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		ofs << "\t" << factorIndeces[factor_i];
	}
	ofs << endl;
	
	// write constraints
	ofs << constraints;
	for (int constraint_i = 0; constraint_i < constraints; constraint_i++) {
		boolConstraints[constraint_i]->writeToStream(ofs);
		ofs << endl;
	}
}

ConstraintGroup::~ConstraintGroup() {
	delete[] weightMin;
	delete[] weightMax;
	
	delete[] factorIndeces;
	
	delete groupLA;
	
	for (int constraint_i = 0; constraint_i < constraints; constraint_i++) {
		delete boolConstraints[constraint_i];
	}
	delete[] boolConstraints;
}

BoolResult::BoolResult(LocatingArray *array) {
	this->array = array;
}

BoolResult *BoolResult::readBoolResult(LocatingArray *array, ifstream &ifs) {
	string type;
	ifs >> type;
	
	if (type == "==") {
		return new EqResult(array, ifs);
	} else if (type == "<=") {
		return new LtEqResult(array, ifs);
	} else if (type == ">") {
		return new GtResult(array, ifs);
	} else if (type == "IF") {
		return new IfResult(array, ifs);
	}
	
}

BoolResult::~BoolResult() {
	// needs to be defined because it must be virtual
}

FloatResult::FloatResult(LocatingArray *array) {
	this->array = array;
}

FloatResult *FloatResult::readFloatResult(LocatingArray *array, ifstream &ifs) {
	string type;
	ifs >> type;
	
	if (type == "+") {
		return new AdditionResult(array, ifs);
	} else if (type == "*") {
		return new MultiplicationResult(array, ifs);
	} else if (type == "/") {
		return new DivisionResult(array, ifs);
	} else if (type == "C") {
		return new ConstantResult(array, ifs);
	} else if (type == "F") {
		return new FactorAssignment(array, ifs);
	}
	
}

FloatResult::~FloatResult() {
	// needs to be defined because it must be virtual
}

EqResult::EqResult(LocatingArray *array, ifstream &ifs): BoolResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

bool EqResult::getResult(int test) {
	return floatResult1->getResult(test) == floatResult2->getResult(test);
}

void EqResult::writeToStream(ofstream &ofs) {
	ofs << "\t==";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

EqResult::~EqResult() {
	delete floatResult1;
	delete floatResult2;
}

LtEqResult::LtEqResult(LocatingArray *array, ifstream &ifs): BoolResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

bool LtEqResult::getResult(int test) {
	return floatResult1->getResult(test) <= floatResult2->getResult(test);
}

void LtEqResult::writeToStream(ofstream &ofs) {
	ofs << "\t<=";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

LtEqResult::~LtEqResult() {
	delete floatResult1;
	delete floatResult2;
}

GtResult::GtResult(LocatingArray *array, ifstream &ifs): BoolResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

bool GtResult::getResult(int test) {
	return floatResult1->getResult(test) > floatResult2->getResult(test);
}

void GtResult::writeToStream(ofstream &ofs) {
	ofs << "\t>";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

GtResult::~GtResult() {
	delete floatResult1;
	delete floatResult2;
}

IfResult::IfResult(LocatingArray *array, ifstream &ifs): BoolResult(array) {
	boolResult1 = BoolResult::readBoolResult(array, ifs);
	boolResult2 = BoolResult::readBoolResult(array, ifs);
}

bool IfResult::getResult(int test) {
	return !boolResult1->getResult(test) || boolResult2->getResult(test);
}

void IfResult::writeToStream(ofstream &ofs) {
	ofs << "\tIF";
	boolResult1->writeToStream(ofs);
	boolResult2->writeToStream(ofs);
}

IfResult::~IfResult() {
	delete boolResult1;
	delete boolResult2;
}

AdditionResult::AdditionResult(LocatingArray *array, ifstream &ifs): FloatResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

float AdditionResult::getResult(int test) {
	return floatResult1->getResult(test) + floatResult2->getResult(test);
}

void AdditionResult::writeToStream(ofstream &ofs) {
	ofs << "\t+";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

AdditionResult::~AdditionResult() {
	delete floatResult1;
	delete floatResult2;
}

MultiplicationResult::MultiplicationResult(LocatingArray *array, ifstream &ifs): FloatResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

float MultiplicationResult::getResult(int test) {
	return floatResult1->getResult(test) * floatResult2->getResult(test);
}

void MultiplicationResult::writeToStream(ofstream &ofs) {
	ofs << "\t*";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

MultiplicationResult::~MultiplicationResult() {
	delete floatResult1;
	delete floatResult2;
}

DivisionResult::DivisionResult(LocatingArray *array, ifstream &ifs): FloatResult(array) {
	floatResult1 = FloatResult::readFloatResult(array, ifs);
	floatResult2 = FloatResult::readFloatResult(array, ifs);
}

float DivisionResult::getResult(int test) {
	return floatResult1->getResult(test) / floatResult2->getResult(test);
}

void DivisionResult::writeToStream(ofstream &ofs) {
	ofs << "\t/";
	floatResult1->writeToStream(ofs);
	floatResult2->writeToStream(ofs);
}

DivisionResult::~DivisionResult() {
	delete floatResult1;
	delete floatResult2;
}

ConstantResult::ConstantResult(LocatingArray *array, ifstream &ifs): FloatResult(array) {
	ifs >> value;
}

float ConstantResult::getResult(int test) {
	return value;
}

void ConstantResult::writeToStream(ofstream &ofs) {
	ofs << "\tC\t" << value;
}

FactorAssignment::FactorAssignment(LocatingArray *array, ifstream &ifs): FloatResult(array) {
	ifs >> factor_i;
}

float FactorAssignment::getResult(int test) {
	return array->getFactorData()->getNumericFactorLevel(factor_i, array->getLevelMatrix()[test][factor_i]);
}

void FactorAssignment::writeToStream(ofstream &ofs) {
	ofs << "\tF\t" << factor_i;
}
