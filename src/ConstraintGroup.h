#ifndef CONSTRAINTGROUP_H
#define CONSTRAINTGROUP_H

#include <cstdlib>
#include <iostream>
#include <string>

#include "LocatingArray.h"
#include "FactorData.h"

using namespace std;

class LocatingArray;

class BoolResult {
protected:
	LocatingArray *array;
public:
	BoolResult(LocatingArray *array);
	
	static BoolResult *readBoolResult(LocatingArray *array, ifstream &ifs);
	
	virtual bool getResult(int test) = 0;
	virtual void writeToStream(ofstream &ofs) = 0;
	
	virtual ~BoolResult();
};

class FloatResult {
protected:
	LocatingArray *array;
public:
	FloatResult(LocatingArray *array);
	
	static FloatResult *readFloatResult(LocatingArray *array, ifstream &ifs);
	
	virtual float getResult(int test) = 0;
	virtual void writeToStream(ofstream &ofs) = 0;
	
	virtual ~FloatResult();
};

class ConstraintGroup {
private:
	LocatingArray *groupLA;
	
	int constraints;
	BoolResult** boolConstraints;
	
	// there exists a weight window for each row (larger weight window means more likely to be chosen)
	int *weightMin;		// bottom of weight window for each row
	int *weightMax;		// top of weight window for each row
	int weightRandMax;
public:
	int factors;
	int *factorIndeces;
	
	ConstraintGroup(LocatingArray *array, ifstream &ifs);
	void populateGroupLA(LocatingArray *array, char *levelRow, char *factorLevels, int fixedFactors, int **settingCount, int &fullFactorialCount);
	bool satisfiableInGroupLA(char *requireLevelRow, char *avoidLevelRow);
	virtual bool getResult(int test);
	void randPopulateLevelRow(char *levelRow);
	void writeToStream(ofstream &ofs);
	
	virtual ~ConstraintGroup();
};

class EqResult: public BoolResult {
private:
	FloatResult *floatResult1; // LHS of Eq equation
	FloatResult *floatResult2; // RHS of Eq equation
public:
	EqResult(LocatingArray *array, ifstream &ifs);
	virtual bool getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~EqResult();
};

class LtEqResult: public BoolResult {
private:
	FloatResult *floatResult1; // LHS of LtEq equation
	FloatResult *floatResult2; // RHS of LtEq equation
public:
	LtEqResult(LocatingArray *array, ifstream &ifs);
	virtual bool getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~LtEqResult();
};

class GtResult: public BoolResult {
private:
	FloatResult *floatResult1; // LHS of Gt equation
	FloatResult *floatResult2; // RHS of Gt equation
public:
	GtResult(LocatingArray *array, ifstream &ifs);
	virtual bool getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~GtResult();
};

class IfResult: public BoolResult {
private:
	BoolResult *boolResult1; // If this is true
	BoolResult *boolResult2; // Then this must also be true
public:
	IfResult(LocatingArray *array, ifstream &ifs);
	virtual bool getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~IfResult();
};

class AdditionResult: public FloatResult {
private:
	FloatResult *floatResult1;
	FloatResult *floatResult2;
public:
	AdditionResult(LocatingArray *array, ifstream &ifs);
	virtual float getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~AdditionResult();
};

class MultiplicationResult: public FloatResult {
private:
	FloatResult *floatResult1;
	FloatResult *floatResult2;
public:
	MultiplicationResult(LocatingArray *array, ifstream &ifs);
	virtual float getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~MultiplicationResult();
};

class DivisionResult: public FloatResult {
private:
	FloatResult *floatResult1;
	FloatResult *floatResult2;
public:
	DivisionResult(LocatingArray *array, ifstream &ifs);
	virtual float getResult(int test);
	virtual void writeToStream(ofstream &ofs);
	virtual ~DivisionResult();
};

class ConstantResult: public FloatResult {
private:
	float value;
public:
	ConstantResult(LocatingArray *array, ifstream &ifs);
	float getResult(int test);
	virtual void writeToStream(ofstream &ofs);
};

class FactorAssignment: public FloatResult {
private:
	int factor_i;
public:
	FactorAssignment(LocatingArray *array, ifstream &ifs);
	float getResult(int test);
	virtual void writeToStream(ofstream &ofs);
};

#endif