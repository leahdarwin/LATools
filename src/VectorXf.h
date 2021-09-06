#ifndef VECTORXF_H
#define VECTORXF_H

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

class VectorXf {
	
	int length;
	float *data;
	
	// indicates whether SStot has been calculated
	bool usesSStot;
	
	// SStot for calculating r-squared
	float SStot;
	
public:
	
	VectorXf(int length);
	
	float *getData();
	int getLength();
	
	void calculateSStot();
	void writeToFile(std::string responseDir, std::string responseCol);
	void loadResponse(string directory, string column, bool performLog);
	
	float getSStot();
	
	~VectorXf();
};

#endif
