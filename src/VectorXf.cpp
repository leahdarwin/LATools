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
#include <Rcpp.h>
#include <sys/types.h>
#include "VectorXf.h"


using namespace Rcpp;
using namespace std;


VectorXf::VectorXf(int length) {
	this->length = length;
	this->usesSStot = false;

	data = new float[length];

}

float *VectorXf::getData() {
	return data;
}

int VectorXf::getLength() {
	return length;
}

void VectorXf::calculateSStot() {
	// calculate SStot for r-squared

	// 1st get yBar
	float yBar, ySum = 0;
	for (int row_i = 0; row_i < length; row_i++) {
		ySum += data[row_i];
	}
	yBar = ySum / length;

	// 2nd get SStot
	SStot = 0;
	for (int row_i = 0; row_i < length; row_i++) {
		SStot += (data[row_i] - yBar) * (data[row_i] - yBar);
	}

	usesSStot = true;
}

float VectorXf::getSStot() {
	// verify SStot has been calculated
	if (!usesSStot) {
		cout << "You have not calculated SStot yet" << endl;
		return 0;
	}

	return SStot;
}


void VectorXf::writeToFile(string responseDir, string responseCol){

	// open a response file for writing
	string file =  "Response_written.tsv";
	ofstream ofs(file.c_str());

	// write response headers
	ofs << length << endl;
	ofs << responseCol << endl;

	float *responseData = this->getData();
	cout << "CHECK Response: " << responseData[0] << endl;

	//write response values
	for (int row_i = 0; row_i < length; row_i++){
		ofs << responseData[row_i] << endl;
	}

}

void VectorXf::loadResponse(string directory, string column, bool performLog){
  struct dirent *dp;
  
  // open the directory
  DIR *dir = opendir(directory.c_str());
  if (dir == NULL) {
    cout << "Cannot open directory \"" << directory << "\"" << endl;
      exit(0);
  }
  
  int cols, col_i;
  string header;
  
  int rows = this->getLength();
  
  string tempString;
  float tempFloat;
  int tempInt;
  
  // the response count vector
  VectorXf *responseCount = new VectorXf(rows);
  
  // grab data vectors
  float *data = this->getData();
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
        data = this->getData();
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
  data = this->getData();
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
  this->calculateSStot();
  
  cout << "Loaded responses" << endl;
  
}

VectorXf::~VectorXf() {
	delete[] data;
}

// [[Rcpp::export]]
SEXP makeVectorXF2(int length, std::string response_dir, std::string response_column, bool log_data){
  VectorXf* vec = new VectorXf(length);
  vec->loadResponse(response_dir, response_column, log_data);
  Rcpp::XPtr<VectorXf> vec_ptr(vec);
  return vec_ptr;
}