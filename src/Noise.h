#ifndef NOISE_H
#define NOISE_H

#include <cstdlib>
#include <iostream>

using namespace std;

class Noise {
private:
	
	float ratio;
	
public:
	
	Noise(float ratio);
	float addNoise(float mean, float range);
	
};

#endif