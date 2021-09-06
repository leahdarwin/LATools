#include "Noise.h"


Noise::Noise(float ratio) {
	this->ratio = ratio;
}

float Noise::addNoise(float mean, float range) {
	float noise = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	float noiseRange = ratio * range;
	float noisy = mean + (noise * noiseRange) - (noiseRange / 2);

	cout << "Range: " << range << ", NoiseRange: " << noiseRange << ", value:noisy: " << mean << "\t" << noisy << endl;

	return noisy;
}
