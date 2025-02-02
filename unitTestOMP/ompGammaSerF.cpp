#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include "c9.h"
#include "templates.hpp"

using namespace std;


inline float ghGamma(float z)
{
	float x = z; 
	float factor = 1.0;
	//shiftAndFactor(&x, &factor);

	float sum = as[0];   
	float xn = 1;

	for(int n = 1; n < 9; n++){
		xn *= x;
		sum += as[n]*xn;
	}
	return factor*sum;
}

inline float ghGammaH(float z)
{
	float x = z; 
	float factor = 1.0;
	shiftAndFactor(&x, &factor);
	
	return  factor * (as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + x*(as[6] + x*(as[7] + as[8]*x)))))));

}

inline float ghGammaH233(float z)
{
	float x = z; 
	float factor = 1.0;
	shiftAndFactor(&x, &factor);
    
 float x2 = x*x ;
 float x3 = x2*x ;

	return ((as[0] + as[1]*x + as[2]*x2) + 
	 (as[3] + as[4]*x + as[5]*x2)*x3 +
		 (as[6] + as[7]*x + as[8]*x2)*x3*x3)*factor ;
}


#include "trapez.hpp"