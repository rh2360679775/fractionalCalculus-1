#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "cxtimers.h"
#include <iostream>
#include "c1024.h"

#include <boost/math/special_functions/gamma.hpp>

using namespace std;

inline FloatType ghGamma(FloatType z)
{

//return  boost::math::tgamma(1+z);
 	FloatType x = z; 

	FloatType val0 =  FloatType(0);
	FloatType val1 =  FloatType(1);
	FloatType factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 

	FloatType sum = as[0];   
	FloatType xn = FloatType(1);
	
	for(int n = 1; n < 1024; n++){
			xn *= x;
			sum += as[n]*xn;
		}
	return factor*sum;
}
#include "trapez.hpp"
