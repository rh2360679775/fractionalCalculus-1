#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "cxtimers.h"
#include <iostream>
#include "c1331.h"
#include "templates.hpp"

#include <boost/math/special_functions/gamma.hpp>

using namespace std;

inline FloatType ghGamma(FloatType z)
{

//return  boost::math::tgamma(1+z);
 	FloatType x = z; 
	FloatType factor = FloatType(1);

	shiftAndFactor(&x, &factor);

	FloatType sum = as[0];   
	FloatType xn = FloatType(1);
	
	for(int n = 1; n < 1331; n++){
			xn *= x;
			sum += as[n]*xn;
		}
	return factor*sum;
}
#include "trapez.hpp"
