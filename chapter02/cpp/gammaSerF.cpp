#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include <numbers>
#include "c9.h"
#include "templates.hpp"

using namespace std;
using namespace numbers;

inline FloatType ghGammaLib(FloatType x)
{
  // Library function
		return tgammaf(FloatType(1) + x);
}

//inline  T ghGammaDyn( T z ){
template <typename T>  
inline  T ghGamma( T z ){
	// normalize
	T x = z; 
	T factor = T(1);
	shiftAndFactor(&x, &factor);
	// run 
	// T result =  factor*polySum(as, x);
	//T result =  factor*polySumUp(as, x);
	T result =  factor*horner(as, x);
	//T result =  factor*horner45(as, x);
	return result;
}


//inline FloatType ghGammaApprox(FloatType z)
inline FloatType ghGammaApprox(FloatType z)
{
  // Library function
	FloatType x = z + FloatType(1);
	return sqrt(2*pi*x)*pow(x/e,x)*pow(x*sinh(1/x),0.5*x)*exp(7./324./(x*x*x*(35.*x*x+33.)))/x;
}



inline FloatType ghGammaH(FloatType z)
{
	FloatType x = z; 
	FloatType factor = 1.0;
	shiftAndFactor(&x, &factor);
	
	return  factor * (as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + x*(as[6] + x*(as[7] + as[8]*x)))))));

}

//inline FloatType ghGammaH33(FloatType z)
inline FloatType ghGamma(FloatType z)
{
	FloatType x = z; 
	FloatType factor = 1.0;
	shiftAndFactor(&x, &factor);
    
 FloatType sum; 
 FloatType x2 = x*x ;
 FloatType x3 = x2*x ;

	return ((as[0] + as[1]*x + as[2]*x2) + 
	        (as[3] + as[4]*x + as[5]*x2)*x3 +
		      (as[6] + as[7]*x + as[8]*x2)*x3*x3)*factor ;
}


#include "trapez.hpp"