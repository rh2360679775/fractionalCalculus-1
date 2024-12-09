#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include <numbers>
#include "c27.h"
#include "templates.hpp"

using namespace std;
using namespace numbers;


inline FloatType ghGammaLib(FloatType x)
{
  // Library function
		return tgammal(FloatType(1) + x);
}


//inline  T ghGammaDyn( T z ){
template <typename T>  
inline  T ghGammaDyn( T z ){
	// normalize
	T x = z; 
	T factor = T(1);
	shiftAndFactor(&x, &factor);
	// run 
	// T result =  factor*polySum(as, x);
	T result =  factor*horner(as, x);
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
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);
	
	return factor*(as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + 
				x*(as[4] + x*  (as[5] + 	x*(as[6] + x*(as[7] + x*(as[8] +x*(as[9] + 	x*	(as[10] +x*	(as[11] +x*(as[12] +x*(as[13] +x*(as[14] +	
		x*	(as[15] +x*	(as[16] +x*		(as[17] +	x*	(as[18] +	x*		(as[19] +		x*
	(as[20] +x*	(as[21] +	x*	(as[22] +	x*	(as[23] +x*	(as[24] +x*		(as[25] +x*	as[26])))))))))))))))))))))))));

}


//inline FloatType ghGammaH74(FloatType z)
inline FloatType ghGamma74(FloatType z)
{
  // infinity-norm 74

	FloatType x = z; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);

   FloatType x7 = x*x*x*x*x*x*x;
			
   FloatType a[4];
			
   a[0] = as[0]  + as[1]*x  + x*x*(as[2]  + x*(as[3]  + x*(as[4]  + x* (as[5]  + x*as[6]))));
   a[1] = as[7]  + as[8]*x  + x*x*(as[9]  + x*(as[10] + x*(as[11] + x* (as[12] + x*as[13]))));
   a[2] = as[14] + as[15]*x + x*x*(as[16] + x*(as[17] + x*(as[18] + x* (as[19] + x*as[20]))));
   a[3] = as[21] + as[22]*x + x*x*(as[23] + x*(as[24] + x*(as[25] + x* as[26] )));
   
   return factor*(a[0]  + a[1]*x7  + x7*x7*(a[2]  + x7*a[3]));

}

//inline FloatType ghGamma333(FloatType q)
inline FloatType ghGamma(FloatType q)
{
  FloatType x = q; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);
			
 std::array<FloatType, 3*3> a;
   a[0] = as[0]  + as[1]*x  + x*x*as[2];
   a[1] = as[3]  + as[4]*x  + x*x*as[5];
   a[2] = as[6]  + as[7]*x  + x*x*as[8];
   a[3] = as[9]  + as[10]*x  + x*x*as[11];
   a[4] = as[12]  + as[13]*x  + x*x*as[14];
   a[5] = as[15]  + as[16]*x  + x*x*as[17];
   a[6] = as[18]  + as[19]*x  + x*x*as[20];
   a[7] = as[21]  + as[22]*x  + x*x*as[23];
   a[8] = as[24]  + as[25]*x  + x*x*as[26];
  
	FloatType  y = x*x*x;
	std::array<FloatType, 3> b;
	b[ 0] = a[ 0] + a[ 1]*y + y*y* a[ 2] ;
	b[ 1] = a[ 3] + a[ 4]*y + y*y* a[ 5] ;
  b[ 2] = a[ 6]  + a[7]*y + y*y* a[ 8];
 	
 FloatType  z = y*y*y;
  return  factor*(b[0] + b[1]*z + z*z* b[2]);

}

#include "trapez.hpp"