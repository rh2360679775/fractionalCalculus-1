#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include "c20.h"
#include "templates.hpp"

using namespace std;

inline double ghGammaFor(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
  double result = as[19];
	for(int i = 18; i >= 0; i--){
		result = result*x+as[i];
	}
  return result*factor;
}

inline double ghGammaFixed(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);

	double result =   as[19];
	result = result*x+as[18];
	result = result*x+as[17];
	result = result*x+as[16];
	result = result*x+as[15];
	result = result*x+as[14];
	result = result*x+as[13];
	result = result*x+as[12];
	result = result*x+as[11];
	result = result*x+as[10];
	result = result*x+as[ 9];
	result = result*x+as[ 8];
	result = result*x+as[ 7];
	result = result*x+as[ 6];
	result = result*x+as[ 5];
	result = result*x+as[ 4];
	result = result*x+as[ 3];
	result = result*x+as[ 2];
	result = result*x+as[ 1];
  result = result*x+as[ 0];
  return result*factor;
}

inline double ghGamma0(double z)
{
  
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
	double sum = as[0];   
	double xn = 1.0;

	for(int n = 1; n < 20; n++){
		xn *= x;
		sum += as[n]*xn;
	}
	return sum*factor;
}

inline double ghGammaH(double z)
{
 	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
	return (as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x* (as[5] + x*(as[6] +  x*(as[7] + x*
	  (as[8] +  x* (as[9] + 	x* (as[10] +  x*(as[11] + x*   (as[12] +  x* (as[13] +  x*(as[14] + x*   (as[15] + 
		x*  (as[16] +  x* (as[17] +  x*(as[18] + as[19]*x))))))))))))) )))))*factor;
}

inline double ghGamma54(double z)
{
 	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	 
	double a[5];
	a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + as[ 3]*x);
	a[1] = as[ 4] + as[ 5]*x + x*x*(as[ 6] + as[ 7]*x);
	a[2] = as[ 8] + as[ 9]*x + x*x*(as[10] + as[11]*x);
	a[3] = as[12] + as[13]*x + x*x*(as[14] + as[15]*x);
	a[4] = as[16] + as[17]*x + x*x*(as[18] + as[19]*x);
	
	return  (a[0] + a[1]*x*x*x*x + x*x*x*x*x*x*x*x*(a[2] + x*x*x*x*(a[3] + a[4]*x*x*x*x)))*factor;
}

inline double ghGamma(double z)
{
	double x = z;
	static const double one{1};  
	double factor = one;
	shiftAndFactor(&x, &factor);
  
	double a[4];
	a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + x*(as[ 3] + as[ 4]*x));
	a[1] = as[ 5] + as[ 6]*x + x*x*(as[ 7] + x*(as[ 8] + as[ 9]*x));
	a[2] = as[10] + as[11]*x + x*x*(as[12] + x*(as[13] + as[14]*x));
	a[3] = as[15] + as[16]*x + x*x*(as[17] + x*(as[18] + as[19]*x));
	
	return (a[0] + a[1]*x*x*x*x*x + x*x*x*x*x*x*x*x*x*x*(a[2] + a[3]*x*x*x*x*x))*factor;
}

inline double ghGamma37(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
	double a[3];
	a[0] = as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + as[6]*x))));
	a[1] = as[7] + as[8]*x + x*x*(as[9] + x*(as[10] + x*(as[11] + x*(as[12] + as[13]*x))));
	a[2] = as[14] + as[15]*x + x*x*(as[16] + x*(as[17] + x*(as[18] + x*as[19] )));
	
	return (a[0] + a[1]*x*x*x*x*x*x*x + x*x*x*x*x*x*x*x*x*x*x*x*x*x*a[2])*factor ;
}
// 45
inline double ghGamma45(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
   
	double a[4];
	a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + x*(as[ 3] + as[ 4]*x));
	a[1] = as[ 5] + as[ 6]*x + x*x*(as[ 7] + x*(as[ 8] + as[ 9]*x));
	a[2] = as[10] + as[11]*x + x*x*(as[12] + x*(as[13] + as[14]*x));
	a[3] = as[15] + as[16]*x + x*x*(as[17] + x*(as[18] + as[19]*x));

	return (a[0] + a[1]*x*x*x*x*x + x*x*x*x*x*x*x*x*x*x*(a[2] + a[3]*x*x*x*x*x))*factor;
}

#include "trapez.hpp"