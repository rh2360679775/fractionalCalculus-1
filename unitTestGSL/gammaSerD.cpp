#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include <numbers>
#include "c20.h"
#include "templates.hpp"

using namespace std;

//inline double ghGammaH120(double z)
inline void printList()
{
 	printList(as, prec);
}


//inline double ghGammaLib(double x)
inline double ghGammaLib(double x)
{
  // Library function
		return tgamma(double(1) + x);
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
	//T result =  factor*polySumUp(as, x);
	//T result =  factor*horner(as, x);
	T result =  factor*horner45(as, x);
	return result;
}

//inline double ghGammaFixed(double z)
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


//inline double ghGammaH120(double z)
inline double ghGammaH120(double z)
{
 	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
	return (as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x* (as[5] + x*(as[6] +  x*(as[7] + x*
	  (as[8] +  x* (as[9] + 	x* (as[10] +  x*(as[11] + x*   (as[12] +  x* (as[13] +  x*(as[14] + x*   (as[15] + 
		x*  (as[16] +  x* (as[17] +  x*(as[18] + as[19]*x))))))))))))) )))))*factor;
}

//inline double ghGamma102(double z)
inline double ghGamma102(double z)
{
 	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	 
	double a[10];
	a[0] = as[ 0] + as[ 1]*x ;
	a[1] = as[ 2] + as[ 3]*x ;
	a[2] = as[ 4] + as[ 5]*x ;
	a[3] = as[ 6] + as[ 7]*x ;
	a[4] = as[ 8] + as[ 9]*x ;
	a[5] = as[10] + as[11]*x ;
	a[6] = as[12] + as[13]*x ;
  a[7] = as[14] + as[15]*x ;
  a[8] = as[16] + as[17]*x ;
  a[9] = as[18] + as[19]*x ;

	double y = x*x;
	
  return (a[0] + a[1]*y + y*y*(a[2] + y*(a[3] + y*(a[4] + y* (a[5] + y*(a[6] +  y*(a[7] + y*
	  (a[8] +  y* (a[9])))))))))*factor;
}



//inline double ghGamma73(double z)
inline double ghGamma73(double z)
{
 	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	 
	double a[7];
	a[0] = as[ 0] + as[ 1]*x + x*x*as[ 2] ;
	a[1] = as[ 3] + as[ 4]*x + x*x*as[ 5] ;
	a[2] = as[ 6] + as[ 7]*x + x*x*as[ 8] ;
	a[3] = as[ 9] + as[10]*x + x*x*as[11] ;
	a[4] = as[12] + as[13]*x + x*x*as[14] ;
	a[5] = as[15] + as[16]*x + x*x*as[17] ;
	a[6] = as[18] + as[19]*x ;

	double y = x*x*x;
	
  return (a[0] + a[1]*y + y*y*(a[2] + y*(a[3] + y*(a[4] + y*(a[5] + y*a[6])))))*factor;
}

//inline double ghGamma54(double z)
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
	double y = x*x*x*x;
	
	return  (a[0] + a[1]*y + y*y*(a[2] + y*(a[3] + a[4]*y)))*factor;
}


// 45
//inline double ghGamma45(double z)
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



//inline double ghGamma37(double z)
inline double ghGamma37(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	double x2 = x*x;
	
	double a[3];
	a[0] = as[ 0] + as[ 1]*x + x2*(as[ 2] + x*(as[ 3] + x*(as[ 4] + x*(as[ 5] + as[ 6]*x))));
	a[1] = as[ 7] + as[ 8]*x + x2*(as[ 9] + x*(as[10] + x*(as[11] + x*(as[12] + as[13]*x))));
	a[2] = as[14] + as[15]*x + x2*(as[16] + x*(as[17] + x*(as[18] + x*as[19] )));
	
	return (a[0] + a[1]*x*x2*x2*x2 + x2*x2*x2*x2*x2*x2*x2*a[2])*factor ;
}


//inline double ghGamma210(double z)
inline double ghGamma210(double z)
{
	double x = z; 
	double factor = double(1);
	shiftAndFactor(&x, &factor);
	
	double a[2];
	a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + x*(as[ 3] + x*(as[ 4] + x*(as[ 5] + x*(as[ 6] + x*(as[ 7] + x*(as[ 8] + as[ 9]*x)))))));
	a[1] = as[10] + as[11]*x + x*x*(as[12] + x*(as[13] + x*(as[14] + x*(as[15] + x*(as[16] + x*(as[17] + x*(as[18] + as[19]*x)))))));
	
//	return (a[0] + a[1]*pow(x,10))*factor ;
	return (a[0] + a[1]*x*x*x*x*x*x*x*x*x*x)*factor ;
}


//inline double ghGammaH3(double q)
inline double ghGamma(double q)
{
  double x = q; 
	double factor = double(1);
 	shiftAndFactor(&x, &factor);
				
	std::array<double, 3*3> a;
	a[ 0] = as[ 0] + as[ 1]*x + x*x* as[ 2] ;
	a[ 1] = as[ 3] + as[ 4]*x + x*x* as[ 5] ;
	a[ 2] = as[ 6] + as[ 7]*x + x*x* as[ 8] ;
	a[ 3] = as[ 9] + as[10]*x + x*x* as[11] ;
	a[ 4] = as[12] + as[13]*x + x*x* as[14] ;
	a[ 5] = as[15] + as[16]*x + x*x* as[17] ;
	a[ 6] = as[18] + as[19]*x ;
	
  double  y = x*x*x;
	std::array<double, 3> b;
	b[ 0] = a[ 0] + a[ 1]*y + y*y* a[ 2] ;
	b[ 1] = a[ 3] + a[ 4]*y + y*y* a[ 5] ;
	b[ 2] = a[ 6]  ;
	
	double  z = y*y*y;
  return  factor*(b[0] + b[1]*z + z*z* b[2]);
}

//inline double ghGammaApprox(double z)
inline double ghGammaApprox(double z)
{
  // Yang, ZH., Tian, JF. An accurate approximation formula for gamma function. 
	// J Inequal Appl 2018, 56 (2018). https://doi.org/10.1186/s13660-018-1646-6
	double x = z + double(1);
	double xm = double(1)/x;
	double x2 = x*x;
	double cn = double(7)/double(324);
	double c35 = double(35);
	double c33 = double(33);
	return sqrt(double(2)*numbers::pi*x)*
		pow(x/numbers::e,x)*
		pow(x*sinh(xm),double(0.5)*x)*
		exp(cn / (x*x2*(c35*x2 + c33)))
		*xm;
}
#include "trapez.hpp"