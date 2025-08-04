



#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <omp.h>
#include "c44.h"
#include <cmath>
#include "templates.hpp"

using namespace std;

	
inline long double ghGamma0(long double x)
{
  // infinity-norm
			
		long double sum = as[0];   
		long double xn = 1.0L;

		for(int n = 1; n < 44; n++){
			xn *= x;
			sum += as[n]*xn;
		}
		return sum;
}

inline long double ghGammaH1(long double x)
{
 				
		return as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + 
         x*(as[4] + x*  (as[5] + 	x*(as[6] + x*(as[7] + x*(as[8] +x*(as[9] + 	x*	(as[10] +x*	(as[11] +x*(as[12] +x*(as[13] +x*(as[14] +	
		 x*	(as[15] +x*	(as[16] +x*		(as[17] +	x*	(as[18] +	x*		(as[19] +		x*
		(as[20] +x*	(as[21] +	x*	(as[22] +	x*	(as[23] +x*	(as[24] +x*		(as[25] +x*	(as[26] +x*	(as[27] +x*
		(as[28] +x*	(as[29] +	x*	(as[30] +	x*	(as[31] +x*	(as[32] +x*		(as[33] +x*	(as[34] +x*	(as[35] +x*
		(as[36] +x*	(as[37] +	x*	(as[38] +	x*	(as[39] +x*	
		(as[40] +x*	(as[41] +x*(as[42] + as[43]*x)))))))))))))))))))))))))))))))))))))))));
}

inline long double ghGamma(long double x)
{
  // infinity-norm
			
   long double a[6];
   long double x2 = x*x;
   long double x8 = x2*x2*x2*x2;

   a[0] = as[0]  + as[1]*x  + x*x*(as[2]  + x*(as[3]  + x*(as[4]  + x* (as[5]  + x*(as[6]  + as[7]*x)))));
   a[1] = as[8]  + as[9]*x  + x*x*(as[10] + x*(as[11] + x*(as[12] + x* (as[13] + x*(as[14] + as[15]*x)))));
   a[2] = as[16] + as[17]*x + x*x*(as[18] + x*(as[19] + x*(as[20] + x* (as[21] + x*(as[22] + as[23]*x)))));
   a[3] = as[24] + as[25]*x + x*x*(as[26] + x*(as[27] + x*(as[28] + x* (as[29] + x*(as[30] + as[31]*x)))));
   a[4] = as[32] + as[33]*x + x*x*(as[34] + x*(as[35] + x*(as[36] + x* (as[37] + x*(as[38] + as[39]*x)))));
   a[5] = as[40] + as[41]*x + x*x*(as[42] + x*as[43] );

	return a[0] + a[1]*x8 + x8*x8*(a[2] + x8*(a[3] + x8*(a[4] + a[5]*x8)));

}


inline long double ghGamma77(long double x)
{
  // infinity-norm
		
   long double x2 = x*x;
   long double x4 = x2*x2;
   long double x7 = x4*x2*x;
			
   long double a[7];
			
   a[0] = as[0]  + as[1]*x  + x2*(as[2]  + x*(as[3]  + x*(as[4]  + x* (as[5]  + x*as[6]))));
   a[1] = as[7]  + as[8]*x  + x2*(as[9]  + x*(as[10] + x*(as[11] + x* (as[12] + x*as[13]))));
   a[2] = as[14] + as[15]*x + x2*(as[16] + x*(as[17] + x*(as[18] + x* (as[19] + x*as[20]))));
   a[3] = as[21] + as[22]*x + x2*(as[23] + x*(as[24] + x*(as[25] + x* (as[26] + x*as[27]))));
   a[4] = as[28] + as[29]*x + x2*(as[30] + x*(as[31] + x*(as[32] + x* (as[33] + x*as[34]))));
   a[5] = as[35] + as[36]*x + x2*(as[37] + x*(as[38] + x*(as[39] + x* (as[40] + x*as[41]))));
   a[6] = as[42] + as[43]*x ;
   
   return a[0]  + a[1]*x7  + x7*x7*(a[2]  + x7*(a[3]  + x7*(a[4]  + x7* (a[5]  + x7*a[6]))));

}

#include "trapez.hpp"