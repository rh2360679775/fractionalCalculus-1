#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include <numbers>
#include "c125.h"
//#include "c1024.h"
//#include "c1331.h"
#include "templates.hpp"
#include <boost/math/special_functions/gamma.hpp>


using namespace std;
using namespace numbers;


//inline FloatType ghGammaLib(FloatType x)
inline FloatType ghGamma(FloatType x)
{
  // Library function
   return  boost::math::tgamma(FloatType(1)+x);
}


//inline  T ghGammaDyn( T z ){
template <typename T>  
inline  T ghGammaDyn( T z ){
	// normalize
	T x = z; 
	T factor = T(1);
	shiftAndFactor(&x, &factor);
	// run 
	 T result =  factor*polySum(as, x);
	//T result =  factor*horner(as, x);
	return result;
}


//inline FloatType ghGammaApprox(FloatType z)
inline FloatType ghGammaApprox(FloatType x)
{
  // Library function
	return sqrt(2*pi*x)*pow(x/e,x)*pow(x*sinh(1/x),0.5*x)*exp(7./324./(x*x*x*(35.*x*x+33.)))/x;
}


//inline FloatType ghGammaFixed(FloatType z)
inline FloatType ghGammaFixed(FloatType z)
{
	FloatType x = z; 
	FloatType factor = FloatType(1);
	shiftAndFactor(&x, &factor);

	FloatType result =    as[124];
	result = result * x + as[123];
	result = result * x + as[122];
	result = result * x + as[121];
	result = result * x + as[120];
	result = result * x + as[119];
	result = result * x + as[118];
	result = result * x + as[117];
	result = result * x + as[116];
	result = result * x + as[115];
	result = result * x + as[114];
	result = result * x + as[113];
	result = result * x + as[112];
	result = result * x + as[111];
	result = result * x + as[110];
	result = result * x + as[109];
	result = result * x + as[108];
	result = result * x + as[107];
	result = result * x + as[106];
	result = result * x + as[105];
	result = result * x + as[104];
	result = result * x + as[103];
	result = result * x + as[102];
	result = result * x + as[101];
	result = result * x + as[100];
	result = result * x + as[99];
	result = result * x + as[98];
	result = result * x + as[97];
	result = result * x + as[96];
	result = result * x + as[95];
	result = result * x + as[94];
	result = result * x + as[93];
	result = result * x + as[92];
	result = result * x + as[91];
	result = result * x + as[90];
	result = result * x + as[89];
	result = result * x + as[88];
	result = result * x + as[87];
	result = result * x + as[86];
	result = result * x + as[85];
	result = result * x + as[84];
	result = result * x + as[83];
	result = result * x + as[82];
	result = result * x + as[81];
	result = result * x + as[80];
	result = result * x + as[79];
	result = result * x + as[78];
	result = result * x + as[77];
	result = result * x + as[76];
	result = result * x + as[75];
	result = result * x + as[74];
	result = result * x + as[73];
	result = result * x + as[72];
	result = result * x + as[71];
	result = result * x + as[70];
	result = result * x + as[69];
	result = result * x + as[68];
	result = result * x + as[67];
	result = result * x + as[66];
	result = result * x + as[65];
	result = result * x + as[64];
	result = result * x + as[63];
	result = result * x + as[62];
	result = result * x + as[61];
	result = result * x + as[60];
	result = result * x + as[59];
	result = result * x + as[58];
	result = result * x + as[57];
	result = result * x + as[56];
	result = result * x + as[55];
	result = result * x + as[54];
	result = result * x + as[53];
	result = result * x + as[52];
	result = result * x + as[51];
	result = result * x + as[50];
	result = result * x + as[49];
	result = result * x + as[48];
	result = result * x + as[47];
	result = result * x + as[46];
	result = result * x + as[45];
	result = result * x + as[44];
	result = result * x + as[43];
	result = result * x + as[42];
	result = result * x + as[41];
	result = result * x + as[40];
	result = result * x + as[39];
	result = result * x + as[38];
	result = result * x + as[37];
	result = result * x + as[36];
	result = result * x + as[35];
	result = result * x + as[34];
	result = result * x + as[33];
	result = result * x + as[32];
	result = result * x + as[31];
	result = result * x + as[30];
	result = result * x + as[29];
	result = result * x + as[28];
	result = result * x + as[27];
	result = result * x + as[26];
	result = result * x + as[25];
	result = result * x + as[24];
	result = result * x + as[23];
	result = result * x + as[22];
	result = result * x + as[21];
	result = result * x + as[20];
	result = result * x + as[19];
	result = result * x + as[18];
	result = result * x + as[17];
	result = result * x + as[16];
	result = result * x + as[15];
	result = result * x + as[14];
	result = result * x + as[13];
	result = result * x + as[12];
	result = result * x + as[11];
	result = result * x + as[10];
	result = result * x + as[9];
	result = result * x + as[8];
	result = result * x + as[7];
	result = result * x + as[6];
	result = result * x + as[5];
	result = result * x + as[4];
	result = result * x + as[3];
	result = result * x + as[2];
	result = result * x + as[1];
	result = result * x + as[0];


 return result*factor;
}



//inline FloatType ghGamma12H12(FloatType z)
inline FloatType ghGamma12H12(FloatType z)
{
  // infinity-norm 74

	FloatType x = z; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);

   		
   FloatType a[12];

	 a[0]= as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + x*(as[6] + x*(as[7] + x*(as[8] + x*(as[9] + x*(as[10] + as[11]*x)))))))));
	 a[1]=as[12] + as[13]*x + x*x*(as[14] + x*(as[15] + x*(as[16] + x*(as[17] + x*(as[18] + x*(as[19] + x*(as[20] + x*(as[21] + x*(as[22] + as[23]*x)))))))));
	 a[2]= as[24] + as[25]*x + x*x*(as[26] + x*(as[27] + x*(as[28] +  x*(as[29] + x*(as[30] + x*(as[31] + x*(as[32] + x*(as[33] + x*(as[34] + as[35]*x)))))))));
   a[3]=as[36] + as[37]*x + x*x*(as[38] + x*(as[39] + x*(as[40] + 
   x*(as[41] + x*(as[42] + x*(as[43] + x*(as[44] + x*(as[45] + x*(as[46] + as[47]*x)))))))));
	 a[4]=	as[48] + as[49]*x + x*x*(as[50] + x*(as[51] + x*(as[52] + x*(as[53] + x*(as[54] + x*(as[55] + x*(as[56] + x*(as[57] + x*(as[58] + as[59]*x)))))))));
	 a[5] = as[60] + as[61]*x + x*x*(as[62] + x*(as[63] + x*(as[64] + x*(as[65] + x*(as[66] + x*(as[67] + x*(as[68] + x*(as[69] + x*(as[70] + as[71]*x)))))))));
	a[6] = as[72] + as[73]*x + x*x*(as[74] + x*(as[75] + x*(as[76] + x*(as[77] + x*(as[78] + x*(as[79] + x*(as[80] + x*(as[81] + x*(as[82] + as[83]*x)))))))));
	a[7] = as[84] + as[85]*x + x*x*(as[86] + x*(as[87] + x*(as[88] + x*(as[89] + x*(as[90] + x*(as[91] + x*(as[92] + x*(as[93] + x*(as[94] + as[95]*x)))))))));
	a[8] = as[96] + as[97]*x + x*x*(as[98] + x*(as[99] + x*(as[100] + x*(as[101] + x*(as[102] + x*(as[103] + x*(as[104] + x*(as[105] + x*(as[106] + as[107]*x)))))))));
	a[9] = as[108] + as[109]*x + x*x*(as[110] + x*(as[111] + x*(as[112] + x*(as[113] + x*(as[114] + x*(as[115] + x*(as[116] + x*(as[117] + x*(as[118] + as[119]*x)))))))));
	a[10] = as[120] + as[121]*x + x*x*(as[122] + x*(as[123] + x*(as[124] )));
	FloatType y = pow(x,12);
			
   
   return factor*( a[0] + a[1]*y + y*y*(a[2] + y*(a[3] + y*(a[4] + y*(a[5] + y*(a[6] + y*(a[7] + y*(a[8] + y*(a[9] + y*(as[10]))))))))));

}

//inline FloatType ghGamma555(FloatType q)
inline FloatType ghGamma555(FloatType q)
{
	FloatType x = q; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);
   			
  std::array<FloatType, 5*5> a;
  
	a[0] = (as[0] + as[1]*x + x*x* (as[2] + x*(as[3] + as[4]*x)));
	a[1] = (as[5] + as[6]*x + x*x* (as[7] + x*(as[8] + as[9]*x)));
	a[2] = (as[10] + as[11]*x + x*x* (as[12] + x*(as[13] + as[14]*x)));
	a[3] = (as[15] + as[16]*x + x*x* (as[17] + x*(as[18] + as[19]*x)));
	a[4] = (as[20] + as[21]*x + x*x* (as[22] + x*(as[23] + as[24]*x)));
	a[5] = (as[25] + as[26]*x + x*x* (as[27] + x*(as[28] + as[29]*x)));
	a[6] = (as[30] + as[31]*x + x*x* (as[32] + x*(as[33] + as[34]*x)));
	a[7] = (as[35] + as[36]*x + x*x* (as[37] + x*(as[38] + as[39]*x)));
	a[8] = (as[40] + as[41]*x + x*x* (as[42] + x*(as[43] + as[44]*x)));
	a[9] = (as[45] + as[46]*x + x*x* (as[47] + x*(as[48] + as[49]*x)));
	a[10] = (as[50] + as[51]*x + x*x* (as[52] + x*(as[53] + as[54]*x)));
	a[11] = (as[55] + as[56]*x + x*x* (as[57] + x*(as[58] + as[59]*x)));
	a[12] = (as[60] + as[61]*x + x*x* (as[62] + x*(as[63] + as[64]*x)));
	a[13] = (as[65] + as[66]*x + x*x* (as[67] + x*(as[68] + as[69]*x)));
	a[14] = (as[70] + as[71]*x + x*x* (as[72] + x*(as[73] + as[74]*x)));
	a[15] = (as[75] + as[76]*x + x*x* (as[77] + x*(as[78] + as[79]*x)));
	a[16] = (as[80] + as[81]*x + x*x* (as[82] + x*(as[83] + as[84]*x)));
	a[17] = (as[85] + as[86]*x + x*x* (as[87] + x*(as[88] + as[89]*x)));
	a[18] = (as[90] + as[91]*x + x*x* (as[92] + x*(as[93] + as[94]*x)));
	a[19] = (as[95] + as[96]*x + x*x* (as[97] + x*(as[98] + as[99]*x)));
	a[20] = (as[100] + as[101]*x + x*x* (as[102] + x*(as[103] + as[104]*x)));
	a[21] = (as[105] + as[106]*x + x*x* (as[107] + x*(as[108] + as[109]*x)));
	a[22] = (as[110] + as[111]*x + x*x* (as[112] + x*(as[113] + as[114]*x)));
	a[23] = (as[115] + as[116]*x + x*x* (as[117] + x*(as[118] + as[119]*x)));
	a[24] = (as[120] + as[121]*x + x*x* (as[122] + x*(as[123] + as[124]*x)));

  std::array<FloatType,  5> b;
  FloatType  y = x*x*x*x*x;
  	 
	b[0] = (a[0] + a[1]*y + y*y* (a[2] + y*(a[3] + a[4]*y)));
	b[1] = (a[5] + a[6]*y + y*y* (a[7] + y*(a[8] + a[9]*y)));
	b[2] = (a[10] + a[11]*y + y*y* (a[12] + y*(a[13] + a[14]*y)));
	b[3] = (a[15] + a[16]*y + y*y* (a[17] + y*(a[18] + a[19]*y)));
	b[4] = (a[20] + a[21]*y + y*y* (a[22] + y*(a[23] + a[24]*y)));

  FloatType  z = y*y*y*y*y;
  
  return  factor*(b[0] + b[1]*z + z*z* (b[2] + z*(b[3] + b[4]*z)));
}


#include "trapez.hpp"