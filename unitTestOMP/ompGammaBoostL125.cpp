#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "c125.h"
#include "templates.hpp"

#include <boost/math/special_functions/gamma.hpp>

using namespace std;

//inline FloatType ghGammaFor(FloatType z)
inline FloatType ghGamma(FloatType z)
{
  //return  boost::math::tgamma(FloatType(1)+z);

  FloatType x = z; 
	FloatType factor = FloatType(1);
	shiftAndFactor(&x, &factor);
	
  FloatType result = as[cSize-1];
	for(int i = cSize-2; i >= 0; i--){
		result = result*x+as[i];
	}
  return result*factor;
}

//inline FloatType ghGammaLoop(FloatType z)
inline FloatType ghGammaLoop(FloatType z)
{

  FloatType x = z; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);
	FloatType sum = as[0];   
	FloatType xn = FloatType(1);
	
	for(int n = 1; n < 125; n++){
			xn *= x;
			sum += as[n]*xn;
		}
	return factor*sum;
}



//inline FloatType ghGammaH1(FloatType z)
inline FloatType ghGammaH1(FloatType z)
{

 	FloatType x = z; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);

return factor*(as[0] + as[1]*x + x*x*
(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + x*(as[6] +x*(as[7] + x*(as[8] +x*(as[9] + 
x*(as[10] +x*(as[11] + x*(as[12] +x*(as[13] + x*(as[14] +x*(as[15] + x*(as[16] +x*(as[17] + 
x*(as[18] +x*(as[19] + x*(as[20] +x*(as[21] + x*(as[22] +x*(as[23] + x*(as[24] +x*(as[25] + 
x*(as[26] +x*(as[27] + x*(as[28] +x*(as[29] + x*(as[30] +x*(as[31] + x*(as[32] +x*(as[33] + 
x*(as[34] +x*(as[35] + x*(as[36] +x*(as[37] + x*(as[38] +x*(as[39] + x*(as[40] +x*(as[41] + 
x*(as[42] +x*(as[43] + x*(as[44] +x*(as[45] + x*(as[46] +x*(as[47] + x*(as[48] +x*(as[49] + 
x*(as[50] +x*(as[51] + x*(as[52] +x*(as[53] + x*(as[54] +x*(as[55] + x*(as[56] +x*(as[57] + 
x*(as[58] +x*(as[59] + x*(as[60] +x*(as[61] + x*(as[62] +x*(as[63] + x*(as[64] +x*(as[65] + 
x*(as[66] +x*(as[67] + x*(as[68] +x*(as[69] + x*(as[70] +x*(as[71] + x*(as[72] +x*(as[73] + 
x*(as[74] +x*(as[75] + x*(as[76] +x*(as[77] + x*(as[78] +x*(as[79] + x*(as[80] +x*(as[81] + 
x*(as[82] +x*(as[83] + x*(as[84] +x*(as[85] + x*(as[86] +x*(as[87] + x*(as[88] +x*(as[89] + 
x*(as[90] +x*(as[91] + x*(as[92] +x*(as[93] + x*(as[94] +x*(as[95] + x*(as[96] +x*(as[97] + 
x*(as[98] +x*(as[99] + x*(as[100] +x*(as[101] + x*(as[102] +x*(as[103] + x*(as[104] +x*(as[105] + 
x*(as[106] +x*(as[107] + x*(as[108] +x*(as[109] + x*(as[110] +x*(as[111] + x*(as[112] +x*(as[113] + 
x*(as[114] +x*(as[115] + x*(as[116] +x*(as[117] + x*(as[118] +x*(as[119] + x*(as[120] +x*(as[121] + 
x*(as[122] +x*(as[123] + as[124]*x)))))))))))))))))))))))))))))))))))))))))))))))))))))))
))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))));
}


//inline FloatType ghGammaH2(FloatType q)
inline FloatType ghGammaH2(FloatType q)
{
 
  FloatType x = q; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);

   FloatType  x2 = x*x;
   		
   std::array<FloatType, 11> a;
   
   a[0]= as[0]+as[1]*x+x2*(as[2]+x*(as[3]+x*(as[4]+x*(as[5]+x*(as[6]+x*(as[7]+x*(as[8]+x*(as[9]+x*(as[10]+as[11]*x)))))))));
   a[1]= as[12]+as[13]*x+x2*(as[14]+x*(as[15]+x*(as[16]+x*(as[17]+x*(as[18]+x*(as[19]+x*(as[20]+x*(as[21]+x*(as[22]+as[23]*x)))))))));
   a[2]= as[24]+as[25]*x+x2*(as[26]+x*(as[27]+x*(as[28]+x*(as[29]+x*(as[30]+x*(as[31]+x*(as[32]+x*(as[33]+x*(as[34]+as[35]*x)))))))));
   a[3]= as[36]+as[37]*x+x2*(as[38]+x*(as[39]+x*(as[40]+x*(as[41]+x*(as[42]+x*(as[43]+x*(as[44]+x*(as[45]+x*(as[46]+as[47]*x)))))))));
   a[4]= as[48]+as[49]*x+x2*(as[50]+x*(as[51]+x*(as[52]+x*(as[53]+x*(as[54]+x*(as[55]+x*(as[56]+x*(as[57]+x*(as[58]+as[59]*x)))))))));
   a[5]= as[60]+as[61]*x+x2*(as[62]+x*(as[63]+x*(as[64]+x*(as[65]+x*(as[66]+x*(as[67]+x*(as[68]+x*(as[69]+x*(as[70]+as[71]*x)))))))));
   			
   a[6]= as[72]+as[73]*x+x2*(as[74]+x*(as[75]+x*(as[76]+x*(as[77]+x*(as[78]+x*(as[79]+x*(as[80]+x*(as[81]+x*(as[82]+as[83]*x)))))))));
   a[7]= as[84]+as[85]*x+x2*(as[86]+x*(as[87]+x*(as[88]+x*(as[89]+x*(as[90]+x*(as[91]+x*(as[92]+x*(as[93]+x*(as[94]+as[95]*x)))))))));
   a[8]= as[96]+as[97]*x+x2*(as[98]+x*(as[99]+x*(as[100]+x*(as[101]+x*(as[102]+x*(as[103]+x*(as[104]+x*(as[105]+x*(as[106]+as[107]*x)))))))));
   a[9]= as[108]+as[109]*x+x2*(as[110]+x*(as[111]+x*(as[112]+x*(as[113]+x*(as[114]+x*(as[115]+x*(as[116]+x*(as[117]+x*(as[118]+as[119]*x)))))))));
   a[10]= as[120]+as[121]*x+x2*(as[122]+x*(as[123]+x*(as[124])));
   
   FloatType  z = x2*x2*x2*x2*x2*x2;

   return factor*(a[0]+a[1]*z+ z*z*(a[2] + z*(a[3]+z*(a[4]+z*(a[5]+z*(a[6]+z*(a[7]+z*(a[8]+z*(a[9]+z*(a[10]))))))))));

}

//inline FloatType ghGammaH3(FloatType q)
inline FloatType ghGammaH3(FloatType q)
{
 
 FloatType x = q; 
	FloatType factor = FloatType(1);
 	shiftAndFactor(&x, &factor);
		
	 
	 FloatType  x2 = x*x;
  			
   std::array<FloatType, 25> a;
   std::array<FloatType,  5> b;

   a[ 0] = (as[ 0] + as[ 1]*x + x2* (as[ 2] + x*(as[ 3] + as[ 4]*x)));
   a[ 1] = (as[ 5] + as[ 6]*x + x2* (as[ 7] + x*(as[ 8] + as[ 9]*x)));
   a[ 2] = (as[10] + as[11]*x + x2* (as[12] + x*(as[13] + as[14]*x)));
   a[ 3] = (as[15] + as[16]*x + x2* (as[17] + x*(as[18] + as[19]*x)));
   a[ 4] = (as[20] + as[21]*x + x2* (as[22] + x*(as[23] + as[24]*x)));
   a[ 5] = (as[25] + as[26]*x + x2* (as[27] + x*(as[28] + as[29]*x)));
   a[ 6] = (as[30] + as[31]*x + x2* (as[32] + x*(as[33] + as[34]*x)));
   a[ 7] = (as[35] + as[36]*x + x2* (as[37] + x*(as[38] + as[39]*x)));
   a[ 8] = (as[40] + as[41]*x + x2* (as[42] + x*(as[43] + as[44]*x)));
   a[ 9] = (as[45] + as[46]*x + x2* (as[47] + x*(as[48] + as[49]*x)));
   a[10] = (as[50] + as[51]*x + x2* (as[52] + x*(as[53] + as[54]*x)));
   a[11] = (as[55] + as[56]*x + x2* (as[57] + x*(as[58] + as[59]*x)));
   a[12] = (as[60] + as[61]*x + x2* (as[62] + x*(as[63] + as[64]*x)));
   a[13] = (as[65] + as[66]*x + x2* (as[67] + x*(as[68] + as[69]*x)));
   a[14] = (as[70] + as[71]*x + x2* (as[72] + x*(as[73] + as[74]*x)));
   a[15] = (as[75] + as[76]*x + x2* (as[77] + x*(as[78] + as[79]*x)));
   a[16] = (as[80] + as[81]*x + x2* (as[82] + x*(as[83] + as[84]*x)));
   a[17] = (as[85] + as[86]*x + x2* (as[87] + x*(as[88] + as[89]*x)));
   a[18] = (as[90] + as[91]*x + x2* (as[92] + x*(as[93] + as[94]*x)));
   a[19] = (as[95] + as[96]*x + x2* (as[97] + x*(as[98] + as[99]*x)));
   a[20] = (as[100] + as[101]*x + x2* (as[102] + x*(as[103] + as[104]*x)));
   a[21] = (as[105] + as[106]*x + x2* (as[107] + x*(as[108] + as[109]*x)));
   a[22] = (as[110] + as[111]*x + x2* (as[112] + x*(as[113] + as[114]*x)));
   a[23] = (as[115] + as[116]*x + x2* (as[117] + x*(as[118] + as[119]*x)));
   a[24] = (as[120] + as[121]*x + x2* (as[122] + x*(as[123] + as[124]*x)));

	FloatType  y = x*x2*x2;
  FloatType  y2 = y*y;
   
   b[0] = (a[ 0] + a[ 1]*y + y2* (a[ 2] + y*(a[ 3] + a[ 4]*y)));
   b[1] = (a[ 5] + a[ 6]*y + y2* (a[ 7] + y*(a[ 8] + a[ 9]*y)));
   b[2] = (a[10] + a[11]*y + y2* (a[12] + y*(a[13] + a[14]*y)));
   b[3] = (a[15] + a[16]*y + y2* (a[17] + y*(a[18] + a[19]*y)));
   b[4] = (a[20] + a[21]*y + y2* (a[22] + y*(a[23] + a[24]*y)));

	FloatType  z = y*y2*y2;
  FloatType  z2 = z*z;
	    
  return  factor*(b[0] + b[1]*z + z2* (b[2] + z*(b[3] + b[4]*z)));
}

#include "trapez.hpp"
