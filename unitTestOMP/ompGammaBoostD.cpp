#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "cxtimers.h"
#include <iostream>

//#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/gamma.hpp>

using boost::multiprecision::cpp_dec_float_50;
using namespace boost::multiprecision;
using namespace std;


   std::array<cpp_dec_float_50, 20>
		as = {1.0E0, 
		-5.7721566490112839913612516377983E-1, 
		9.8905599526004546488561968743339E-1, 
		-9.0747907205154988315825417112244E-1, 
		9.8172796446769040219142487649730E-1, 
		-9.8199282456119311153012614684595E-1, 
		9.9312179432071013302064800481304E-1, 
		-9.9576722393326028337686469083995E-1, 
		9.9662498968300581987064380393876E-1, 
		-9.9193816665161287149311120508649E-1, 
		9.7319991296155473144438879723251E-1, 
		-9.2254446683703155299124831330242E-1, 
		8.1800576403405562896572614793714E-1, 
		-6.5063455752880833774388636414585E-1, 
		4.4314628622691644741140621352038E-1, 
		-2.4602120185468968949442836220557E-1, 
		1.0527746780281570888539551947701E-1, 
		-3.2239556919065630956147031443838E-2, 
		6.2446860911293176575928532278774E-3, 
		-5.72125609583894452654455305E-4};	

inline cpp_dec_float_50 ghGamma(cpp_dec_float_50 z)
{
 	cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
		
	cpp_dec_float_50 sum = as[0];   
	cpp_dec_float_50 xn = cpp_dec_float_50(1);

		for(int n = 1; n < 20; n++){
			xn *= x;
			sum += as[n]*xn;
		}
		return factor*sum;
}
inline cpp_dec_float_50 ghGammaH1(cpp_dec_float_50 z)
{
  cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 

	return factor*(as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x* (as[5] + x*(as[6] +  x*(as[7] + x*  (as[8] +  x* (as[9] + 
                    x* (as[10] +  x*(as[11] + x*   (as[12] +  x* (as[13] +  x*(as[14] + x*   (as[15] + 
                    x*  (as[16] +  x* (as[17] +  x*(as[18] + as[19]*x))))))))))))) )))));
}


inline cpp_dec_float_50 ghGamma54(cpp_dec_float_50 z)
{
   
	cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
	 
	 
	 cpp_dec_float_50 a[5];
   a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + as[ 3]*x);
   a[1] = as[ 4] + as[ 5]*x + x*x*(as[ 6] + as[ 7]*x);
   a[2] = as[ 8] + as[ 9]*x + x*x*(as[10] + as[11]*x);
   a[3] = as[12] + as[13]*x + x*x*(as[14] + as[15]*x);
   a[4] = as[16] + as[17]*x + x*x*(as[18] + as[19]*x);
   
   return  factor*(a[0] + a[1]*x*x*x*x + x*x*x*x*x*x*x*x*(a[2] + x*x*x*x*(a[3] + a[4]*x*x*x*x)));
   
}

inline cpp_dec_float_50 ghGamma45(cpp_dec_float_50 z)
{
  cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
	

   cpp_dec_float_50 a[4];
   a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + x*(as[ 3] + as[ 4]*x));
   a[1] = as[ 5] + as[ 6]*x + x*x*(as[ 7] + x*(as[ 8] + as[ 9]*x));
   a[2] = as[10] + as[11]*x + x*x*(as[12] + x*(as[13] + as[14]*x));
   a[3] = as[15] + as[16]*x + x*x*(as[17] + x*(as[18] + as[19]*x));
   
   return factor*(a[0] + a[1]*x*x*x*x*x + x*x*x*x*x*x*x*x*x*x*(a[2] + a[3]*x*x*x*x*x));
   

}

inline cpp_dec_float_50 ghGamma37(cpp_dec_float_50 z)
{
  cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
	
   cpp_dec_float_50 a[3];
   
   a[0] = as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + x*(as[4] + x*(as[5] + as[6]*x))));
   a[1] = as[7] + as[8]*x + x*x*(as[9] + x*(as[10] + x*(as[11] + x*(as[12] + as[13]*x))));
   a[2] = as[14] + as[15]*x + x*x*(as[16] + x*(as[17] + x*(as[18] + x*as[19] )));
   
   return factor*(a[0] + a[1]*x*x*x*x*x*x*x + x*x*x*x*x*x*x*x*x*x*x*x*x*x*a[2] );
   

}

int main(int argc,char *argv[])
{
	int intervals = (argc > 1) ? pow(10,atoi(argv[1])) : pow(10,7);                // number of intervals 
	int threads   = (argc > 2) ? atoi(argv[2]) : omp_get_max_threads();   // omp threads

	;
		
	cpp_dec_float_50 a =  cpp_dec_float_50(0);               // lower limit of integral
	cpp_dec_float_50 b =  cpp_dec_float_50(1);               // upper limit of integral
	cpp_dec_float_50 h = (b-a)/intervals;                    // intervalsize

  auto start = chrono::steady_clock::now();
 
 	// trapez rule begins
	
	// start with the endpoints , so we can treat all values simultaneously
	cpp_dec_float_50 total = -cpp_dec_float_50(0.5)*(ghGamma(a) - ghGamma(b));

  omp_set_num_threads(threads); // openMP                            
  #pragma omp parallel for reduction(+: total) // openMP

  for (int i = 0; i < intervals; i++) {
    cpp_dec_float_50 x = a + h * i; // threads !
    total += ghGamma(x); // get function values f(x), a <= f(x) <= b
  }
	
	total *= h ;
	
	// trapez end
    auto end = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_ms = end - start;
    cout << "omN booD integral = " << fixed << setprecision(36) << total
	     <<",intervals "<< intervals
		 <<",threads "  << threads
		 <<" time "     << fixed << setprecision(3) << elapsed_ms.count() << " ms"<< endl;
	//cout << "omM serD g(1+2.1) = " << ghGamma(2.1)<< endl;
	//cout << "omM serD g(1-0.1) = " << ghGamma(-0.1)<< endl;
	//cout << "omM serD g(1-0.9) = " << ghGamma(-0.9)<< endl;
	return 0;
}
