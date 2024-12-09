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

 // infinity-norm
 std::array<cpp_dec_float_50, 44>
		c = {cpp_dec_float_50("1.0E0"), 
		cpp_dec_float_50("-5.772156649015328606065120900811147869971364976889090534598610754E-1"), 
		cpp_dec_float_50("9.890559953279725553953956506240423671202674552699000416815617590E-1"), 
		cpp_dec_float_50("-9.074790760808862890165599403326630695172499443566639141228350735E-1"), 
		cpp_dec_float_50("9.817280868344001873363489528342861049688256032190769915794711079E-1"), 
		cpp_dec_float_50("-9.819950689031452021020117320612718583442918449183142167585530299E-1"), 
		cpp_dec_float_50("9.931491146212761929967542720767229875051768187033365212431423127E-1"), 
		cpp_dec_float_50("-9.960017604424315273348899829193040711937972812246485390362523714E-1"), 
		cpp_dec_float_50("9.981056937831287104598267998312028488211738882849988309136402569E-1"), 
		cpp_dec_float_50("-9.990252676219496108548168120879641529493498958990440051922734814E-1"), 
		cpp_dec_float_50("9.995156560726729950352311527499018806138949805332372083099874680E-1"), 
		cpp_dec_float_50("-9.997565975069089326147533648042511793930733109784959623289806169E-1"), 
		cpp_dec_float_50("9.998782712924048749693272288003650387791672041394025956092421232E-1"), 
		cpp_dec_float_50("-9.999390639500361838376302390162654659981997834270604580313626238E-1"), 
		cpp_dec_float_50("9.999695153063227517196779222021399032838314435617313279370236657E-1"), 
		cpp_dec_float_50("-9.999847325077349843270781730947969152973048166685092926261672694E-1"), 
		cpp_dec_float_50("9.999922310501464132856505045954559614279785743407558980033752476E-1"), 
		cpp_dec_float_50("-9.999952997975882015030905231681311778395598749375787960906202725E-1"), 
		cpp_dec_float_50("9.999932923329043915545626014247801736970172024923877560396011216E-1"), 
		cpp_dec_float_50("-9.999761726410551623318833114330639033103128152258282916461652035E-1"), 
		cpp_dec_float_50("9.999031699825277916545840475064771436502247872354775478033814950E-1"), 
		cpp_dec_float_50("-9.996395079911130131293345845236534154532310050646656557410541415E-1"), 
		cpp_dec_float_50("9.988001399465882671230072300724713571671982836349981339942501656E-1"), 
		cpp_dec_float_50("-9.964299086378231916970242038424730705692064898468098244450242840E-1"), 
		cpp_dec_float_50("9.904803721945194828051480771783097059273366309925329755674642882E-1"), 
		cpp_dec_float_50("-9.771903959668231660001650344209014376844610541065140171146703776E-1"), 
		cpp_dec_float_50("9.507566285290634152273554068635265441805657963383035505439166131E-1"), 
		cpp_dec_float_50("-9.039394503087400748422098261186427985594552404103133133431830190E-1"), 
		cpp_dec_float_50("8.301408601281595784518412615615398433984511210128695422756952569E-1"), 
		cpp_dec_float_50("-7.267156509987080673453942827462472722157318814528795473944927099E-1"), 
		cpp_dec_float_50("5.980656283176606038535497894497100212984050287593212555580102086E-1"), 
		cpp_dec_float_50("-4.563668566500621862462690466158450409222802664765250655445463874E-1"), 
		cpp_dec_float_50("3.186094143974604670922360433556823785326295951609043933507911320E-1"), 
		cpp_dec_float_50("-2.008851614419559675511242425558453969449382800634259878370712186E-1"), 
		cpp_dec_float_50("1.129163100266353906085243110807774649562870582986292181201189752E-1"), 
		cpp_dec_float_50("-5.582222178976752061541203950456824258556878363518235246017163945E-2"),
		cpp_dec_float_50("2.391142886993550611650414375338213728346144145032559139356178414E-2"), 
		cpp_dec_float_50("-8.720473300827586927801133083738729174457855062108678086878403895E-3"),
	    cpp_dec_float_50("2.649319597817444957091198369908709540187229812679815061870225761E-3"), 
		cpp_dec_float_50("-6.513625508929841713853194217456528556405311677773251464098718441E-4"),
		cpp_dec_float_50("1.243622417930514061732771456020937579791277533668802765005137035E-4"), 
		cpp_dec_float_50("-1.728310866146560966393780006667255436587181342713110952039735320E-5"), 
		cpp_dec_float_50("1.554063789826756803020042506421929437083429409306232969482125683E-6"), 
		cpp_dec_float_50("-6.781853572014058307188623677747988302999896395397948173E-8")};


inline cpp_dec_float_50 ghGamma(cpp_dec_float_50 z)
{
 	cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
		
		return factgor*(as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + 
         x*(as[4] + x*  (as[5] + 	x*(as[6] + x*(as[7] + 		x*		(as[8] + 		x*		(as[9] + 		x*		(as[10] +		x*		(as[11] +		x*
				(as[12] +		x*		(as[13] +		x*		(as[14] +		x*		(as[15] +		x*		(as[16] +		x*		(as[17] +		x*		(as[18] +		x*		(as[19] +		x*
				(as[20] +		x*		(as[21] +		x*		(as[22] +		x*		(as[23] +		x*		(as[24] +		x*		(as[25] +		x*		(as[26] +		x*		(as[27] +		x*
				(as[28] +		x*		(as[29] +		x*		(as[30] +		x*		(as[31] +		x*		(as[32] +		x*		(as[33] +		x*		(as[34] +		x*		(as[35] +		x*
				(as[36] +		x*		(as[37] +		x*		(as[38] +		x*		(as[39] +		x*		
				(as[40] +		x*		(as[41] +	x*(as[42] + as[43]*x)))))))))))))))))))))))))))))))))))))))));
}


inline cpp_dec_float_50 ghGammaH2(cpp_dec_float_50 z)
{
 
  cpp_dec_float_50 x = z; 
	cpp_dec_float_50 val0 =  cpp_dec_float_50(0);
	cpp_dec_float_50 val1 =  cpp_dec_float_50(1);
	cpp_dec_float_50 factor = val1;

  while(x > val1){ factor *= x; x -= val1; } 
  while(x < val0){ factor /= val1+x; x += val1; } 
   cpp_dec_float_50  x2 = x*x;
   cpp_dec_float_50  x4 = x2*x2;
   cpp_dec_float_50  x7 = x4*x2*x;
			
   std::array<cpp_dec_float_50, 7> a;
			
   a[0] = as[0]  + as[1]*x  + x2*(as[2]  + x*(as[3]  + x*(as[4]  + x* (as[5]  + x*as[6]))));
   a[1] = as[7]  + as[8]*x  + x2*(as[9]  + x*(as[10] + x*(as[11] + x* (as[12] + x*as[13]))));
   a[2] = as[14] + as[15]*x + x2*(as[16] + x*(as[17] + x*(as[18] + x* (as[19] + x*as[20]))));
   a[3] = as[21] + as[22]*x + x2*(as[23] + x*(as[24] + x*(as[25] + x* (as[26] + x*as[27]))));
   a[4] = as[28] + as[29]*x + x2*(as[30] + x*(as[31] + x*(as[32] + x* (as[33] + x*as[34]))));
   a[5] = as[35] + as[36]*x + x2*(as[37] + x*(as[38] + x*(as[39] + x* (as[40] + x*as[41]))));
   a[6] = as[42] + as[43]*x ;
   
   return factor*(a[0]  + a[1]*x7  + x7*x7*(a[2]  + x7*(a[3]  + x7*(a[4]  + x7* (a[5]  + x7*a[6])))));

}



int main(int argc,char *argv[])
{
	int steps   = (argc > 1) ? atoi(argv[1]) : 100000000; // get command
	int threads = (argc > 2) ? atoi(argv[2]) : 4;        // omp threads
	cx::timer tim;
	
	cpp_dec_float_50  og = cpp_dec_float_50(1.0);
	cpp_dec_float_50  step_size = og / (steps-1); // NB n-1 steps between n points
	cpp_dec_float_50  omp_sum = cpp_dec_float_50(0.0);

	omp_set_num_threads(threads);                   // OpenMP 
    #pragma omp parallel for reduction (+:omp_sum)  // OpenMP
	for(int step = 0; step < steps; step++){
		cpp_dec_float_50 x = step_size*step;
		omp_sum += ghGamma(x);   // get sum of Taylor series
	}
	double  cpu_time = tim.lap_ms(); // get elapsed time
	// Trapezoidal Rule correction for end points
	omp_sum -= cpp_dec_float_50(0.5)*(ghGamma(cpp_dec_float_50(0.0))+ghGamma(cpp_dec_float_50(og)));
	omp_sum *= step_size;
	
	std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);
    std::cout <<"omp booL integral = "<< omp_sum <<" steps "<< steps << " time "<< cpu_time<< " ms" << std::endl;

	return 0;
}