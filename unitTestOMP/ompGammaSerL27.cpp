#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <omp.h>
#include <cmath>
using namespace std;
	long double as[27] = {
		1.0E0, 
		-5.772156649015328574847107526236417396399473414593683294694448354E-1, 
		9.890559953279716055879006187112349085181194245439605371493621847E-1, 
		-9.074790760807840913849637364128458684942765833469234123768037827E-1, 
		9.817280868287493872829494368343906581423463320746975472327806233E-1, 
		-9.819950687133308024036797799970825430041145353647966988985601196E-1, 
		9.931491103560720422277583602583809840640935608838069770575025500E-1, 
		-9.960016921875578703578265696608380067250364412371040786498723143E-1, 
		9.981048811182184740194832873848235437462323425587398485676418276E-1, 
		-9.990178339256502246393676805933536852541867471072467946582140434E-1, 
		9.994621279751838577149190275163130262747245024804575366693054501E-1, 
		-9.994473777770612291622853078667857569155792808053497759251338579E-1, 
		9.984234865818553609719362575102820190846558249735448507547378456E-1, 
		-9.942964224792548885499600606423283534205252993190213223445278409E-1, 
		9.817415529373300718490015920863342801408285485623247052344540861E-1, 
		-9.505087356260045072671857543877141073379640221543396483474575002E-1, 
		8.862340777038458017327321534931532223720878671817192171976377175E-1, 
		-7.766020439022703527799336629113315125364278308070082251553428261E-1, 
		6.218261781120564563123419874557456000659195749648797901211510819E-1, 
		-4.417060348412027926063218282522105999297358388447257554947061329E-1, 
		2.701475353497292098269379690209491963786354304362498578566182867E-1, 
		-1.378768248000485686016129493050214125669796790493931341733992288E-1, 
		5.666741300908340721257290134619610032425452294036514887478396272E-2, 
		-1.791721124454662588127948572017495010218971969715187449742012969E-2, 
		4.073484259812306386189879368456543383087672500179665265686023159E-3, 
		-5.909388467936236862167480758102751061875702154896091292021547692E-4, 
		4.099576613045368062084546287872853816528530699297713842313E-5};
			
	


inline long double ghGamma0(long double z)
{
  
	long double x = z; 
	long double factor = 1.0;

  while(x > 1.0L){ factor *= x; x -= 1.0; } 
  while(x < 0.0L){ factor /= 1.0L+x; x += 1.0L; } 
		long double sum = as[0];   
		long double xn = 1.0L;

		for(int n = 1; n < 27; n++){
			xn *= x;
			sum += as[n]*xn;
		}
		return factor*sum;
}

inline long double ghGamma(long double z)
{

	long double x = z; 
	long double factor = 1.0;

  while(x > 1.0L){ factor *= x; x -= 1.0; } 
  while(x < 0.0L){ factor /= 1.0L+x; x += 1.0L; } 
			
		return factor*(as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + 
         x*(as[4] + x*  (as[5] + 	x*(as[6] + x*(as[7] + x*(as[8] +x*(as[9] + 	x*	(as[10] +x*	(as[11] +x*(as[12] +x*(as[13] +x*(as[14] +	
		 x*	(as[15] +x*	(as[16] +x*		(as[17] +	x*	(as[18] +	x*		(as[19] +		x*
		(as[20] +x*	(as[21] +	x*	(as[22] +	x*	(as[23] +x*	(as[24] +x*		(as[25] +x*	as[26])))))))))))))))))))))))));
}

inline long double ghGammaH2(long double z)
{
  // infinity-norm 74

	long double x = z; 
	long double factor = 1.0;

  while(x > 1.0L){ factor *= x; x -= 1.0; } 
  while(x < 0.0L){ factor /= 1.0L+x; x += 1.0L; } 

   long double x7 = x*x*x*x*x*x*x;
			
   long double a[4];
			
   a[0] = as[0]  + as[1]*x  + x*x*(as[2]  + x*(as[3]  + x*(as[4]  + x* (as[5]  + x*as[6]))));
   a[1] = as[7]  + as[8]*x  + x*x*(as[9]  + x*(as[10] + x*(as[11] + x* (as[12] + x*as[13]))));
   a[2] = as[14] + as[15]*x + x*x*(as[16] + x*(as[17] + x*(as[18] + x* (as[19] + x*as[20]))));
   a[3] = as[21] + as[22]*x + x*x*(as[23] + x*(as[24] + x*(as[25] + x* as[26] )));
   
   return factor*(a[0]  + a[1]*x7  + x7*x7*(a[2]  + x7*a[3]));

}



inline long double ghGammaH3(long double z)
{
  // infinity-norm
  // 333

	long double x = z; 
	long double factor = 1.0;

  while(x > 1.0L){ factor *= x; x -= 1.0; } 
  while(x < 0.0L){ factor /= 1.0L+x; x += 1.0L; } 

			
   long double a[9];
   a[0] = as[0]  + as[1]*x  + x*x*as[2];
   a[1] = as[3]  + as[4]*x  + x*x*as[5];
   a[2] = as[6]  + as[7]*x  + x*x*as[8];
   a[3] = as[9]  + as[10]*x  + x*x*as[11];
   a[4] = as[12]  + as[13]*x  + x*x*as[14];
   a[5] = as[15]  + as[16]*x  + x*x*as[17];
   a[6] = as[18]  + as[19]*x  + x*x*as[20];
   a[7] = as[21]  + as[22]*x  + x*x*as[23];
   a[8] = as[24]  + as[25]*x  + x*x*as[26];
  
   long double b[3];
   long double x3 = x*x*x;
   b[0] = a[0]  + a[1]*x3  + x3*x3*a[2];
   b[1] = a[3]  + a[4]*x3  + x3*x3*a[5];
   b[2] = a[6]  + a[7]*x3  + x3*x3*a[8];
   
   long double x9 = x3*x3*x3;
   return factor*(b[0]  + b[1]*x9  + x9*x9*b[2]);

}



inline long double ghGammaLib(long double x)
{
  // ibrary function
		return tgammal( 1.0L + x );
}

int main(int argc,char *argv[])
{
	int intervals = (argc > 1) ? pow(10,atoi(argv[1])) : pow(10,7);                // number of intervals 
	int threads   = (argc > 2) ? atoi(argv[2]) : omp_get_max_threads();   // omp threads
		
	long double a = 0.0L;               // lower limit of integral
	long double b = 1.0L;               // upper limit of integral
	long double h = (b-a)/intervals;    // intervalsize

    auto start = chrono::steady_clock::now();
 
 	// trapez rule begins
	
	// start with the endpoints , so we can treat all values simultaneously
	long double total = -0.5L*( ghGamma(a) - ghGamma(b) );
   	
	omp_set_num_threads(threads);                                             
    
	#pragma omp parallel for reduction (+:total)  // OpenMP
	for( int i = 0; i < intervals; i++ )
	{
		long double x = a + h*i;     // thread safe
		total += ghGamma(x);   // get function values f(x), a <= f(x) <= b
	}
	
	total *= h ;
	
	// trapez end
    auto end = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_ms = end - start;
    cout << "omp s27L integral = " << fixed << setprecision(36) << total
	     <<",intervals "<< intervals
		 <<",threads "  << threads
		 <<" time "     << fixed << setprecision(3) << elapsed_ms.count() << " ms"<< endl;
	//cout << "omM serD g(1+2.1) = " << ghGamma(2.1)<< endl;
	//cout << "omM serD g(1-0.1) = " << ghGamma(-0.1)<< endl;
	//cout << "omM serD g(1-0.9) = " << ghGamma(-0.9)<< endl;
	return 0;
}
