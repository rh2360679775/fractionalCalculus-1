#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <omp.h>
#include <cmath>
using namespace std;
	long double as[30] = {
		1.0E0L, 
		-5.772156649015328605867504180524819743117047335253698707946215782E-1L, 
		9.890559953279725480567736016011074551731970227871841326147670037E-1L, 
		-9.074790760808853204034041123531852470048241450183333319688654403E-1L, 
		9.817280868343342993508170961593701086734340526645559305830260364E-1L, 
		-9.819950689004158144339628200322199063804287933962962701660364352E-1L, 
		9.931491145454546596782858789496350328652763778299599764441614315E-1L, 
		-9.960017589382732748089175162811376314697147473176727548253925349E-1L, 
		9.981056715129189680990573784917422247740403929981223999668834399E-1L, 
		-9.990250134118399800372439103878567997541294491709751397892570239E-1L, 
		9.995133628162019845482052880889581723670991136329077245238140199E-1L, 
		-9.997399284276842800847569140354098696960804838510334107125689736E-1L, 
		9.997791255148989102928786523422115508693106406801813139004551520E-1L, 
		-9.994504397726759377378756316403384848980389881027808702884505848E-1L, 
		9.979535245458726816208470913414942086863398570568581090328391650E-1L, 
		-9.929609973542025129822480482007778526218927585141017688579208531E-1L, 
		9.791739016563532773539533100165369613035793704671064296668045075E-1L, 
		-9.471479549842206850430750762894205529764109064566593394834907327E-1L, 
		8.843618782148606539970266389247349247667405257184289291251064020E-1L, 
		-7.804429582293780035462100559354943723915093847420921135949438237E-1L, 
		6.354645773560346047205100084882105703634791106498237509281422778E-1L, 
		-4.655773241310680445401866356009338719839776846509169710182875371E-1L, 
		2.993245854990955480240411051361101341114304008683617635575278364E-1L, 
		-1.646081873675282259931192131477486890184996839736595319623786902E-1L, 
		7.531660607553232952371244924139683428299053155028415770684736347E-2L, 
		-2.774103935426569153750679132953790742164427656251619171270142394E-2L, 
		7.871475255988880238358932061028742971085369173140485201668656701E-3L, 
		-1.610055659291006362895556315972035309013456479231051353987482957E-3L, 
		2.108086307668984701877376536205400787900905232817074763827402896E-4, 
		-1.324627302460587650246889364226604892336483855819419952292E-5L};
			
	


inline long double ghGammaR(long double x)
{
  // infinity-norm
	
	double result =   as[29];
	result = result*x+as[28];
	result = result*x+as[27];
	result = result*x+as[26];
	result = result*x+as[25];
	result = result*x+as[24];
	result = result*x+as[23];
	result = result*x+as[22];
	result = result*x+as[21];
	result = result*x+as[20];
	result = result*x+as[19];
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
	
    return result;
}

inline long double ghGammaH1(long double x)
{
 				
		return as[0] + as[1]*x + x*x*(as[2] + x*(as[3] + 
         x*(as[4] + x*  (as[5] + 	x*(as[6] + x*(as[7] + x*(as[8] +x*(as[9] + 	x*	(as[10] +x*	(as[11] +x*(as[12] +x*(as[13] +x*(as[14] +	
		 x*	(as[15] +x*	(as[16] +x*		(as[17] +	x*	(as[18] +	x*		(as[19] +		x*
		(as[20] +x*	(as[21] +	x*	(as[22] +	x*	(as[23] +x*	(as[24] +x*		(as[25] +x*	(as[26] +x*	(as[27] +x*
		(as[28] +x*	as[29] )))))))))))))))))))))))))));
}

inline long double ghGamma(long double x)
{
  // infinity-norm
		
   long double x2 = x*x;
   long double x6 = x2*x2*x2;
  			
   long double a[5];
			
   a[0] = as[ 0]  + as[ 1]*x  + x2*(as[ 2]  + x*(as[ 3]  + x*(as[ 4]  + x* as[ 5])));
   a[1] = as[ 6]  + as[ 7]*x  + x2*(as[ 8]  + x*(as[ 9]  + x*(as[10]  + x* as[11])));
   a[2] = as[12]  + as[13]*x  + x2*(as[14]  + x*(as[15]  + x*(as[16]  + x* as[17])));
   a[3] = as[18]  + as[19]*x  + x2*(as[20]  + x*(as[21]  + x*(as[22]  + x* as[23])));
   a[4] = as[24]  + as[25]*x  + x2*(as[26]  + x*(as[27]  + x*(as[28]  + x* as[29])));
   
   return a[0]  + a[1]*x6  + x6*x6*(a[2]  + x6*(a[3]  + x6*a[4] ));

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
    cout << "omp s30L integral = " << fixed << setprecision(36) << total
	     <<",intervals "<< intervals
		 <<",threads "  << threads
		 <<" time "     << fixed << setprecision(3) << elapsed_ms.count() << " ms"<< endl;
	//cout << "omM serD g(1+2.1) = " << ghGamma(2.1)<< endl;
	//cout << "omM serD g(1-0.1) = " << ghGamma(-0.1)<< endl;
	//cout << "omM serD g(1-0.9) = " << ghGamma(-0.9)<< endl;
	return 0;
}
