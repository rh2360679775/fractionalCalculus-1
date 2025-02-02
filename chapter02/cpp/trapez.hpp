#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <omp.h>
using namespace std;
using namespace chrono;


int main(int argc, char * argv[]) {
  int pwr = (argc > 1) ? atoi(argv[1]) : 7;  
  int intervals = pow(10, pwr) ;      // number of intervals 10**pwr
    
  auto start = steady_clock::now();
  
  // trapez rule : integrate gamma from 0 to 1  
  //FloatType total = trapez( ghGamma, FloatType(18.5), FloatType(19.5), intervals);

  if( omp_get_max_threads() < 4 ){
    cout << "maximum  available threads <4 , stopping this calculation"<<endl;
    exit(-1);
  }

// Parallel region with 4 threads
   
  FloatType total1, total2, total3, total4; 

    #pragma omp parallel num_threads(4)
    {
        int thread_id = omp_get_thread_num();
        if (thread_id == 0) {
              total1 = trapez( ghGamma, FloatType(0), FloatType(0.5), intervals/4);
        } 
        else if (thread_id == 1) {
              total2 = trapez( ghGamma, FloatType(0.5), FloatType(1.0), intervals/4);
        }
        else if (thread_id == 2) {
              total3 = trapez( ghGamma, FloatType(0.0), FloatType(0.5), intervals/4, 1);
        }
        else if (thread_id == 3) {
              total4 = trapez( ghGamma, FloatType(0.5), FloatType(1.0), intervals/4, 1);
        }
    }

  FloatType total = 2*(total1 + total2)  + total3 + total4;
  
  auto stop = steady_clock::now();
  duration < double, milli > elapsed_ms = stop - start;
  
  cout << " trapez " << "integral = " << fixed << setprecision(prec) << total <<
    ",intervals 10**" << pwr  <<
    " time " << fixed << setprecision(3) << elapsed_ms.count() << " ms" << endl;
  return 0;
}