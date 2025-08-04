#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;

int main(int argc, char * argv[]) {
  int pwr = (argc > 1) ? atoi(argv[1]) : 7;  
  int intervals = pow(10, pwr) ;      // number of intervals 10**pwr
    
  auto start = steady_clock::now();
  
  // trapezoidal rule : integrate gamma from 0 to 2  
  FloatType a = double(0);
  FloatType b = double(2);
  FloatType total = trapez( ghGamma, a, b, intervals );
  
  auto stop = steady_clock::now();
  duration < double, milli > elapsed_ms = stop - start;

   // write results 
  cout << ",intervals 10**" << pwr  << endl;
  cout << " time " << elapsed_ms.count() << " ms" << endl;
  cout << " total = " << setprecision(prec) << total << endl;
  return 0;
}