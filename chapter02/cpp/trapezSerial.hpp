#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
using namespace std;
using namespace chrono;


int main(int argc, char * argv[]) {
  int pwr = (argc > 1) ? atoi(argv[1]) : 7;  
  int intervals = pow(10, pwr) ;      // number of intervals 10**pwr
    
  auto start = steady_clock::now();
  
  // trapez rule : integrate gamma from 0 to 1  
  //FloatType total = trapez( ghGamma, FloatType(18.5), FloatType(19.5), intervals);


// Serial region    
  FloatType total1, total2; 
  int i2 = intervals/2;
    total1 = trapez( ghGamma, FloatType(0.0), FloatType(1.0), i2);
    total2 = trapez( ghGamma, FloatType(0.0), FloatType(1.0), i2, 1);
  
  FloatType total = 2*total1 + total2;
  
  auto stop = steady_clock::now();
  duration < double, milli > elapsed_ms = stop - start;
  
  cout << " trapez " << "integral = " << fixed << setprecision(prec) << total <<
    ",intervals 10**" << pwr  <<
    " time " << fixed << setprecision(3) << elapsed_ms.count() << " ms" << endl;
  return 0;
}