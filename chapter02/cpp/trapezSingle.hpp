#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace std;
using namespace chrono;


int main(int argc, char * argv[]) {
  int pwr = (argc > 1) ? atoi(argv[1]) : 7;  
  int up = (argc > 2) ? atoi(argv[2]) : 1;  
  int intervals = pow(10, pwr) ;      // number of intervals 10**pwr
    
  auto start = steady_clock::now();
  
  // trapez rule : integrate gamma from 0 to 1  
  FloatType og = FloatType(up);
  FloatType total = trapez( ghGamma, FloatType(0), og, intervals);
  auto stop = steady_clock::now();
  duration < double, milli > elapsed_ms = stop - start;
  
  cout << " startSer " << "integral = " << fixed << setprecision(prec) << total <<
    ",intervals 10**" << pwr  <<
    ",og " << setprecision(3) << og  <<
    " time " << fixed << setprecision(3) << elapsed_ms.count() << " ms" << endl;
  return 0;
}