#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>                                                        // openMP
using namespace std;

int main(int argc, char * argv[]) {
  int intervals = (argc > 1) ? pow(10, atoi(argv[1])) : pow(10, 7); // number of intervals 
  int threads = (argc > 2) ? atoi(argv[2]) : omp_get_max_threads(); // openMP threads

  FloatType a = FloatType(0); // lower limit of integral
  FloatType b = FloatType(1); // upper limit of integral
  FloatType h = (b - a) / intervals; // intervalsize

  auto start = chrono::steady_clock::now();

  // trapez rule begins
  // total = sum of all
  // start with the endpoints , so we can treat all values simultaneously
  FloatType total = -FloatType(0.5) * (ghGamma(a) - ghGamma(b));

  omp_set_num_threads(threads); // openMP                            
  #pragma omp parallel for reduction(+: total) // openMP

  for (int i = 0; i < intervals; i++) {
    FloatType x = a + h * i; // threads !
    total += ghGamma(x); // get function values f(x), a <= f(x) <= b
  }

  total *= h; // sum times interval size

  // trapez rule end
  auto stop = chrono::steady_clock::now();
  chrono::duration < double, milli > elapsed_ms = stop - start;
  cout << name << "  integral = " << fixed << setprecision(prec) << total <<
    ",intervals " << intervals << ",threads " << threads <<
    " time " << fixed << setprecision(3) << elapsed_ms.count() << " ms" << endl;
  return 0;
}