#include <iostream>
#include "globalPadeMLF.hpp"

using namespace std;

int main() {
    // create an instance of the class GlobalPadeMLF
  GlobalPadeMLF mlfPadeApprox;  // constructor

  double alpha{0.9};
  double beta{1.5};
  double x{-1};
    // exact value of E_{9/10, 3/2}(-1) from CAS
  double exact{0.59595802527072791093339988837073}; 

    // print title and column headers
  cout << "E("<<alpha<<","<<beta<<","<<x<<")" << endl;
  cout << "n log10(abs(err))" << endl;

    // print error for increasing orders n=2...15
  for (int n = 2; n < 15; n++) {
    mlfPadeApprox.setOrder(n);  
    double approx = mlfPadeApprox.mittagLefflerE(alpha, beta, -1.0); 
    double relErr = ( approx - exact ) / exact; 
    cout << n << " " << log10(abs(relErr)) << endl; 
  }

  return 0;
}