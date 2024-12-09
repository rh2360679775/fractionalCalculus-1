#include <iomanip>
#include <iostream>
#include "globalPadeMLF.hpp"

using namespace std;

int main() {
  GlobalPadeMLF mlfPadeApprox;

  double closeTo1[5] = {   // alpha = 1.0 - 10^k, beta = 1, x = -0.5
    0.60340549869586096761551528315510,
    0.60608995263141647835498381023141,
    0.60648529133691131557613292598209,
    0.60652610988754118255503652642507,
	  0.60653020460023866289227390850454
  };
  double diff;
  // try 1.0 and then 0.3
  mlfPadeApprox.setThresh(1.0);
  
  for (int i = 2; i < 20; i++) {
    diff = 0.1;
    mlfPadeApprox.setOrder(i);
    for (int j = 0; j < 5;j++){
      cout << i << " " << 1.0-diff << " " << mlfPadeApprox.mittagLefflerE(1.0-diff, 1.0, -0.5) - closeTo1[j] <<  endl;
      diff *= 0.1;
    }
  }

  
  return 0;
}