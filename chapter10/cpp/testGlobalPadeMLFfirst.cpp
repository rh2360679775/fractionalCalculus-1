#include <iostream>
#include "globalPadeMLF.hpp"
#include <chrono>
#include <array>  // For counter (threadsafe)



using namespace std;

int main() {
    // create an instance of the class GlobalPadeMLF
  GlobalPadeMLF mlfPadeApprox;  

   mlfPadeApprox.setOrder(10);  
   
  
  cout << setprecision(32);
  

  double aa = 0.4;
  int cnt = 500;
  double delta = 1.0/(double )cnt;
  int cstart, cstop;

  auto start = std::chrono::steady_clock::now();


 mlfPadeApprox.setOrder(5);  
 
 mlfPadeApprox.setThresh(0.89);  

double p = 4;  
double steps = pow(10,p);
double total = 0;
double del = 1.0/steps;
const std::array<double, 8> lst0510 = {  // a = 0.5 , b = 1.0 
   0.61935048353641165864651948054595, // p=1
   0.64440426938301900671649421778058, // p=2
   0.64697308821739005202848649860626, // p=3
   0.64723060505683548115477251311871, // p=4
   0.647256363090395569743291226399,   // p=5
   0.647258938957247739605256642325,   // p=6
   0.647259196544567918202038954940,   // p=7
   0.6472593   // p=8
};

const std::array<double, 8> lst = {  // a = 0.9 , b = 1.0 
   0.60066246256512242540038297471695, // p=1
   0.62816430526792030755735528843546, // p=2
   0.63096625491132908348317583495221, // p=3
   0.63124696766894645199375982012512, // p=4
   0.6312750441226549137086949299756,  // p=5
   0.6,   // p=6
   0.6,   // p=7
   0.6 // p=8
};

 for (int i = 1; i <= steps; i++) {
  double x =  i * del;
 // total += MittagLefflerZAB(-x, 0.5,1.0);
  total += mlfPadeApprox.mittagLefflerE(0.9, 1.0, -x);
 }
 total /= steps;  

auto stop = std::chrono::steady_clock::now();
std::chrono::duration <double, milli> elapsed_ms = stop - start;
 
  cout << " steps 10**" << p << endl;
  cout << " total = " << total << endl;
  cout << " error = " << (total-lst[p-1])/total  << endl;
  cout << " time    " << fixed << setprecision(3) << elapsed_ms.count() << " ms" << endl;

  return 0;
//////////////////////////////////////////////////////////////////////////////////



  double alpha{0.9};
  double beta{1.5};
  double x{-1};
    // exact value of E_{9/10, 3/2}(-1) from CAS
  double exact = x < 0 ? 0.59595802527072791093339988837073 : 2.4811589951245268838032857244302 ; 

    // print title and column headers
  cout << "E("<<alpha<<","<<beta<<","<<x<<")" << endl;
  cout << "n log10(abs(err))" << endl;

    // print error for increasing orders n=2...15
  for (int n = 2; n < 10; n++) {
    mlfPadeApprox.setOrder(n);  
    double approx = mlfPadeApprox.mittagLefflerE(alpha, beta, x); 
    double relErr = ( approx - exact ) / exact; 
    cout << n << " " << log10(abs(relErr)) << endl; 
  }

  return 0;
}    