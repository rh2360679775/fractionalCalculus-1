#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

const double PI = numbers::pi; // or std::numbers::pi in C++20

class GlobalPadeMLF {

public:
    // constructor
  GlobalPadeMLF(): n(0), n2(0), alpha(-1), beta(-1) {
    setOrder(2); // Initialize n to 2 by default
  }

  void setThresh(double new_tresh){
    alphaThresh = new_tresh;
  }
  
  inline void setOrder(int new_n) {
    // initialize values when n changes
    if (new_n != n) {
      n = new_n;
      n2 = 2 * n;

      matM.resize(n2, n2);
      matM.setZero();

      vecB.resize(n2);
      vecB.setZero();

      vecPsQs.resize(n2);

      lstGammasM.resize(n + 1);
      lstGammasP.resize(n + 1);

      alpha = -1;  // triggers a change
      beta = -1;   // triggers a change
    }
  }
  
  inline double mittagLefflerE( double alpha_new, double beta_new, double x ) {
    double result = 0;
 
    // special cases: 
    if ( alpha_new > 2.0 ) {
      cout<< "case : alpha > 2, not yet implemented , returning NaN " <<endl;
      return numeric_limits<double>::quiet_NaN();   // not yet implemented
    }

    // purpose: evaluate f[alpha, beta, x] = MittagLefflerE[alpha, beta, -x] 
    
    x = -x;  // from now on we calculate with f(...,x) = MittagLefflerE(..., -x)
    
    // from numerical experiments: for significantly smaller errors    
    // use reduction formula if alpha > alphaThresh << 1.0
    // see e.g. eq. (2.1) in:
    // Seybold and Hilfer (2008), 
    // NUMERICAL ALGORITHM FOR CALCULATING THE GENERALIZED MITTAG-LEFFLER FUNCTION
    // SIAM Vol. 47, No.1, pp 69-88 ,  doi. 10.1137/070700280
    
    if (alpha_new > alphaThresh) {              // reduction start
      int m = floor(alpha_new / alphaThresh);
      int mm = 2 * m + 1; 
      double mp = 1.0 / mm;
      complex<double> z;
      complex<double> cResult;

      SetAlphaBetaSolveEqsGetPsQs(alpha_new * mp, beta_new);

      // now we are ready for the reduction formula
      for (int k = 0; k < mm; k++) {
        z = pow(x, mp) * exp(2i * PI * mp * (double) k);
        cResult = getPadeApproximant(z,n,n);
        result += cResult.real();                // all imags cancel as long as x /in real
      }
      result *= mp;
    }                                            // reduction end 
    else {                                       // no reduction necessary
      SetAlphaBetaSolveEqsGetPsQs(alpha_new, beta_new);
      result = getPadeApproximant(x,n,n);
    }

    return result;
  }

  // add on  
  // incomplete gamma for x /in R (not regularized) 
  // extending boost::math::tgamma to negative argument tgamma(a,-x)

  inline complex <double> tgammaX(double a, double x, int n_new = 5) {

    if (x > 0) return boost::math::tgamma( max(a, 0.0), x ); // works only for a>0, x>0
    // x<=0, a> -1
    setOrder(n_new);                                      // setOrder
    a = max(a, -1.0);                                     // set limits for beta = 1+a 
    SetAlphaBetaSolveEqsGetPsQs(1.0, 1.0 + a);            // set alpha, beta and solve corresponding linear system
    double fct = getPadeApproximant(-x, n, n);  // MittagLefflerE[1,1+a,-x]
    complex<double> resC = ( 1.0 - exp(-x)*fct*pow(x + 0i, a))*tgamma(a);
    return resC ;

  }

  private: 
    double alpha;                // current alpha
    double beta;                 // current beta
    double alphaThresh{1.0};     // splitting threshold
    
    int n ;                      // order of Pade Approx pade = P(n)/Q(n)                       
    int n2;                      // total number of coefficients = 2*n
                                  // Linear Solve using  Eigen library 
    Eigen::VectorXd lstGammasP;  // 1/Gamma[beta + i*alpha] elements for M,B
    Eigen::VectorXd lstGammasM;  // 1/Gamma[beta - i*alpha] elements for M,B
    Eigen::MatrixXd matM;        // matM * vecPsQs = vecB
    Eigen::VectorXd vecB;        // inhomogeneity b    
    Eigen::VectorXd vecPsQs;     // solution vector {Pn, Qn} 

  // create list of gammas
  inline void setLstGammas() {
    double value;
    for (int i = 0; i <= n; i++) {
      value = i * alpha - beta; // test for negative integer , in that case
      lstGammasM(i) =
        (almostEqual(floor(value), value)) ? 0 : 1 / tgamma(-value); // -i
      lstGammasP(i) = 1 / tgamma(beta + i * alpha); // +i
    }
  }

  // setup vecB
  inline void setVectorB() {
    for (int i = 0; i < n; i++)
      vecB(n + i) = (((n + i) % 2) ? +1 : -1) * lstGammasM(i + 1);
    vecB(n) *= -1;
  }

  // set up the matrix of size n
  inline void setMatrixM() {
    // upper left quadrant
    for (int i = 0; i < n; i++) matM(i, i) = lstGammasM(1);
    // lower left quadrant
    for (int i = 1; i < n; i++) matM(i + n, n - i) = lstGammasM(1);
    // upper right quadrant
    for (int i = 1; i <= n; i++)
      for (int j = 0; j < n - i + 1; j++)
        matM(i + j, n - 1 + i) = (((j) % 2) ? +1 : -1) * lstGammasP(j);
    // lower right quadrant
    for (int i = 1; i < n; i++)
      for (int j = 0; j < n - i; j++)
        matM(n + i + j, n2 - i) = (((j) % 2) ? +1 : -1) * lstGammasM(j + 1);
  }

  // set alpha/beta dependent quantities and 
  // solve linear equation for the coefficients of Pade approximation
  inline void SetAlphaBetaSolveEqsGetPsQs(double alpha_new, double beta_new = 1) {
    // if alpha or beta changed, there is need for preparation
    if (!almostEqual(alpha_new, alpha) || !almostEqual(beta_new, beta)) {
      alpha = alpha_new;
      beta = beta_new;
      setLstGammas();                            // refresh gammas
      setMatrixM();                              // now populate matrix
      setVectorB();                              // now populate vector
      vecPsQs = matM.fullPivLu().solve(vecB);    // solve the linear system of eqs
      vecPsQs = (((n) % 2) ? +1 : -1) * vecPsQs; // correct sign for given odd/even n
    }
  }

  // Evaluate the Pade Approximant using Horner's scheme
  // template, because it works for float and complex
  template < typename T >
  inline T getPadeApproximant(T x, int n1 , int n2 ) {
    T nominator = 1.0;
    T deNominator = 1.0;
    int i;
    for (i = n1 - 1; i > 0; i--)              // for first n1 -> P(n1)
      nominator = nominator * x + vecPsQs[i];
    for (i =  n1 + n2  - 1; i >= n1 ; i--)    // for next  n2 -> Q(n2)
      deNominator = deNominator * x + vecPsQs[i];

    return ( nominator / deNominator * lstGammasM[1]);  // P/Q
  }

  // a==b for double / float / complex
  template < typename T >
  inline bool almostEqual(T a, T b) {
    return (a == 0.0 || b == 0.0) ?
      (abs(a - b) <= std::numeric_limits<T>::epsilon()) :
      ((abs(a - b) <= std::numeric_limits<T>::epsilon() * abs(a)) || (abs(a - b) <= std::numeric_limits<T>::epsilon() * abs(b)));
  }


}; // fin

/*
!
!    This program is free software; you can redistribute it and/or modify it under the terms of
!    the GNU General Public License version 2 as published by the Free Software Foundation.
!
!    This  program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.
!    See the GNU General Public License for more details.
!
!    Author  : R.Herrmann
!    Email   : r.herrmann@gigahedron.de
!    Version : 0.99
!    Date    : June,5 2024
!
*/

