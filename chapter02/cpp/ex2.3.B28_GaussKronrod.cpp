#include <atomic>  // For counter (thread-safe)
#include <numbers> // C++20: pi
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

atomic<long int> counter{0};  // thread save counter 

inline auto K = [](double z, double a, double b,
                   atomic<long int> &counter) {
   // create the integrand as a lambda function
 return [ = , &counter](double t) {
  double x = t^2; // coordinate transform x = t^2

  double result = pow(x,(1-b)/a)*exp(-pow(x, 1/a))*
   (x*sin(numbers::pi*(1-b))-z*sin(numbers::pi*(1-b+a)))/
   (x*x - 2*x*z*cos(numbers::pi*a) + z*z);

  counter++;         // count function calls
  return 2*t*result; // volume element dx->dt=2*t*dt
 };
};

inline double MittagLefflerZAB(double z,double a,double b)
{  
   // ToDo: specials, only one example implemented
 if( almostEqual(z, 0.0) ){ 
   return 1.0 / boost::math::tgamma(b);
 }  
    // prepare input for quadrature
 auto integrand = K(z,a,b,counter);
    // define the interval for the transformed variable t
 double lowerLimit = 0.0; 
 double upperLimit = numeric_limits<double>::infinity();
    // prepare precision attributes 
 int numNodes = 15;             // number of nodes
 double errorTolerance = 1.e-6; // error tolerance 
 int adaptLevel = 10;           // iterations 
    // central function call
 double Q = boost::math::quadrature::
    gauss_kronrod<double, numNodes>::integrate(integrand,
    lowerLimit, upperLimit, adaptLevel, errorTolerance);
  
 return Q / (a * numbers::pi);
}
