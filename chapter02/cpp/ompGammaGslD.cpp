
#include <gsl/gsl_sf_gamma.h>
#include "templates.hpp"


using FloatType = double;         // Use float
FloatType prec = 16;
// id 
const std::string name = "GSLDouble"; 


inline FloatType ghGamma(FloatType x)
{
  // Library function
		return gsl_sf_gamma(1.0 + x);
}

#include "trapez.hpp"

