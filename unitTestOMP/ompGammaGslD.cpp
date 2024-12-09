#include <cmath>
#include "c20.h"
#include "templates.hpp"   // home of almostEqual


#include <gsl/gsl_sf_gamma.h>

inline double ghGamma(double x)
{
  // Library function
		return gsl_sf_gamma(1.0 + x);
}

#include "trapez.hpp"