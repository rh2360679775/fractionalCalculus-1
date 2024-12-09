#include <cmath>
#include "templates.hpp"


using FloatType = double;         // Use float
FloatType prec = 16;

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgamma(1.0 + x);
}

#include "trapez.hpp"

