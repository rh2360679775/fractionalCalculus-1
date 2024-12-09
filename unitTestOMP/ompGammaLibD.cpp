#include <cmath>

using FloatType = double;         // Use float

FloatType prec = 33;

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgamma(1.0 + x);
}

#include "trapez.hpp"

