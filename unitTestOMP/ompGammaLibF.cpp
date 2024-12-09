#include <cmath>

using FloatType = float;         // Use float

FloatType prec = 33;

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgammaf(1.0 + x);
}

#include "trapez.hpp"