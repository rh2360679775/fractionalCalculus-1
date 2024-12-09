#include <cmath>

using FloatType = long double;         // Use float

FloatType prec = 33;

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgammal(1.0L + x);
}

#include "trapez.hpp"