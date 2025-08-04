#include <cmath>
#include "templates.hpp"


using FloatType = long double;         // Use float
FloatType prec = 44;
// id 
const std::string name = "LibLongDouble"; 

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgamma(FloatType(1) + x);
}

#include "trapez.hpp"

