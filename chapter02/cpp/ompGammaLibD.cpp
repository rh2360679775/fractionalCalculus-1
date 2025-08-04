#include <cmath>
#include "templates.hpp"


using FloatType = double;         // Use float
FloatType prec = 16;
// id 
const std::string name = "LibDouble"; 

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgamma(FloatType(1) + x);
}

#include "trapez.hpp"

