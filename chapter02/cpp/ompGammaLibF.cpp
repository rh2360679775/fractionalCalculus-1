#include <cmath>
#include "templates.hpp"


using FloatType = float;         // Use float
FloatType prec = 6;
// id 
const std::string name = "LibFloat"; 

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgamma(FloatType(1) + x);
}

#include "trapez.hpp"

