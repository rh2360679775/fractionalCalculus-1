// c++ 20
// 
#include <concepts>
#include <boost/multiprecision/cpp_bin_float.hpp>

// Define the Boost multiprecision types
namespace mp = boost::multiprecision;

template <unsigned Digits>
using HighPrecision = mp::number<mp::cpp_bin_float<Digits>>;

// Define the concept
template<typename T>
concept FloatTypesAllowed = 
	std::is_same_v<T, float> || 
	std::is_same_v<T, double> || 
	std::is_same_v<T, long double> ||
// boost multiprecision types
	std::is_same_v<T, HighPrecision<99>> ||
	std::is_same_v<T, HighPrecision<778>> ||
	std::is_same_v<T, HighPrecision<1008>>;


