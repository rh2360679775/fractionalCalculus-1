#include <iostream>
#include <iomanip> // for setprecision
#include <chrono>
#include <cmath>
#include <numbers>
#include "c9.h"
#include "templates.hpp"

using namespace std;
using namespace numbers;

inline FloatType ghGamma(FloatType x)
{
  // Library function
		return tgammaf(FloatType(1) + x);
}

#include "trapez.hpp"