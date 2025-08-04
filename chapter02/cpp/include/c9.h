// c++20 version
#include <iostream>
#include <iomanip> // for setprecision
#include <array>
#include <string>

using FloatType = float;
// precision achieved e-16
const int prec = 6;
// order of the polynome = 9
const std::string name = "seriesF"; 

const std::array<FloatType, 9> as = {
   1.0E0, 
	-5.7719156460972292235905772951355E-1, 
	 9.8820481604618728448450314968256E-1, 
	-8.9705315901264346229071766875510E-1, 
	 9.1820617820812772320771798684085E-1, 
	-7.5672630257745595578426061862168E-1, 
	 4.8224860961030939451976064440311E-1, 
	-1.9357009831920499786845025692432E-1, 
	 3.58815206544029360905044928881E-2};



