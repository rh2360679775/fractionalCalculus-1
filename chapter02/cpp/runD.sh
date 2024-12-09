
# single cpu test
#cp trapezSingle.hpp trapez.hpp  
# single cpu test
#cp trapezSerial.hpp trapez.hpp  
# 4 threads test
cp trapezOMP4.hpp trapez.hpp  

# include directory
incl="/home/r"

g++ -std=c++20 -Ofast -fopenmp -march=native -I./include gammaLibD.cpp -o gammaLibD
g++ -std=c++20 -Ofast -fopenmp -march=native -I./include gammaSerD.cpp -o gammaSerD
g++ -std=c++20 -Ofast -fopenmp -march=native -I./include gammaLanD.cpp -o gammaLanD
g++ -std=c++20 -Ofast -fopenmp -march=native -I$incl/gsl/include -I./include -c gammaGslD.cpp
g++ -std=c++20 -Ofast -fopenmp -march=native -I$incl/gsl/include -I./include -L $incl/gsl/lib gammaGslD.o -lgsl -lgslcblas -lm -o gammaGslD

steps="8"

echo "LibD"
./gammaLibD $steps
echo "SerD"
./gammaSerD $steps
echo "LanD"
./gammaLanD $steps
echo "GslD"
./gammaGslD $steps
