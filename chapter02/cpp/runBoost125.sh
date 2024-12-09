# single cpu test
#cp trapezSingle.hpp trapez.hpp  
# serial cpu test
#cp trapezSerial.hpp trapez.hpp  
# 4 threads test
cp trapezOMP4.hpp trapez.hpp  


# clear
rm gammaSerB125.o
rm gammaSerB125
g++ -Ofast -std=c++20  -fopenmp gammaSerB125.cpp -c -I/home/r/boost_1_85_0 -I./include
g++ -Ofast -std=c++20  -fopenmp gammaSerB125.o  -lm -o gammaSerB125

# steps: number of intervals
steps="3"
./gammaSerB125 $steps 2
