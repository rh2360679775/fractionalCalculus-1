rm ompGammaGslD

g++ -std=c++20 -ffast-math -Ofast -march=native -fopenmp -I/home/r/boost_1.86.0  -I./include -c ompGammaGslD.cpp
g++ -std=c++20 -ffast-math -Ofast -march=native -fopenmp  ompGammaGslD.o -lgsl -lgslcblas -lm -o ompGammaGslD

steps="4"
tasks="4"

./ompGammaGslD $steps $tasks
