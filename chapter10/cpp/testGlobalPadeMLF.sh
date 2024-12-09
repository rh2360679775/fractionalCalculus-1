#!/bin/bash

# remove old stuff
rm testGlobalPadeMLF
g++ -std=c++20 -Ofast -I/home/r/eigen-3.4.0 -I/home/r/boost_1_85_0 testGlobalPadeMLF.cpp  -o testGlobalPadeMLF
./testGlobalPadeMLF




 