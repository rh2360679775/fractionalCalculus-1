#!/bin/bash

# remove old stuff
rm testGlobalPadeMLFfirst
g++ -std=c++20 -Ofast -I/home/r/eigen-3.4.0 -I/home/r/boost_1_85_0 testGlobalPadeMLFfirst.cpp  -o testGlobalPadeMLFfirst
./testGlobalPadeMLFfirst




 