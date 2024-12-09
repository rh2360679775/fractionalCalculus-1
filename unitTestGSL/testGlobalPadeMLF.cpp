#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include "globalPadeMLF.hpp"

using namespace std;

// Function to print usage information
void printUsage() {
    cout << "Usage: unitTestInstallationML [options]" << endl;
    cout << "Options:" << endl;
    cout << "  --help            Display this help message" << endl;
    cout << "  --thresh <value>  Set the threshold value (double, default is 0.3)" << endl;
    cout << endl;
    cout << "Example usage:" << endl;
    cout << "  unitTestInstallationML --thresh 0.5 " << endl;
    cout << "  unitTestInstallationML --help" << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program calculates the global Pade approximation for the Mittag-Leffler function, E(alpha, beta=1, x = -0.5)." << endl;
    cout << "  It uses the summation formula with a user-defined threshold " << endl;
    cout << "  The program then compares the computed values with a set of known values for the Mittag-Leffler function." << endl;
    cout << endl;
}

// Main function
int main(int argc, char* argv[]) {
  double thresh = 0.3;   // Default threshold value
  int order = 2;         // Default order
  double alpha = 1.0;    // Default alpha value
 // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            printUsage();
            return 0;
        } else if (strcmp(argv[i], "--thresh") == 0 && i + 1 < argc) {
            thresh = atof(argv[++i]);
        } else {
            cerr << "Unknown option: " << argv[i] << endl;
            printUsage();
            return 1;
        }
    }
  GlobalPadeMLF mlfPadeApprox;
    // Set the threshold based on the command line argument
  mlfPadeApprox.setThresh(thresh);
  mlfPadeApprox.setOrder(order);

    // alpha values to check E(alpha, beta=1, x = -0.5)
  constexpr double closeTo1[5] = {   // alpha = 1.0 - 10^k, beta = 1, x = -0.5
    0.60340549869586096761551528315510,
    0.60608995263141647835498381023141,
    0.60648529133691131557613292598209,
    0.60652610988754118255503652642507,
    0.60653020460023866289227390850454
  };
    
  cout << "Testing the Mittag-Leffler function for orders 2 .. 25" <<  endl;
 
  for (int i = 2; i < 25; i++) {
    double product = 0.1;
    mlfPadeApprox.setOrder(i);
    cout << "MittagLefflerE(alpha, beta=1, x=-0.5)," << endl;
    cout << "Using default threshold value: " << thresh << endl;
    cout << "order\talpha\terror" << endl;
    for (int j = 0; j < 5; j++) {
        double alpha = 1.0 - product;  // alpha value
        double error = mlfPadeApprox.mittagLefflerE(alpha, 1.0, -0.5) - closeTo1[j]; // error value
        
        // Print values with titles: "order", "alpha", and "error"
        cout << i << "\t" << alpha << "\t" << error << endl;
        
        product *= 0.1;  // Update product for the next iteration
    }
  }

  return 0;
}

/*
!
!    This program is free software; you can redistribute it and/or modify it under the terms of
!    the GNU General Public License version 2 as published by the Free Software Foundation.
!
!    This  program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.
!    See the GNU General Public License for more details.
!
!    Author  : R.Herrmann
!    Email   : r.herrmann@gigahedron.de
!    Version : 0.99
!    Date    : June,5 2024
!
*/

