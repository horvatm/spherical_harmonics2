/*
  Compare the values of spherical harmonics and their derivatives wrt to theta
  and phi generated with C++ with the ones computered using Wolfram mathematica.

  Author: Martin Horvat, December 2020
*/
#include <iostream>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include <fstream>
#include <chrono>

#include "spher_harm.h"
#include "constants.h"

#include "print.h"

int main(int argc, char **argv){

  if (argc != 6) {
    std::cerr
      << "./test_spher_harm <values file> <full report> <short report> <summary> <type>\n"
      << "notes:\n"
      << "  type - double, long_double\n";

    return EXIT_FAILURE;
  }

  char *type = argv[5];

  if (cmp_str("double", type))
    check_values<double>(argv[1], argv[2], argv[3], argv[4]);
  else if (cmp_str("long_double", type))
    check_values<long double>(argv[1], argv[2], argv[3], argv[4]);
  else {
    std::cerr << "This type is not supported.\n";
  }

  return EXIT_FAILURE;
}
