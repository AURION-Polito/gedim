// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __TEST_BrentMinimization_H
#define __TEST_BrentMinimization_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "BrentMinimization.hpp"
#include "GeometryUtilities.hpp"

namespace UnitTesting
{
  TEST(TestBrentMinimization, TestBrentMinimization_Functions)
  {
    struct TestCase final
    {
        std::string name;
        std::function<double(double)> f;
        double a;
        double b;
        double expected_x;
        double tolerance;
        unsigned int max_iterations;
    };

    const double tol = 1.0e-14;
    const unsigned int max_it = 500;

    const std::vector<TestCase> tests = {
      {"Quadratic",
       [](double x) { return std::pow(x - 2.0, 2); },
       0.0, 5.0, 2.0, tol, max_it},

      {"Quartic",
       [](double x) { return std::pow(x - 1.0, 4); },
       -2.0, 4.0, 1.0, tol, max_it},

      {"Absolute",
       [](double x) { return std::abs(x - 1.0); },
       -2.0, 4.0, 1.0, tol, max_it},

      {"Exponential",
       [](double x) { return std::exp(x) - 3.0*x; },
       -2.0, 3.0, std::log(3.0), tol, max_it},

      {"Near zero",
       [](double x) {
         return std::pow(x - 1e-10, 2);
       },
       -1.0, 1.0, 1.0e-10, tol, max_it},

      {"Near left boundary",
       [](double x) {
         return std::pow(x - 1e-8, 2);
       },
       0.0, 1.0, 1.0e-8, tol, max_it},

      {"Near right boundary",
       [](double x) {
         return std::pow(x - (1.0 - 1e-8), 2);
       },
       0.0, 1.0, 1.0 - 1.0e-8, tol, max_it},

      {"Cosine",
       [](double x) { return -std::cos(x); },
       -1.0, 1.0, 0.0, tol, max_it}
    };

    for (const auto& test : tests)
    {
      const auto min_result = Gedim::Algorithms::Minimization::brent_minimize(
                                test.f, test.a, test.b, test.tolerance, test.max_iterations
                                );

//      std::cout.precision(2);
//      std::cout<< std::scientific
//               << test.name
//               << " | x_nm " << min_result.min
//               << " | x_ex " << test.expected_x
//               << " | f_x_nm " << min_result.f_min
//               << " | f_x_ex " << counted_f(test.expected_x)
//               << " | it " << min_result.iterations
//               << " | err_x " << std::abs(min_result.min - test.expected_x) / (std::abs(test.expected_x) < 1.0 ? 1.0 : std::abs(test.expected_x))
//               << " | err_ " << std::abs(min_result.f_min - counted_f(test.expected_x)) / (std::abs(counted_f(test.expected_x)) < 1.0 ? 1.0 : std::abs(counted_f(test.expected_x)))
//               << std::endl;

      ASSERT_TRUE(min_result.converged);
      ASSERT_TRUE(std::abs(min_result.f_min - test.f(test.expected_x)) < test.tolerance * (std::abs(test.f(test.expected_x)) < 1.0 ? 1.0 : std::abs(test.f(test.expected_x))));
    }
  }
} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
