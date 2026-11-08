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

#include "BrentSegmentMinimization.hpp"
#include "BrentTriangleMinimization.hpp"
#include "GeometryUtilities.hpp"

namespace UnitTesting
{
  TEST(TestBrentMinimization, TestBrentSegmentMinimization)
  {
    struct TestCase final
    {
        std::string name;
        std::function<double(double)> f;
        double a;
        double b;
        std::vector<double> expected_x;
        double tolerance;
        unsigned int max_iterations;
    };

    const double tol = 1.0e-14;
    const unsigned int max_it = 500;

    const std::vector<TestCase> tests = {
      {"Quadratic",
       [](double x) { return (x - 2.0) * (x - 2.0); },
       0.0, 5.0, { 2.0 }, tol, max_it},

      {"Quartic",
       [](double x) { return std::pow(x - 1.0, 4); },
       -2.0, 4.0, { 1.0 }, tol, max_it},

      {"Absolute",
       [](double x) { return std::abs(x - 1.0); },
       -2.0, 4.0, { 1.0 }, tol, max_it},

      {"Exponential",
       [](double x) { return std::exp(x) - 3.0*x; },
       -2.0, 3.0, { std::log(3.0) }, tol, max_it},

      {"Near zero",
       [](double x) {
         return (x - 1.0e-10) * (x - 1.0e-10);
       },
       -1.0, 1.0, { 1.0e-10 }, tol, max_it},

      {"Near left boundary",
       [](double x) {
         return (x - 1.0e-8) * (x - 1.0e-8);
       },
       0.0, 1.0, { 1.0e-8 }, tol, max_it},

      {"Near right boundary",
       [](double x) {
         return (x - (1.0 - 1.0e-8)) * (x - (1.0 - 1.0e-8));
       },
       0.0, 1.0, { 1.0 - 1.0e-8 }, tol, max_it},

      {"Cosine",
       [](double x) { return -std::cos(x); },
       -1.0, 1.0, { 0.0 }, tol, max_it}
    };

    for (const auto& test : tests)
    {
      const auto min_result = Gedim::Algorithms::Minimization::brent_segment_minimization(
                                test.f, test.a, test.b, test.tolerance, test.max_iterations
                                );

//      std::cout<< std::scientific
//               << test.name
//               << " | x_nm " << min_result.min.at(0)
//               << " | x_ex " << test.expected_x.at(0)
//               << " | f_x_nm " << min_result.f_min
//               << " | f_x_ex " << test.f(test.expected_x.at(0))
//               << " | it " << min_result.iterations
//               << " | err_x " << std::abs(min_result.min.at(0) - test.expected_x.at(0)) / (std::abs(test.expected_x.at(0)) < 1.0 ? 1.0 : std::abs(test.expected_x.at(0)))
//               << " | err_f " << std::abs(min_result.f_min - test.f(test.expected_x.at(0))) / (std::abs(test.f(test.expected_x.at(0))) < 1.0 ? 1.0 : std::abs(test.f(test.expected_x.at(0))))
//               << std::endl;

      ASSERT_TRUE(min_result.converged);
      ASSERT_TRUE(std::abs(min_result.f_min - test.f(test.expected_x.at(0))) < test.tolerance * (std::abs(test.f(test.expected_x.at(0))) < 1.0 ? 1.0 : std::abs(test.f(test.expected_x.at(0)))));
    }
  }

  TEST(TestBrentMinimization, TestBrentTriangleMinimization)
  {
    struct TestCase final
    {
        std::string name;
        std::function<double(double,double)> f;
        std::array<double, 2> a;
        std::array<double, 2> b;
        std::array<double, 2> c;
        std::vector<double> expected_x;
        double tolerance;
        unsigned int max_iterations;
    };

    const double tol = 1.0e-12;
    const unsigned int max_it = 500;

    const std::vector<TestCase> tests = {

      {
        "Quadratic - interior",
        [](double x, double y) {
          return (x - 0.3) * (x - 0.3)
          + (y - 0.4) * (y - 0.4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.4},
        tol, max_it
      },
      {
        "Quadratic - edge BC",
        [](double x, double y) {
          return (x - 0.6) * (x - 0.6)
          + (y - 0.8) * (y - 0.8);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.4, 0.6},
        tol, max_it
      },
      {
        "Quadratic - edge AC",
        [](double x, double y) {
          return (x + 0.5) * (x + 0.5)
          + (y - 0.4) * (y - 0.4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.0, 0.4},
        tol, max_it
      },
      {
        "Quadratic - edge AB",
        [](double x, double y) {
          return (x - 0.4) * (x - 0.4)
          + (y + 0.5) * (y + 0.5);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.4, 0.0},
        tol, max_it
      },
      {
        "Linear - vertex A",
        [](double x, double y) {
          return x + y;
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.0, 0.0},
        tol, max_it
      },
      {
        "Linear - vertex B",
        [](double x, double y) {
          return -x;
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {1.0, 0.0},
        tol, max_it
      },
      {
        "Linear - vertex C",
        [](double x, double y) {
          return -y;
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.0, 1.0},
        tol, max_it
      },
      {
        "Absolute value - interior",
        [](double x, double y) {
          return std::abs(x - 0.3)
          + std::abs(y - 0.4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.4},
        tol, max_it
      },
      {
        "Absolute value - edge BC",
        [](double x, double y) {
          return std::abs(x - 0.4)
          + std::abs(y - 0.6);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.4, 0.6},
        tol, max_it
      },
      {
        "Maximum of absolute values",
        [](double x, double y) {
          return std::max(
          std::abs(x - 0.2),
          std::abs(y - 0.3)
          );
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.2, 0.3},
        tol, max_it
      },
      {
        "Maximum of affine functions",
        [](double x, double y) {
          return std::max({
            x - 0.3,
            y - 0.3,
            0.6 - x - y
          });
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.3},
        tol, max_it
      },
      {
        "Diagonal kink",
        [](double x, double y) {
          return std::abs(x - y)
          + (x + y - 0.6) * (x + y - 0.6);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.3},
        tol, max_it
      },
      {
        "Ill-conditioned quadratic",
        [](double x, double y) {
          return 1000.0 * (x - 0.2) * (x - 0.2)
          + 0.01 * (y - 0.3) * (y - 0.3);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.2, 0.3},
        tol, 10000
      },
      {
        "Quartic - flat minimum",
        [](double x, double y) {
          return std::pow(x - 0.3, 4)
          + std::pow(y - 0.4, 4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.4},
        tol, max_it
      },
      {
        "Weighted L1",
        [](double x, double y) {
          return 10.0 * std::abs(x - 0.2)
          + std::abs(y - 0.4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.2, 0.4},
        tol, max_it
      },
      {
        "Asymmetric piecewise linear",
        [](double x, double y) {
          return std::max(
          2.0 * (x - 0.3),
          0.3 - x
          ) + std::abs(y - 0.2);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.2},
        tol, max_it
      },
      {
        "Clockwise triangle - interior",
        [](double x, double y) {
          return std::abs(x - 0.2)
          + std::abs(y - 0.3);
        },
        {0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0},
        {0.2, 0.3},
        tol, max_it
      },
      {
        "Permuted vertices - interior",
        [](double x, double y) {
          return (x - 0.3) * (x - 0.3)
          + (y - 0.2) * (y - 0.2);
        },
        {1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0},
        {0.3, 0.2},
        tol, max_it
      },
      {
        "Subtriangle - interior",
        [](double x, double y) {
          return (x - 0.3) * (x - 0.3)
          + (y - 0.3) * (y - 0.3);
        },
        {0.1, 0.1}, {0.6, 0.1}, {0.1, 0.6},
        {0.3, 0.3},
        tol, max_it
      },
      {
        "Subtriangle - boundary",
        [](double x, double y) {
          return (x - 0.6) * (x - 0.6)
          + (y - 0.6) * (y - 0.6);
        },
        {0.1, 0.1}, {0.6, 0.1}, {0.1, 0.6},
        {0.35, 0.35},
        tol, max_it
      },
      {
        "Coupled quadratic",
        [](double x, double y) {
          return (x + y - 0.6) * (x + y - 0.6)
          + (x - y - 0.2) * (x - y - 0.2);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.4, 0.2},
        tol, max_it
      },
      {
        "Maximum of quadratics",
        [](double x, double y) {
          return std::max(
          (x - 0.3) * (x - 0.3),
          (y - 0.4) * (y - 0.4)
          );
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.4},
        tol, max_it
      },
      {
        "Exponential - interior",
        [](double x, double y) {
          return std::exp(x - 0.3)
          + std::exp(0.3 - x)
          + std::exp(y - 0.2)
          + std::exp(0.2 - y);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.2},
        tol, max_it
      },
      {
        "Absolute ridge with quadratic regularization",
        [](double x, double y) {
          return std::abs(x + y - 0.7)
          + (x - 0.3) * (x - 0.3)
          + (y - 0.4) * (y - 0.4);
        },
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
        {0.3, 0.4},
        tol, max_it
      }
    };

    for (const auto& test : tests)
    {
      const auto min_result = Gedim::Algorithms::Minimization::brent_triangle_minimization(
                                test.f, test.a, test.b, test.c, test.tolerance, test.max_iterations
                                );

//      std::cout.precision(2);
//      std::cout<< std::scientific
//               << test.name
//               << " | x_nm " << min_result.min.at(0)
//               << " | y_nm " << min_result.min.at(1)
//               << " | x_ex " << test.expected_x.at(0)
//               << " | y_ex " << test.expected_x.at(1)
//               << " | f_x_nm " << min_result.f_min
//               << " | f_x_ex " << test.f(test.expected_x.at(0), test.expected_x.at(1))
//               << " | it " << min_result.iterations
//               << " | err_x " << std::abs(min_result.min.at(0) - test.expected_x.at(0)) / (std::abs(test.expected_x.at(0)) < 1.0 ? 1.0 : std::abs(test.expected_x.at(0)))
//               << " | err_y " << std::abs(min_result.min.at(1) - test.expected_x.at(1)) / (std::abs(test.expected_x.at(1)) < 1.0 ? 1.0 : std::abs(test.expected_x.at(1)))
//               << " | err_f " << std::abs(min_result.f_min - test.f(test.expected_x.at(0), test.expected_x.at(1))) / (std::abs(test.f(test.expected_x.at(0), test.expected_x.at(1))) < 1.0 ? 1.0 : std::abs(test.f(test.expected_x.at(0), test.expected_x.at(1))))
//               << std::endl;

      ASSERT_TRUE(min_result.converged);
      ASSERT_TRUE(std::abs(min_result.f_min - test.f(test.expected_x.at(0), test.expected_x.at(1))) < test.tolerance * 1e1 * (std::abs(test.f(test.expected_x.at(0), test.expected_x.at(1))) < 1.0 ? 1.0 : std::abs(test.f(test.expected_x.at(0), test.expected_x.at(1)))));
    }
  }
} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
