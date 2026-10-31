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
  TEST(TestBrentMinimization, TestBrentMinimization_Parabola)
  {
    std::string exportFolder = "./Export/TestBrentMinimization/TestBrentMinimization_Parabola";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    auto f = [](double x) {
      return (x - 2.0) * (x - 2.0) + 1.0;
    };

    const double tol = 1.0e-10;
    const double exact_min = 2.0;
    const auto min_result = Gedim::Algorithms::Minimization::brent_minimize(
            f,
            -5.0,
            10.0,
            tol,
            500
        );

    ASSERT_TRUE(min_result.converged);
    ASSERT_TRUE(std::abs(min_result.min - exact_min) < tol * std::abs(exact_min));
    ASSERT_TRUE(std::abs(min_result.f_min - f(exact_min)) < tol * std::abs(f(exact_min)));
  }
} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
