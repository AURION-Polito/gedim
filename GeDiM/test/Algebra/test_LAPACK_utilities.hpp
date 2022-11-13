#ifndef __TEST_LAPACK_UTILITIES_H
#define __TEST_LAPACK_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "LAPACK_utilities.hpp"
#include "Eigen_Utilities.hpp"

namespace UnitTesting
{
  TEST(TestLAPACK_utilities, TestSVD)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-15;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      const Eigen::MatrixXd A = (Eigen::MatrixXd(3, 3)<< 8.147236863931789e-01, 9.133758561390194e-01, 2.784982188670484e-01, 9.057919370756192e-01, 6.323592462254095e-01, 5.468815192049838e-01, 1.269868162935061e-01, 9.754040499940952e-02, 9.575068354342976e-01).finished();

      Eigen::VectorXd sv = LAPACK_utilities::svd(A);
      Eigen::VectorXd exptected_sv = (Eigen::VectorXd(3)<< 1.816794219434571e+00, 8.389148459479050e-01, 1.815153182691790e-01).finished();

      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(exptected_sv[0], sv[0]));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(exptected_sv[1], sv[1]));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(exptected_sv[2], sv[2]));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.000904076172981e+01, LAPACK_utilities::cond(sv)));
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }
}

#endif // __TEST_LAPACK_UTILITIES_H
