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

#ifndef __TEST_LAPACK_UTILITIES_H
#define __TEST_LAPACK_UTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Eigen_Utilities.hpp"
#include "GeometryUtilities.hpp"
#include "LAPACK_utilities.hpp"

namespace UnitTesting
{
TEST(TestLAPACK_utilities, TestSVD)
{
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-15;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        const Eigen::MatrixXd A = (Eigen::MatrixXd(3, 3) << 8.147236863931789e-01,
                                   9.133758561390194e-01,
                                   2.784982188670484e-01,
                                   9.057919370756192e-01,
                                   6.323592462254095e-01,
                                   5.468815192049838e-01,
                                   1.269868162935061e-01,
                                   9.754040499940952e-02,
                                   9.575068354342976e-01)
                                      .finished();

        Eigen::VectorXd sv = LAPACK_utilities::svd(A);
        Eigen::VectorXd exptected_sv =
            (Eigen::VectorXd(3) << 1.816794219434571e+00, 8.389148459479050e-01, 1.815153182691790e-01).finished();

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(exptected_sv[0], sv[0], geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(exptected_sv[1], sv[1], geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(exptected_sv[2], sv[2], geometryUtilities.Tolerance1D()));
        ASSERT_TRUE(
            geometryUtilities.AreValuesEqual(1.000904076172981e+01, LAPACK_utilities::cond(sv), geometryUtilities.Tolerance1D()));
}

TEST(TestLAPACK_utilities, TestQR_3x2)
{
        const double tolerance = 1.0e-15;

        const Eigen::MatrixXd A = (Eigen::MatrixXd(3, 2) <<
                                   1.0, 2.0, 4.0,
                                   5.0, 7.0, 8.0)
                                      .finished();

        const auto qr_data = LAPACK_utilities::QR(A,
                                                  tolerance);

        ASSERT_TRUE((qr_data.Q * qr_data.Q.transpose() -
                     Eigen::MatrixXd::Identity(qr_data.Q.cols(),
                                                                                   qr_data.Q.rows())).norm() <=
                    tolerance * Eigen::MatrixXd::Identity(qr_data.Q.cols(),
                                                                                                                                                      qr_data.Q.rows()).norm());
        ASSERT_TRUE((qr_data.Q * qr_data.R - A).norm() <= tolerance * A.norm());
        ASSERT_EQ(2,
                  qr_data.Space_Dimension);

}

TEST(TestLAPACK_utilities, TestQR_2x3)
{
        const double tolerance = 1.0e-15;

        const Eigen::MatrixXd A = (Eigen::MatrixXd(3, 2) <<
                                   1.0, 2.0, 4.0,
                                   5.0, 7.0, 8.0)
                                      .finished();

        const auto qr_data = LAPACK_utilities::QR(A.transpose(),
                                                  tolerance);

        ASSERT_TRUE((qr_data.Q * qr_data.Q.transpose() -
                     Eigen::MatrixXd::Identity(qr_data.Q.cols(),
                                                                                   qr_data.Q.rows())).norm() <=
                    tolerance * Eigen::MatrixXd::Identity(qr_data.Q.cols(),
                                                                                                                                                      qr_data.Q.rows()).norm());
        ASSERT_TRUE((qr_data.Q * qr_data.R - A.transpose()).norm() <=
                    tolerance * A.transpose().norm());
        ASSERT_EQ(2,
                  qr_data.Space_Dimension);

}

} // namespace UnitTesting

#endif // __TEST_LAPACK_UTILITIES_H
