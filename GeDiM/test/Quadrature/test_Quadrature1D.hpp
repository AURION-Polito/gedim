#ifndef __TEST_QUADRATURE1D_H
#define __TEST_QUADRATURE1D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_GaussLobatto1D.hpp"

namespace UnitTesting
{
TEST(TestQuadrature1D, TestQuadrature_GaussLobatto1D)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 43;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points =
                Gedim::Quadrature::Quadrature_GaussLobatto1D::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, weights.sum(), geometryUtilities.Tolerance1D()));

            for (unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsX = (points.row(0));
                Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
                double result = pointsXPow.dot(weights);
                double expectedResult = 1.0 / (ord + 1.0);
                ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
            }
        }
    }
    catch (const std::exception &exception)
    {
        std::cerr << exception.what() << std::endl;
        FAIL();
    }
}

TEST(TestQuadrature1D, TestQuadrature_Gauss1D)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 31;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, weights.sum(), geometryUtilities.Tolerance1D()));

            for (unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsX = (points.row(0));
                Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
                double result = pointsXPow.dot(weights);
                double expectedResult = 1.0 / (ord + 1.0);

                ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
            }
        }
    }
    catch (const std::exception &exception)
    {
        std::cerr << exception.what() << std::endl;
        FAIL();
    }
}
} // namespace UnitTesting

#endif // __TEST_QUADRATURE1D_H
