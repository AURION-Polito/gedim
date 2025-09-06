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

#ifndef __TEST_QUADRATURE2D_H
#define __TEST_QUADRATURE2D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "Quadrature_Gauss2D_Square.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"

namespace UnitTesting
{
TEST(TestQuadrature2D, TestQuadrature_Gauss2D_Triangle_XToN)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
        geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 30;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points =
                Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

            for (unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsX = (points.row(0));
                Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
                double result = pointsXPow.dot(weights);
                double expectedResult = 1.0 / ((ord + 1) * (ord + 2));

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

TEST(TestQuadrature2D, TestQuadrature_Gauss2D_Triangle_YToN)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
        geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 30;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points =
                Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

            for (unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsY = (points.row(1));
                Eigen::VectorXd pointsYPow = (pointsY.array()).pow(ord);
                double result = pointsYPow.dot(weights);
                double expectedResult = 1.0 / ((ord + 1) * (ord + 2));

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

TEST(TestQuadrature2D, TestQuadrature_Gauss2D_Triangle_XToNY)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
        geometryUtilitiesConfig.Tolerance2D = 1.0e-11;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 30;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points =
                Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

            for (unsigned int ord = 0; ord < orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsX = points.row(0).array().pow(ord);
                Eigen::VectorXd pointsY = points.row(1);
                Eigen::VectorXd cwiseProd = pointsX.cwiseProduct(pointsY);
                double result = cwiseProd.dot(weights);
                double expectedResult = 1.0 / ((ord + 1) * (ord + 2) * (ord + 3));

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

TEST(TestQuadrature2D, TestQuadrature_Gauss2D_Square_XToNYToN)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
        geometryUtilitiesConfig.Tolerance2D = 1.0e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        unsigned int minOrder = 0;
        unsigned int maxOrder = 30;

        Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
        Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

        for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
        {
            const auto quadrature_points =
                Gedim::Quadrature::Quadrature_Gauss2D_Square::FillPointsAndWeights(quadratureOrders[numOrd]);

            const Eigen::MatrixXd &points = quadrature_points.Points;
            const Eigen::VectorXd &weights = quadrature_points.Weights;

            ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, weights.sum(), geometryUtilities.Tolerance2D()));

            for (unsigned int ord = 0; ord < orderMax[numOrd]; ord++)
            {
                Eigen::VectorXd pointsX = points.row(0).array().pow(ord);
                Eigen::VectorXd pointsY = points.row(1).array().pow(ord);
                Eigen::VectorXd cwiseProd = pointsX.cwiseProduct(pointsY);
                double result = cwiseProd.dot(weights);

                double expectedResult = 1.0 / ((ord + 1) * (ord + 1));

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

#endif // __TEST_QUADRATURE2D_H
