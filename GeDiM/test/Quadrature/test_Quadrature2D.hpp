#ifndef __TEST_QUADRATURE2D_H
#define __TEST_QUADRATURE2D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "Quadrature_Gauss2D_Square.hpp"

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
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                 points,
                                                                 weights);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = (points.row(0));
          Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
          double result = pointsXPow.dot(weights);
          double expectedResult = 1.0 / ((ord + 1) * (ord + 2));

          ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
        }
      }
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
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
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                 points,
                                                                 weights);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsY = (points.row(1));
          Eigen::VectorXd pointsYPow = (pointsY.array()).pow(ord);
          double result = pointsYPow.dot(weights);
          double expectedResult = 1.0 / ((ord + 1) * (ord + 2));

          ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
        }
      }
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
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
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                 points,
                                                                 weights);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(0.5, weights.sum(), geometryUtilities.Tolerance2D()));

        for(unsigned int ord = 0; ord < orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = points.row(0).array().pow(ord);
          Eigen::VectorXd pointsY = points.row(1);
          Eigen::VectorXd cwiseProd = pointsX.cwiseProduct(pointsY);
          double result = cwiseProd.dot(weights);
          double expectedResult	= 1.0 / ((ord + 1) * (ord + 2) * (ord + 3));

          ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
        }
      }
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
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
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss2D_Square::FillPointsAndWeights(quadratureOrders[numOrd],
                                                               points,
                                                               weights);

        ASSERT_TRUE(geometryUtilities.AreValuesEqual(1.0, weights.sum(), geometryUtilities.Tolerance2D()));

        for(unsigned int ord = 0; ord < orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = points.row(0).array().pow(ord);
          Eigen::VectorXd pointsY = points.row(1).array().pow(ord);
          Eigen::VectorXd cwiseProd = pointsX.cwiseProduct(pointsY);
          double result = cwiseProd.dot(weights);

          double expectedResult	= 1.0/ ((ord+1)*(ord+1));

          ASSERT_TRUE(geometryUtilities.AreValuesEqual(expectedResult, result, geometryUtilities.Tolerance1D()));
        }
      }
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }
}

#endif // __TEST_QUADRATURE2D_H
