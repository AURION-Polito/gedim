#ifndef __TEST_QUADRATURE1D_H
#define __TEST_QUADRATURE1D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "Quadrature_GaussLobatto1D.hpp"
#include "Quadrature_Gauss1D.hpp"

namespace UnitTesting
{
  TEST(TestQuadrature1D, TestQuadrature_GaussLobatto1D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-15;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 15;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1,minOrder,maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1,minOrder,maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_GaussLobatto1D::FillPointsAndWeights(quadratureOrders[numOrd],
                                                               points,
                                                               weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = (points.row(0));
          Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
          double result = pointsXPow.dot(weights);
          double expectedResult = 1.0 / (ord + 1.0);

          ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(expectedResult, result));
        }
      }
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }

  TEST(TestQuadrature1D, TestQuadrature_Gauss1D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-14;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 31;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1,minOrder,maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1,minOrder,maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss1D::FillPointsAndWeights(quadratureOrders[numOrd],
                                                        points,
                                                        weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = (points.row(0));
          Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
          double result = pointsXPow.dot(weights);
          double expectedResult = 1.0 / (ord + 1.0);

          ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(expectedResult, result));
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

#endif // __TEST_QUADRATURE1D_H
