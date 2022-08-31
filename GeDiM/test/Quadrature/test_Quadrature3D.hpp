#ifndef __TEST_QUADRATURE3D_H
#define __TEST_QUADRATURE3D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "Quadrature_Gauss3D_Hexahedron.hpp"

namespace UnitTesting
{
  TEST(TestQuadrature3D, TestQuadrature_Gauss3D_Tetrahedron_XToN)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 15;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                    points,
                                                                    weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0 / 6.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = (points.row(0));
          Eigen::VectorXd pointsXPow = (pointsX.array()).pow(ord);
          double result = pointsXPow.dot(weights);

          double expectedResult =  1.0 / ((ord+1) * (ord+2) * (ord+3));

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

  TEST(TestQuadrature3D, TestQuadrature_Gauss3D_Tetrahedron_YToN)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 15;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                    points,
                                                                    weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0 / 6.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsY = (points.row(1));
          Eigen::VectorXd pointsYPow = (pointsY.array()).pow(ord);
          double result = pointsYPow.dot(weights);

          double expectedResult =  1.0/((ord+1) * (ord+2) * (ord+3));

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

  TEST(TestQuadrature3D, TestQuadrature_Gauss3D_Tetrahedron_ZToN)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 15;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                    points,
                                                                    weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0 / 6.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsZ = (points.row(2));
          Eigen::VectorXd pointsZPow = (pointsZ.array()).pow(ord);
          double result = pointsZPow.dot(weights);

          double expectedResult =  1.0/((ord+1) * (ord+2) * (ord+3));

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

  TEST(TestQuadrature3D, TestQuadrature_Gauss3D_Hexahedron_XToN_YToN)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      unsigned int minOrder = 0;
      unsigned int maxOrder = 15;

      Eigen::VectorXi quadratureOrders = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);
      Eigen::VectorXi orderMax = Eigen::VectorXi::LinSpaced(maxOrder + 1, minOrder, maxOrder);

      for (unsigned int numOrd = 0; numOrd < orderMax.size(); numOrd++)
      {
        Eigen::MatrixXd points;
        Eigen::VectorXd weights;

        Gedim::Quadrature_Gauss3D_Hexahedron::FillPointsAndWeights(quadratureOrders[numOrd],
                                                                   points,
                                                                   weights);

        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0, weights.sum()));

        for(unsigned int ord = 0; ord <= orderMax[numOrd]; ord++)
        {
          Eigen::VectorXd pointsX = points.row(0).array().pow(ord);
          Eigen::VectorXd pointsY = points.row(1).array().pow(ord);
          Eigen::VectorXd pointsZ = points.row(2).array().pow(ord);
          Eigen::VectorXd cwiseProd = pointsX.cwiseProduct(pointsY);
          cwiseProd = cwiseProd.cwiseProduct(pointsZ);
          double result = cwiseProd.dot(weights);

          double expectedResult	= 1.0 / ((ord+1) * (ord+1) * (ord+1));

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

#endif // __TEST_QUADRATURE3D_H
