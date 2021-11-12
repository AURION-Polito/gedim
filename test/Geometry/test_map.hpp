#ifndef __TEST_QUADRATUREMAP_H
#define __TEST_QUADRATUREMAP_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MapQuadrilateral.hpp"
#include "MapTriangle.hpp"
#include "Eigen"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestQuadratureMap, TestMappingQuadrilateral)
  {
    Eigen::MatrixXd referencePoints;
    referencePoints.setZero(3, 9);
    referencePoints.row(0)<< 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01, 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01, 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01;
    referencePoints.row(1)<< 1.127016653792583e-01, 1.127016653792583e-01, 1.127016653792583e-01, 5.000000000000000e-01, 5.000000000000000e-01, 5.000000000000000e-01, 8.872983346207417e-01, 8.872983346207417e-01, 8.872983346207417e-01;

    Eigen::VectorXd referenceWeights(9);
    referenceWeights<< 7.716049382716050e-02, 1.234567901234568e-01, 7.716049382716050e-02, 1.234567901234568e-01, 1.975308641975309e-01, 1.234567901234568e-01, 7.716049382716050e-02, 1.234567901234568e-01, 7.716049382716050e-02;


    {
      Eigen::MatrixXd vertices;
      vertices.setZero(3, 4);
      vertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
      vertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

      Gedim::MapQuadrilateral mapping;

      Eigen::MatrixXd mappedPoints = mapping.F(vertices,
                                               referencePoints);

      Eigen::MatrixXd expectedPoints;
      expectedPoints.setZero(3, 9);
      expectedPoints.row(0)<< 1.551915495751552e+00, 3.293649167310371e+00, 5.035382838869189e+00, 1.628266328441182e+00, 3.100000000000000e+00, 4.571733671558819e+00, 1.704617161130811e+00, 2.906350832689629e+00, 4.108084504248447e+00;
      expectedPoints.row(1)<< 2.421754163448146e+00, 1.223346994592885e+00, 2.493982573762388e-02, 3.508407168855261e+00, 2.850000000000000e+00, 2.191592831144739e+00, 4.595060174262376e+00, 4.476653005407115e+00, 4.358245836551855e+00;

      ASSERT_TRUE((expectedPoints - mappedPoints).norm() < 1e-14);

      Eigen::VectorXd mappedWeights = referenceWeights.array() * mapping.DetJ(vertices,
                                                                              referencePoints).array().abs();

      Eigen::VectorXd expectedWeights(9);
      expectedWeights<< 1.020658186245617e+00, 2.140844247829071e+00, 1.655397123540722e+00, 1.357640948929348e+00, 2.984691358024691e+00, 2.373223248601515e+00, 6.763929999160688e-01, 1.590019949701793e+00, 1.311131937211173e+00;

      ASSERT_TRUE((expectedWeights - mappedWeights).norm() < 1e-14);
      ASSERT_DOUBLE_EQ(mappedWeights.sum(), 1.511000000000000e+01);
    }

    {
      Eigen::MatrixXd vertices;
      vertices.setZero(3, 4);
      vertices.row(0) << -7.7312156205493587e-01,  6.3425787364898834e-01,  7.7312156205493587e-01, -6.3425787364898834e-01;
      vertices.row(1) <<  6.3425787364898834e-01, -7.7312156205493587e-01, -6.3425787364898834e-01,  7.7312156205493587e-01;

      Gedim::MapQuadrilateral mapping;
      Eigen::MatrixXd mappedPoints = mapping.F(vertices,
                                               referencePoints);

      Eigen::MatrixXd expectedPoints;
      expectedPoints.setZero(3, 9);
      expectedPoints.row(0)<< -5.988573868865261e-01, -5.378167525891708e-02, 4.912940363686920e-01, -5.450757116276090e-01, 0, 5.450757116276090e-01, -4.912940363686920e-01, 5.378167525891706e-02, 5.988573868865261e-01;
      expectedPoints.row(1)<< 4.912940363686920e-01, -5.378167525891708e-02, -5.988573868865261e-01, 5.450757116276090e-01, 0, -5.450757116276090e-01, 5.988573868865261e-01, 5.378167525891706e-02, -4.912940363686920e-01;

      ASSERT_TRUE((expectedPoints - mappedPoints).norm() < 1e-14);

      Eigen::VectorXd mappedWeights = referenceWeights.array() * mapping.DetJ(vertices,
                                                                              referencePoints).array().abs();

      Eigen::VectorXd expectedWeights(9);
      expectedWeights<< 3.015955238094568e-02, 4.825528380951308e-02, 3.015955238094568e-02, 4.825528380951308e-02, 7.720845409522094e-02, 4.825528380951308e-02, 3.015955238094568e-02, 4.825528380951308e-02, 3.015955238094568e-02;

      ASSERT_TRUE((expectedWeights - mappedWeights).norm() < 1e-14);
      ASSERT_DOUBLE_EQ(mappedWeights.sum(), 3.908677988570559e-01);
    }

    {
      Eigen::MatrixXd vertices;
      vertices.setZero(3, 4);
      vertices.row(0) << -1.0630237354616701e+00,  1.1148026974627239e+00,  1.1494634009478966e+00, -2.2942118276715412e-01;
      vertices.row(1) <<  1.5565368774740850e-03, -1.1354686296150518e+00, -9.6872630175294960e-01,  8.5220794112921472e-01;

      Gedim::MapQuadrilateral mapping;
      Eigen::MatrixXd mappedPoints = mapping.F(vertices,
                                               referencePoints);

      Eigen::MatrixXd expectedPoints;
      expectedPoints.setZero(3, 9);
      expectedPoints.row(0)<< -7.337785656371913e-01, 7.481683848006121e-02, 8.834122425973138e-01, -4.457988316904019e-01, 2.429552950454490e-01, 9.317094217813001e-01, -1.578190977436126e-01, 4.110937516108370e-01, 9.800066009652866e-01;
      expectedPoints.row(1)<< -3.940504728612496e-02, -5.096250623909643e-01, -9.798450774958036e-01, 2.601987632229570e-01, -3.126076133403282e-01, -8.854139899036132e-01, 5.598025737320390e-01, -1.155901642896920e-01, -7.909829023114229e-01;

      ASSERT_TRUE((expectedPoints - mappedPoints).norm() < 1e-14);

      Eigen::VectorXd mappedWeights = referenceWeights.array() * mapping.DetJ(vertices,
                                                                              referencePoints).array().abs();

      Eigen::VectorXd expectedWeights(9);
      expectedWeights<< 1.942757711316908e-01, 1.961888653655394e-01, 5.096031057523357e-02, 3.056049290531979e-01, 3.055240969728514e-01, 7.630019216286643e-02, 1.877303901848067e-01, 1.857162558505248e-01, 4.441492962834947e-02;

      ASSERT_TRUE((expectedWeights - mappedWeights).norm() < 1e-14);
      ASSERT_DOUBLE_EQ(mappedWeights.sum(), 1.546715740925060e+00);
    }
  }

  TEST(TestQuadratureMap, TestMapTriangle)
  {
    Eigen::MatrixXd vertices;
    vertices.setZero(3, 3);
    vertices.row(0) << -1.0, +5.0, +4.0;
    vertices.row(1) << -2.0, -1.0, +5.0;

    Gedim::MapTriangle mapping;

    Eigen::MatrixXd points;
    points.setZero(3, 7);
    points.row(0)<< 3.333333333333330e-01, 4.701420641051151e-01, 4.701420641051151e-01, 5.971587178976985e-02, 1.012865073234563e-01, 1.012865073234563e-01, 7.974269853530873e-01;
    points.row(1)<< 3.333333333333330e-01, 4.701420641051151e-01, 5.971587178976985e-02, 4.701420641051151e-01, 1.012865073234563e-01, 7.974269853530873e-01, 1.012865073234563e-01;

    Eigen::MatrixXd mappedPoints = mapping.F(vertices,
                                             points);

    Eigen::MatrixXd expectedPoints;
    expectedPoints.setZero(3, 7);
    expectedPoints.row(0)<< 2.666666666666663e+00, 4.171562705156266e+00, 2.119431743579540e+00, 1.709005551264195e+00, 1.141515805580195e-01, 3.594853970706175e+00, 4.290994448735805e+00;
    expectedPoints.row(1)<< 6.666666666666639e-01, 1.761136512840920e+00, -1.111846833366496e+00, 1.350710320525575e+00, -1.189707941412349e+00, 3.683275404795067e+00, -4.935674633827185e-01;

    ASSERT_TRUE((expectedPoints - mappedPoints).norm() < 1e-14);

    Eigen::VectorXd weights(7);
    weights<< 1.125000000000000e-01, 6.619707639425310e-02, 6.619707639425310e-02, 6.619707639425310e-02, 6.296959027241358e-02, 6.296959027241358e-02, 6.296959027241358e-02;

    Eigen::VectorXd mappedWeights = weights.array() * mapping.DetJ(vertices,
                                                                   points).array().abs();

    Eigen::VectorXd expectedWeights(7);
    expectedWeights<< 4.162500000000001e+00, 2.449291826587364e+00, 2.449291826587364e+00, 2.449291826587364e+00, 2.329874840079303e+00, 2.329874840079303e+00, 2.329874840079303e+00;

    ASSERT_TRUE((expectedWeights - mappedWeights).norm() < 1e-14);
    ASSERT_TRUE(abs((mappedWeights).sum() - 1.850000000000000e+01) < 1e-14);
  }
}

#endif // __TEST_QUADRATUREMAP_H
