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
    Eigen::MatrixXd vertices;
    vertices.setZero(3, 4);
    vertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
    vertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

    Gedim::MapQuadrilateral mapping;

    Eigen::MatrixXd points;
    points.setZero(3, 9);
    points.row(0)<< 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01, 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01, 1.127016653792583e-01, 5.000000000000000e-01, 8.872983346207417e-01;
    points.row(1)<< 1.127016653792583e-01, 1.127016653792583e-01, 1.127016653792583e-01, 5.000000000000000e-01, 5.000000000000000e-01, 5.000000000000000e-01, 8.872983346207417e-01, 8.872983346207417e-01, 8.872983346207417e-01;

    Eigen::MatrixXd mappedPoints = mapping.F(vertices,
                                             points);

    Eigen::MatrixXd expectedPoints;
    expectedPoints.setZero(3, 9);
    expectedPoints.row(0)<< 1.551915495751552e+00, 3.293649167310371e+00, 5.035382838869189e+00, 1.628266328441182e+00, 3.100000000000000e+00, 4.571733671558819e+00, 1.704617161130811e+00, 2.906350832689629e+00, 4.108084504248447e+00;
    expectedPoints.row(1)<< 2.421754163448146e+00, 1.223346994592885e+00, 2.493982573762388e-02, 3.508407168855261e+00, 2.850000000000000e+00, 2.191592831144739e+00, 4.595060174262376e+00, 4.476653005407115e+00, 4.358245836551855e+00;

    ASSERT_TRUE((expectedPoints - mappedPoints).norm() < 1e-14);

    Eigen::VectorXd weights(9);
    weights<< 7.716049382716050e-02, 1.234567901234568e-01, 7.716049382716050e-02, 1.234567901234568e-01, 1.975308641975309e-01, 1.234567901234568e-01, 7.716049382716050e-02, 1.234567901234568e-01, 7.716049382716050e-02;

    Eigen::VectorXd mappedWeights = weights.array() * mapping.DetJ(vertices,
                                                                   points).array().abs();

    Eigen::VectorXd expectedWeights(9);
    expectedWeights<< 1.020658186245617e+00, 2.140844247829071e+00, 1.655397123540722e+00, 1.357640948929348e+00, 2.984691358024691e+00, 2.373223248601515e+00, 6.763929999160688e-01, 1.590019949701793e+00, 1.311131937211173e+00;

    ASSERT_TRUE((expectedWeights - mappedWeights).norm() < 1e-14);
    ASSERT_TRUE(abs((mappedWeights).sum() - 1.511000000000000e+01) < 1e-14);
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
