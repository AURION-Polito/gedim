#include "MapHexahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  bool MapHexahedron::TestMapConfiguration(const Eigen::MatrixXd& vertices,
                                           const vector<unsigned int>& coordinateSystem,
                                           const Eigen::MatrixXd& referencePoints,
                                           const unsigned int& secondVertexIndex,
                                           const unsigned int& thirdVertexIndex,
                                           const unsigned int& fourthVertexIndex,
                                           MapHexahedron::MapHexahedronData& result) const
  {
    result.Q = Q(vertices.col(coordinateSystem[0]),
        vertices.col(coordinateSystem[secondVertexIndex]),
        vertices.col(coordinateSystem[thirdVertexIndex]),
        vertices.col(coordinateSystem[fourthVertexIndex]));
    result.b = b(vertices.col(coordinateSystem[0]));

    return (geometryUtilities.IsValueZero((vertices -
                                           F(result,
                                             referencePoints)).norm(),
                                          geometryUtilities.Tolerance1D()));
  }
  // ***************************************************************************
  MapHexahedron::MapHexahedronData MapHexahedron::Compute(const Eigen::MatrixXd& vertices,
                                                          const vector<unsigned int>& coordinateSystem) const
  {
    MapHexahedronData result;

    MatrixXd referencePoints;
    referencePoints.resize(3, 8);
    referencePoints.col(0)<< 0.0, 0.0, 0.0;
    referencePoints.col(1)<< 1.0, 0.0, 0.0;
    referencePoints.col(2)<< 1.0, 1.0, 0.0;
    referencePoints.col(3)<< 0.0, 1.0, 0.0;
    referencePoints.col(4)<< 0.0, 0.0, 1.0;
    referencePoints.col(5)<< 1.0, 0.0, 1.0;
    referencePoints.col(6)<< 1.0, 1.0, 1.0;
    referencePoints.col(7)<< 0.0, 1.0, 1.0;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             1, 2, 3,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             1, 3, 2,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             2, 1, 3,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             2, 3, 1,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             3, 1, 2,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             coordinateSystem,
                             referencePoints,
                             3, 2, 1,
                             result))
      return result;

    throw runtime_error("Hexahedron cannot be mapped");
  }
  // ***************************************************************************
  MatrixXd MapHexahedron::J(const MapHexahedronData& mapData,
                            const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);
    Matrix3d q = mapData.Q;

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = q;

    return jacb;
  }
  // ***************************************************************************
  VectorXd MapHexahedron::DetJ(const MapHexahedronData& mapData,
                               const MatrixXd& x) const
  {
    MatrixXd jacb = J(mapData,
                      x);

    const unsigned int numPoints = x.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
      detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
  }
  // ***************************************************************************
}
