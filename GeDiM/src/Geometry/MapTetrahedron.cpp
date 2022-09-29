#include "MapTetrahedron.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  bool MapTetrahedron::TestMapConfiguration(const Eigen::MatrixXd& vertices,
                                            const Eigen::MatrixXd& referencePoints,
                                            const unsigned int& secondVertexIndex,
                                            const unsigned int& thirdVertexIndex,
                                            const unsigned int& fourthVertexIndex,
                                            MapTetrahedron::MapTetrahedronData& result) const
  {
    result.Q = Q(vertices.col(0),
                 vertices.col(secondVertexIndex),
                 vertices.col(thirdVertexIndex),
                 vertices.col(fourthVertexIndex));
    result.b = b(vertices.col(0));

    return (geometryUtilities.IsValue1DZero((vertices -
                                           F(result,
                                             referencePoints)).norm()));
  }
  // ***************************************************************************
  MapTetrahedron::MapTetrahedronData MapTetrahedron::Compute(const Eigen::MatrixXd& vertices) const
  {
    MapTetrahedronData result;

    MatrixXd referencePoints;
    referencePoints.resize(3, 4);
    referencePoints.col(0)<< 0.0, 0.0, 0.0;
    referencePoints.col(1)<< 1.0, 0.0, 0.0;
    referencePoints.col(2)<< 0.0, 1.0, 0.0;
    referencePoints.col(3)<< 0.0, 0.0, 1.0;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             1, 2, 3,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             1, 3, 2,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             2, 1, 3,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             2, 3, 1,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             3, 1, 2,
                             result))
      return result;

    if (TestMapConfiguration(vertices,
                             referencePoints,
                             3, 2, 1,
                             result))
      return result;

    throw runtime_error("Tetrahedron cannot be mapped");
  }
  // ***************************************************************************
  MatrixXd MapTetrahedron::J(const MapTetrahedronData& mapData,
                             const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = mapData.Q;

    return jacb;
  }
  // ***************************************************************************
  VectorXd MapTetrahedron::DetJ(const MapTetrahedronData& mapData,
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
