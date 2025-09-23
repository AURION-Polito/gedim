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

#include "MapParallelepiped.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
// ***************************************************************************
MapParallelepiped::MapParallelepiped(const Gedim::GeometryUtilities &geometryUtilities) : geometryUtilities(geometryUtilities)
{
    ReferenceVertices.resize(3, 8);
    ReferenceVertices.col(0) << 0.0, 0.0, 0.0;
    ReferenceVertices.col(1) << 1.0, 0.0, 0.0;
    ReferenceVertices.col(2) << 1.0, 1.0, 0.0;
    ReferenceVertices.col(3) << 0.0, 1.0, 0.0;
    ReferenceVertices.col(4) << 0.0, 0.0, 1.0;
    ReferenceVertices.col(5) << 1.0, 0.0, 1.0;
    ReferenceVertices.col(6) << 1.0, 1.0, 1.0;
    ReferenceVertices.col(7) << 0.0, 1.0, 1.0;
}
// ***************************************************************************
bool MapParallelepiped::TestMapConfiguration(const Eigen::MatrixXd &vertices,
                                             const vector<unsigned int> &coordinateSystem,
                                             const Eigen::MatrixXd &referencePoints,
                                             const unsigned int &secondVertexIndex,
                                             const unsigned int &thirdVertexIndex,
                                             const unsigned int &fourthVertexIndex,
                                             MapParallelepiped::MapParallelepipedData &result) const
{
    result.Q = Q(vertices.col(coordinateSystem[0]),
                 vertices.col(coordinateSystem[secondVertexIndex]),
                 vertices.col(coordinateSystem[thirdVertexIndex]),
                 vertices.col(coordinateSystem[fourthVertexIndex]));
    result.b = b(vertices.col(coordinateSystem[0]));

    if (!geometryUtilities.IsValueZero((vertices - F(result, referencePoints)).norm(), geometryUtilities.Tolerance1D()))
        return false;

    result.QInv = result.Q.inverse();
    result.DetQ = result.Q.determinant();
    result.DetQInv = result.QInv.determinant();

    return true;
}
// ***************************************************************************
MapParallelepiped::MapParallelepipedData MapParallelepiped::Compute(const Eigen::MatrixXd &vertices,
                                                                    const vector<unsigned int> &coordinateSystem) const
{
    MapParallelepipedData result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 1, 2, 3, result))
        return result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 1, 3, 2, result))
        return result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 2, 1, 3, result))
        return result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 2, 3, 1, result))
        return result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 3, 1, 2, result))
        return result;

    if (TestMapConfiguration(vertices, coordinateSystem, ReferenceVertices, 3, 2, 1, result))
        return result;

    throw runtime_error("Hexahedron cannot be mapped");
}
// ***************************************************************************
MatrixXd MapParallelepiped::J(const MapParallelepipedData &mapData, const MatrixXd &x)
{
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);
    Matrix3d q = mapData.Q;

    for (unsigned int p = 0; p < numPoints; p++)
        jacb.block(0, 3 * p, 3, 3) = q;

    return jacb;
}
// ***************************************************************************
VectorXd MapParallelepiped::DetJ(const MapParallelepipedData &mapData, const MatrixXd &x)
{
    MatrixXd jacb = J(mapData, x);

    const unsigned int numPoints = x.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
        detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
}
// ***************************************************************************
} // namespace Gedim
