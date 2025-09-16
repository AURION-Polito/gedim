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

#include "MapHexahedron.hpp"
#include "CommonUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
// ***************************************************************************
MapHexahedron::MapHexahedronData MapHexahedron::Compute(const Eigen::MatrixXd &vertices, const std::vector<Eigen::MatrixXi> &faces)
{
    MapHexahedron::MapHexahedronData result;

    result.ReferenceVertices.resize(3, 8);
    result.ReferenceVertices.col(0) << 0.0, 0.0, 0.0;
    result.ReferenceVertices.col(1) << 1.0, 0.0, 0.0;
    result.ReferenceVertices.col(2) << 1.0, 1.0, 0.0;
    result.ReferenceVertices.col(3) << 0.0, 1.0, 0.0;
    result.ReferenceVertices.col(4) << 0.0, 0.0, 1.0;
    result.ReferenceVertices.col(5) << 1.0, 0.0, 1.0;
    result.ReferenceVertices.col(6) << 1.0, 1.0, 1.0;
    result.ReferenceVertices.col(7) << 0.0, 1.0, 1.0;

    result.A.resize(8, 8);
    result.A.row(0) << 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0;
    result.A.row(1) << 0.0, 1.0, 0.0, 0.0, -1.0, -1.0, 0.0, 1.0;
    result.A.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0;
    result.A.row(3) << 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 1.0;
    result.A.row(4) << 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, -1.0, 1.0;
    result.A.row(5) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0;
    result.A.row(6) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    result.A.row(7) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0;

    result.Vertices.resize(3, 8);

    for (unsigned int v = 0; v < 4; v++)
        result.VertexOrder.insert({faces[0](0, v), v});

    unsigned int vertex_index = 4;
    for (unsigned int v = 0; v < 8; v++)
    {
        const auto it = result.VertexOrder.find(v);
        if (it == result.VertexOrder.end())
        {
            result.VertexOrder.insert({v, vertex_index});
            vertex_index++;
        }
    }

    for (auto idx : result.VertexOrder)
        result.Vertices.col(idx.second) = vertices.col(idx.first);

    result.Coefficients = result.Vertices * result.A;

    return result;
}
// ***************************************************************************
MatrixXd MapHexahedron::F(const MapHexahedron::MapHexahedronData &mapData, const MatrixXd &referencePoints)
{
    const MatrixXd psi_values = Psi(referencePoints);
    return mapData.Vertices * psi_values;
}
// ***************************************************************************
MatrixXd MapHexahedron::J(const MapHexahedron::MapHexahedronData &mapData, const MatrixXd &referencePoints)
{
    const unsigned int numPoints = referencePoints.cols();
    MatrixXd jacb = Eigen::MatrixXd::Zero(3, 3 * numPoints);
    std::vector<Eigen::MatrixXd> dpsi_values = dPsi(referencePoints);

    MatrixXd pointJacb;
    for (unsigned int p = 0; p < numPoints; p++)
    {
        pointJacb.setZero(3, 8);

        for (unsigned int d = 0; d < 3; d++)
            pointJacb.row(d) << dpsi_values[d].col(p).transpose();

        jacb.block(0, 3 * p, 3, 3) = mapData.Vertices * pointJacb.transpose();
    }

    return jacb;
}
// ***************************************************************************
MatrixXd MapHexahedron::FInv(const MapHexahedron::MapHexahedronData &mapData, const MatrixXd &points)
{
    Gedim::Utilities::Unused(mapData);
    Gedim::Utilities::Unused(points);

    throw std::runtime_error("not implemented method");
}
// ***************************************************************************
MatrixXd MapHexahedron::JInv(const MapHexahedron::MapHexahedronData &mapData, const MatrixXd &referencePoints)
{
    const unsigned int numPoints = referencePoints.cols();
    const MatrixXd jacb = J(mapData, referencePoints);

    MatrixXd jacbInv = Eigen::MatrixXd::Zero(3, 3 * numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
        jacbInv.block(0, 3 * p, 3, 3) = jacb.block(0, 3 * p, 3, 3).inverse();

    return jacbInv;
}
// ***************************************************************************
VectorXd MapHexahedron::DetJ(const MapHexahedron::MapHexahedronData &mapData, const MatrixXd &referencePoints)
{
    MatrixXd jacb = J(mapData, referencePoints);

    const unsigned int numPoints = referencePoints.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
        detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
}
// ***************************************************************************
} // namespace Gedim
