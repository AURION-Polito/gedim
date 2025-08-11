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

#include "MapQuadrilateral.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Gedim
{
// ***************************************************************************
MatrixXd MapQuadrilateral::F(const MatrixXd &vertices, const MatrixXd &referencePoints) const
{

    const MatrixXd psi_values = Psi(referencePoints);

    return vertices * psi_values;
}
// ***************************************************************************
MatrixXd MapQuadrilateral::J(const MatrixXd &vertices, const MatrixXd &referencePoints) const
{
    const unsigned int numPoints = referencePoints.cols();
    MatrixXd jacb = Eigen::MatrixXd::Zero(3, 3 * numPoints);
    std::vector<Eigen::MatrixXd> dpsi_values = dPsi(referencePoints);

    MatrixXd pointJacb;
    for (unsigned int p = 0; p < numPoints; p++)
    {
        pointJacb.setZero(2, 4);
        pointJacb.row(0) << dpsi_values[0](0, p), dpsi_values[0](1, p), dpsi_values[0](2, p), dpsi_values[0](3, p);
        pointJacb.row(1) << dpsi_values[1](0, p), dpsi_values[1](1, p), dpsi_values[1](2, p), dpsi_values[1](3, p);

        jacb.block(0, 3 * p, 2, 2) = vertices.block(0, 0, 2, 4) * pointJacb.transpose();
        jacb.block(2, 3 * p, 1, 3) << 0.0, 0.0, 1.0;
    }

    return jacb;
}
// ***************************************************************************
MatrixXd MapQuadrilateral::FInv(const MatrixXd &vertices, const MatrixXd &points) const
{
    throw std::runtime_error("not implemented method");
}
// ***************************************************************************
MatrixXd MapQuadrilateral::JInv(const MatrixXd &vertices, const MatrixXd &referencePoints) const
{
    const unsigned int numPoints = referencePoints.cols();
    const MatrixXd jacb = J(vertices, referencePoints);

    MatrixXd jacbInv = Eigen::MatrixXd::Zero(3, 3 * numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
    {
        jacbInv.block(0, 3 * p, 2, 2) = jacb.block(0, 3 * p, 2, 2).inverse();
        jacbInv.block(2, 3 * p, 1, 3) << 0.0, 0.0, 1.0;
    }

    return jacbInv;
}
// ***************************************************************************
VectorXd MapQuadrilateral::DetJ(const MatrixXd &vertices, const MatrixXd &referencePoints) const
{
    MatrixXd jacb = J(vertices, referencePoints);

    const unsigned int numPoints = referencePoints.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
        detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
}
// ***************************************************************************
} // namespace Gedim
