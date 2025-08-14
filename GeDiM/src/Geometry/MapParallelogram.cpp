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

#include "MapParallelogram.hpp"

namespace Gedim
{
// ***************************************************************************
Eigen::MatrixXd MapParallelogram::J(const MapParallelogramData &mapData, const Eigen::MatrixXd &x)
{
    const unsigned int numPoints = x.cols();
    Eigen::MatrixXd jacb(3, 3 * numPoints);

    for (unsigned int p = 0; p < numPoints; p++)
        jacb.block(0, 3 * p, 3, 3) = mapData.B;

    return jacb;
}
// ***************************************************************************
MapParallelogram::MapParallelogramData MapParallelogram::Compute(const Eigen::MatrixXd &vertices)
{
    MapParallelogram::MapParallelogramData result;

    result.B = B(vertices);
    result.BInv = result.B.inverse();
    result.b = b(vertices);
    result.DetB = result.B.determinant();
    result.DetBInv = result.BInv.determinant();

    return result;
}
// ***************************************************************************
} // namespace Gedim
