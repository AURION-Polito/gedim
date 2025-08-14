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

#ifndef __MapParallelogram_H
#define __MapParallelogram_H

#include "Eigen/Eigen"

namespace Gedim
{
class MapParallelogram
{
  public:
    struct MapParallelogramData final
    {
        Eigen::Matrix3d B;
        Eigen::Matrix3d BInv;
        Eigen::Vector3d b;
        double DetB;
        double DetBInv;
    };
    Eigen::MatrixXd ReferenceVertices;

  private:
    static inline Eigen::Matrix3d B(const Eigen::MatrixXd &vertices)
    {
        Eigen::Matrix3d B;
        B.row(0) << vertices(0, 1) - vertices(0, 0), vertices(0, 3) - vertices(0, 0), 0.0;
        B.row(1) << vertices(1, 1) - vertices(1, 0), vertices(1, 3) - vertices(1, 0), 0.0;
        B.row(2) << 0.0, 0.0, 1.0;
        return B;
    }

    /// translation b for linear map x = B * x_r + b from reference triangle [0,1]x[0,1]/2 to triangle with x points
    /// vertices the triangle to map vertices, size 3 x 3
    /// return the resulting value, size 3 x 3
    static inline Eigen::Vector3d b(const Eigen::MatrixXd &vertices)
    {
        return vertices.col(0);
    }

  public:
    MapParallelogram()
    {
        ReferenceVertices.resize(3, 4);
        ReferenceVertices.col(0) << 0.0, 0.0, 0.0;
        ReferenceVertices.col(1) << 1.0, 0.0, 0.0;
        ReferenceVertices.col(2) << 1.0, 1.0, 0.0;
        ReferenceVertices.col(3) << 0.0, 1.0, 0.0;
    }
    ~MapParallelogram()
    {
    }

    MapParallelogramData Compute(const Eigen::MatrixXd &vertices);

    static Eigen::MatrixXd F(const MapParallelogramData &mapData, const Eigen::MatrixXd &x)
    {
        return (mapData.B * x).colwise() + mapData.b;
    }

    static Eigen::MatrixXd FInv(const MapParallelogramData &mapData, const Eigen::MatrixXd &x)
    {
        return mapData.BInv * (x.colwise() - mapData.b);
    }

    /// Compute the jacobian matrix of the transformation F
    /// \param mapData the map data
    /// \param x points in reference triangle, size 3 x numPoints
    /// \return the B matrix for each points, size 3 x (3 * numPoints)
    static Eigen::MatrixXd J(const MapParallelogramData &mapData, const Eigen::MatrixXd &x);

    static Eigen::VectorXd DetJ(const MapParallelogramData &mapData, const Eigen::MatrixXd &x)
    {
        return Eigen::VectorXd::Constant(x.cols(), mapData.DetB);
    }

    static double DetJ(const MapParallelogramData &mapData)
    {
        return mapData.DetB;
    }
};
} // namespace Gedim

#endif // __MappingSquare_H
