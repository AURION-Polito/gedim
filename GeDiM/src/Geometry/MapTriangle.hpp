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

#ifndef __MapTriangle_H
#define __MapTriangle_H

#include "Eigen/Eigen"

namespace Gedim
{
class MapTriangle
{
  public:
    struct MapTriangleData final
    {
        Eigen::Matrix3d B;
        Eigen::Matrix3d BInv;
        Eigen::Vector3d b;
        double DetB;
        double DetBInv;
    };

  private:
    /// Matrix B for linear map x = B * x_r + b from reference triangle [0,1]x[0,1]/2 to triangle with x points
    /// vertices the triangle to map vertices, size 3 x 3
    /// return the resulting value, size 3 x 3
    static inline Eigen::Matrix3d B(const Eigen::Matrix3d &vertices)
    {
        Eigen::Matrix3d B;
        B.row(0) << vertices(0, 1) - vertices(0, 0), vertices(0, 2) - vertices(0, 0), 0.0;
        B.row(1) << vertices(1, 1) - vertices(1, 0), vertices(1, 2) - vertices(1, 0), 0.0;
        B.row(2) << 0.0, 0.0, 1.0;
        return B;
    }

    /// translation b for linear map x = B * x_r + b from reference triangle [0,1]x[0,1]/2 to triangle with x points
    /// vertices the triangle to map vertices, size 3 x 3
    /// return the resulting value, size 3 x 3
    static inline Eigen::Vector3d b(const Eigen::Matrix3d &vertices)
    {
        return vertices.col(0);
    }

  public:
    MapTriangle()
    {
    }
    ~MapTriangle()
    {
    }

    /// Map from the triangle reference element [0,1]x[0,1]/2 to the polygon x = F(x_r) = B * x_r + b
    /// \param vertices the triangle 2D to map vertices, size 3 x 3
    /// \return the map data
    MapTriangleData Compute(const Eigen::Matrix3d &vertices) const;

    /// Map from the triangle reference element [0,1]x[0,1] to the polygon x = F(x_r) = B * x_r + b
    /// \param mapData the map data
    /// \param x points in reference triangle, size 3 x numPoints
    /// \return the mapped polygon points, size 3 x numPoints
    static inline Eigen::MatrixXd F(const MapTriangleData &mapData, const Eigen::MatrixXd &x)
    {
        return (mapData.B * x).colwise() + mapData.b;
    }
    /// Map from the polygon x to the triangle reference element [0,1]x[0,1] x_r = F^-1(x_r) = B^-1 * (x - b)
    /// \param mapData the map data
    /// \param x points in polygon, size 3 x numPoints
    /// \return the mapped reference points, size 3 x numPoints
    static inline Eigen::MatrixXd FInv(const MapTriangleData &mapData, const Eigen::MatrixXd &x)
    {
        return mapData.BInv * (x.colwise() - mapData.b);
    }
    /// Compute the jacobian matrix of the transformation F
    /// \param mapData the map data
    /// \param x points in reference triangle, size 3 x numPoints
    /// \return the B matrix for each points, size 3 x (3 * numPoints)
    static Eigen::MatrixXd J(const MapTriangleData &mapData, const Eigen::MatrixXd &x);
    /// Compute the determinant of the jacobian matrix of the trasformation
    /// \param mapData the map data
    /// \param x points in reference triangle, size 3 x numPoints
    /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
    static inline Eigen::VectorXd DetJ(const MapTriangleData &mapData, const Eigen::MatrixXd &x)
    {
        return Eigen::VectorXd::Constant(x.cols(), mapData.DetB);
    }

    /// Compute the determinant of the jacobian matrix of the trasformation
    /// \param mapData the map data
    /// \param x points in reference triangle, size 3 x numPoints
    /// \return the determinant of Jacobian matrix
    static inline double DetJ(const MapTriangleData &mapData)
    {
        return mapData.DetB;
    };
};
} // namespace Gedim

#endif // __MapTriangle_H
