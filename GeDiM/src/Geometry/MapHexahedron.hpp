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

#ifndef __MapHexahedron_H
#define __MapHexahedron_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"

namespace Gedim
{
class MapHexahedron
{
  public:
    struct MapHexahedronData final
    {
        Eigen::Matrix3d Q;
        Eigen::Vector3d b;
        Eigen::Matrix3d QInv;
        double DetQ;
        double DetQInv;
    };

  private:
    const GeometryUtilities &geometryUtilities;

    bool TestMapConfiguration(const Eigen::MatrixXd &vertices,
                              const std::vector<unsigned int> &coordinateSystem,
                              const Eigen::MatrixXd &referencePoints,
                              const unsigned int &secondVertexIndex,
                              const unsigned int &thirdVertexIndex,
                              const unsigned int &fourthVertexIndex,
                              MapHexahedron::MapHexahedronData &result) const;

  public:
    MapHexahedron(const GeometryUtilities &geometryUtilities) : geometryUtilities(geometryUtilities)
    {
    }
    ~MapHexahedron()
    {
    }

    /// Matrix Q for linear map x = Q * x_r + b from reference Hexahedron [0,1]x[0,1]x[0,1] to Hexahedron with x points
    /// vertices the Hexahedron to map vertices, size 3 x 4
    /// return the resulting value, size 3 x 3
    static Eigen::Matrix3d Q(const Eigen::Vector3d &firstVertex,
                             const Eigen::Vector3d &secondVertex,
                             const Eigen::Vector3d &thirdVertex,
                             const Eigen::Vector3d &fourthVertex)
    {
        Eigen::Matrix3d Q;
        Q.row(0) << secondVertex.x() - firstVertex.x(), thirdVertex.x() - firstVertex.x(),
            fourthVertex.x() - firstVertex.x();
        Q.row(1) << secondVertex.y() - firstVertex.y(), thirdVertex.y() - firstVertex.y(),
            fourthVertex.y() - firstVertex.y();
        Q.row(2) << secondVertex.z() - firstVertex.z(), thirdVertex.z() - firstVertex.z(),
            fourthVertex.z() - firstVertex.z();

        return Q;
    }

    /// Matrix Q for linear map x = Q * x_r + b from reference Hexahedron [0,1]x[0,1]x[0,1] to Hexahedron with x points
    /// vertices the Hexahedron to map vertices, size 3 x 4
    /// return the resulting value, size 3 x 3
    static inline Eigen::Vector3d b(const Eigen::Vector3d &firstVertex)
    {
        return firstVertex;
    }

    /// Map from the hexahedron reference element [0,1]x[0,1]x[0,1]/2 to the polygon x = F(x_r) = Q * x_r + b
    /// \param vertices the hexahedron to map vertices, size 3 x 4
    /// \param edges the hexahedron edges
    /// \return the map data
    MapHexahedronData Compute(const Eigen::MatrixXd &vertices, const std::vector<unsigned int> &coordinateSystem) const;

    /// Map from the Hexahedron reference element [0,1]x[0,1]x[0,1] to the polygon x = F(x_r) = Q * x_r + b
    /// \param mapData the map data computed
    /// \param x points in reference Hexahedron, size 3 x numPoints
    /// \return the mapped points, size 3 x numPoints
    static inline Eigen::MatrixXd F(const MapHexahedronData &mapData, const Eigen::MatrixXd &x)
    {
        return (mapData.Q * x).colwise() + mapData.b;
    }

    static inline Eigen::MatrixXd FInv(const MapHexahedronData &mapData, const Eigen::MatrixXd &x)
    {
        return mapData.QInv * (x.colwise() - mapData.b);
    }

    /// Compute the jacobian matrix of the transformation F
    /// \param mapData the map data computed
    /// \param x points in reference Hexahedron, size 3 x numPoints
    /// \return the Q matrix for each points, size 2 x (2 * numPoints)
    static Eigen::MatrixXd J(const MapHexahedronData &mapData, const Eigen::MatrixXd &x);

    /// Compute the determinant of the jacobian matrix of the trasformation
    /// \param mapData the map data computed
    /// \param x points in reference Hexahedron, size 3 x numPoints
    /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
    static Eigen::VectorXd DetJ(const MapHexahedronData &mapData, const Eigen::MatrixXd &x);
};
} // namespace Gedim

#endif // __MapHexahedron_H
