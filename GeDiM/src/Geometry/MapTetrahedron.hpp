#ifndef __MapTetrahedron_H
#define __MapTetrahedron_H

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"

namespace Gedim
{
  class MapTetrahedron
  {
    public:
      struct MapTetrahedronData final
      {
          Eigen::Matrix3d Q;
          Eigen::Vector3d b;
      };

    private:
      const GeometryUtilities& geometryUtilities;

      bool TestMapConfiguration(const Eigen::MatrixXd& vertices,
                                const Eigen::MatrixXd& referencePoints,
                                const unsigned int& secondVertexIndex,
                                const unsigned int& thirdVertexIndex,
                                const unsigned int& fourthVertexIndex,
                                MapTetrahedron::MapTetrahedronData& result) const;

    public:
      MapTetrahedron(const GeometryUtilities& geometryUtilities) :
        geometryUtilities(geometryUtilities) {}
      ~MapTetrahedron() {}

      /// Matrix Q for linear map x = Q * x_r + b from reference tetrahedron [0,1]x[0,1]x[0,1]/2 to tetrahedron with x points
      /// vertices the tetrahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Matrix3d Q(const Eigen::Vector3d& firstVertex,
                               const Eigen::Vector3d& secondVertex,
                               const Eigen::Vector3d& thirdVertex,
                               const Eigen::Vector3d& fourthVertex) const
      {
        Eigen::Matrix3d Q;
        Q.row(0)<< secondVertex.x() - firstVertex.x(), thirdVertex.x() - firstVertex.x(), fourthVertex.x() - firstVertex.x();
        Q.row(1)<< secondVertex.y() - firstVertex.y(), thirdVertex.y() - firstVertex.y(), fourthVertex.y() - firstVertex.y();
        Q.row(2)<< secondVertex.z() - firstVertex.z(), thirdVertex.z() - firstVertex.z(), fourthVertex.z() - firstVertex.z();

        return Q;
      }

      /// Matrix Q for linear map x = Q * x_r + b from reference tetrahedron [0,1]x[0,1]x[0,1]/2 to tetrahedron with x points
      /// vertices the tetrahedron to map vertices, size 3 x 4
      /// return the resulting value, size 3 x 3
      inline Eigen::Vector3d b(const Eigen::Vector3d& firstVertex) const
      {
        return firstVertex;
      }

      /// Map from the tetrahedron reference element [0,1]x[0,1]x[0,1]/2 to the polygon x = F(x_r) = Q * x_r + b
      /// \param vertices the tetrahedron to map vertices, size 3 x 4
      /// \return the map data
      MapTetrahedronData Compute(const Eigen::MatrixXd& vertices) const;

      /// Map from the tetrahedron reference element [0,1]x[0,1]x[0,1]/2 to the polygon x = F(x_r) = Q * x_r + b
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the mapped points, size 3 x numPoints
      inline Eigen::MatrixXd F(const MapTetrahedronData& mapData,
                               const Eigen::MatrixXd& x) const
      { return (mapData.Q * x).colwise() + mapData.b; }

      inline Eigen::MatrixXd FInv(const MapTetrahedronData& mapData,
                                  const Eigen::MatrixXd& x) const
      { return mapData.Q.inverse() * (x.colwise() - mapData.b); }

      /// Compute the jacobian matrix of the transformation F
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the Q matrix for each points, size 2 x (2 * numPoints)
      Eigen::MatrixXd J(const MapTetrahedronData& mapData,
                        const Eigen::MatrixXd& x) const;
      /// Compute the determinant of the jacobian matrix of the trasformation
      /// \param mapData the map data computed
      /// \param x points in reference tetrahedron, size 3 x numPoints
      /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
      Eigen::VectorXd DetJ(const MapTetrahedronData& mapData,
                           const Eigen::MatrixXd& x) const;

      /// Compute the determinant of the jacobian matrix of the trasformation
      /// \param mapData the map data computed
      /// \return the determinant of Jacobian matrix
      inline double DetJ(const MapTetrahedronData& mapData) const
      {
        return mapData.Q.determinant();
      };
  };
}

#endif // __MapTetrahedron_H
