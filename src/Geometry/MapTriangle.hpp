#ifndef __MapTriangle_H
#define __MapTriangle_H

#include "Eigen"

using namespace std;

namespace Gedim
{
  class MapTriangle
  {
    private:
      /// Matrix Q for linear map x = Q * x_r + b from reference triangle [0,1]x[0,1]/2 to triangle with x points
      /// vertices the triangle to map vertices, size 3 x 3
      /// return the resulting value, size 3 x 3
      inline Eigen::Matrix3d Q(const Eigen::Matrix3d& vertices) const
      {
        Eigen::Matrix3d Q;
        Q.row(0)<< vertices(0, 1) - vertices(0, 0), vertices(0, 2) - vertices(0, 0), 0.0;
        Q.row(1)<< vertices(1, 1) - vertices(1, 0), vertices(1, 2) - vertices(1, 0), 0.0;
        Q.row(2)<< 0.0, 0.0, 1.0;
        return Q;
      }

      /// Matrix Q for linear map x = Q * x_r + b from reference triangle [0,1]x[0,1]/2 to triangle with x points
      /// vertices the triangle to map vertices, size 3 x 3
      /// return the resulting value, size 3 x 3
      inline Eigen::Vector3d b(const Eigen::Matrix3d& vertices) const
      {
        return vertices.col(0);
      }

    public:
      MapTriangle() {}
      ~MapTriangle() {}

      /// Map from the triangle reference element [0,1]x[0,1] to the polygon x = F(x_r) = Q * x_r + b
      /// vertices the triangle to map vertices, size 3 x 3
      /// x points in reference triangle, size 3 x numPoints
      /// \return the mapped points, size 3 x numPoints
      Eigen::MatrixXd F(const Eigen::Matrix3d& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the jacobian matrix of the transformation F
      /// x matrix of points between 0.0 and 1.0, size 2 x numPoints
      /// \return the Q matrix for each points, size 2 x (2 * numPoints)
      Eigen::MatrixXd J(const Eigen::Matrix3d& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the determinant of the jacobian matrix of the trasformation
      /// x matrix of points between 0.0 and 1.0, size 2 x numPoints
      /// \return the determinant of Jacobian matrix for each points, size 1 x numPoints
      Eigen::VectorXd DetJ(const Eigen::Matrix3d& vertices,
                           const Eigen::MatrixXd& x) const;
  };
}

#endif // __MapTriangle_H
