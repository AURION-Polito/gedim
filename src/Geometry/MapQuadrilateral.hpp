#ifndef __MappingQuadrilateral_H
#define __MappingQuadrilateral_H

#include "Eigen"

using namespace std;

namespace Gedim
{
  class MapQuadrilateral
  {
    private:
      /// map square quadrature points in [-1,1]
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting shifted points, size 3 x numPoints
      inline void ShiftSquareQuadraturePoints(const Eigen::MatrixXd& x,
                                              Eigen::MatrixXd& result) const
      { result = 2.0 * x.array() - 1.0; }
      /// Shape functions 1
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// return the resulting value, size numPoints
      inline Eigen::VectorXd Psi1(const Eigen::MatrixXd& x) const
      { return (1.0 - x.row(0).array()) * (1.0 - x.row(1).array()) / 4.0; }
      /// Shape functions 2
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// return the resulting value, size numPoints
      inline Eigen::VectorXd Psi2(const Eigen::MatrixXd& x) const
      { return (1.0 + x.row(0).array())*(1.0 - x.row(1).array()) / 4.0; }
      /// Shape functions 3
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// return the resulting value, size numPoints
      inline Eigen::VectorXd Psi3(const Eigen::MatrixXd& x) const
      { return (1.0 + x.row(0).array())*(1.0 + x.row(1).array()) / 4.0; }
      /// Shape functions 4
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// return the resulting value, size numPoints
      inline Eigen::VectorXd Psi4(const Eigen::MatrixXd& x) const
      { return (1.0 - x.row(0).array())*(1.0 + x.row(1).array()) / 4.0; }
      /// Shape functions derivatives 1
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi11(const Eigen::MatrixXd& x) const
      { return -(1.0 - x.row(1).array()) / 4.0; }
      /// Shape functions x derivatives 2
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi21(const Eigen::MatrixXd& x) const
      { return (1.0 - x.row(1).array()) / 4.0; }
      /// Shape functions x derivatives 3
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi31(const Eigen::MatrixXd& x) const
      { return (1.0 + x.row(1).array()) / 4.0; }
      /// Shape functions x derivatives 4
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi41(const Eigen::MatrixXd& x) const
      { return -(1.0 + x.row(1).array()) / 4.0; }
      /// Shape functions x.row(1).array() derivatives 1
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi12(const Eigen::MatrixXd& x) const
      { return -(1.0 - x.row(0).array()) / 4.0; }
      /// Shape functions x.row(1).array() derivatives 2
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi22(const Eigen::MatrixXd& x) const
      { return -(1.0 + x.row(0).array()) / 4.0; }
      /// Shape functions x.row(1).array() derivatives 3
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi32(const Eigen::MatrixXd& x) const
      { return (1.0 + x.row(0).array()) / 4.0; }
      /// Shape functions x.row(1).array() derivatives 4
      /// x matrix of quadrature points between 0.0 and 1.0, size 3 x numPoints
      /// result the resulting value, size numPoints
      inline Eigen::VectorXd dPsi42(const Eigen::MatrixXd& x) const
      { return (1.0 - x.row(0).array()) / 4.0; }

    public:
      MapQuadrilateral() {}
      ~MapQuadrilateral() {}

      /// Map from the square reference element [0,1]x[0,1] to the polygon
      /// x matrix of points between 0.0 and 1.0, size 3 x numPoints
      /// \return the mapped points, size 3 x numPoints
      Eigen::MatrixXd F(const Eigen::MatrixXd& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the jacobian matrix of the transformation
      /// x matrix of points between 0.0 and 1.0, size 3 x numPoints
      /// \return the Jacobian matrix for each points, size 2 x (3 x numPoints)
      Eigen::MatrixXd J(const Eigen::MatrixXd& vertices,
                        const Eigen::MatrixXd& x) const;
      /// Compute the determinant of the jacobian matrix of the trasformation
      /// x matrix of points between 0.0 and 1.0, size 3 x numPoints
      /// \return the Jacobian matrix for each points, size 2 x (3 x numPoints)
      Eigen::VectorXd DetJ(const Eigen::MatrixXd& vertices,
                           const Eigen::MatrixXd& x) const;
  };
}

#endif // __MappingQuadrilateral_H
