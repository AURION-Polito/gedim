#ifndef __QuadratureData_H
#define __QuadratureData_H

#include "Eigen/Eigen"

namespace Gedim
{
  namespace Quadrature
  {
    struct QuadratureData final
    {
        Eigen::MatrixXd Points;
        Eigen::VectorXd Weights;
    };
  }
}

#endif // __QuadratureData_H
