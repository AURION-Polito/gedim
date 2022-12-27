#include "MapTriangle.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  MapTriangle::MapTriangleData MapTriangle::Compute(const Eigen::Matrix3d& vertices) const
  {
    MapTriangleData result;

    result.B = B(vertices);
    result.BInv = result.B.inverse();
    result.b = b(vertices);
    result.DetB = result.B.determinant();
    result.DetBInv = result.BInv.determinant();

    return result;
  }
  // ***************************************************************************
  MatrixXd MapTriangle::J(const MapTriangleData& mapData,
                          const MatrixXd& x) const
  {
    const unsigned int numPoints = x.cols();
    MatrixXd jacb(3, 3 * numPoints);

    for (unsigned int p = 0; p < numPoints; p++)
      jacb.block(0, 3 * p, 3, 3) = mapData.B;

    return jacb;
  }
  // ***************************************************************************
}
