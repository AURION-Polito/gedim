#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  Eigen::Vector3d GeometryUtilities::PolygonNormal(const MatrixXd& polygonVertices) const
  {
    Vector3d normal;

    normal.setZero();
    unsigned int numVertices =  polygonVertices.cols();

    for (unsigned int i = 0; i < numVertices; i++)
    {
      Vector3d edge = polygonVertices.col((i + 1) % numVertices) - polygonVertices.col(i);
      Vector3d edgePrevious = polygonVertices.col((i - 1) % numVertices) - polygonVertices.col(i);
      normal.noalias() += edge.cross(edgePrevious);
    }

    return normal.normalized();
  }
  // ***************************************************************************
  void GeometryUtilities::PolygonRotation(const MatrixXd& polygonVertices,
                                          const Vector3d& normal,
                                          Matrix3d& rotationMatrix,
                                          Vector3d& translation) const
  {
    Output::Assert(Compare1DValues(normal.norm(), 1.0) == CompareTypes::Coincident);

    unsigned int numVertices = polygonVertices.cols();
    MatrixXd Z(3, numVertices);
    MatrixXd W(3, numVertices);
    Matrix3d H;
    Vector3d V1mV0 = polygonVertices.col(1) -  polygonVertices.col(0);
    double normVectorOne = V1mV0.norm();
    Z.col(0) = V1mV0;
    W.col(0) << normVectorOne, 0.0, 0.0;
    for (unsigned int i = 2; i < numVertices; i++)
    {
      Vector3d VimV0 = polygonVertices.col(i) - polygonVertices.col(0);
      Z.col(i - 1) = VimV0;

      double normVectorI = VimV0.norm();
      double cosTheta = VimV0.dot(V1mV0) / (normVectorOne * normVectorI);

      if (Compare1DValues(cosTheta, 1.0) == CompareTypes::SecondBeforeFirst)
        W.col(i - 1) << normVectorI, 0.0, 0.0;
      else if (Compare1DValues(cosTheta, -1.0) == CompareTypes::FirstBeforeSecond)
        W.col(i - 1) << -normVectorI, 0.0, 0.0;
      else if (Compare1DValues(cosTheta, 0.0) == CompareTypes::Coincident)
        W.col(i - 1) << 0.0, normVectorI, 0.0;
      else
        W.col(i - 1) << normVectorI * cosTheta, normVectorI * sqrt(1.0 - cosTheta*cosTheta), 0;
    }
    Z.col(numVertices - 1) = normal;
    W.col(numVertices - 1)<< 0.0, 0.0, 1.0;
    H = W * Z.transpose();
    JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

    rotationMatrix =  svd.matrixV() * (svd.matrixU()).transpose();
    translation = polygonVertices.col(0);
  }
  // ***************************************************************************
  bool GeometryUtilities::PolygonIsConvex(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.row(2).isZero(_configuration.Tolerance));

    const unsigned int numVertices = polygonVertices.cols();
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d edgeOrigin = polygonVertices.col(v);
      const Eigen::Vector3d edgeEnd = polygonVertices.col((v + 1) % numVertices);
      const Eigen::Vector3d vertexNextEdge = polygonVertices.col((v + 2) % numVertices);

      if (PointSegmentPosition(vertexNextEdge,
                               edgeOrigin,
                               edgeEnd) == PointSegmentPositionTypes::RightTheSegment)
        return false;
    }

    return true;
  }
  // ***************************************************************************
}
