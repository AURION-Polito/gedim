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
    Output::Assert(PointsAre2D(polygonVertices));

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
  vector<unsigned int> GeometryUtilities::PolygonTriangulation(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    list<unsigned int> triangleList;

    const unsigned int numPolygonVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      const unsigned int nextVertex = (v + 1) % numPolygonVertices;
      const unsigned int nextNextVertex = (v + 2) % numPolygonVertices;

      if (nextNextVertex == 0)
        break;

      triangleList.push_back(0);
      triangleList.push_back(nextVertex);
      triangleList.push_back(nextNextVertex);
    }

    Output::Assert(triangleList.size() % 3 == 0);

    return vector<unsigned int>(triangleList.begin(), triangleList.end());
  }
  // ***************************************************************************
  GeometryUtilities::PolygonCirclePositionTypes GeometryUtilities::PolygonCirclePosition(const Eigen::MatrixXd& polygonVertices,
                                                                                         const Eigen::Vector3d& circleCenter,
                                                                                         const double& circleRadius,
                                                                                         const IntersectionPolygonCircleResult& polygonCircleIntersections) const
  {
    GeometryUtilities::PointPolygonPositionResult centerPosition = PointPolygonPosition(circleCenter,
                                                                                        polygonVertices);
    Output::Assert(centerPosition.Type != Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Unknown);

    const unsigned int numVertices = polygonVertices.cols();
    vector<PointCirclePositionResult> vertexPositions(numVertices);
    bool oneVertexOutsideCircle = false;
    bool oneVertexOnCircleBorder = false;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      vertexPositions[v] = PointCirclePosition(polygonVertices.col(v),
                                               circleCenter,
                                               circleRadius);
      if (vertexPositions[v] == PointCirclePositionResult::Outside)
        oneVertexOutsideCircle = true;
      if (vertexPositions[v] == PointCirclePositionResult::OnBorder)
        oneVertexOnCircleBorder = true;
    }

    switch (centerPosition.Type)
    {
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside:
      {
        if (oneVertexOutsideCircle &&
            polygonCircleIntersections.Intersections.size() == 0)
          return PolygonCirclePositionTypes::CircleNotIntersectPolygon;
        else if (!oneVertexOutsideCircle &&
                 polygonCircleIntersections.Intersections.size() == 0)
          return PolygonCirclePositionTypes::PolygonInsideCircle;
        else if (polygonCircleIntersections.Intersections.size() == 1)
          return PolygonCirclePositionTypes::CircleIntersectsPolygonWithOnePoint;
      }
      break;
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
      {
        if (oneVertexOutsideCircle &&
            polygonCircleIntersections.Intersections.size() == 0)
          return PolygonCirclePositionTypes::CircleInsidePolygon;
        else if (!oneVertexOutsideCircle &&
                 polygonCircleIntersections.Intersections.size() == 0)
          return PolygonCirclePositionTypes::PolygonInsideCircle;
        else if (oneVertexOnCircleBorder &&
                 polygonCircleIntersections.Intersections.size() > 0)
          return PolygonCirclePositionTypes::PolygonIntersectsCircleOnBorder;
      }
      break;
      default:
      break;
    }

    return PolygonCirclePositionTypes::Unknown;
  }
  // ***************************************************************************
}
