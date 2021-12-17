#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  Eigen::Vector3d GeometryUtilities::PolygonNormal(const MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    Vector3d normal;

    normal.setZero();
    const unsigned int& numVertices = polygonVertices.cols();

    for (unsigned int i = 0; i < numVertices; i++)
    {
      Vector3d edge = polygonVertices.col((i + 1) % numVertices) - polygonVertices.col(i);
      Vector3d edgePrevious = polygonVertices.col((i - 1) % numVertices) - polygonVertices.col(i);
      normal.noalias() += edge.cross(edgePrevious);
    }

    return normal.normalized();
  }
  // ***************************************************************************
  VectorXd GeometryUtilities::PolygonEdgeLengths(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    VectorXd edgeLengths(numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeLengths[v] = SegmentLength(polygonVertices.col(v),
                                     polygonVertices.col((v + 1) % numVertices));
    }

    return edgeLengths;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolygonEdgeTangents(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    MatrixXd edgeTangents(3, numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeTangents.col(v) = SegmentTangent(polygonVertices.col(v),
                                           polygonVertices.col((v + 1) % numVertices));
    }

    return edgeTangents;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolygonEdgeNormals(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    MatrixXd edgeNormals(3, numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeNormals.col(v) = SegmentNormal(polygonVertices.col(v),
                                         polygonVertices.col((v + 1) % numVertices));
    }

    return edgeNormals;
  }
  // ***************************************************************************
  Vector3d GeometryUtilities::PolygonCentroid(const Eigen::MatrixXd& polygonVertices,
                                              const double& polygonArea) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    Eigen::Vector3d centroid;
    centroid.setZero();

    const unsigned int& numVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d& vertex = polygonVertices.col(v);
      const Eigen::Vector3d& nextVertex = polygonVertices.col((v + 1) % numVertices);
      centroid.x() += (vertex.x() + nextVertex.x()) *
                      (vertex.x() * nextVertex.y() - nextVertex.x() * vertex.y());

      centroid.y() += (vertex.y() + nextVertex.y()) *
                      (vertex.x() * nextVertex.y() - nextVertex.x() * vertex.y());
    }

    centroid /= 6.0 * polygonArea;

    return centroid;
  }
  // ***************************************************************************
  double GeometryUtilities::PolygonDiameter(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    double diameter = 0.0;

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d& vertexOne = polygonVertices.col(v);
      for (unsigned int w = v + 1; w < numVertices; w++)
      {
        const Eigen::Vector3d& vertexTwo = polygonVertices.col(w);
        double distance = PointDistance(vertexOne, vertexTwo);
        if (Compare1DValues(diameter, distance) == CompareTypes::FirstBeforeSecond)
          diameter = distance;
      }
    }

    return diameter;
  }
  // ***************************************************************************
  Matrix3d GeometryUtilities::PolygonRotationMatrix(const Eigen::MatrixXd& polygonVertices,
                                                    const Eigen::Vector3d& polygonNormal,
                                                    const Eigen::Vector3d& polygonTranslation) const
  {
    Output::Assert(Compare1DValues(polygonNormal.norm(), 1.0) == CompareTypes::Coincident);

    const unsigned int& numVertices = polygonVertices.cols();
    MatrixXd Z(3, numVertices);
    MatrixXd W(3, numVertices);
    Matrix3d H;
    Vector3d V1mV0 = polygonVertices.col(1) - polygonTranslation;
    double normVectorOne = V1mV0.norm();
    Z.col(0) = V1mV0;
    W.col(0) << normVectorOne, 0.0, 0.0;
    for (unsigned int i = 2; i < numVertices; i++)
    {
      Vector3d VimV0 = polygonVertices.col(i) - polygonTranslation;
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
    Z.col(numVertices - 1) = polygonNormal;
    W.col(numVertices - 1)<< 0.0, 0.0, 1.0;
    H = W * Z.transpose();
    JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

    return svd.matrixV() * (svd.matrixU()).transpose();
  }
  // ***************************************************************************
  bool GeometryUtilities::PolygonIsConvex(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
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
  GeometryUtilities::PolygonTypes GeometryUtilities::PolygonType(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
    if (numVertices == 3)
      return PolygonTypes::Triangle;
    else if (numVertices == 4)
      return PolygonTypes::Quadrilateral;
    else
      return PolygonTypes::Generic;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolygonTriangulationByFirstVertex(const Eigen::MatrixXd& polygonVertices) const
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
  vector<unsigned int> GeometryUtilities::PolygonTriangulationByInternalPoint(const Eigen::MatrixXd& polygonVertices,
                                                                              const Eigen::Vector3d& point) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);
    Output::Assert(PointPolygonPosition(point, polygonVertices).Type == PointPolygonPositionResult::Types::Inside);

    const unsigned int numPolygonVertices = polygonVertices.cols();
    vector<unsigned int> triangles(3 * numPolygonVertices);

    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      triangles[3 * v] = numPolygonVertices;
      triangles[3 * v + 1] = v;
      triangles[3 * v + 2] = (v + 1) % numPolygonVertices;
    }

    Output::Assert(triangles.size() % 3 == 0);

    return triangles;
  }
  // ***************************************************************************
  double GeometryUtilities::PolygonArea(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
    double area = 0.0;
    Eigen::VectorXd xPoints(numVertices + 1);
    Eigen::VectorXd yPoints(numVertices + 1);
    xPoints<< polygonVertices.row(0).transpose(), polygonVertices(0, 0);
    yPoints<< polygonVertices.row(1).transpose(), polygonVertices(1, 0);

    return 0.5 * (xPoints.segment(0, numVertices).dot(yPoints.segment(1, numVertices)) -
                  xPoints.segment(1, numVertices).dot(yPoints.segment(0, numVertices)));
  }
  // ***************************************************************************
  double GeometryUtilities::PolygonArea(const Eigen::MatrixXd& polygonVertices,
                                        const vector<unsigned int>& polygonTriangulation) const
  {
    Output::Assert(PointsAre2D(polygonVertices) &&
                   polygonVertices.cols() > 2 &&
                   polygonTriangulation.size() > 0 &&
                   polygonTriangulation.size() % 3 == 0);

    const unsigned int numTriangles = polygonTriangulation.size() / 3;
    double area = 0.0;
    for (unsigned int t = 0; t < numTriangles; t++)
    {
      Eigen::Matrix3d triangleVertices;
      triangleVertices.col(0)<< polygonVertices.col(polygonTriangulation[3 * t]);
      triangleVertices.col(1)<< polygonVertices.col(polygonTriangulation[3 * t + 1]);
      triangleVertices.col(2)<< polygonVertices.col(polygonTriangulation[3 * t + 2]);
      area += PolygonArea(triangleVertices);
    }

    return area;
  }
  // ***************************************************************************
  GeometryUtilities::PolygonCirclePositionTypes GeometryUtilities::PolygonCirclePosition(const Eigen::MatrixXd& polygonVertices,
                                                                                         const Eigen::Vector3d& circleCenter,
                                                                                         const double& circleRadius,
                                                                                         const vector<PointCirclePositionResult>& vertexPositions,
                                                                                         const IntersectionPolygonCircleResult& polygonCircleIntersections) const
  {
    // check circle center position respect the polygon
    GeometryUtilities::PointPolygonPositionResult centerPosition = PointPolygonPosition(circleCenter,
                                                                                        polygonVertices);
    Output::Assert(centerPosition.Type != Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Unknown);

    const unsigned int& numVertices = polygonVertices.cols();
    const unsigned int& numIntersections = polygonCircleIntersections.Intersections.size();

    // compute polygon vertices position respect the circle
    bool oneVertexOutsideCircle = false;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      if (vertexPositions[v] == PointCirclePositionResult::Outside)
        oneVertexOutsideCircle = true;
    }

    switch (centerPosition.Type)
    {
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside:
      {
        // the circle is outside the polygon
        if (oneVertexOutsideCircle &&
            numIntersections == 0) // vertices are far from the circle, no intersection found
          return PolygonCirclePositionTypes::PolygonOutsideCircleNoIntersection;
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 0) // all vertices are inside the circle and there are no intersections then the polygon is inside the circle
          return PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection;
        else if (oneVertexOutsideCircle &&
                 numIntersections == 1) // at least one vertex is outside the circle and only one intersection found, the circle is outside and touches the polygon in one vertex
        {
          if (polygonCircleIntersections.Intersections[0].IndexType == IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex)
            return PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionOnVertex;
          else if (polygonCircleIntersections.Intersections[0].IndexType == IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge &&
                   polygonCircleIntersections.Intersections[0].Type == IntersectionPolygonCircleResult::Intersection::Types::Tangent)
            return PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionTangentOnEdge;
        }
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 1) // all vertices are inside the circle and one intersection found
          return PolygonCirclePositionTypes::PolygonInsideCircleOneVertexIntersection;
        else if (numIntersections > 1)
          return PolygonCirclePositionTypes::CirclePolygonMultipleIntersections;
      }
      break;
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
      {
        // the circle center is inside the polygon
        if (oneVertexOutsideCircle &&
            numIntersections == 0) // if at least one vertex is ouside the circle and there are no intersections than the circle is inside the polygon
          return PolygonCirclePositionTypes::CircleInsidePolygonNoIntersection;
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 0) // if at least a vertex is not outside the circle and there are no intersections then the polygon is inside the circle
          return PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection;
        else if (numIntersections == 1) // only one intersection found, the circle touch the polygon tangent in one edge
          return PolygonCirclePositionTypes::CircleInsidePolygonOneIntersectionTangentOnEdge;
        else if (numIntersections > 1)
        {
          // check if the circle intersects only the vertices (a sort of circumscribed circle)
          bool onlyVerticesIntersections = true;
          for (unsigned int i = 0; i < numIntersections; i++)
          {
            if (polygonCircleIntersections.Intersections[i].IndexType !=
                IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex)
            {
              onlyVerticesIntersections = false;
              break;
            }
          }

          if (onlyVerticesIntersections)
            return PolygonCirclePositionTypes::PolygonInsideCircleIntersectionOnlyOnVertices;

          return PolygonCirclePositionTypes::CirclePolygonMultipleIntersections;
        }
      }
      break;
      default:
      break;
    }

    throw runtime_error("PolygonCirclePosition failed");
  }
  // ***************************************************************************
}
