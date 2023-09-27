#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  Eigen::VectorXd GeometryUtilities::PointDistances(const Eigen::MatrixXd& points,
                                                    const Eigen::Vector3d& point) const
  {
    return (points.colwise() - point).colwise().norm();
  }
  // ***************************************************************************
  Eigen::MatrixXd GeometryUtilities::PointsDistance(const Eigen::MatrixXd& points) const
  {
    Output::Assert(points.rows() == 3);

    const unsigned int& numPoints = points.cols();

    Eigen::MatrixXd distances = Eigen::MatrixXd::Constant(numPoints,
                                                          numPoints,
                                                          -1.0);

    for (unsigned int v = 0; v < numPoints; v++)
    {
      const Eigen::Vector3d& vertexOne = points.col(v);
      for (unsigned int w = v + 1; w < numPoints; w++)
      {
        const Eigen::Vector3d& vertexTwo = points.col(w);
        distances(v, w) = PointDistance(vertexOne, vertexTwo);
      }
    }

    return distances;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PointsBoundingBox(const Eigen::MatrixXd& points) const
  {
    Eigen::MatrixXd boundingBox(3, 2);

    boundingBox(0, 0) = points.row(0).minCoeff();
    boundingBox(0, 1) = points.row(0).maxCoeff();
    boundingBox(1, 0) = points.row(1).minCoeff();
    boundingBox(1, 1) = points.row(1).maxCoeff();
    boundingBox(2, 0) = points.row(2).minCoeff();
    boundingBox(2, 1) = points.row(2).maxCoeff();

    return boundingBox;
  }
  // ***************************************************************************
  double GeometryUtilities::PointsMaxDistance(const Eigen::MatrixXd& points) const
  {
    Output::Assert(points.rows() == 3);

    const unsigned int& numVertices = points.cols();

    double maxDistance = 0.0;

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d& vertexOne = points.col(v);
      for (unsigned int w = v + 1; w < numVertices; w++)
      {
        const Eigen::Vector3d& vertexTwo = points.col(w);
        double distance = PointDistance(vertexOne, vertexTwo);
        if (Compare1DValues(maxDistance, distance) == CompareTypes::FirstBeforeSecond)
          maxDistance = distance;
      }
    }

    return maxDistance;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::FindPointInPoints(const Eigen::MatrixXd& points,
                                                            const Eigen::Vector3d& point) const
  {
    VectorXd pointDistances = PointDistances(points,
                                             point);
    list<unsigned int> indices;
    for (unsigned int p = 0; p < pointDistances.size(); p++)
    {
      if (IsValue1DZero(pointDistances[p]))
        indices.push_back(p);
    }

    return vector<unsigned int>(indices.begin(), indices.end());
  }
  // ***************************************************************************
  bool GeometryUtilities::IsPointOnLine(const Eigen::Vector3d& point,
                                        const Eigen::Vector3d& lineOrigin,
                                        const Eigen::Vector3d& lineTangent,
                                        const double& lineTangentSquaredLength) const
  {
    const Eigen::Vector3d pointDirection = (point - lineOrigin);

    if (IsValue2DZero(pointDirection.squaredNorm()) ||
        IsValue2DZero(pointDirection.cross(lineTangent).squaredNorm() / lineTangentSquaredLength))
      return true;

    return false;
  }
  // ***************************************************************************
  GeometryUtilities::PointSegmentPositionTypes GeometryUtilities::PointSegmentPosition(const Vector3d& point,
                                                                                       const Vector3d& segmentOrigin,
                                                                                       const Vector3d& segmentEnd) const
  {
    const Vector3d segmentTangent = SegmentTangent(segmentOrigin, segmentEnd).normalized();
    Vector3d pointDirection = (point - segmentOrigin).normalized();
    double pointDirectionSquaredNorm = (point - segmentOrigin).squaredNorm();

    // check if point is on line
    if (IsValue2DZero(pointDirectionSquaredNorm) ||
        IsValue2DZero(pointDirection.cross(segmentTangent).squaredNorm()))
    {
      // compute curvilinear coordinate of point (segment is between 0.0 and 1.0)
      return PointSegmentPosition(PointCurvilinearCoordinate(point, segmentOrigin, segmentEnd));
    }

    // check point position out of line, supported only in 2D plane
    Output::Assert(IsValue1DZero(segmentTangent.z()) && IsValue1DZero(pointDirection.z()));
    // rotate the 2D tangent by 90 degrees
    Vector3d normalTangent(segmentTangent.y(), -segmentTangent.x(), 0.0);

    if (IsValue1DPositive(pointDirection.dot(normalTangent)))
      return PointSegmentPositionTypes::RightTheSegment;
    else
      return PointSegmentPositionTypes::LeftTheSegment;
  }
  // ***************************************************************************
  GeometryUtilities::PointSegmentPositionTypes GeometryUtilities::PointSegmentPosition(const double& curvilinearCoordinate) const
  {
    if (IsValue1DZero(curvilinearCoordinate))
      return PointSegmentPositionTypes::OnSegmentOrigin;

    if (Compare1DValues(1.0, curvilinearCoordinate) == CompareTypes::Coincident)
      return PointSegmentPositionTypes::OnSegmentEnd;

    if (IsValue1DPositive(curvilinearCoordinate) &&
        Compare1DValues(1.0, curvilinearCoordinate) == CompareTypes::SecondBeforeFirst)
      return PointSegmentPositionTypes::InsideSegment;

    if (IsValue1DNegative(curvilinearCoordinate))
      return PointSegmentPositionTypes::OnSegmentLineBeforeOrigin;

    if (Compare1DValues(1.0, curvilinearCoordinate) == CompareTypes::FirstBeforeSecond)
      return PointSegmentPositionTypes::OnSegmentLineAfterEnd;

    throw runtime_error("PointSegmentPosition failed");
  }
  // ***************************************************************************
  double GeometryUtilities::PointPlaneDistance(const Eigen::Vector3d& point,
                                               const std::array<Eigen::Vector3d, 3>& planePoints) const
  {
    const Eigen::Vector3d pad = planePoints[0] - point;
    const double padNorm = pad.norm();
    if (IsValue1DZero(padNorm))
      return 0.0;

    const Eigen::Vector3d pbd = planePoints[1] - point;
    const double pbdNorm = pbd.norm();
    if (IsValue1DZero(pbdNorm))
      return 0.0;

    const Eigen::Vector3d pcd = planePoints[2] - point;
    const double pcdNorm = pcd.norm();
    if (IsValue1DZero(pcdNorm))
      return 0.0;

    return (pad[0] * (pbd[1] * pcd[2] - pbd[2] * pcd[1])
        + pbd[0] * (pcd[1] * pad[2] - pcd[2] * pad[1])
        + pcd[0] * (pad[1] * pbd[2] - pad[2] * pbd[1])) /
        (padNorm * pbdNorm * pcdNorm);
  }
  // ***************************************************************************
  double GeometryUtilities::PointPlaneDistance(const Eigen::Vector3d& point,
                                               const Eigen::Vector3d& planeNormal,
                                               const Eigen::Vector3d& planeOrigin) const
  {
    const Eigen::Vector3d pod = point - planeOrigin;
    const double podNorm = pod.norm();
    if (IsValue1DZero(podNorm))
      return 0.0;

    return planeNormal.dot(pod) / podNorm;
  }
  // ***************************************************************************
  GeometryUtilities::PointPlanePositionTypes GeometryUtilities::PointPlanePosition(const double& pointPlaneDistance) const
  {
    if (IsValue1DZero(pointPlaneDistance))
      return PointPlanePositionTypes::OnPlane;
    else if (IsValue1DNegative(pointPlaneDistance))
      return PointPlanePositionTypes::Negative;
    else
      return PointPlanePositionTypes::Positive;
  }
  // ***************************************************************************
  GeometryUtilities::PointPolygonPositionResult GeometryUtilities::PointPolygonPosition(const Vector3d& point,
                                                                                        const MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.cols() > 2 && PointsAre2D(point));

    GeometryUtilities::PointPolygonPositionResult result;

    unsigned int numVertices =  polygonVertices.cols();
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Vector3d& vertex = polygonVertices.col(v);
      const Vector3d& nextVertex = polygonVertices.col((v + 1) % numVertices);

      GeometryUtilities::PointSegmentPositionTypes resultSegment = PointSegmentPosition(point, vertex, nextVertex);

      if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::RightTheSegment)
      {
        result.Type = GeometryUtilities::PointPolygonPositionResult::Types::Outside;
        return result;
      }
      else if (resultSegment != GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment)
      {
        if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
        {
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge;
          result.BorderIndex = v;
          return result;
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
          result.BorderIndex = v;
          return result;
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
          result.BorderIndex = (v + 1) % numVertices;
          return result;
        }
        else
          continue; // check for future aligned points on edges
      }
    }

    result.Type = GeometryUtilities::PointPolygonPositionResult::Types::Inside;
    return result;
  }
  // ***************************************************************************
  GeometryUtilities::PointPolygonPositionResult GeometryUtilities::PointPolygonPosition_RayCasting(const Eigen::Vector3d& point,
                                                                                                   const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.cols() > 2 && PointsAre2D(point));

    GeometryUtilities::PointPolygonPositionResult result;

    const unsigned int numVertices =  polygonVertices.cols();
    unsigned int numIntersections = 0;

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const unsigned int v_next = (v + 1) % numVertices;
      const Vector3d& edgeOrigin = polygonVertices.col(v);
      const Vector3d& edgeEnd = polygonVertices.col(v_next);
      const Vector3d edgeTangent = SegmentTangent(edgeOrigin,
                                                  edgeEnd);

      if (IsValue1DZero(edgeTangent.y()))
      {
        if (IsValue1DZero(point.y() - edgeOrigin.y()))
        {
          const double intersection_point_curvilinearCoordinate = PointLineCurvilinearCoordinate(point,
                                                                                                 edgeOrigin,
                                                                                                 edgeTangent,
                                                                                                 edgeTangent.squaredNorm());
          const PointSegmentPositionTypes intersection_point_position = PointSegmentPosition(intersection_point_curvilinearCoordinate);

          switch (intersection_point_position)
          {
            case PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
              continue; // intersection parallel on edge
            case PointSegmentPositionTypes::OnSegmentLineAfterEnd:
              continue; // intersection not interesting
            case PointSegmentPositionTypes::OnSegmentOrigin:
            {
              result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
              result.BorderIndex = v;
              return result;
            }
            case PointSegmentPositionTypes::InsideSegment:
            {
              result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge;
              result.BorderIndex = v;
              return result;
            }
            case PointSegmentPositionTypes::OnSegmentEnd:
            {
              result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
              result.BorderIndex = v_next;
              return result;
            }
            default:
              throw runtime_error("Intersection point not expected");
          }
        }
        else
          continue;
      }

      const double intersection_edge_curvilinearCoordinate = (point.y() - edgeOrigin.y()) / edgeTangent.y();

      const PointSegmentPositionTypes intersection_edge_position = PointSegmentPosition(intersection_edge_curvilinearCoordinate);

      switch (intersection_edge_position)
      {
        case PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
        case PointSegmentPositionTypes::OnSegmentLineAfterEnd:
          continue; // intersection outside the edge
        case PointSegmentPositionTypes::OnSegmentOrigin:
        case PointSegmentPositionTypes::InsideSegment:
        case PointSegmentPositionTypes::OnSegmentEnd:
        {
          const double x_intersection = edgeOrigin.x() +
                                        intersection_edge_curvilinearCoordinate * edgeTangent.x();

          const double intersection_point_curvilinearCoordinate = (x_intersection - point.x());
          const PointSegmentPositionTypes intersection_point_position = PointSegmentPosition(intersection_point_curvilinearCoordinate);

          switch (intersection_point_position)
          {
            case PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
              continue; // intersection not interesting
            case PointSegmentPositionTypes::OnSegmentOrigin:
            {
              if (intersection_edge_position == PointSegmentPositionTypes::OnSegmentOrigin)
              {
                result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
                result.BorderIndex = v;
                return result;
              }

              if (intersection_edge_position == PointSegmentPositionTypes::OnSegmentEnd)
              {
                result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
                result.BorderIndex = v_next;
                return result;
              }

              result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge;
              result.BorderIndex = v;
              return result;
            }
            case PointSegmentPositionTypes::InsideSegment:
            case PointSegmentPositionTypes::OnSegmentEnd:
            case PointSegmentPositionTypes::OnSegmentLineAfterEnd:
            {
              if (intersection_edge_position == PointSegmentPositionTypes::OnSegmentOrigin)
              {
                // the inserction count only if the other edge vertex is uder the ray
                if (IsValue1DGreater(point.y(), edgeEnd.y()))
                  numIntersections++;
                continue;
              }

              if (intersection_edge_position == PointSegmentPositionTypes::OnSegmentEnd)
              {
                // the inserction count only if the other edge vertex is uder the ray
                if (IsValue1DGreater(point.y(), edgeOrigin.y()))
                  numIntersections++;
                continue;
              }

              numIntersections++;
              continue;
            }
            default:
              throw runtime_error("Intersection point not expected");
          }

          break;
        }
        default:
          throw runtime_error("Intersection edge not expected");
      }
    }

    result.Type = (numIntersections % 2) == 1 ? GeometryUtilities::PointPolygonPositionResult::Types::Inside :
                                                GeometryUtilities::PointPolygonPositionResult::Types::Outside;
    return result;
  }
  // ***************************************************************************
  GeometryUtilities::PointPolyhedronPositionResult GeometryUtilities::PointPolyhedronPosition(const Eigen::Vector3d& point,
                                                                                              const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                              const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                              const vector<Eigen::MatrixXd>& polyhedronFaceRotatedVertices,
                                                                                              const vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                                              const vector<bool>& polyhedronFaceNormalDirections,
                                                                                              const vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                                              const vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const
  {
    PointPolyhedronPositionResult result;

    unsigned int negativeFacePositions = 0;
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::MatrixXd& faceVertices3D = polyhedronFaceVertices[f];
      const Eigen::MatrixXd& faceVertices2D = polyhedronFaceRotatedVertices[f];
      const double faceOutgoingNormal = polyhedronFaceNormalDirections[f] ? 1.0 : -1.0;

      const PointPlanePositionTypes pointFacePlanePosition = PointPlanePosition(
                                                               PointPlaneDistance(point,
                                                                                  faceOutgoingNormal * polyhedronFaceNormals[f],
                                                                                  faceVertices3D.col(0)));

      switch (pointFacePlanePosition)
      {
        case PointPlanePositionTypes::Positive:
          continue;
        case PointPlanePositionTypes::Negative:
          negativeFacePositions++;
          continue;
        case PointPlanePositionTypes::Unknown:
          throw runtime_error("Not managed pointFacePlanePosition");
        default:
          break;
      }

      // Point is on face
      const Eigen::Vector3d point2D = RotatePointsFrom3DTo2D(point,
                                                             polyhedronFaceRotationMatrices[f].transpose(),
                                                             polyhedronFaceTranslations[f]);

      PointPolygonPositionResult pointFacePosition = PointPolygonPosition(point2D,
                                                                          faceVertices2D);

      switch (pointFacePosition.Type)
      {
        case PointPolygonPositionResult::Types::BorderVertex:
        {
          result.Type = PointPolyhedronPositionResult::Types::BorderVertex;
          result.BorderIndex = polyhedronFaces[f](0, pointFacePosition.BorderIndex);
          return result;
        }
        case PointPolygonPositionResult::Types::BorderEdge:
        {
          result.Type = PointPolyhedronPositionResult::Types::BorderEdge;
          result.BorderIndex = polyhedronFaces[f](1, pointFacePosition.BorderIndex);
          return result;
        }
        case PointPolygonPositionResult::Types::Inside:
        {
          result.Type = PointPolyhedronPositionResult::Types::BorderFace;
          result.BorderIndex = f;
          return result;
        }
        case PointPolygonPositionResult::Types::Unknown:
          throw runtime_error("Not managed pointFacePosition");
        default:
          continue;
      }
    }

    result.Type = (negativeFacePositions == polyhedronFaces.size()) ?
                    GeometryUtilities::PointPolyhedronPositionResult::Types::Inside :
                    GeometryUtilities::PointPolyhedronPositionResult::Types::Outside;
    return result;
  }
  // ***************************************************************************
  GeometryUtilities::PointCirclePositionResult GeometryUtilities::PointCirclePosition(const Eigen::Vector3d& point,
                                                                                      const Eigen::Vector3d& circleCenter,
                                                                                      const double& circleRadius) const
  {
    CompareTypes comparison = Compare1DValues(circleRadius, PointDistance(circleCenter, point));
    if (comparison == CompareTypes::FirstBeforeSecond)
      return PointCirclePositionResult::Outside;
    else if (comparison == CompareTypes::Coincident)
      return PointCirclePositionResult::OnBorder;
    else
      return PointCirclePositionResult::Inside;
  }
  // ***************************************************************************
  vector<GeometryUtilities::PointCirclePositionResult> GeometryUtilities::PointCirclePositions(const Eigen::MatrixXd& points,
                                                                                               const Eigen::Vector3d& circleCenter,
                                                                                               const double& circleRadius) const
  {
    const unsigned int numVertices = points.cols();
    vector<PointCirclePositionResult> positions(numVertices);

    for (unsigned int v = 0; v < numVertices; v++)
    {
      positions[v] = PointCirclePosition(points.col(v),
                                         circleCenter,
                                         circleRadius);
    }

    return positions;
  }
  // ***************************************************************************
  vector<bool> GeometryUtilities::PointsAreOnLine(const Eigen::MatrixXd& points,
                                                  const Eigen::Vector3d& lineOrigin,
                                                  const Eigen::Vector3d& lineTangent) const
  {
    Output::Assert(points.rows() == 3 && points.cols() > 0);
    const unsigned int numPoints = points.cols();

    const Eigen::Vector3d t = lineTangent.normalized();

    vector<bool> aligned(numPoints, false);
    for (unsigned int p = 0; p < numPoints; p++)
    {
      const Eigen::Vector3d s = (points.col(p) - lineOrigin).normalized();
      aligned[p] = IsValue1DZero(t.cross(s).norm());
    }
    return aligned;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::UnalignedPoints(const Eigen::MatrixXd& points,
                                                          const unsigned int numDesiredUnalignedPoints) const
  {
    Output::Assert(points.rows() == 3 && points.cols() > 1);
    const unsigned int& numPoints = points.cols();

    if (numPoints == 2)
      return { 0 , 1 };

    std::list<unsigned int> unalignedPoints;

    for (unsigned int p = 0; p < numPoints; p++)
    {
      const Eigen::Vector3d& segmentOrigin = points.col(p == 0 ? numPoints - 1 : p - 1);
      const Eigen::Vector3d& segmentEnd = points.col(p);
      const Eigen::Vector3d& nextPoint = points.col((p + 1) % numPoints);

      if (!PointIsAligned(segmentOrigin,
                          segmentEnd,
                          nextPoint))
        unalignedPoints.push_back(p);

      if (numDesiredUnalignedPoints > 0 &&
          unalignedPoints.size() == numDesiredUnalignedPoints)
        break;
    }

    return std::vector<unsigned int>(unalignedPoints.begin(), unalignedPoints.end());
  }
  // ***************************************************************************
}
