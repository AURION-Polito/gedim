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

    PointSegmentPositionTypes result;

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
  GeometryUtilities::PointPlanePositionTypes GeometryUtilities::PointPlanePosition(const Eigen::Vector3d& point,
                                                                                   const Eigen::Vector3d& planeNormal,
                                                                                   const Eigen::Vector3d& planeOrigin) const
  {

    const double position = planeNormal.dot(point - planeOrigin) / (point - planeOrigin).norm();

    if (IsValue1DZero(position))
      return PointPlanePositionTypes::OnPlane;
    else if (IsValue1DNegative(position))
      return PointPlanePositionTypes::Negative;
    else
      return PointPlanePositionTypes::Positive;
  }
  // ***************************************************************************
  GeometryUtilities::PointPolygonPositionResult GeometryUtilities::PointPolygonPosition(const Vector3d& point,
                                                                                        const MatrixXd& polygonVertices) const
  {
    Output::Assert(point.rows() == 3 && point.cols() > 0 && PointsAre2D(point));

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
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
          result.BorderIndex = v;
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex;
          result.BorderIndex = (v + 1) % numVertices;
        }
        else
          result.Type = GeometryUtilities::PointPolygonPositionResult::Types::Outside;

        return result;
      }
    }

    result.Type = GeometryUtilities::PointPolygonPositionResult::Types::Inside;
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
  vector<bool> GeometryUtilities::PointsAreAligned(const Eigen::Vector3d& segmentOrigin,
                                                   const Eigen::Vector3d& segmentEnd,
                                                   const Eigen::MatrixXd& points) const
  {
    Output::Assert(points.rows() == 3 && points.cols() > 0);
    const unsigned int numPoints = points.cols();

    const Eigen::Vector3d t = SegmentTangent(segmentOrigin, segmentEnd).normalized();

    vector<bool> aligned(numPoints, false);
    for (unsigned int p = 0; p < numPoints; p++)
    {
      const Eigen::Vector3d s = (points.col(p) - segmentOrigin).normalized();
      aligned[p] = IsValue1DZero(t.cross(s).norm());
    }
    return aligned;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::UnalignedPoints(const Eigen::MatrixXd& points) const
  {
    Output::Assert(points.rows() == 3 && points.cols() > 1);
    const unsigned int& numPoints = points.cols();

    list<unsigned int> unalignedPoints;

    unalignedPoints.push_back(0);
    unalignedPoints.push_back(1);
    for (unsigned int p = 0; p < numPoints - 2; p++)
    {
      const unsigned int pointIndexToTest = (p + 2) % numPoints;
      const Eigen::Vector3d& segmentOrigin = points.col(p);
      const Eigen::Vector3d& segmentEnd = points.col((p + 1) % numPoints);
      const Eigen::Vector3d& nextPoint = points.col((p + 2) % numPoints);

      if (PointsAreAligned(segmentOrigin,
                           segmentEnd,
                           nextPoint)[0])
        unalignedPoints.pop_back();

      unalignedPoints.push_back(pointIndexToTest);
    }

    return vector<unsigned int>(unalignedPoints.begin(), unalignedPoints.end());
  }
  // ***************************************************************************
}
