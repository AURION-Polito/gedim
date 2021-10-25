#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::PointSegmentPositionTypes GeometryUtilities::PointSegmentPosition(const Vector3d& point,
                                                                                       const Vector3d& segmentOrigin,
                                                                                       const Vector3d& segmentEnd) const
  {
    Vector3d segmentTangent = (segmentEnd - segmentOrigin).normalized();
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

    if (Compare1DValues(curvilinearCoordinate, 1.0) == CompareTypes::Coincident)
      return PointSegmentPositionTypes::OnSegmentEnd;

    if (IsValue1DPositive(curvilinearCoordinate) &&
        Compare1DValues(curvilinearCoordinate, 1.0) == CompareTypes::FirstBeforeSecond)
      return PointSegmentPositionTypes::InsideSegment;

    if (IsValue1DNegative(curvilinearCoordinate))
      return PointSegmentPositionTypes::OnSegmentLineBeforeOrigin;

    if (Compare1DValues(curvilinearCoordinate, 1.0) == CompareTypes::SecondBeforeFirst)
      return PointSegmentPositionTypes::OnSegmentLineAfterEnd;

    return PointSegmentPositionTypes::Unknown;
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
        result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::Outside;
        return result;
      }
      else if (resultSegment != GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment)
      {
        if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
        {
          result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderEdge;
          result.BorderIndex = v;
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderVertex;
          result.BorderIndex = v;
        }
        else if (resultSegment == GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderVertex;
          result.BorderIndex = (v + 1) % numVertices;
        }
        else
          result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::Outside;

        return result;
      }
    }

    result.PositionType = GeometryUtilities::PointPolygonPositionResult::PositionTypes::Inside;
    return result;
  }
  // ***************************************************************************
}
