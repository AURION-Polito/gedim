#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::GeometryUtilities(const GeometryUtilitiesConfig& configuration) :
    _configuration(configuration)
  {
  }
  GeometryUtilities::~GeometryUtilities()
  {
  }
  // ***************************************************************************
  GeometryUtilities::CompareTypes GeometryUtilities::CompareValues(const double& first,
                                                                   const double& second,
                                                                   const double& tolerance) const
  {
    double maxFirstSecondOne = std::max( { 1.0, abs(first) , abs(second) } ) ;
    double difference = second - first;

    if (abs(difference) <= tolerance * maxFirstSecondOne)
      return CompareTypes::Coincident;
    else if (difference < -tolerance * maxFirstSecondOne)
      return CompareTypes::SecondBeforeFirst;
    else
      return CompareTypes::FirstBeforeSecond;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionSegmentSegmentResult GeometryUtilities::IntersectionSegmentSegment(const Vector3d& firstSegmentOrigin,
                                                                                                    const Vector3d& firstSegmentEnd,
                                                                                                    const Vector3d& secondSegmentOrigin,
                                                                                                    const Vector3d& secondSegmentEnd) const
  {
    GeometryUtilities::IntersectionSegmentSegmentResult result;

    // segments are x2->x1 and x4->x3
    // intersect two lines: r1 = s*t1+x1 and r2 = q*t2+x3
    // with t1 = x2 - x1 and t2 = x4 - x3

    const Vector3d t1 = firstSegmentEnd - firstSegmentOrigin;
    const Vector3d t2 = secondSegmentEnd - secondSegmentOrigin;

    // check if t1 and t2 are not zero
    Gedim::Output::Assert(IsValue2DPositive(t1.squaredNorm()));
    Gedim::Output::Assert(IsValue2DPositive(t2.squaredNorm()));

    // coplanarity check: (x3-x1).dot((x2-x1) x (x4-x1)) = det([x3-x1; x2-x1; x4-x1]) = 0
    // see https://en.wikipedia.org/wiki/Coplanarity
    Matrix3d coplanarMatrix;
    coplanarMatrix.col(0)<< (secondSegmentOrigin - firstSegmentOrigin).normalized();
    coplanarMatrix.col(1)<< t1.normalized();
    coplanarMatrix.col(2)<< (secondSegmentEnd - firstSegmentOrigin).normalized();

    // Check non-coplanarity of segments: (x3-x1) != 0 && (x4-x1) != 0 && (x3-x2) != 0 && (x4-x2) != 0 && nDotRhs != 0
    if (IsValue2DPositive((secondSegmentOrigin - firstSegmentOrigin).squaredNorm()) &&
        IsValue2DPositive((secondSegmentEnd - firstSegmentOrigin).squaredNorm()) &&
        IsValue2DPositive((secondSegmentOrigin - firstSegmentEnd).squaredNorm()) &&
        IsValue2DPositive((secondSegmentEnd - firstSegmentEnd).squaredNorm()) &&
        IsValue1DPositive(abs(coplanarMatrix.determinant())))
    {
      // segments are not on the same plane
      result.IntersectionLinesType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::OnDifferentPlanes;
      result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection;
      return result;
    }

    const double l1 = t1.norm();
    const double l2 = t2.norm();

    Matrix<double, 3, 2> matrixTangentVector;
    // t1 = x2 - x1
    matrixTangentVector.col(0) = t1.normalized();
    // t1 = x4 - x3
    matrixTangentVector.col(1) = t2.normalized();
    // rhs = x3 - x1
    Vector3d rightHandSide = secondSegmentOrigin - firstSegmentOrigin;

    // tangentsDot = t1.dot(t2) / (||t1|| * ||t2||)
    double tangentsDot = matrixTangentVector.col(0).dot(matrixTangentVector.col(1));
    // Check parallelism of segments
    // segments are not parallel if squared norm of cross product ||t1 x t2||^2 / (||t1||^2 * ||t2||^2) = 1.0 - (t1.dot(t2) / (||t1|| * ||t2||))^2 > 0
    // which means to check 1.0 - (t1.dot(t2) / (||t1|| * ||t2||))^2 > tol^2 => (1.0 - t1.dot(t2) / (||t1|| * ||t2||)) * (1.0 + t1.dot(t2) / (||t1|| * ||t2||)) > tol^2
    // NB: check on abs(t1.dot(t2) / (||t1|| * ||t2||)) != 1.0 is done to avoid numerical loss of significance
    if (IsValue1DPositive(abs(abs(tangentsDot) - 1.0)) &&
        IsValue2DPositive((1.0 - tangentsDot) * (1.0 + tangentsDot)))
    {
      // no parallel segments
      result.IntersectionLinesType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting;

      // intersecting lines, only one intersection
      result.SecondIntersectionRelation.resize(1);
      result.SecondIntersectionRelation[0] = 0;

      result.FirstSegmentIntersections.resize(1);
      result.SecondSegmentIntersections.resize(1);

      // Having r1 = s*t1+x1 and r2 = q*t2+x3 solve system (t1, t2) * (s, -q) = rsh to find intersection point
      Vector2d resultParametricCoordinates = matrixTangentVector.fullPivHouseholderQr().solve(rightHandSide);
      //Vector2d resultParametricCoordinates = matrixTangentVector.bdcSvd(ComputeFullU | ComputeFullV).solve(rightHandSide);
      result.FirstSegmentIntersections[0].CurvilinearCoordinate = resultParametricCoordinates[0] / l1;
      result.SecondSegmentIntersections[0].CurvilinearCoordinate = -resultParametricCoordinates[1] / l2;

      // Check intersection position
      result.FirstSegmentIntersections[0].Type = PointSegmentPosition(result.FirstSegmentIntersections[0].CurvilinearCoordinate);
      result.SecondSegmentIntersections[0].Type = PointSegmentPosition(result.SecondSegmentIntersections[0].CurvilinearCoordinate);

      if ((result.FirstSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentOrigin ||
           result.FirstSegmentIntersections[0].Type == PointSegmentPositionTypes::InsideSegment ||
           result.FirstSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentEnd) &&
          (result.SecondSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentOrigin ||
           result.SecondSegmentIntersections[0].Type == PointSegmentPositionTypes::InsideSegment ||
           result.SecondSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentEnd))
        result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection;
      else
        result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection;
    }
    else
    {
      // segments are parallel
      result.IntersectionLinesType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarParallel;

      double checkCurvilinearCoordinate = PointCurvilinearCoordinate(secondSegmentOrigin, firstSegmentOrigin, firstSegmentEnd);
      if (PointDistance(firstSegmentOrigin + checkCurvilinearCoordinate * t1,
                        secondSegmentOrigin) > _configuration.Tolerance * secondSegmentOrigin.norm())
      {
        // segments are parallel on different lines
        result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection;
        return result;
      }
      else
      {
        // segments are on the same line, check multiple intersections
        const double r1 = l1 * 0.5;
        const double r2 = l2 * 0.5;

        Vector3d centroid1 = 0.5 * (firstSegmentEnd + firstSegmentOrigin);
        Vector3d centroid2 = 0.5 * (secondSegmentEnd + secondSegmentOrigin);
        double distance = (centroid2 - centroid1).norm();

        // Check distance of spheres to exclude intersections
        if (IsValue1DPositive((distance - (r1 + r2))))
        {
          result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection;
          return result;
        }

        // There are multiple intersections
        result.IntersectionSegmentsType = GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections;
        result.SecondIntersectionRelation.resize(2);
        result.FirstSegmentIntersections.resize(2);
        result.SecondSegmentIntersections.resize(2);

        result.FirstSegmentIntersections[0].CurvilinearCoordinate = PointCurvilinearCoordinate(secondSegmentOrigin, firstSegmentOrigin, firstSegmentEnd);
        result.FirstSegmentIntersections[1].CurvilinearCoordinate = PointCurvilinearCoordinate(secondSegmentEnd, firstSegmentOrigin, firstSegmentEnd);

        result.SecondSegmentIntersections[0].CurvilinearCoordinate = PointCurvilinearCoordinate(firstSegmentOrigin, secondSegmentOrigin, secondSegmentEnd);
        result.SecondSegmentIntersections[1].CurvilinearCoordinate = PointCurvilinearCoordinate(firstSegmentEnd, secondSegmentOrigin, secondSegmentEnd);

        if (result.FirstSegmentIntersections[0].CurvilinearCoordinate > result.FirstSegmentIntersections[1].CurvilinearCoordinate)
        {
          double temp = result.FirstSegmentIntersections[0].CurvilinearCoordinate;
          result.FirstSegmentIntersections[0].CurvilinearCoordinate = result.FirstSegmentIntersections[1].CurvilinearCoordinate;
          result.FirstSegmentIntersections[1].CurvilinearCoordinate = temp;
        }

        if (result.SecondSegmentIntersections[0].CurvilinearCoordinate > result.SecondSegmentIntersections[1].CurvilinearCoordinate)
        {
          double temp = result.SecondSegmentIntersections[0].CurvilinearCoordinate;
          result.SecondSegmentIntersections[0].CurvilinearCoordinate = result.SecondSegmentIntersections[1].CurvilinearCoordinate;
          result.SecondSegmentIntersections[1].CurvilinearCoordinate = temp;
        }

        // Check intersection position
        result.FirstSegmentIntersections[0].Type = PointSegmentPosition(result.FirstSegmentIntersections[0].CurvilinearCoordinate);
        result.FirstSegmentIntersections[1].Type = PointSegmentPosition(result.FirstSegmentIntersections[1].CurvilinearCoordinate);

        if (result.FirstSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
        {
          result.FirstSegmentIntersections[0].CurvilinearCoordinate = 0.0;
          result.FirstSegmentIntersections[0].Type = PointSegmentPositionTypes::OnSegmentOrigin;
        }
        else if (result.FirstSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentLineAfterEnd)
        {
          result.FirstSegmentIntersections[0].CurvilinearCoordinate = 0.0;
          result.FirstSegmentIntersections[0].Type = PointSegmentPositionTypes::OnSegmentOrigin;
        }

        if (result.FirstSegmentIntersections[1].Type == PointSegmentPositionTypes::OnSegmentLineAfterEnd)
        {
          result.FirstSegmentIntersections[1].CurvilinearCoordinate = 1.0;
          result.FirstSegmentIntersections[1].Type = PointSegmentPositionTypes::OnSegmentEnd;
        }

        result.SecondSegmentIntersections[0].Type = PointSegmentPosition(result.SecondSegmentIntersections[0].CurvilinearCoordinate);
        result.SecondSegmentIntersections[1].Type = PointSegmentPosition(result.SecondSegmentIntersections[1].CurvilinearCoordinate);

        if (result.SecondSegmentIntersections[0].Type == PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
        {
          result.SecondSegmentIntersections[0].CurvilinearCoordinate = 0.0;
          result.SecondSegmentIntersections[0].Type = PointSegmentPositionTypes::OnSegmentOrigin;
        }

        if (result.SecondSegmentIntersections[1].Type == PointSegmentPositionTypes::OnSegmentLineAfterEnd)
        {
          result.SecondSegmentIntersections[1].CurvilinearCoordinate = 1.0;
          result.SecondSegmentIntersections[1].Type = PointSegmentPositionTypes::OnSegmentEnd;
        }

        // check relation between intersections
        result.SecondIntersectionRelation[0] = 0;
        result.SecondIntersectionRelation[1] = 1;

        Vector3d firstPoint = firstSegmentOrigin + result.FirstSegmentIntersections[0].CurvilinearCoordinate * t1;
        if (PointDistance(firstPoint,
                          secondSegmentOrigin + result.SecondSegmentIntersections[0].CurvilinearCoordinate * t2) >
            _configuration.Tolerance * firstPoint.norm())
        {
          result.SecondIntersectionRelation[0] = 1;
          result.SecondIntersectionRelation[1] = 0;
        }
      }
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionSegmentPlaneResult GeometryUtilities::IntersectionSegmentPlane(const Eigen::Vector3d& segmentOrigin,
                                                                                                const Eigen::Vector3d& segmentEnd,
                                                                                                const Eigen::Vector3d& planeNormal,
                                                                                                const Eigen::Vector3d& planeOrigin) const
  {
    GeometryUtilities::IntersectionSegmentPlaneResult result;

    const Vector3d t = segmentEnd - segmentOrigin;
    const Vector3d n = planeNormal.normalized();

    // check if t is not zero
    Gedim::Output::Assert(IsValue2DPositive(t.squaredNorm()));

    // check if the plane normal n is perpendicular to segment tangent t
    if (IsValue1DZero(n.dot(t.normalized())))
    {
      // compare if n * segmentOrigin = n * planeOrigin
      if (Compare1DValues(n.dot(segmentOrigin), n.dot(planeOrigin)) == GeometryUtilities::CompareTypes::Coincident)
      {
        // multiple intersection, the segment is coplanar to plane
        result.Type = GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections;
      }
      else
      {
        // no intersection, the segment is on a parallel plane
        result.Type = GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection;
      }
    }
    else
    {
      // plane and segment have a single intersection
      result.Type = GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection;
      result.SingleIntersection.CurvilinearCoordinate = n.dot(planeOrigin - segmentOrigin) / n.dot(t);
      result.SingleIntersection.Type = PointSegmentPosition(result.SingleIntersection.CurvilinearCoordinate);
    }

    return result;
  }
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

    // check point position out of line
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
  GeometryUtilities::SplitPolygonResult GeometryUtilities::SplitPolygon(const GeometryUtilities::SplitPolygonInput& input) const
  {
    GeometryUtilities::SplitPolygonResult result;

    // check if segment is on polygon vertices in contigous edges, no split needed
    if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex &&
        input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex)
    {
      const unsigned int& originIndex = input.Segment.Origin.Index;
      const unsigned int& endIndex = input.Segment.End.Index;

      if ((originIndex + 1) % input.NumberPolygonVertices == endIndex ||
          (endIndex + 1) % input.NumberPolygonVertices == originIndex)
      {
        // No split needed
        result.Type = SplitPolygonResult::Types::NoAction;
        return result;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // No split needed
          result.Type = SplitPolygonResult::Types::NoAction;
          return result;
        }
      }
    }

    // check if segment is on polygon vertex and polygon edge in contigous edges, only update needed
    if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex &&
        input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge)
    {
      const unsigned int& originIndex = input.Segment.Origin.Index;
      const unsigned int& endIndex = input.Segment.End.Index;

      if (originIndex == endIndex)
      {
        // Only update needed
        result.Type = SplitPolygonResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == endIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonResult::NewVertex> newVertexTypes;
        vector<unsigned int> newVertices(updatedPolygon.Vertices.begin(), updatedPolygon.Vertices.end());
        for (unsigned int v = 0; v < newVertices.size(); v++)
        {
          unsigned int origin = newVertices[v];
          unsigned int end = newVertices[(v + 1) % newVertices.size()];

          if (origin < input.NumberPolygonVertices &&
              end < input.NumberPolygonVertices)
          {
            // old edge
            updatedPolygon.Edges.push_back(origin);
            continue;
          }
          else
          {
            // new edge
            updatedPolygon.Edges.push_back(newEdgeNumber++);
            result.NewEdges.push_back(SplitPolygonResult::NewEdge());
            SplitPolygonResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin new end
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(end,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end new origin
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(origin,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else
              throw runtime_error("This case is not possibile because origin is vertex type");
          }
        }

        for (map<unsigned int, SplitPolygonResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonResult::NewVertex());
          SplitPolygonResult::NewVertex& newVertex = result.NewVertices.back();
          newVertex = it->second;
        }
        return result;
      }
    }

    // check if segment is on polygon vertex and polygon edge in contigous edges, only update needed
    if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge &&
        input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex)
    {
      const unsigned int& originIndex = input.Segment.Origin.Index;
      const unsigned int& endIndex = input.Segment.End.Index;

      if (originIndex == endIndex)
      {
        // Only update needed
        result.Type = SplitPolygonResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == originIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonResult::NewVertex> newVertexTypes;
        vector<unsigned int> newVertices(updatedPolygon.Vertices.begin(), updatedPolygon.Vertices.end());
        for (unsigned int v = 0; v < newVertices.size(); v++)
        {
          unsigned int origin = newVertices[v];
          unsigned int end = newVertices[(v + 1) % newVertices.size()];

          if (origin < input.NumberPolygonVertices &&
              end < input.NumberPolygonVertices)
          {
            // old edge
            updatedPolygon.Edges.push_back(origin);
            continue;
          }
          else
          {
            // new edge
            updatedPolygon.Edges.push_back(newEdgeNumber++);
            result.NewEdges.push_back(SplitPolygonResult::NewEdge());
            SplitPolygonResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(end,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;
              continue;

              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(origin,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;
              continue;
            }
            else
              throw runtime_error("This case is not possibile because end is vertex type");
          }
        }

        for (map<unsigned int, SplitPolygonResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonResult::NewVertex());
          SplitPolygonResult::NewVertex& newVertex = result.NewVertices.back();
          newVertex = it->second;
        }
        return result;
      }
    }

    // check if segment is on polygon edge in contigous edges, only update needed
    if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge &&
        input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge)
    {
      const unsigned int& originIndex = input.Segment.Origin.Index;
      const unsigned int& endIndex = input.Segment.End.Index;

      if (originIndex == endIndex)
      {
        // Only update needed
        result.Type = SplitPolygonResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];

        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == originIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);

          if (v == endIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices + 1);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonResult::NewVertex> newVertexTypes;
        vector<unsigned int> newVertices(updatedPolygon.Vertices.begin(), updatedPolygon.Vertices.end());
        for (unsigned int v = 0; v < newVertices.size(); v++)
        {
          unsigned int origin = newVertices[v];
          unsigned int end = newVertices[(v + 1) % newVertices.size()];

          if (origin < input.NumberPolygonVertices &&
              end < input.NumberPolygonVertices)
          {
            // old edge
            updatedPolygon.Edges.push_back(origin);
            continue;
          }
          else
          {
            // new edge
            updatedPolygon.Edges.push_back(newEdgeNumber++);
            result.NewEdges.push_back(SplitPolygonResult::NewEdge());
            SplitPolygonResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(end,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = (end == 4) ? SplitPolygonResult::NewVertex::Types::SegmentOrigin :
                                            SplitPolygonResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(origin,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = (origin == 4) ? SplitPolygonResult::NewVertex::Types::SegmentOrigin :
                                               SplitPolygonResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else
            {
              // new origin and new end
              newEdge.OldEdgeId = input.Segment.Origin.Index;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(origin,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newOrigin = newVertexTypes[origin];
              newOrigin.Type = (origin == 4) ? SplitPolygonResult::NewVertex::Types::SegmentOrigin :
                                               SplitPolygonResult::NewVertex::Types::SegmentEnd;

              newVertexTypes.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(end,
                                                                                     SplitPolygonResult::NewVertex()));

              SplitPolygonResult::NewVertex& newEnd = newVertexTypes[end];
              newEnd.Type = (end == 4) ? SplitPolygonResult::NewVertex::Types::SegmentOrigin :
                                         SplitPolygonResult::NewVertex::Types::SegmentEnd;
              continue;
            }
          }
        }

        for (map<unsigned int, SplitPolygonResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonResult::NewVertex());
          SplitPolygonResult::NewVertex& newVertex = result.NewVertices.back();
          newVertex = it->second;
        }
        return result;
      }
    }

    // segment is on not contigous edges and generates two new polygons
    list<SplitPolygonResult::NewPolygon> polygons;

    map<unsigned int, SplitPolygonResult::NewVertex> newVertices;
    map<unsigned int, SplitPolygonResult::NewEdge> newEdges;
    vector<bool> visited(input.NumberPolygonVertices, false);
    bool allVerticesVisited = false;
    unsigned int v = 0;
    do
    {
      // starting new polygon creation
      unsigned int startingVertex = v;
      polygons.push_back(SplitPolygonResult::NewPolygon());
      SplitPolygonResult::NewPolygon& newPolygon = polygons.back();

      do
      {
        visited[v] = true;

        if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex &&
            input.Segment.Origin.Index == v)
        {
          // origin segment is a vertex of polygon

          newPolygon.Vertices.push_back(input.Segment.Origin.Index);
          if (input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex)
          {
            // end segment is on vertex too, add no new vertices and one edges in polygon
            newPolygon.Vertices.push_back(input.Segment.End.Index);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.Segment.End.Index);

            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdge = newEdges[input.NumberPolygonVertices];
            newEdge.OriginId = input.Segment.Origin.Index;
            newEdge.EndId = input.Segment.End.Index;
            newEdge.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdge.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
          }
          else
          {
            // end segment is on an edge, add one new vertex and three edges in polygon
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 1);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices];
            newEnd.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;


            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.Segment.Origin.Index;
            newEdgeOne.EndId = input.NumberPolygonVertices;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.NumberPolygonVertices;
            newEdgeTwo.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.End.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.Segment.End.Index;
            newEdgeThree.EndId = input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.End.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
          }

          v = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
          continue;
        }
        else if (input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex &&
                 input.Segment.End.Index == v)
        {
          // end segment is a vertex of polygon

          newPolygon.Vertices.push_back(input.Segment.End.Index);
          if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex)
          {
            // origin segment is on vertex too, add no new vertices and one new edges in polygon
            newPolygon.Vertices.push_back(input.Segment.Origin.Index);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.Segment.Origin.Index);

            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdge = newEdges[input.NumberPolygonVertices];
            newEdge.OriginId = input.Segment.Origin.Index;
            newEdge.EndId = input.Segment.End.Index;
            newEdge.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdge.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
          }
          else
          {
            // origin segment is on edge, add one new vertex and three edges in polygon
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 2);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;

            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.Segment.End.Index;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
          }

          v = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
          continue;
        }

        // add current edge
        newPolygon.Vertices.push_back(v);

        // check if there is an intersection in edge
        if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge &&
            input.Segment.Origin.Index == v)
        {
          // origin segment is on a edge of polygon
          if (input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge)
          {
            // end of segment is on edge too, add two vertices and five edges
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Vertices.push_back(input.NumberPolygonVertices + 1);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 1);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 4);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;
            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices + 1,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices + 1];
            newEnd.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;


            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 3,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 4,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.NumberPolygonVertices + 1;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonResult::NewEdge& newEdgeFour = newEdges[input.NumberPolygonVertices + 3];
            newEdgeFour.OriginId = input.Segment.End.Index;
            newEdgeFour.EndId = input.NumberPolygonVertices + 1;
            newEdgeFour.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeFour.OldEdgeId = input.Segment.End.Index;
            newEdgeFour.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonResult::NewEdge& newEdgeFive = newEdges[input.NumberPolygonVertices + 4];
            newEdgeFive.OriginId = input.NumberPolygonVertices + 1;
            newEdgeFive.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeFive.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeFive.OldEdgeId = input.Segment.End.Index;
            newEdgeFive.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
          }
          else
          {
            // end is on polygon vertex, add one new vertex and three new edges
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Vertices.push_back(input.Segment.End.Index);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 1);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.Segment.End.Index);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;


            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.Segment.End.Index;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
          }

          v = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
          continue;
        }
        else if (input.Segment.End.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge &&
                 input.Segment.End.Index == v)
        {
          // end segment is on a edge of polygon

          if (input.Segment.Origin.Type == GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge)
          {
            // origin segment is also on edge polygon, add two new vertices and five new edges
            newPolygon.Vertices.push_back(input.NumberPolygonVertices + 1);
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 3);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 2);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonResult::NewVertex::Types::SegmentOrigin;
            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices + 1,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices + 1];
            newEnd.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;

            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 3,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 4,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.NumberPolygonVertices + 1;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeFour = newEdges[input.NumberPolygonVertices + 3];
            newEdgeFour.OriginId = input.Segment.End.Index;
            newEdgeFour.EndId = input.NumberPolygonVertices + 1;
            newEdgeFour.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeFour.OldEdgeId = input.Segment.End.Index;
            newEdgeFour.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonResult::NewEdge& newEdgeFive = newEdges[input.NumberPolygonVertices + 4];
            newEdgeFive.OriginId = input.NumberPolygonVertices + 1;
            newEdgeFive.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeFive.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeFive.OldEdgeId = input.Segment.End.Index;
            newEdgeFive.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
          }
          else
          {
            // origin is on polygon vertex, add one new vertex and three new edges
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Vertices.push_back(input.Segment.Origin.Index);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 2);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.Segment.Origin.Index);

            newVertices.insert(pair<unsigned int,SplitPolygonResult::NewVertex>(input.NumberPolygonVertices,
                                                                                SplitPolygonResult::NewVertex()));

            SplitPolygonResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices];
            newEnd.Type = SplitPolygonResult::NewVertex::Types::SegmentEnd;

            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                            SplitPolygonResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                            SplitPolygonResult::NewEdge()));
            SplitPolygonResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.Segment.Origin.Index;
            newEdgeOne.EndId = input.NumberPolygonVertices;
            newEdgeOne.Type = SplitPolygonResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.NumberPolygonVertices;
            newEdgeTwo.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.End.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.Segment.End.Index;
            newEdgeThree.EndId = input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.End.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
          }

          v = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
          continue;
        }

        // no intersection found, add current edge
        newPolygon.Edges.push_back(v);
        v = (v + 1) % input.NumberPolygonVertices;
      }
      while (v != startingVertex);

      // check new starting vertex
      v = input.NumberPolygonVertices + 1;
      for (unsigned int k = 0; k < input.NumberPolygonVertices; k++)
      {
        if (!visited[k])
        {
          v = k;
          break;
        }
      }

      allVerticesVisited = (v == input.NumberPolygonVertices + 1);
    }
    while (!allVerticesVisited);

    result.Type = SplitPolygonResult::Types::PolygonCreation;
    for (map<unsigned int, SplitPolygonResult::NewVertex>::const_iterator it = newVertices.begin();
         it != newVertices.end(); it++)
    {
      result.NewVertices.push_back(SplitPolygonResult::NewVertex());
      SplitPolygonResult::NewVertex& newVertex = result.NewVertices.back();
      newVertex = it->second;
    }
    for (map<unsigned int, SplitPolygonResult::NewEdge>::const_iterator it = newEdges.begin();
         it != newEdges.end(); it++)
    {
      result.NewEdges.push_back(SplitPolygonResult::NewEdge());
      SplitPolygonResult::NewEdge& newEdge = result.NewEdges.back();
      newEdge = it->second;
    }

    result.NewPolygons.resize(polygons.size());
    copy(polygons.begin(), polygons.end(), result.NewPolygons.begin());

    return result;
  }
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
  Matrix3d GeometryUtilities::PlaneRotation(const Eigen::Vector3d& planeNormal) const
  {
    Matrix3d Q;
    Q.setIdentity();

    Output::Assert(Compare1DValues(planeNormal.norm(), 1.0) == CompareTypes::Coincident);

    // if planeNormal is already oriented as z-axis the return the identity
    const double n_xy = sqrt(planeNormal.x() * planeNormal.x() + planeNormal.y() * planeNormal.y());
    if (IsValue1DZero(n_xy))
      return Q;

    Q.col(0)<< -planeNormal.y() / n_xy, planeNormal.x() / n_xy, 0.0;
    Q.col(1)<< -planeNormal.x() * planeNormal.z() / n_xy, -planeNormal.y() * planeNormal.z() / n_xy, n_xy;
    Q.col(2)<< planeNormal;

    return Q;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::ConvexHull(const Eigen::MatrixXd& points)
  {
    // pseudocode
    // // S is the set of points
    // // P will be the set of points which form the convex hull. Final set size is i.
    // pointOnHull = leftmost point in S // which is guaranteed to be part of the CH(S)
    // i := 0
    // repeat
    //     P[i] := pointOnHull
    //     endpoint := S[0]      // initial endpoint for a candidate edge on the hull
    //     for j from 0 to |S| do
    //         // endpoint == pointOnHull is a rare case and can happen only when j == 1 and a better endpoint has not yet been set for the loop
    //         if (endpoint == pointOnHull) or (S[j] is on right of line from P[i] to endpoint) then
    //             endpoint := S[j]   // found greater right turn, update endpoint
    //     i := i + 1
    //     pointOnHull = endpoint
    // until endpoint = P[0]      // wrapped around to first hull point

    Output::Assert(points.rows() == 3 && points.cols() > 0);
    const unsigned int numPoints = points.cols();
    unsigned int leftMost = 0;
    for (unsigned int p = 1; p < numPoints; p++)
    {
      if (points(0, p) < points(0, leftMost))
        leftMost = p;
    }
    list<unsigned int> convexHull;

    unsigned int pointOnHull = leftMost;
    do
    {
      convexHull.push_back(pointOnHull);
      unsigned int endpoint = 0;
      for (unsigned int j = 0; j < numPoints; j++)
      {
        if ((endpoint == pointOnHull) ||
            PointSegmentPosition(points.col(j),
                                 points.col(pointOnHull),
                                 points.col(endpoint)) == GeometryUtilities::PointSegmentPositionTypes::RightTheSegment)
          endpoint = j;
      }
      pointOnHull = endpoint;
    }
    while (pointOnHull != convexHull.front());

    return vector<unsigned int>(convexHull.begin(), convexHull.end());
  }
  // ***************************************************************************
}
