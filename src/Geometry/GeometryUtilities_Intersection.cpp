#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
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

      // check if are on the same line checking if (x2->x1)x(x3->x1)=(x2->x1)x(x4->x1)=0.0
      double checkOne = t1.cross(secondSegmentOrigin - firstSegmentOrigin).norm();
      double checkTwo = t1.cross(secondSegmentEnd - firstSegmentOrigin).norm();

      if (!(IsValue1DZero(checkOne) && IsValue1DZero(checkTwo)))
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

    const Vector3d t = SegmentTangent(segmentOrigin, segmentEnd);

    // check if t is not zero and plane normal is normalized
    Gedim::Output::Assert(Compare1DValues(planeNormal.norm(), 1.0) == CompareTypes::Coincident);
    Gedim::Output::Assert(IsValue2DPositive(t.squaredNorm()));

    // check if the plane normal n is perpendicular to segment tangent t
    if (IsValue1DZero(planeNormal.dot(t.normalized())))
    {
      // compare if n * segmentOrigin = n * planeOrigin
      if (Compare1DValues(planeNormal.dot(segmentOrigin),
                          planeNormal.dot(planeOrigin)) == GeometryUtilities::CompareTypes::Coincident)
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
      result.SingleIntersection.CurvilinearCoordinate = planeNormal.dot(planeOrigin - segmentOrigin) /
                                                        planeNormal.dot(t);
      result.SingleIntersection.Type = PointSegmentPosition(result.SingleIntersection.CurvilinearCoordinate);
      if (result.SingleIntersection.Type == PointSegmentPositionTypes::OnSegmentOrigin)
        result.SingleIntersection.CurvilinearCoordinate = 0.0;
      else if (result.SingleIntersection.Type == PointSegmentPositionTypes::OnSegmentEnd)
        result.SingleIntersection.CurvilinearCoordinate = 1.0;
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionPolyhedronPlaneResult GeometryUtilities::IntersectionPolyhedronPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                                                                      const Eigen::MatrixXi& polyhedronEdges,
                                                                                                      const vector<Eigen::MatrixXi> polyhedronFaces,
                                                                                                      const Eigen::Vector3d& planeNormal,
                                                                                                      const Eigen::Vector3d& planeOrigin,
                                                                                                      const Eigen::Matrix3d& planeRotationMatrix,
                                                                                                      const Eigen::Vector3d& planeTranslation) const
  {
    GeometryUtilities::IntersectionPolyhedronPlaneResult result;

    // check if plane normal is normalized
    Gedim::Output::Assert(Compare1DValues(planeNormal.norm(), 1.0) == CompareTypes::Coincident);

    unsigned int numberOfIntersections = 0;
    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    const unsigned int numPolyhedronEdges = polyhedronEdges.cols();
    const unsigned int numPolyhedronFaces = polyhedronFaces.size();

    list<GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection> intersectionsList;
    list<Eigen::Vector3d> intersectionCoordinates;

    result.VertexIntersections.resize(numPolyhedronVertices);
    result.EdgeIntersections.resize(numPolyhedronEdges);
    result.FaceIntersections.resize(numPolyhedronFaces);

    for (auto& vertexIntersection : result.VertexIntersections)
      vertexIntersection.Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection;
    for (auto& faceIntersection : result.FaceIntersections)
      faceIntersection.Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection;

    for (unsigned int e = 0; e < numPolyhedronEdges; e++)
    {
      const unsigned int edgeOriginId = polyhedronEdges(0, e);
      const unsigned int edgeEndId = polyhedronEdges(1, e);

      const Vector3d edgeOrigin = polyhedronVertices.col(edgeOriginId);
      const Vector3d edgeEnd = polyhedronVertices.col(edgeEndId);

      result.EdgeIntersections[e].Intersection =  GeometryUtilities::IntersectionSegmentPlane(edgeOrigin,
                                                                                              edgeEnd,
                                                                                              planeNormal,
                                                                                              planeOrigin);
      const IntersectionSegmentPlaneResult& intersectionEdge = result.EdgeIntersections[e].Intersection;

      if (intersectionEdge.Type == Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection)
        continue;

      if (intersectionEdge.Type == Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections)
      {
        // edge intersection
        if (result.VertexIntersections[edgeOriginId].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
        {
          result.VertexIntersections[edgeOriginId].Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection;
          numberOfIntersections++;
        }

        if (result.VertexIntersections[edgeEndId].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
        {
          result.VertexIntersections[edgeEndId].Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection;
          numberOfIntersections++;
        }

        continue;
      }

      switch (intersectionEdge.Type)
      {
        case Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection:
        {
          switch (intersectionEdge.SingleIntersection.Type)
          {
            case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
            {
              // edge origin intersection (vertex)
              if (result.VertexIntersections[edgeOriginId].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
              {
                result.VertexIntersections[edgeOriginId].Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection;
                numberOfIntersections++;

                intersectionsList.push_back(GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection());
                intersectionCoordinates.push_back(polyhedronVertices.col(edgeOriginId));
                GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection& intersection = intersectionsList.back();
                intersection.Type = GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Vertex;
                intersection.VertexId = edgeOriginId;
              }
            }
              break;
            case Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
            {
              // inside edge intersection
              numberOfIntersections++;

              intersectionsList.push_back(GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection());
              intersectionCoordinates.push_back(edgeOrigin + intersectionEdge.SingleIntersection.CurvilinearCoordinate * (edgeEnd - edgeOrigin));
              GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection& intersection = intersectionsList.back();
              intersection.Type = GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge;
              intersection.EdgeId = e;
            }
              break;
            case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd:
            {
              // edge end intersection (vertex)
              if (result.VertexIntersections[edgeEndId].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
              {
                result.VertexIntersections[edgeEndId].Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection;
                numberOfIntersections++;

                intersectionsList.push_back(GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection());
                intersectionCoordinates.push_back(polyhedronVertices.col(edgeEndId));
                GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection& intersection = intersectionsList.back();
                intersection.Type = GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Vertex;
                intersection.VertexId = edgeEndId;
              }
            }
              break;
            default:
              continue;
          }
        }
          break;
        default:
          throw runtime_error("Unknwon intersection edge type");
      }
    }

    switch (numberOfIntersections)
    {
      case 0:
      {
        // no intersection found
        result.Type = IntersectionPolyhedronPlaneResult::Types::None;
      }
        break;
      case 1:
      {
        // one intersection found, single vertex intersection
        result.Type = IntersectionPolyhedronPlaneResult::Types::OnVertex;
        for (unsigned int v = 0; v < numPolyhedronVertices; v++)
        {
          if (result.VertexIntersections[v].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
            continue;

          result.IntersectionId = v;
        }
      }
        break;
      case 2:
      {
        // two intersections found, edge intersection
        result.Type = IntersectionPolyhedronPlaneResult::Types::OnEdge;
        for (unsigned int e = 0; e < numPolyhedronEdges; e++)
        {
          if (result.EdgeIntersections[e].Intersection.Type != IntersectionSegmentPlaneResult::Types::MultipleIntersections)
            continue;

          result.IntersectionId = e;
        }
      }
        break;
      default:
      {
        // more than two intersections found, of polyhedron face intersection, or new polygon intersection
        // check intersection on face
        int faceIntersection = -1;
        for (unsigned int f = 0; f < numPolyhedronFaces; f++)
        {
          bool faceVerticesIntersection = true;
          for (unsigned int v = 0; v < polyhedronFaces[f].cols(); v++)
          {
            if (result.VertexIntersections[polyhedronFaces[f](0, v)].Type != Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection)
            {
              faceVerticesIntersection = false;
              break;
            }
          }

          if (!faceVerticesIntersection)
            continue;

          bool faceEdgesIntersection = true;
          for (unsigned int e = 0; e < polyhedronFaces[f].cols(); e++)
          {
            if (result.EdgeIntersections[polyhedronFaces[f](1, e)].Intersection.Type != Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections)
            {
              faceEdgesIntersection = false;
              break;
            }
          }

          if (!faceEdgesIntersection)
            continue;

          result.FaceIntersections[f].Type = Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::Intersection;
          faceIntersection = f;
          break;
        }

        if (faceIntersection >= 0)
        {
          // face intersection
          result.Type = IntersectionPolyhedronPlaneResult::Types::OnFace;
          result.IntersectionId = faceIntersection;
        }
        else
        {
          // inside polyhedron intersection
          result.Type = IntersectionPolyhedronPlaneResult::Types::NewPolygon;

          // create new polygon
          const unsigned int numIntersions = intersectionCoordinates.size();
          Eigen::MatrixXd convexHull3DPoints(3, numIntersions);
          unsigned int numIntersection = 0;
          for (const auto& intersection : intersectionCoordinates)
            convexHull3DPoints.col(numIntersection++)<< intersection;

          Eigen::MatrixXd convexHull2DPoints = RotatePointsFrom3DTo2D(convexHull3DPoints,
                                                                      planeRotationMatrix,
                                                                      planeTranslation);

          vector<unsigned int> convexHull = ConvexHull(convexHull2DPoints);
          Output::Assert(convexHull.size() == numIntersions);
          vector<GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection> intersections(intersectionsList.begin(), intersectionsList.end());

          result.Intersections.resize(numIntersions);
          result.IntersectionCoordinates.resize(3, numIntersions);

          for (unsigned int c = 0; c < numIntersions; c++)
          {
            result.Intersections[c] = intersections[convexHull[c]];
            result.IntersectionCoordinates.col(c) << convexHull3DPoints.col(convexHull[c]);
          }
        }
      }
        break;
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionPolyhedronLineResult GeometryUtilities::IntersectionPolyhedronLine(const Eigen::MatrixXd& polyhedronVertices,
                                                                                                    const Eigen::MatrixXi& polyhedronEdges,
                                                                                                    const vector<Eigen::MatrixXi> polyhedronFaces,
                                                                                                    const Eigen::Vector3d& lineTangent,
                                                                                                    const Eigen::Vector3d& lineOrigin) const
  {
    IntersectionPolyhedronLineResult result;

    unsigned int numberOfIntersections = 0;
    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    const unsigned int numPolyhedronEdges = polyhedronEdges.cols();
    const unsigned int numPolyhedronFaces = polyhedronFaces.size();

    list<GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection> intersectionsList;
    list<Eigen::Vector3d> intersectionCoordinates;

    result.PolyhedronVertexIntersections.resize(numPolyhedronVertices);
    result.PolyhedronEdgeIntersections.resize(numPolyhedronEdges);
    result.PolyhedronFaceIntersections.resize(numPolyhedronFaces);
    result.LineIntersections.resize(2);

    list<int> vertIntersection; //contiene l'indice dei vertici su cui c'è intersezione
    list<int> edgeIntersection;
    list<int> faceIntersection;

    bool flag = false;

    for (auto& vertexIntersection : result.PolyhedronVertexIntersections)
          vertexIntersection.Type = Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronVertexIntersection::Types::NoIntersection;
    for (auto& faceIntersection : result.PolyhedronFaceIntersections)
          faceIntersection.Type = Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::NoIntersection;

    //scorro i vertici e controllo se ci sono intersezioni
    for (unsigned int i=0; i<numPolyhedronVertices; i++)
    {
        //controllo se il vertice sta nella retta e metto flag = true
        //casi due componenti della tangente nulle
        if (lineTangent(1)==0 && lineTangent(2)==0 &&
            polyhedronVertices(1,i)==lineOrigin(1) && polyhedronVertices(2,i)==lineOrigin(2))
        {
            flag = true;
        }
        else if (lineTangent(0)==0 && lineTangent(1)==0 &&
            polyhedronVertices(1,i)==lineOrigin(1) && polyhedronVertices(0,i)==lineOrigin(0))
        {
            flag = true;
        }
        else if (lineTangent(0)==0 && lineTangent(2)==0 &&
            polyhedronVertices(0,i)==lineOrigin(0) && polyhedronVertices(2,i)==lineOrigin(2))
        {
            flag = true;
        }
        // una componente nulla
        else if (lineTangent(0)==0 && polyhedronVertices(0,i)==lineOrigin(0) &&
                 lineTangent(2)*(polyhedronVertices(1,i)-lineOrigin(1))==lineTangent(1)*(polyhedronVertices(2,i)-lineOrigin(2)))
        {
            flag = true;
        }
        else if (lineTangent(1)==0 && polyhedronVertices(1,i)==lineOrigin(1) &&
                 lineTangent(2)*(polyhedronVertices(0,i)-lineOrigin(0))==lineTangent(0)*(polyhedronVertices(2,i)-lineOrigin(2)))
        {
            flag = true;
        }
        else if (lineTangent(2)==0 && polyhedronVertices(2,i)==lineOrigin(2) &&
                 lineTangent(1)*(polyhedronVertices(0,i)-lineOrigin(0))==lineTangent(0)*(polyhedronVertices(1,i)-lineOrigin(1)))
        {
            flag = true;
        }
        // se sono in uno dei casi precedenti oppure caso nessuna componente nulla
        if (flag == true || (lineTangent(0)!=0 && lineTangent(1)!=0 && lineTangent(2)!=0))
        {
            if (flag == true || (lineTangent(1)*(polyhedronVertices(0,i)-lineOrigin(0))==lineTangent(0)*(polyhedronVertices(1,i)-lineOrigin(1))
                 && lineTangent(2)*(polyhedronVertices(1,i)-lineOrigin(1))==lineTangent(1)*(polyhedronVertices(2,i)-lineOrigin(2))))
             {
                 // aggiorno la lista dei vertici segnando che nel vertice che sto considerando c'è intersezione e scrivo il numero dell'intersezione
                 result.PolyhedronVertexIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronVertexIntersection::Types::Intersection;
                 result.PolyhedronVertexIntersections[i].LineIntersectionIndex = numberOfIntersections;

                 // aggiorno le informazioni dal punto di vista della retta
                 result.LineIntersections[numberOfIntersections].CurvilinearCoordinate = 0.0;
                 result.LineIntersections[numberOfIntersections].PolyhedronType = GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnVertex;
                 result.LineIntersections[numberOfIntersections].PolyhedronIndex = i;

                 vertIntersection.push_back(i);

                 numberOfIntersections++;
                 flag = false;
             }
             else
             {
                 result.PolyhedronVertexIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronVertexIntersection::Types::NoIntersection;
                 result.PolyhedronVertexIntersections[i].LineIntersectionIndex = 0;
             }
        }
    }

    //scorro i lati e controllo se ci sono intersezioni
    //distanza tra tutti i vertici e l'origine retta, prendo la più grande
    Eigen::VectorXd distance = GeometryUtilities::PointDistances(polyhedronVertices, lineOrigin);
    double max = distance(0);
    for (unsigned int i=1; i<numPolyhedronVertices; i++)
    {
         if (distance(i)>max)
             max = distance(i);
    }

    //calcolo l'estremo della retta, essendo sicura che è fuori dal poliedro
    const Eigen::Vector3d& s2 = lineOrigin + max*lineTangent;
    // per ogni lato salvo origine e fine
    for (unsigned int i=0; i<numPolyhedronEdges; i++)
    {
         const Eigen::Vector3d& edgeOrigin = polyhedronVertices.col(polyhedronEdges(0,i));
         const Eigen::Vector3d& edgeEnd = polyhedronVertices.col(polyhedronEdges(1,i));

         GeometryUtilitiesConfig geometryUtilityConfig;
         GeometryUtilities geometryUtility(geometryUtilityConfig);
         GeometryUtilities::IntersectionSegmentSegmentResult r = geometryUtility.IntersectionSegmentSegment(edgeOrigin, edgeEnd,
                                                                                                            lineOrigin, s2);
         // se ho un'intersezione
         if (r.IntersectionSegmentsType == GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
         {
             // se l'intersezione è interna al segmento (se è sui vertici l'ho già considerata)
             if (r.SecondSegmentIntersections[0].Type == GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
             {
                 double c = r.SecondSegmentIntersections[0].CurvilinearCoordinate;
                 result.LineIntersections[numberOfIntersections].CurvilinearCoordinate = c;
                 result.LineIntersections[numberOfIntersections].PolyhedronType = GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnEdge;
                 result.LineIntersections[numberOfIntersections].PolyhedronIndex = i;

                 result.PolyhedronEdgeIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronEdgeIntersection::Types::Intersection;
                 result.PolyhedronEdgeIntersections[i].LineIntersectionIndex = numberOfIntersections;

                 edgeIntersection.push_back(i);
                 numberOfIntersections++;
             }
             // se non è all'interno del segmento
             else
             {
                 result.PolyhedronEdgeIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronEdgeIntersection::Types::NoIntersection;
                 result.PolyhedronEdgeIntersections[i].LineIntersectionIndex = 0;
             }
         }
         // se non è intersezione singola
         else
         {
             result.PolyhedronEdgeIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronEdgeIntersection::Types::NoIntersection;
             result.PolyhedronEdgeIntersections[i].LineIntersectionIndex = 0;
         }
    }

    //scorro le facce e controllo se ci sono intersezioni
    for (unsigned int i=0; i<numPolyhedronFaces; i++)
    {
         //prendo il primo vertice della faccia come origine del piano e salvo gli altri indici
         const int a = polyhedronFaces[i](0,0); //contiene indice 0
         const int b = polyhedronFaces[i](0,1); //1
         const int c = polyhedronFaces[i](0,2);
         const int d = polyhedronFaces[i](0,3);
         const Eigen::Vector3d& planeOrigin = polyhedronVertices.col(a);
         Eigen::Vector3d x(+1.0, +0.0, +0.0);
         Eigen::Vector3d y(+0.0, +1.0, +0.0);
         Eigen::Vector3d z(+0.0, +0.0, +1.0);
         Eigen::Vector3d planeNormal;

         //calcolo la normale alla faccia
         //se i vertici della faccia hanno stessa coordinata x
         if (polyhedronVertices(0,a) == polyhedronVertices(0,b) && polyhedronVertices(0,a) == polyhedronVertices(0,c) && polyhedronVertices(0,a) == polyhedronVertices(0,d))
             planeNormal = x;
         //se i vertici della faccia hanno stessa coordinata y
         else if (polyhedronVertices(1,a) == polyhedronVertices(1,b) && polyhedronVertices(1,a) == polyhedronVertices(1,c) && polyhedronVertices(1,a) == polyhedronVertices(1,d))
                  planeNormal = y;
         //se i vertici della faccia hanno stessa coordinata z
         else if (polyhedronVertices(2,a) == polyhedronVertices(2,b) && polyhedronVertices(2,a) == polyhedronVertices(2,c) && polyhedronVertices(2,a) == polyhedronVertices(2,d))
                  planeNormal = z;

         GeometryUtilitiesConfig geometryUtilityConfig;
         GeometryUtilities geometryUtility(geometryUtilityConfig);
         GeometryUtilities::IntersectionSegmentPlaneResult r = geometryUtility.IntersectionSegmentPlane(lineOrigin, s2,
                                                                                                        planeNormal, planeOrigin);
         // se c'è intersezione con il piano contenente la faccia
         if (r.Type == GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection)
         {
             // salvo la coord curvilinea
             double coord_curv = r.SingleIntersection.CurvilinearCoordinate;
             // calcolo le coordinate di intersezione
             Eigen::Vector3d inters;
             inters = lineOrigin + coord_curv*(s2-lineOrigin);
             //Eigen::MatrixXd polygonVertices(3,4);
             //polygonVertices.col(0) = polyhedronVertices.col(a);
             //polygonVertices.col(1) = polyhedronVertices.col(b);
             //polygonVertices.col(2) = polyhedronVertices.col(c);
             //polygonVertices.col(3) = polyhedronVertices.col(d);

             //GeometryUtilities::PointPolygonPositionResult s = geometryUtility.PointPolygonPosition(inters,
             //                                                                                       polygonVertices);

             //if (s.Type == GeometryUtilities::PointPolygonPositionResult::Types::Inside)


             // controllo se l'intersezione è interna alla faccia in base alla posizione della faccia (in base alla normale)
             flag = false; //se è interna diventa true
             if (planeNormal==x)
             {
                 if (inters(1)>polyhedronVertices(1,0)&&inters(1)<polyhedronVertices(1,3) &&
                     inters(2)>polyhedronVertices(2,0)&&inters(2)<polyhedronVertices(2,4))
                     flag = true;
             }
             else if (planeNormal==y)
             {
                 if (inters(0)>polyhedronVertices(0,0)&&inters(0)<polyhedronVertices(0,1) &&
                     inters(2)>polyhedronVertices(2,0)&&inters(2)<polyhedronVertices(2,4))
                     flag = true;
             }
             else if (planeNormal==z)
             {
                 if (inters(0)>polyhedronVertices(0,0)&&inters(0)<polyhedronVertices(0,1) &&
                     inters(1)>polyhedronVertices(1,0)&&inters(1)<polyhedronVertices(1,3))
                     flag = true;
             }
             // se è interna alla faccia, aggiungo le informazioni relative all'intersezione
             if (flag==true)
             {
                 result.LineIntersections[numberOfIntersections].CurvilinearCoordinate = coord_curv;
                 result.LineIntersections[numberOfIntersections].PolyhedronType = GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnFace;
                 result.LineIntersections[numberOfIntersections].PolyhedronIndex = i;

                 result.PolyhedronFaceIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::Intersection;
                 result.PolyhedronFaceIntersections[i].LineIntersectionIndex = numberOfIntersections;

                 faceIntersection.push_back(i);
                 numberOfIntersections++;
             }
             else
             {
                 result.PolyhedronFaceIntersections[i].Type = GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::NoIntersection;
                 result.PolyhedronFaceIntersections[i].LineIntersectionIndex = 0;
             }
         }
     }

    // in base al numero di intersezioni trovate, aggiorno il tipo
    if (numberOfIntersections==0)
        result.Type = GeometryUtilities::IntersectionPolyhedronLineResult::Types::None;
    else if (numberOfIntersections==1)
        result.Type = GeometryUtilities::IntersectionPolyhedronLineResult::Types::OneIntersection;
    else if (numberOfIntersections==2)
        result.Type = GeometryUtilities::IntersectionPolyhedronLineResult::Types::TwoIntersections;

    result.LineIntersections.resize(numberOfIntersections);

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionSegmentCircleResult GeometryUtilities::IntersectionSegmentCircle(const Eigen::Vector3d& segmentOrigin,
                                                                                                  const Eigen::Vector3d& segmentEnd,
                                                                                                  const Eigen::Vector3d& circleCenter,
                                                                                                  const double& circleRadius) const
  {
    GeometryUtilities::IntersectionSegmentCircleResult result;

    const Vector3d d = SegmentTangent(segmentOrigin, segmentEnd);
    Vector3d f = segmentOrigin - circleCenter;

    double a = d.dot(d);
    double b = 2.0 * f.dot(d) ;
    double c = f.dot(f) - circleRadius*circleRadius;

    Output::Assert(IsValue2DPositive(a));

    double discriminant = b * b - 4.0 * a * c;
    if (IsValue2DNegative(discriminant))
    {
      // no intersection found
      result.Type = GeometryUtilities::IntersectionSegmentCircleResult::Types::NoIntersection;
    }
    else if (IsValue2DZero(discriminant))
    {
      // one intersection found
      double intersection = -b / (2.0 * a);
      result.Type = GeometryUtilities::IntersectionSegmentCircleResult::Types::TangentIntersection;
      result.SegmentIntersections.resize(1);
      result.SegmentIntersections[0].CurvilinearCoordinate = intersection;
      result.SegmentIntersections[0].Type = PointSegmentPosition(intersection);
    }
    else
    {
      // two intersections found
      discriminant = sqrt(discriminant);

      // either solution may be on or off the ray so need to test both
      // t1 is always the smaller value, because BOTH discriminant and
      // a are nonnegative.
      double t1 = (-b - discriminant)/(2.0 * a);
      double t2 = (-b + discriminant)/(2.0 * a);

      result.Type = GeometryUtilities::IntersectionSegmentCircleResult::Types::TwoIntersections;
      result.SegmentIntersections.resize(2);
      result.SegmentIntersections[0].CurvilinearCoordinate = t1;
      result.SegmentIntersections[0].Type = PointSegmentPosition(t1);
      result.SegmentIntersections[1].CurvilinearCoordinate = t2;
      result.SegmentIntersections[1].Type = PointSegmentPosition(t2);
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::IntersectionPolygonCircleResult GeometryUtilities::IntersectionPolygonCircle(const Eigen::MatrixXd& polygonVertices,
                                                                                                  const Eigen::Vector3d& circleCenter,
                                                                                                  const double& circleRadius) const
  {
    GeometryUtilities::IntersectionPolygonCircleResult result;

    list<IntersectionPolygonCircleResult::Intersection> intersections;
    set<unsigned int> vertexIntersections;
    const unsigned int numEdges = polygonVertices.cols();
    for (unsigned int e = 0; e < numEdges; e++)
    {
      const unsigned int vertexOrigin = e;
      const unsigned int vertexEnd = (e + 1) % numEdges;
      const Vector3d& edgeOrigin = polygonVertices.col(vertexOrigin);
      const Vector3d& edgeEnd = polygonVertices.col(vertexEnd);

      IntersectionSegmentCircleResult intersection = IntersectionSegmentCircle(edgeOrigin,
                                                                               edgeEnd,
                                                                               circleCenter,
                                                                               circleRadius);

      if (intersection.Type == IntersectionSegmentCircleResult::Types::NoIntersection)
        continue;

      for (unsigned int i = 0; i < intersection.SegmentIntersections.size(); i++)
      {
        IntersectionSegmentCircleResult::IntersectionPosition& position = intersection.SegmentIntersections[i];
        Output::Assert(position.Type != PointSegmentPositionTypes::Unknown);
        switch (position.Type)
        {
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin:
          {
            if (vertexIntersections.find(vertexOrigin) == vertexIntersections.end())
            {
              vertexIntersections.insert(vertexOrigin);
              intersections.push_back(IntersectionPolygonCircleResult::Intersection());
              IntersectionPolygonCircleResult::Intersection& vertexIntersection = intersections.back();
              vertexIntersection.Type = (intersection.Type == IntersectionSegmentCircleResult::Types::TangentIntersection) ?
                                          IntersectionPolygonCircleResult::Intersection::Types::Tangent :
                                          IntersectionPolygonCircleResult::Intersection::Types::Secant;
              vertexIntersection.Index = vertexOrigin;
              vertexIntersection.IndexType = IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex;
            }
          }
            break;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment:
          {
            intersections.push_back(IntersectionPolygonCircleResult::Intersection());
            IntersectionPolygonCircleResult::Intersection& edgeIntersection = intersections.back();
            edgeIntersection.Type = (intersection.Type == IntersectionSegmentCircleResult::Types::TangentIntersection) ?
                                      IntersectionPolygonCircleResult::Intersection::Types::Tangent :
                                      IntersectionPolygonCircleResult::Intersection::Types::Secant;
            edgeIntersection.Index = e;
            edgeIntersection.IndexType = IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge;
            edgeIntersection.CurvilinearCoordinate = position.CurvilinearCoordinate;
          }
            break;
          case Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd:
          {
            if (vertexIntersections.find(vertexEnd) == vertexIntersections.end())
            {
              vertexIntersections.insert(vertexEnd);
              intersections.push_back(IntersectionPolygonCircleResult::Intersection());
              IntersectionPolygonCircleResult::Intersection& vertexIntersection = intersections.back();
              vertexIntersection.Type = (intersection.Type == IntersectionSegmentCircleResult::Types::TangentIntersection) ?
                                          IntersectionPolygonCircleResult::Intersection::Types::Tangent :
                                          IntersectionPolygonCircleResult::Intersection::Types::Secant;
              vertexIntersection.Index = vertexEnd;
              vertexIntersection.IndexType = IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex;
            }
          }
            break;
          default:
            continue;
        }
      }
    }

    result.Intersections = vector<IntersectionPolygonCircleResult::Intersection>(intersections.begin(),
                                                                                 intersections.end());
    return result;
  }
  // ***************************************************************************
}
