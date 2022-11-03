#include "IntersectorMesh2DSegment.hpp"
#include "IOUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  IntersectorMesh2DSegment::IntersectorMesh2DSegment(const Gedim::IMeshDAO& mesh,
                                                     const Gedim::GeometryUtilities& geometryUtilities) :
    _mesh(mesh),
    _geometryUtilities(geometryUtilities)
  {
  }
  IntersectorMesh2DSegment::~IntersectorMesh2DSegment()
  {
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::ToCurvilinearCoordinates(const IntersectorMesh2DSegment::IntersectionMesh& intersectingMesh,
                                                          vector<double>& curvilinearCoordinates)
  {
    curvilinearCoordinates.reserve(intersectingMesh.Points.size());
    for (std::map<double,
         IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>::const_iterator it = intersectingMesh.Points.begin();
         it != intersectingMesh.Points.end(); it++)
    {
      curvilinearCoordinates.push_back(it->first);
    }
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::ToString(const IntersectorMesh2DSegment::IntersectionMesh& intersectingMesh)
  {
    cerr.precision(16);
    cerr<< "Points:"<< endl;
    for_each(intersectingMesh.Points.begin(), intersectingMesh.Points.end(), [](const std::pair<double,
             IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>& p)
    { cerr<< scientific << "{ Key: " << p.first<< "; Value: V: "<< p.second.Vertex2DIds<< " E: "<< p.second.Edge2DIds<< " C: "<< p.second.Cell2DIds<< " }\n"; });
    cerr<< "Segments:"<< endl;
    for_each(intersectingMesh.Segments.begin(), intersectingMesh.Segments.end(), [](
             const IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshSegment& p)
    { cerr<< scientific << "{ E: "<< p.Edge2DIds<< " C: "<< p.Cell2DIds<< " }\n"; });
  }
  // ***************************************************************************
  IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint& IntersectorMesh2DSegment::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                                                     IntersectorMesh2DSegment::IntersectionMesh& result,
                                                                                                                     bool& found)
  {
    double foundCoordinate = -1.0;
    for (std::map<double,
         IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>::const_iterator it = result.Points.begin();
         it != result.Points.end(); it++)
    {
      if (!_geometryUtilities.IsValue1DPositive(abs(it->first - curvilinearCoordinate)))
      {
        foundCoordinate = it->first;
        break;
      }
    }

    if (foundCoordinate != -1.0)
    {
      found = true;
      return result.Points[foundCoordinate];
    }

    result.Points.insert(pair<double,
                         IntersectionMesh::IntersectionMeshPoint>(curvilinearCoordinate,
                                                                  IntersectionMesh::IntersectionMeshPoint()));
    found = false;
    return result.Points[curvilinearCoordinate];
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::CheckOriginAndEndSegmentPosition(const Vector3d& segmentOrigin,
                                                                  const Vector3d& segmentEnd,
                                                                  IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    for (unsigned int c = 0; c < _mesh.Cell2DTotalNumber(); c++)
    {
      if (!_mesh.Cell2DIsActive(c))
        continue;

      const MatrixXd cell2DVertices = _mesh.Cell2DVerticesCoordinates(c);

      Gedim::GeometryUtilities::PointPolygonPositionResult pointPolygonPositionResult = _geometryUtilities.PointPolygonPosition(segmentOrigin,
                                                                                                                                cell2DVertices);

      bool cellFound = false;
      switch (pointPolygonPositionResult.Type) {
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(0.0,
                                                                                        result,
                                                                                        found);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(0.0,
                                                                                        result,
                                                                                        found);

          unsigned int edgeId = _mesh.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
          intersection.Edge2DIds.insert(edgeId);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(0.0,
                                                                                        result,
                                                                                        found);
          unsigned int vertexId = _mesh.Cell2DVertex(c, pointPolygonPositionResult.BorderIndex);
          intersection.Vertex2DIds.insert(vertexId);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        default:
          continue;
      }

      if (cellFound)
        break;
    }

    // check if segment end is inside a single cell
    for (unsigned int c = 0; c < _mesh.Cell2DTotalNumber(); c++)
    {
      if (!_mesh.Cell2DIsActive(c))
        continue;

      const MatrixXd cell2DVertices = _mesh.Cell2DVerticesCoordinates(c);

      // check end position
      Gedim::GeometryUtilities::PointPolygonPositionResult pointPolygonPositionResult = _geometryUtilities.PointPolygonPosition(segmentEnd,
                                                                                                                                cell2DVertices);

      bool cellFound = false;
      switch (pointPolygonPositionResult.Type) {
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(1.0,
                                                                                        result,
                                                                                        found);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(1.0,
                                                                                        result,
                                                                                        found);
          unsigned int edgeId = _mesh.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
          intersection.Edge2DIds.insert(edgeId);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(1.0,
                                                                                        result,
                                                                                        found);
          unsigned int vertexId = _mesh.Cell2DVertex(c, pointPolygonPositionResult.BorderIndex);
          intersection.Vertex2DIds.insert(vertexId);
          intersection.Cell2DIds.insert(c);
          cellFound = true;
          break;
        }
        default:
          continue;
      }

      if (cellFound)
        break;
    }
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::CreateIntersectionPoints(const Vector3d& segmentOrigin,
                                                          const Vector3d& segmentEnd,
                                                          const Vector3d& ,
                                                          const Vector3d& segmentBarycenter,
                                                          const double& segmentLength,
                                                          IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    for (unsigned int e = 0; e < _mesh.Cell1DTotalNumber(); e++)
    {
      if (!_mesh.Cell1DIsActive(e))
        continue;

      const unsigned int edgeOriginId = _mesh.Cell1DOrigin(e);
      const unsigned int edgeEndId = _mesh.Cell1DEnd(e);
      const Vector3d edgeOrigin = _mesh.Cell1DOriginCoordinates(e);
      const Vector3d edgeEnd = _mesh.Cell1DEndCoordinates(e);

      const Eigen::Vector3d edgeBarycenter = _geometryUtilities.SegmentBarycenter(edgeOrigin,
                                                                                  edgeEnd);
      const double edgeLength = _geometryUtilities.SegmentLength(edgeOrigin,
                                                                 edgeEnd);

      if (_geometryUtilities.CheckNoSpheresIntersection(edgeBarycenter,
                                                        segmentBarycenter,
                                                        edgeLength,
                                                        segmentLength))
        continue;

      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult intersectionSegmentSegmentResult = _geometryUtilities.IntersectionSegmentSegment(edgeOrigin,
                                                                                                                                                  edgeEnd,
                                                                                                                                                  segmentOrigin,
                                                                                                                                                  segmentEnd);
      // no intersection found
      if (intersectionSegmentSegmentResult.IntersectionSegmentsType ==
          Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection)
        continue;

      // check insersection type
      if (intersectionSegmentSegmentResult.IntersectionSegmentsType ==
          Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
      {
        // intersection not considered if outside the segment
        const Gedim::GeometryUtilities::PointSegmentPositionTypes& segmentIntersectionType = intersectionSegmentSegmentResult.SecondSegmentIntersections[0].Type;
        if (segmentIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin &&
            segmentIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment &&
            segmentIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
          continue;

        // intersection not considered if outside the edge
        const Gedim::GeometryUtilities::PointSegmentPositionTypes& edgeIntersectionType = intersectionSegmentSegmentResult.FirstSegmentIntersections[0].Type;
        if (edgeIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin &&
            edgeIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment &&
            edgeIntersectionType != Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
          continue;

        // only one intersection found, the segment intersects the edge
        const double& curvilinearCoordinate = intersectionSegmentSegmentResult.SecondSegmentIntersections[0].CurvilinearCoordinate;

        // check the position on the edge
        if (edgeIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeOriginId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeEndId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
      }
      else if (intersectionSegmentSegmentResult.IntersectionSegmentsType ==
               Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections)
      {
        // multiple intersections found, segment is coincident with edge
        const double& segmentFirstCurvilinearCoordinate = intersectionSegmentSegmentResult.SecondSegmentIntersections[0].CurvilinearCoordinate;
        const double& segmentSecondCurvilinearCoordinate = intersectionSegmentSegmentResult.SecondSegmentIntersections[1].CurvilinearCoordinate;
        const unsigned int& indexEdgeFirstIntersection = intersectionSegmentSegmentResult.SecondIntersectionRelation[0];
        const unsigned int& indexEdgeSecondIntersection = intersectionSegmentSegmentResult.SecondIntersectionRelation[1];

        const Gedim::GeometryUtilities::PointSegmentPositionTypes& edgeFirstIntersectionType = intersectionSegmentSegmentResult.FirstSegmentIntersections[indexEdgeFirstIntersection].Type;
        const Gedim::GeometryUtilities::PointSegmentPositionTypes& edgeSecondIntersectionType = intersectionSegmentSegmentResult.FirstSegmentIntersections[indexEdgeSecondIntersection].Type;

        if (edgeFirstIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentFirstCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeOriginId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeFirstIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentFirstCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeEndId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeFirstIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentFirstCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }

        if (edgeSecondIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentSecondCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeOriginId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeSecondIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentSecondCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Vertex2DIds.insert(edgeEndId);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
        else if (edgeSecondIntersectionType == Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment)
        {
          bool found;
          IntersectionMesh::IntersectionMeshPoint& intersection = InsertNewIntersection(segmentSecondCurvilinearCoordinate,
                                                                                        result,
                                                                                        found);
          intersection.Edge2DIds.insert(e);
          for (unsigned int n = 0; n < _mesh.Cell1DNumberNeighbourCell2D(e); n++)
          {
            if (!_mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            intersection.Cell2DIds.insert(_mesh.Cell1DNeighbourCell2D(e, n));
          }
        }
      }
      else
        throw runtime_error("Intersection type not managed");
    }
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::CreateIntersectionSegments(IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    result.Segments.resize(result.Points.size() - 1);
    map<double, IntersectionMesh::IntersectionMeshPoint>::const_iterator itPoint = result.Points.begin();
    map<double, IntersectionMesh::IntersectionMeshPoint>::const_iterator itPointNext = result.Points.begin();
    itPointNext++;
    for (unsigned int p = 0; p < result.Segments.size(); p++)
    {
      const double& curvilinearCoordinatePoint = itPoint->first;
      const double& curvilinearCoordinatePointNext = itPointNext->first;
      const IntersectionMesh::IntersectionMeshPoint& intersectionPoint = itPoint->second;
      const IntersectionMesh::IntersectionMeshPoint& intersectionPointNext = itPointNext->second;

      // fill origin and end of segment
      IntersectionMesh::IntersectionMeshSegment& meshSegment = result.Segments[p];
      meshSegment.Points.resize(2);
      meshSegment.Points[0] = curvilinearCoordinatePoint;
      meshSegment.Points[1] = curvilinearCoordinatePointNext;

      // fill the mesh 2D edges
      set_intersection(intersectionPoint.Edge2DIds.begin(),
                       intersectionPoint.Edge2DIds.end(),
                       intersectionPointNext.Edge2DIds.begin(),
                       intersectionPointNext.Edge2DIds.end(),
                       std::inserter(meshSegment.Edge2DIds,
                                     meshSegment.Edge2DIds.begin()));
      // fill the mesh 2D cells
      set_intersection(intersectionPoint.Cell2DIds.begin(),
                       intersectionPoint.Cell2DIds.end(),
                       intersectionPointNext.Cell2DIds.begin(),
                       intersectionPointNext.Cell2DIds.end(),
                       std::inserter(meshSegment.Cell2DIds,
                                     meshSegment.Cell2DIds.begin()));
      itPoint++;
      itPointNext++;
    }
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::CheckOriginAndEndPointExistence(IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    bool found = false;
    const pair<double, IntersectionMesh::IntersectionMeshPoint>& firstElement = *result.Points.begin();
    IntersectionMesh::IntersectionMeshPoint& firstIntersection = InsertNewIntersection(0.0, result, found);

    if (!found)
    {
      // the origin of the segment is outside the mesh
      // replace the first intersection found with 0.0
      firstIntersection = firstElement.second;
      result.Points.erase(firstElement.first);
    }

    found = false;
    const pair<double, IntersectionMesh::IntersectionMeshPoint>& lastElement = *result.Points.rbegin();
    IntersectionMesh::IntersectionMeshPoint& lastIntersection = InsertNewIntersection(1.0, result, found);

    if (!found)
    {
      // the end of the segment is outside the mesh
      // replace the last intersection found with 1.0
      lastIntersection = lastElement.second;
      result.Points.erase(lastElement.first);
    }
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::SmoothIntersections(const Vector3d& ,
                                                     const Vector3d& ,
                                                     IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    //const Vector3d segmentTangent = _geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);

    for (map<double, IntersectionMesh::IntersectionMeshPoint>::iterator itPoint = result.Points.begin();
         itPoint != result.Points.end();
         itPoint++)
    {
      IntersectionMesh::IntersectionMeshPoint& intersectionPoint = itPoint->second;

      // Extract active edges in the mesh
      list<unsigned int> activeEdges;
      for (unsigned int edgeId : intersectionPoint.Edge2DIds)
      {
        if (!_mesh.Cell1DIsActive(edgeId))
          continue;

        activeEdges.push_back(edgeId);
      }

      // if an intersection has more then one edge and no vertices then the intersection shall be
      // near the vertex in common of the edges
      // error: the program should be started with greater tolerance
      if (intersectionPoint.Vertex2DIds.size() == 0 &&
          activeEdges.size() > 1 &&
          intersectionPoint.Cell2DIds.size() > 1)
      {
        unsigned int vertexInCommon = _mesh.Cell0DTotalNumber();
        unsigned int numEdges = activeEdges.size();
        list<unsigned int>::const_iterator itEdge = activeEdges.begin();
        list<unsigned int>::const_iterator itEdgeNext = activeEdges.begin();
        itEdgeNext++;

        for (unsigned int e = 0; e < numEdges - 1; e++)
        {
          if (e == 0)
          {
            Gedim::Output::Assert(_mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                                  _mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext) ||
                                  _mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                                  _mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext));
            vertexInCommon = (_mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                              _mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext)) ? _mesh.Cell1DOrigin(*itEdge) :
                                                                                             _mesh.Cell1DEnd(*itEdge);
          }
          else
          {
            Gedim::Output::Assert(_mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                                  _mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext) ||
                                  _mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                                  _mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext));
            if (_mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                _mesh.Cell1DOrigin(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext))
              Gedim::Output::Assert(vertexInCommon == _mesh.Cell1DOrigin(*itEdge));
            if (_mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DOrigin(*itEdgeNext) ||
                _mesh.Cell1DEnd(*itEdge) == _mesh.Cell1DEnd(*itEdgeNext))
              Gedim::Output::Assert(vertexInCommon == _mesh.Cell1DEnd(*itEdge));
          }
        }
        Gedim::Output::Assert(vertexInCommon != _mesh.Cell0DTotalNumber());

        throw runtime_error("Find problem with tolerance. Try to restart computations with a greater tolerance.");
      }
    }

    map<double, IntersectionMesh::IntersectionMeshPoint> smoothedPoints;
    smoothedPoints.insert(*result.Points.begin());

    const unsigned int numSegments = result.Points.size() - 1;
    map<double, IntersectionMesh::IntersectionMeshPoint>::const_iterator itPoint = result.Points.begin();
    map<double, IntersectionMesh::IntersectionMeshPoint>::const_iterator itPointNext = result.Points.begin();
    itPointNext++;
    for (unsigned int p = 0; p < numSegments; p++)
    {
      const double& curvilinearCoordinatePoint = itPoint->first;
      const double& curvilinearCoordinatePointNext = itPointNext->first;
      const IntersectionMesh::IntersectionMeshPoint& intersectionPoint = itPoint->second;
      const IntersectionMesh::IntersectionMeshPoint& intersectionPointNext = itPointNext->second;

      // If the two intersection points have the same set of vertices then the point is the same
      if (intersectionPoint.Vertex2DIds.size() > 0 &&
          intersectionPointNext.Vertex2DIds.size() > 0 &&
          intersectionPoint.Vertex2DIds == intersectionPointNext.Vertex2DIds)
      {
        const double avgCoordinate = (curvilinearCoordinatePoint + curvilinearCoordinatePointNext) / 2.0;
        Output::Assert(!_geometryUtilities.IsValue1DPositive(abs(curvilinearCoordinatePoint - avgCoordinate)));
        Output::Assert(!_geometryUtilities.IsValue1DPositive(abs(curvilinearCoordinatePointNext - avgCoordinate)));

        IntersectionMesh::IntersectionMeshPoint avgPoint;
        avgPoint.Vertex2DIds = intersectionPoint.Vertex2DIds;
        avgPoint.Edge2DIds = intersectionPoint.Edge2DIds;
        avgPoint.Cell2DIds = intersectionPoint.Cell2DIds;

        // unify the Cell1Ds and Cell2Ds
        std::merge(intersectionPoint.Edge2DIds.begin(),
                   intersectionPoint.Edge2DIds.end(),
                   intersectionPointNext.Edge2DIds.begin(),
                   intersectionPointNext.Edge2DIds.end(),
                   std::inserter(avgPoint.Edge2DIds,
                                 avgPoint.Edge2DIds.begin()));

        std::merge(intersectionPoint.Cell2DIds.begin(),
                   intersectionPoint.Cell2DIds.end(),
                   intersectionPointNext.Cell2DIds.begin(),
                   intersectionPointNext.Cell2DIds.end(),
                   std::inserter(avgPoint.Cell2DIds,
                                 avgPoint.Cell2DIds.begin()));

        smoothedPoints.erase(curvilinearCoordinatePoint);
        result.Points.erase(curvilinearCoordinatePoint);
        result.Points.erase(curvilinearCoordinatePointNext);
        result.Points.insert(make_pair(avgCoordinate, avgPoint));
        smoothedPoints.insert(make_pair(avgCoordinate, avgPoint));

        itPoint = result.Points.begin();
        itPointNext = result.Points.begin();
        itPointNext++;
        for (unsigned int s = 0; s <= p; s++)
        {
          itPoint++;
          itPointNext++;
        }


        continue;
      }

      // insert the next point on the smoothed list
      smoothedPoints.insert(*itPointNext);

      itPoint++;
      itPointNext++;
    }

    result.Points = smoothedPoints;
  }
  // ***************************************************************************
  void IntersectorMesh2DSegment::CreateIntersectionMesh(const Vector3d& segmentOrigin,
                                                        const Vector3d& segmentEnd,
                                                        const Vector3d& segmentTangent,
                                                        const Vector3d& segmentBarycenter,
                                                        const double& segmentLength,
                                                        IntersectorMesh2DSegment::IntersectionMesh& result)
  {
    // check if segment origin is inside a single cell
    CheckOriginAndEndSegmentPosition(segmentOrigin,
                                     segmentEnd,
                                     result);

    // check all segment intersection points with mesh edges
    CreateIntersectionPoints(segmentOrigin,
                             segmentEnd,
                             segmentTangent,
                             segmentBarycenter,
                             segmentLength,
                             result);

    Gedim::Output::Assert(!result.Points.empty() && result.Points.size() > 1);

    // check origin and end of segment existence in intersection mesh
    CheckOriginAndEndPointExistence(result);

    // smooth the intersections
    SmoothIntersections(segmentOrigin, segmentEnd, result);

    // create 1D intersection segments
    CreateIntersectionSegments(result);
  }
  // ***************************************************************************
}
