#include "ConformerMeshSegment.hpp"

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  ConformerMeshSegment::ConformerMeshSegment(const GeometryUtilities& geometryUtilities) :
    _geometryUtilities(geometryUtilities)
  {
  }
  ConformerMeshSegment::~ConformerMeshSegment()
  {
  }
  // ***************************************************************************
  void ConformerMeshSegment::ToCurvilinearCoordinates(const ConformerMeshSegment::ConformMesh& conformMesh,
                                                      vector<double>& curvilinearCoordinates)
  {
    curvilinearCoordinates.reserve(conformMesh.Points.size());
    for (std::map<double,
         ConformerMeshSegment::ConformMesh::ConformMeshPoint>::const_iterator it = conformMesh.Points.begin();
         it != conformMesh.Points.end(); it++)
    {
      curvilinearCoordinates.push_back(it->first);
    }
  }
  // ***************************************************************************
  void ConformerMeshSegment::ToString(const ConformerMeshSegment::ConformMesh& conformMesh)
  {
    cerr.precision(16);
    cerr<< "Points:"<< endl;
    for_each(conformMesh.Points.begin(), conformMesh.Points.end(), [](const std::pair<double,
             ConformerMeshSegment::ConformMesh::ConformMeshPoint>& p)
    { cerr<< scientific << "{ Key: " << p.first<< "; Value: V: "<< p.second.Vertex2DIds<< " E: "<< p.second.Edge2DIds<< " C: "<< p.second.Cell2DIds<< " }\n"; });
    cerr<< "Segments:"<< endl;
    for_each(conformMesh.Segments.begin(), conformMesh.Segments.end(), [](
             const ConformerMeshSegment::ConformMesh::ConformMeshSegment& p)
    { cerr<< scientific << "{ Key: "<< p.Points[0]<< "->"<< p.Points[1] << "; Value: E: "<< p.Edge2DIds<< " C: "<< p.Cell2DIds<< " }\n"; });
  }
  // ***************************************************************************
  ConformerMeshSegment::ConformMesh::ConformMeshPoint& ConformerMeshSegment::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                                   ConformerMeshSegment::ConformMesh& result,
                                                                                                   bool& found)
  {
    double foundCoordinate = -1.0;
    for (std::map<double,
         ConformMesh::ConformMeshPoint>::const_iterator it = result.Points.begin();
         it != result.Points.end(); it++)
    {
      GeometryUtilities::CompareTypes result = _geometryUtilities.Compare1DValues(it->first,
                                                                                  curvilinearCoordinate);

      if (result == GeometryUtilities::CompareTypes::Coincident)
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
                         ConformMesh::ConformMeshPoint>(curvilinearCoordinate,
                                                        ConformMesh::ConformMeshPoint()));
    found = false;
    return result.Points[curvilinearCoordinate];
  }
  // ***************************************************************************
  void ConformerMeshSegment::CreateConformPoints(const Gedim::IntersectorMesh2DSegment::IntersectionMesh& meshIntersection,
                                                 const Gedim::UnionMeshSegment::UnionMesh& meshUnion,
                                                 const unsigned int& meshIntersectionPosition,
                                                 ConformerMeshSegment::ConformMesh& result)
  {
    map<double,
        Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>::const_iterator itIntersectionPoint = meshIntersection.Points.begin();
    unsigned int intersectionSegmentIndex = 0;
    map<double,
        UnionMeshSegment::UnionMesh::UnionMeshPoint>::const_iterator itUnionPoint = meshUnion.Points.begin();

    UnionMeshSegment::UnionMesh::UnionMeshPoint::Types otherIntersectionMeshUnionType = (meshIntersectionPosition == 0) ?
                                                                                          UnionMeshSegment::UnionMesh::UnionMeshPoint::Second :
                                                                                          UnionMeshSegment::UnionMesh::UnionMeshPoint::First;
    for (unsigned int p = 0; p < meshUnion.Points.size(); p++)
    {
      const double& unionCurvilinearCoordinate = itUnionPoint->first;

      const UnionMeshSegment::UnionMesh::UnionMeshPoint& unionPoint = itUnionPoint->second;

      bool found;
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& conformPoint = InsertNewIntersection(unionCurvilinearCoordinate,
                                                                                                result,
                                                                                                found);
      // Point already added, not possible
      Output::Assert(!found);

      // check if union point is in intersection mesh
      if (unionPoint.Type != otherIntersectionMeshUnionType)
      {
        // union point belongs to mine intersection mesh
        // inherit properties from intersection mesh point
        const IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint& intersectionPoint = itIntersectionPoint->second;
        conformPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;

        conformPoint.Vertex2DIds.resize(intersectionPoint.Vertex2DIds.size());
        conformPoint.Edge2DIds.resize(intersectionPoint.Edge2DIds.size());
        conformPoint.Cell2DIds.resize(intersectionPoint.Cell2DIds.size());
        copy(intersectionPoint.Vertex2DIds.begin(), intersectionPoint.Vertex2DIds.end(), conformPoint.Vertex2DIds.begin());
        copy(intersectionPoint.Edge2DIds.begin(), intersectionPoint.Edge2DIds.end(), conformPoint.Edge2DIds.begin());
        copy(intersectionPoint.Cell2DIds.begin(), intersectionPoint.Cell2DIds.end(), conformPoint.Cell2DIds.begin());

        itIntersectionPoint++;
        intersectionSegmentIndex++;
      }
      else
      {
        Output::Assert(intersectionSegmentIndex > 0 && intersectionSegmentIndex <= meshIntersection.Segments.size());

        // union point belongs to the other intersection mesh
        // add the new point to conform mesh and inherit properties from intersection mesh segment
        const IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshSegment& intersectionSegment = meshIntersection.Segments[intersectionSegmentIndex - 1];
        conformPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::Inherited;

        conformPoint.Edge2DIds.resize(intersectionSegment.Edge2DIds.size());
        conformPoint.Cell2DIds.resize(intersectionSegment.Cell2DIds.size());
        copy(intersectionSegment.Edge2DIds.begin(), intersectionSegment.Edge2DIds.end(), conformPoint.Edge2DIds.begin());
        copy(intersectionSegment.Cell2DIds.begin(), intersectionSegment.Cell2DIds.end(), conformPoint.Cell2DIds.begin());
      }

      itUnionPoint++;
    }
  }
  // ***************************************************************************
  void ConformerMeshSegment::CreateConformSegments(ConformerMeshSegment::ConformMesh& result)
  {
    result.Segments.resize(result.Points.size() - 1);
    map<double, ConformMesh::ConformMeshPoint>::const_iterator itPoint = result.Points.begin();
    map<double, ConformMesh::ConformMeshPoint>::const_iterator itPointNext = result.Points.begin();
    itPointNext++;
    for (unsigned int p = 0; p < result.Segments.size(); p++)
    {
      const double& curvilinearCoordinatePoint = itPoint->first;
      const double& curvilinearCoordinatePointNext = itPointNext->first;
      const ConformMesh::ConformMeshPoint& intersectionPoint = itPoint->second;
      const ConformMesh::ConformMeshPoint& intersectionPointNext = itPointNext->second;

      // fill origin and end of segment
      ConformMesh::ConformMeshSegment& meshSegment = result.Segments[p];
      meshSegment.Points.resize(2);
      meshSegment.Points[0] = curvilinearCoordinatePoint;
      meshSegment.Points[1] = curvilinearCoordinatePointNext;

      // fill the mesh 2D edges
      set_intersection(intersectionPoint.Edge2DIds.begin(),
                       intersectionPoint.Edge2DIds.end(),
                       intersectionPointNext.Edge2DIds.begin(),
                       intersectionPointNext.Edge2DIds.end(),
                       back_inserter(meshSegment.Edge2DIds));
      // fill the mesh 2D cells
      set_intersection(intersectionPoint.Cell2DIds.begin(),
                       intersectionPoint.Cell2DIds.end(),
                       intersectionPointNext.Cell2DIds.begin(),
                       intersectionPointNext.Cell2DIds.end(),
                       back_inserter(meshSegment.Cell2DIds));
      itPoint++;
      itPointNext++;
    }
  }
  // ***************************************************************************
  void ConformerMeshSegment::CreateConformMesh(const Gedim::IntersectorMesh2DSegment::IntersectionMesh& meshIntersection,
                                               const UnionMeshSegment::UnionMesh& meshUnion,
                                               const unsigned int& meshIntersectionPosition,
                                               ConformMesh& result)
  {
    CreateConformPoints(meshIntersection,
                        meshUnion,
                        meshIntersectionPosition,
                        result);
    CreateConformSegments(result);
  }
  // ***************************************************************************
  void ConformerMeshSegment::InsertExternalPoint(const Vector3d& segmentOrigin,
                                                 const Vector3d& segmentEnd,
                                                 const Gedim::IMeshDAO& mesh2D,
                                                 const double& curvilinearCoordinate,
                                                 ConformerMeshSegment::ConformMesh& result)
  {
    bool found;
    ConformerMeshSegment::ConformMesh::ConformMeshPoint& newPoint = InsertNewIntersection(curvilinearCoordinate,
                                                                                          result,
                                                                                          found);
    // Point already added
    if (found)
      return;

    newPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
    const Vector3d newPoint2D = segmentOrigin + curvilinearCoordinate * (segmentEnd - segmentOrigin);

    // look for 2D cell which contains the new point
    for (unsigned int c = 0; c < mesh2D.Cell2DTotalNumber(); c++)
    {
      if (!mesh2D.Cell2DIsActive(c))
        continue;

      const MatrixXd cell2DVertices = mesh2D.Cell2DVerticesCoordinates(c);

      GeometryUtilities::PointPolygonPositionResult pointPolygonPositionResult;
      _geometryUtilities.PointPolygonPosition(newPoint2D,
                                              cell2DVertices,
                                              pointPolygonPositionResult);

      if (pointPolygonPositionResult.PositionType == GeometryUtilities::PointPolygonPositionResult::Outside)
        continue;

      newPoint.Cell2DIds.push_back(c);
      if (pointPolygonPositionResult.PositionType == GeometryUtilities::PointPolygonPositionResult::Inside)
        break;

      else if (pointPolygonPositionResult.PositionType == GeometryUtilities::PointPolygonPositionResult::BorderEdge)
      {
        const unsigned int edgeToInsert = mesh2D.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
        if (find(newPoint.Edge2DIds.begin(), newPoint.Edge2DIds.end(), edgeToInsert) == newPoint.Edge2DIds.end())
          newPoint.Edge2DIds.push_back(edgeToInsert);
      }
      else if (pointPolygonPositionResult.PositionType == GeometryUtilities::PointPolygonPositionResult::BorderVertex)
      {
        const unsigned int vertexToInsert = mesh2D.Cell2DVertex(c, pointPolygonPositionResult.BorderIndex);
        int previousVertex = pointPolygonPositionResult.BorderIndex - 1;
        const unsigned int edgeToInsert = mesh2D.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
        const unsigned int previousEdgeToInsert = mesh2D.Cell2DEdge(c, (previousVertex < 0) ? mesh2D.Cell2DNumberVertices(c) - 1 :
                                                                                              previousVertex);
        if (find(newPoint.Vertex2DIds.begin(), newPoint.Vertex2DIds.end(), vertexToInsert) == newPoint.Vertex2DIds.end())
          newPoint.Vertex2DIds.push_back(vertexToInsert);
        if (find(newPoint.Edge2DIds.begin(), newPoint.Edge2DIds.end(), edgeToInsert) == newPoint.Edge2DIds.end())
          newPoint.Vertex2DIds.push_back(edgeToInsert);
        if (find(newPoint.Edge2DIds.begin(), newPoint.Edge2DIds.end(), previousEdgeToInsert) == newPoint.Edge2DIds.end())
          newPoint.Vertex2DIds.push_back(previousEdgeToInsert);
      }
    }

    // recreate conform mesh segments
    result.Segments.clear();
    CreateConformSegments(result);
  }
  // ***************************************************************************
}
