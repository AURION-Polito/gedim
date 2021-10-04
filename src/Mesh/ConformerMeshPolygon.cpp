#include "ConformerMeshPolygon.hpp"

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  ConformerMeshPolygon::ConformerMeshPolygon(const GeometryUtilities& geometryUtilities) :
    _geometryUtilities(geometryUtilities)
  {
  }
  ConformerMeshPolygon::~ConformerMeshPolygon()
  {
  }
  // ***************************************************************************
  ConformerMeshSegment::ConformMesh::ConformMeshPoint& ConformerMeshPolygon::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                                   ConformerMeshSegment::ConformMesh& result,
                                                                                                   bool& found)
  {
    double foundCoordinate = -1.0;
    for (std::map<double,
         ConformerMeshSegment::ConformMesh::ConformMeshPoint>::const_iterator it = result.Points.begin();
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
                         ConformerMeshSegment::ConformMesh::ConformMeshPoint>(curvilinearCoordinate,
                                                                              ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
    found = false;
    return result.Points[curvilinearCoordinate];
  }
  // ***************************************************************************
  void ConformerMeshPolygon::CheckSegmentOriginAndEndIntersections(const Vector3d& segmentOrigin,
                                                                   const Vector3d& segmentEnd,
                                                                   ConformerMeshSegment::ConformMesh& mesh1D,
                                                                   const Gedim::IMeshDAO& mesh2DReader)
  {
    const double originCurvilinearCoordinate = mesh1D.Segments[0].Points[0];
    const double endCurvilinearCoordinate = mesh1D.Segments[mesh1D.Segments.size() - 1].Points[1];

    const ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points.at(originCurvilinearCoordinate);
    const ConformerMeshSegment::ConformMesh::ConformMeshPoint& end = mesh1D.Points.at(endCurvilinearCoordinate);

    // check if nothing to do
    if ((origin.Vertex2DIds.size() != 0 || origin.Edge2DIds.size() != 0) &&
        (end.Vertex2DIds.size() != 0 || end.Edge2DIds.size() != 0))
      return;

    // check if origin is inside cell
    if (origin.Vertex2DIds.size() == 0 && origin.Edge2DIds.size() == 0)
    {
      if (origin.Cell2DIds.size() != 1)
        throw runtime_error("Error on Cell0DMesh1D origin Cell2DMesh2D");

      const unsigned int cell2DMesh2DId =  origin.Cell2DIds.front();

      // check intersection between Cell1DMesh2D and segment
      for (unsigned int e = 0; e < mesh2DReader.Cell2DNumberEdges(cell2DMesh2DId); e++)
      {
        const unsigned int edgeId = mesh2DReader.Cell2DEdge(cell2DMesh2DId, e);
        const unsigned int edgeOriginId = mesh2DReader.Cell1DOrigin(edgeId);
        const unsigned int edgeEndId = mesh2DReader.Cell1DEnd(edgeId);
        Vector3d edgeOrigin = mesh2DReader.Cell0DCoordinates(edgeOriginId);
        Vector3d edgeEnd = mesh2DReader.Cell0DCoordinates(edgeEndId);

        GeometryUtilities::IntersectionSegmentSegmentResult result = _geometryUtilities.IntersectionSegmentSegment(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   edgeOrigin,
                                                                                                                   edgeEnd);

        if (result.IntersectionLinesType != GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting ||
            result.IntersectionSegmentsType ==  GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
          continue;

        if (result.FirstSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
          continue;

        if (result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::InsideSegment &&
            result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin &&
            result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
          continue;

        // insert the new intersection in mesh1D
        const double intersectionCurvilinearCoordinate = result.FirstSegmentIntersections[0].CurvilinearCoordinate;

        bool found = false;
        ConformerMeshSegment::ConformMesh::ConformMeshPoint& newPoint = InsertNewIntersection(intersectionCurvilinearCoordinate,
                                                                                              mesh1D,
                                                                                              found);

        newPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        newPoint.Edge2DIds.push_back(edgeId);

        if (find(newPoint.Cell2DIds.begin(), newPoint.Cell2DIds.end(), cell2DMesh2DId) == newPoint.Cell2DIds.end())
          newPoint.Cell2DIds.push_back(cell2DMesh2DId);

        for (unsigned int en = 0; en < mesh2DReader.Cell1DNumberNeighbourCell2D(edgeId); en++)
        {
          if (!mesh2DReader.Cell1DHasNeighbourCell2D(edgeId, en))
            continue;

          if (mesh2DReader.Cell1DNeighbourCell2D(edgeId, en) != cell2DMesh2DId)
            newPoint.Cell2DIds.push_back(mesh2DReader.Cell1DNeighbourCell2D(edgeId, en));
        }

        if (result.SecondSegmentIntersections[0].Type == GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          if (find(newPoint.Vertex2DIds.begin(), newPoint.Vertex2DIds.end(), edgeOriginId) == newPoint.Vertex2DIds.end())
            newPoint.Vertex2DIds.push_back(edgeOriginId);
        }
        else if (result.SecondSegmentIntersections[0].Type == GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          if (find(newPoint.Vertex2DIds.begin(), newPoint.Vertex2DIds.end(), edgeEndId) == newPoint.Vertex2DIds.end())
            newPoint.Vertex2DIds.push_back(edgeEndId);
        }
      }
    }

    // check if origin is inside cell
    if (end.Vertex2DIds.size() == 0 && end.Edge2DIds.size() == 0)
    {
      if (end.Cell2DIds.size() != 1)
        throw runtime_error("Error on Cell0DMesh1D end Cell2DMesh2D");

      const unsigned int cell2DMesh2DId =  end.Cell2DIds.front();

      // check intersection between Cell1DMesh2D and segment
      for (unsigned int e = 0; e < mesh2DReader.Cell2DNumberEdges(cell2DMesh2DId); e++)
      {
        const unsigned int edgeId = mesh2DReader.Cell2DEdge(cell2DMesh2DId, e);
        const unsigned int edgeOriginId = mesh2DReader.Cell1DOrigin(edgeId);
        const unsigned int edgeEndId = mesh2DReader.Cell1DEnd(edgeId);
        Vector3d edgeOrigin = mesh2DReader.Cell0DCoordinates(edgeOriginId);
        Vector3d edgeEnd = mesh2DReader.Cell0DCoordinates(edgeEndId);

        GeometryUtilities::IntersectionSegmentSegmentResult result = _geometryUtilities.IntersectionSegmentSegment(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   edgeOrigin,
                                                                                                                   edgeEnd);

        if (result.IntersectionLinesType != GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting ||
            result.IntersectionSegmentsType ==  GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
          continue;

        if (result.FirstSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd)
          continue;

        if (result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::InsideSegment &&
            result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin &&
            result.SecondSegmentIntersections[0].Type != GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
          continue;

        // insert the new intersection in mesh1D
        const double intersectionCurvilinearCoordinate = result.FirstSegmentIntersections[0].CurvilinearCoordinate;

        bool found = false;
        ConformerMeshSegment::ConformMesh::ConformMeshPoint& newPoint = InsertNewIntersection(intersectionCurvilinearCoordinate,
                                                                                              mesh1D,
                                                                                              found);

        newPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
        newPoint.Edge2DIds.push_back(edgeId);
        newPoint.Cell2DIds.push_back(cell2DMesh2DId);
        for (unsigned int en = 0; en < mesh2DReader.Cell1DNumberNeighbourCell2D(edgeId); en++)
        {
          if (!mesh2DReader.Cell1DHasNeighbourCell2D(edgeId, en))
            continue;

          if (mesh2DReader.Cell1DNeighbourCell2D(edgeId, en) != cell2DMesh2DId)
            newPoint.Cell2DIds.push_back(mesh2DReader.Cell1DNeighbourCell2D(edgeId, en));
        }

        if (result.SecondSegmentIntersections[0].Type == GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin)
        {
          if (find(newPoint.Vertex2DIds.begin(), newPoint.Vertex2DIds.end(), edgeOriginId) == newPoint.Vertex2DIds.end())
            newPoint.Vertex2DIds.push_back(edgeOriginId);
        }
        else if (result.SecondSegmentIntersections[0].Type == GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd)
        {
          if (find(newPoint.Vertex2DIds.begin(), newPoint.Vertex2DIds.end(), edgeEndId) == newPoint.Vertex2DIds.end())
            newPoint.Vertex2DIds.push_back(edgeEndId);
        }
      }
    }

    // recreate conform mesh segments
    mesh1D.Segments.clear();
    ConformerMeshSegment::CreateConformSegments(mesh1D);
  }
  // ***************************************************************************
  void ConformerMeshPolygon::Cell2DMesh2DListCell1DToListCell0D(const Gedim::IMeshDAO& mesh2DReader,
                                                                const unsigned int& cell1DMesh2DId,
                                                                const list<unsigned int>& cell1DMesh2DUpdated,
                                                                vector<unsigned int>& cell0DMesh2Ds,
                                                                vector<unsigned int>& cell1DMesh2Ds)
  {
    cell0DMesh2Ds.resize(cell1DMesh2DUpdated.size() + 1);
    cell1DMesh2Ds.resize(cell1DMesh2DUpdated.size());

    // compute new vertices and edges order
    cell0DMesh2Ds[0] = mesh2DReader.Cell1DOrigin(cell1DMesh2DId);
    cell0DMesh2Ds[cell1DMesh2DUpdated.size()] = mesh2DReader.Cell1DEnd(cell1DMesh2DId);

    bool edgeListDirection = (mesh2DReader.Cell1DOrigin(cell1DMesh2DUpdated.front()) == cell0DMesh2Ds[0] ||
                             mesh2DReader.Cell1DEnd(cell1DMesh2DUpdated.front()) == cell0DMesh2Ds[0]); // true if forward, false if backward

    unsigned int counter = 1;
    if (edgeListDirection)
    {
      for (list<unsigned int>::const_iterator iter = cell1DMesh2DUpdated.begin();
           iter != cell1DMesh2DUpdated.end();
           iter++)
      {
        cell0DMesh2Ds[counter] = (mesh2DReader.Cell1DOrigin(*iter) == cell0DMesh2Ds[counter - 1]) ?
              mesh2DReader.Cell1DEnd(*iter) : mesh2DReader.Cell1DOrigin(*iter);
        cell1DMesh2Ds[counter - 1] = *iter;
        counter++;
      }
    }
    else
    {
      for (list<unsigned int>::const_reverse_iterator iter = cell1DMesh2DUpdated.rbegin();
           iter != cell1DMesh2DUpdated.rend();
           iter++)
      {
        cell0DMesh2Ds[counter] = (mesh2DReader.Cell1DOrigin(*iter) == cell0DMesh2Ds[counter - 1]) ?
              mesh2DReader.Cell1DEnd(*iter) : mesh2DReader.Cell1DOrigin(*iter);
        cell1DMesh2Ds[counter - 1] = *iter;
        counter++;
      }
    }
  }
  // ***************************************************************************
  void ConformerMeshPolygon::Cell2DMesh2DUpdatedPointsAndEdges(const Gedim::IMeshDAO& mesh2DReader,
                                                               const unsigned int& cell2DMesh2DId,
                                                               const unsigned int& oldCell1DMesh2DId,
                                                               const vector<unsigned int>& cell0DMesh2DUpdated,
                                                               const vector<unsigned int>& cell1DMesh2DUpdated,
                                                               vector<unsigned int>& cell2DMesh2DNewVertices,
                                                               vector<unsigned int>& cell2DMesh2DNewEdges)
  {
    cell2DMesh2DNewVertices.resize(mesh2DReader.Cell2DNumberVertices(cell2DMesh2DId) + cell1DMesh2DUpdated.size() - 1);
    cell2DMesh2DNewEdges.resize(mesh2DReader.Cell2DNumberEdges(cell2DMesh2DId) + cell1DMesh2DUpdated.size() - 1);

    bool edgeListDirection = true; // true if forward, false if backward

    unsigned int vertexCounter = 0;
    unsigned int oldVertexCounter = 0;
    while (oldVertexCounter < mesh2DReader.Cell2DNumberVertices(cell2DMesh2DId))
    {
      const unsigned int oldVertex = mesh2DReader.Cell2DVertex(cell2DMesh2DId,
                                                               oldVertexCounter);
      const unsigned int oldNextVertex = mesh2DReader.Cell2DVertex(cell2DMesh2DId,
                                                                   (oldVertexCounter + 1) % mesh2DReader.Cell2DNumberVertices(cell2DMesh2DId));
      if (oldVertex == cell0DMesh2DUpdated[0] &&
          oldNextVertex == cell0DMesh2DUpdated[cell0DMesh2DUpdated.size() - 1])
      {
        edgeListDirection = true;
        for (unsigned int nv = 0; nv < cell0DMesh2DUpdated.size() - 1; nv++)
          cell2DMesh2DNewVertices[vertexCounter++] = cell0DMesh2DUpdated[nv];
      }
      else if (oldVertex == cell0DMesh2DUpdated[cell0DMesh2DUpdated.size() - 1] &&
               oldNextVertex == cell0DMesh2DUpdated[0])
      {
        edgeListDirection = false;
        for (unsigned int nv = cell0DMesh2DUpdated.size() - 1; nv > 0; nv--)
          cell2DMesh2DNewVertices[vertexCounter++] = cell0DMesh2DUpdated[nv];
      }
      else
        cell2DMesh2DNewVertices[vertexCounter++] = oldVertex;

      oldVertexCounter++;
    }

    unsigned int edgeCounter = 0;
    unsigned int oldEdgeCounter = 0;
    while (oldEdgeCounter < mesh2DReader.Cell2DNumberEdges(cell2DMesh2DId))
    {
      const unsigned int oldEdge = mesh2DReader.Cell2DEdge(cell2DMesh2DId,
                                                           oldEdgeCounter);
      if (oldEdge == oldCell1DMesh2DId)
      {
        if (edgeListDirection)
        {
          for (unsigned int ne = 0; ne < cell1DMesh2DUpdated.size(); ne++)
            cell2DMesh2DNewEdges[edgeCounter++] = cell1DMesh2DUpdated[ne];
        }
        else
        {
          for (unsigned int ne = cell1DMesh2DUpdated.size(); ne > 0; ne--)
            cell2DMesh2DNewEdges[edgeCounter++] = cell1DMesh2DUpdated[ne - 1];
        }
      }
      else
        cell2DMesh2DNewEdges[edgeCounter++] = oldEdge;

      oldEdgeCounter++;
    }
  }
  // ***************************************************************************
  void ConformerMeshPolygon::Cell2DMesh2DToMaps(const Gedim::IMeshDAO& mesh2DReader,
                                                const unsigned int& cell2DMesh2DId,
                                                map<unsigned int, unsigned int>& verticesMap,
                                                map<unsigned int, unsigned int>& edgesMap)
  {
    for (unsigned int v = 0; v < mesh2DReader.Cell2DNumberVertices(cell2DMesh2DId); v++)
    {
      verticesMap.insert(pair<unsigned int, unsigned int>(mesh2DReader.Cell2DVertex(cell2DMesh2DId, v), v));
      edgesMap.insert(pair<unsigned int, unsigned int>(mesh2DReader.Cell2DEdge(cell2DMesh2DId, v), v));
    }

    // if (cell2DMesh2DId == 10)
    // {
    //   cerr<< "INPUT "<< cell2DMesh2DId<< ": "<< endl;
    //   for (const auto& vertex : verticesMap)
    //     cerr<< vertex.first<< " - "<< vertex.second<< endl;
    //   for (const auto& edge : edgesMap)
    //     cerr<< edge.first<< " ("<< mesh2DReader.Cell1DOrigin(edge.first)<< ", "<< mesh2DReader.Cell1DEnd(edge.first)<< ")"<< " - "<< edge.second<< endl;
    // }
  }
  // ***************************************************************************
  void ConformerMeshPolygon::Cell2DMesh2DToSplitInput(const list<unsigned int> cell1DMesh1DIds,
                                                      const ConformerMeshSegment::ConformMesh& mesh1D,
                                                      const Gedim::IMeshDAO& mesh2D,
                                                      const unsigned int& cell2DMesh2DId,
                                                      const map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                                                      const map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                                                      GeometryUtilities::SplitPolygonInput& splitInput)
  {
    const ConformerMeshSegment::ConformMesh::ConformMeshSegment& firstCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.front()];
    const ConformerMeshSegment::ConformMesh::ConformMeshSegment& lastCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.back()];

    const double originCurvilinearCoordinate = firstCell1DMesh1D.Points[0];
    const double endCurvilinearCoordinate = lastCell1DMesh1D.Points[1];
    const ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points.at(originCurvilinearCoordinate);
    const ConformerMeshSegment::ConformMesh::ConformMeshPoint& end = mesh1D.Points.at(endCurvilinearCoordinate);

    splitInput.NumberPolygonVertices = mesh2D.Cell2DNumberVertices(cell2DMesh2DId);
    if (origin.Vertex2DIds.size() == 1)
    {
      splitInput.Segment.Origin.Type = GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
      splitInput.Segment.Origin.Index = cell2DMesh2DVerticesMap.at(origin.Vertex2DIds.front());
    }
    else if (origin.Vertex2DIds.size() == 0)
    {
      Output::Assert(origin.Edge2DIds.size() > 0);
      splitInput.Segment.Origin.Type = GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
      splitInput.Segment.Origin.Index = cell2DMesh2DEdgesMap.at(origin.Edge2DIds.front());
    }
    else
      throw runtime_error("cell1Ds origin of mesh1D are not correct");

    if (end.Vertex2DIds.size() == 1)
    {
      splitInput.Segment.End.Type = GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
      splitInput.Segment.End.Index = cell2DMesh2DVerticesMap.at(end.Vertex2DIds.front());
    }
    else if (end.Vertex2DIds.size() == 0)
    {
      Output::Assert(end.Edge2DIds.size() > 0);
      splitInput.Segment.End.Type = GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
      splitInput.Segment.End.Index = cell2DMesh2DEdgesMap.at(end.Edge2DIds.front());
    }
    else
      throw runtime_error("cell1Ds end of mesh1D are not correct");
  }
  // ***************************************************************************
  void ConformerMeshPolygon::SplitCell2DMesh2D(const Vector3d& segmentOrigin,
                                               const Vector3d& segmentTangent,
                                               const list<unsigned int> cell1DMesh1DIds,
                                               ConformerMeshSegment::ConformMesh& mesh1D,
                                               Gedim::IMeshDAO& mesh2D,
                                               const unsigned int& cell2DMesh2DId,
                                               map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                                               map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                                               const GeometryUtilities::SplitPolygonResult& splitResult)
  {
    ConformerMeshSegment::ConformMesh::ConformMeshSegment& firstCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.front()];
    ConformerMeshSegment::ConformMesh::ConformMeshSegment& lastCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.back()];

    const double originCurvilinearCoordinate = firstCell1DMesh1D.Points[0];
    const double endCurvilinearCoordinate = lastCell1DMesh1D.Points[1];

    ConformerMeshSegment::ConformMesh::ConformMeshPoint& originCell0DMesh1D = mesh1D.Points[originCurvilinearCoordinate];
    ConformerMeshSegment::ConformMesh::ConformMeshPoint& endCell0DMesh1D = mesh1D.Points[endCurvilinearCoordinate];

    GeometryUtilities::CompareTypes segment1DLenght = _geometryUtilities.Compare1DValues(originCurvilinearCoordinate,
                                                                                         endCurvilinearCoordinate);

    if (segment1DLenght == GeometryUtilities::CompareTypes::Coincident)
      throw runtime_error("Mesh1D cell1d length under tolerance");

    const unsigned int numOriginalVerticesCell2DMesh2D = mesh2D.Cell2DNumberVertices(cell2DMesh2DId);

    // if (cell2DMesh2DId == 10)
    // {
    //   cerr<< "Splitting cell "<< cell2DMesh2DId<< " in mesh "<< mesh2D.Cell0DTotalNumber()<< "x"<< mesh2D.Cell1DTotalNumber()<< "x"<< mesh2D.Cell2DTotalNumber()<< endl;
    //   cerr<< "Mesh1D segments "<< cell1DMesh1DIds<< endl;
    // }

    // add new cell0D to mesh2D
    unsigned int v = mesh2D.Cell0DAppend(splitResult.NewVertices.size());

    unsigned int nv = numOriginalVerticesCell2DMesh2D;
    map<unsigned int, unsigned int> newVertexMap;
    for (const GeometryUtilities::SplitPolygonResult::NewVertex& newVertex : splitResult.NewVertices)
    {
      if (newVertex.Type == GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin)
      {
        const Vector3d newCell0DCoordinates = segmentOrigin + originCurvilinearCoordinate * segmentTangent;
        mesh2D.Cell0DInsertCoordinates(v, newCell0DCoordinates);
        mesh2D.Cell0DSetMarker(v, 0);
        mesh2D.Cell0DSetState(v, true);

        originCell0DMesh1D.Vertex2DIds.push_back(v);
      }
      else if (newVertex.Type == GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentEnd)
      {
        const Vector3d newCell0DCoordinates = segmentOrigin + endCurvilinearCoordinate * segmentTangent;
        mesh2D.Cell0DInsertCoordinates(v, newCell0DCoordinates);
        mesh2D.Cell0DSetMarker(v, 0);
        mesh2D.Cell0DSetState(v, true);

        endCell0DMesh1D.Vertex2DIds.push_back(v);
      }
      else
        throw runtime_error("New cell0d of mesh2D not recognized");

      cell2DMesh2DVerticesMap.insert(pair<unsigned int, unsigned int>(v, nv));
      newVertexMap.insert(pair<unsigned int, unsigned int>(nv, v));
      nv++;
      v++;
    }

    // add new cell1D and cell2D to mesh2D
    unsigned int e = mesh2D.Cell1DAppend(splitResult.NewEdges.size());
    unsigned int c = mesh2D.Cell2DAppend(splitResult.NewPolygons.size());

    // if (cell2DMesh2DId == 10)
    // {
    //   cerr<< "Splitted in "<< splitResult.NewPolygons.size()<< " polygons: "<< endl;
    //   for (unsigned int pt = 0; pt < splitResult.NewPolygons.size(); pt++)
    //   {
    //     cerr<< "Vertices: "<< splitResult.NewPolygons[pt].Vertices<< endl;
    //     cerr<< "Edges: "<< splitResult.NewPolygons[pt].Edges<< endl;
    //   }
    //
    //   cerr<< "New edges: "<< splitResult.NewEdges.size()<< endl;
    //   for (const GeometryUtilities::SplitPolygonResult::NewEdge& edget : splitResult.NewEdges)
    //   {
    //     cerr<< " Type: "<< edget.Type<< " OldEdge: "<< edget.OldEdgeId<< " edge: "<< edget.OriginId<< "-"<< edget.EndId<< endl;
    //   }
    //
    //   cerr<< "New vertices: "<< splitResult.NewVertices.size()<< endl;
    //   for (const GeometryUtilities::SplitPolygonResult::NewVertex& vertext : splitResult.NewVertices)
    //   {
    //     cerr<< "Type: "<< vertext.Type<< endl;
    //   }
    // }

    unsigned int ne = numOriginalVerticesCell2DMesh2D;
    map<unsigned int, unsigned int> newEdgeMap;
    for (const GeometryUtilities::SplitPolygonResult::NewEdge& newEdge : splitResult.NewEdges)
    {
      const unsigned int cell1DOrigin = (newEdge.OriginId < numOriginalVerticesCell2DMesh2D) ? mesh2D.Cell2DVertex(cell2DMesh2DId,
                                                                                                                   newEdge.OriginId) :
                                                                                               newVertexMap.at(newEdge.OriginId);
      const unsigned int cell1DEnd = (newEdge.EndId < numOriginalVerticesCell2DMesh2D) ? mesh2D.Cell2DVertex(cell2DMesh2DId,
                                                                                                             newEdge.EndId) :
                                                                                         newVertexMap.at(newEdge.EndId);

      mesh2D.Cell1DInsertExtremes(e, cell1DOrigin, cell1DEnd);

      mesh2D.Cell1DSetMarker(e, 0);
      mesh2D.Cell1DSetState(e, true);

      // update conform mesh 1D
      if (originCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DOrigin(e) ||
          originCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DEnd(e))
      {
        originCell0DMesh1D.Edge2DIds.push_back(e);
      }

      if (endCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DOrigin(e) ||
          endCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DEnd(e))
      {
        endCell0DMesh1D.Edge2DIds.push_back(e);
      }

      // update mesh2D
      if (newEdge.Type == GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate)
      {
        const unsigned int edgeIdToUpdate = mesh2D.Cell2DEdge(cell2DMesh2DId,
                                                              newEdge.OldEdgeId);

        mesh2D.Cell1DInsertUpdatedCell1D(edgeIdToUpdate, e);
        mesh2D.Cell1DSetState(edgeIdToUpdate, 0);

        mesh2D.Cell1DSetMarker(e, mesh2D.Cell1DMarker(edgeIdToUpdate));

        if (newEdge.OriginId >= numOriginalVerticesCell2DMesh2D)
          mesh2D.Cell0DSetMarker(newVertexMap.at(newEdge.OriginId),
                                 mesh2D.Cell1DMarker(edgeIdToUpdate));
        if (newEdge.EndId >= numOriginalVerticesCell2DMesh2D)
          mesh2D.Cell0DSetMarker(newVertexMap.at(newEdge.EndId),
                                 mesh2D.Cell1DMarker(edgeIdToUpdate));

        mesh2D.Cell1DInitializeNeighbourCell2Ds(e,
                                                mesh2D.Cell1DNumberNeighbourCell2D(edgeIdToUpdate));

        for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(edgeIdToUpdate); n++)
        {
          if (!mesh2D.Cell1DHasNeighbourCell2D(edgeIdToUpdate, n))
            continue;

          if (mesh2D.Cell1DNeighbourCell2D(edgeIdToUpdate, n) == cell2DMesh2DId)
            mesh2D.Cell1DInsertNeighbourCell2D(e,
                                               n,
                                               c + newEdge.Cell2DNeighbours[0]);
          else
            mesh2D.Cell1DInsertNeighbourCell2D(e,
                                               n,
                                               mesh2D.Cell1DNeighbourCell2D(edgeIdToUpdate, n));
        }

        if (splitResult.Type == GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate)
        {
          if ((originCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DOrigin(e) ||
               originCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DEnd(e)) &&
              (endCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DOrigin(e) ||
               endCell0DMesh1D.Vertex2DIds.front() == mesh2D.Cell1DEnd(e)))
          {
            for (const unsigned int& mesh1DCell1DId : cell1DMesh1DIds)
              mesh1D.Segments[mesh1DCell1DId].Edge2DIds.push_back(e);
          }
        }
      }
      else if (newEdge.Type == GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew)
      {
        mesh2D.Cell1DInitializeNeighbourCell2Ds(e,
                                                2);

        for (unsigned int n = 0; n < 2; n++)
        {
          mesh2D.Cell1DInsertNeighbourCell2D(e,
                                             newEdge.Cell2DNeighbours[n],
                                             c + newEdge.Cell2DNeighbours[n]);
        }

        for (const unsigned int& mesh1DCell1DId : cell1DMesh1DIds)
          mesh1D.Segments[mesh1DCell1DId].Edge2DIds.push_back(e);
      }

      cell2DMesh2DEdgesMap.insert(pair<unsigned int, unsigned int>(e, ne));
      newEdgeMap.insert(pair<unsigned int, unsigned int>(ne, e));
      ne++;
      e++;
    }

    // if (cell2DMesh2DId == 10)
    // {
    //   cerr<< "Splitting cell "<< cell2DMesh2DId<< " with "<< originCurvilinearCoordinate<< "->"<<endCurvilinearCoordinate<< endl;
    //
    //   for (const auto& newVertex : newVertexMap)
    //     cerr<< "New vertex "<< newVertex.second<< endl;
    //   for (const auto& newEdge : newEdgeMap)
    //     cerr<< "New edge "<< newEdge.second<< " "<< mesh2D.Cell1DOrigin(newEdge.second)<< "x"<< mesh2D.Cell1DEnd(newEdge.second)<< endl;
    //
    //   unsigned int cTemp = mesh2D.Cell2DTotalNumber();
    //   for (const auto& newPolygon : splitResult.NewPolygons)
    //   {
    //     cerr<< "Polygon "<< cTemp<< " with "<< newPolygon.Vertices.size()<< " vertices and "<< newPolygon.Edges.size()<< " edges"<< endl;
    //     cerr<< "{";
    //     for (const unsigned int& vId : newPolygon.Vertices)
    //       cerr<< ((vId < numOriginalVerticesCell2DMesh2D) ? mesh2D.Cell2DVertex(cell2DMesh2DId, vId) : newVertexMap.at(vId))<< ", ";
    //     cerr<< "}"<< endl;
    //
    //     cerr<< "{";
    //     for (const unsigned int& eId : newPolygon.Edges)
    //       cerr<< ((eId < numOriginalVerticesCell2DMesh2D) ? mesh2D.Cell2DEdge(cell2DMesh2DId, eId) : newEdgeMap.at(eId))<< ", ";
    //     cerr<< "}"<< endl;
    //
    //     cTemp++;
    //   }
    // }

    // insert new cell2D data in mesh2D
    unsigned int nc = 0;
    for (const GeometryUtilities::SplitPolygonResult::NewPolygon& newPolygon : splitResult.NewPolygons)
    {
      mesh2D.Cell2DInitializeVertices(c, newPolygon.Vertices.size());
      mesh2D.Cell2DInitializeEdges(c, newPolygon.Edges.size());

      unsigned int cv = 0;
      for (const unsigned int& vId : newPolygon.Vertices)
      {
        // update mesh2D
        mesh2D.Cell2DInsertVertex(c,
                                  cv,
                                  (vId < numOriginalVerticesCell2DMesh2D) ? mesh2D.Cell2DVertex(cell2DMesh2DId,
                                                                                                vId) :
                                                                            newVertexMap.at(vId));
        cv++;
      }

      unsigned int ce = 0;
      for (const unsigned int& eId : newPolygon.Edges)
      {
        if (eId < numOriginalVerticesCell2DMesh2D)
        {
          const unsigned int edgeId = mesh2D.Cell2DEdge(cell2DMesh2DId, eId);
          mesh2D.Cell2DInsertEdge(c,
                                  ce,
                                  edgeId);

          for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(edgeId); n++)
          {
            if (!mesh2D.Cell1DHasNeighbourCell2D(edgeId, n))
              continue;

            if (mesh2D.Cell1DNeighbourCell2D(edgeId, n) == cell2DMesh2DId)
              mesh2D.Cell1DInsertNeighbourCell2D(edgeId,
                                                 n,
                                                 c);
          }
        }
        else
        {
          // update mesh2D
          mesh2D.Cell2DInsertEdge(c,
                                  ce,
                                  newEdgeMap.at(eId));
        }

        ce++;
      }

      // update conform mesh 1D
      for (const unsigned int& mesh1DCell1DId : cell1DMesh1DIds)
      {
        ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[mesh1DCell1DId];
        ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points.at(segment.Points[0]);

        origin.Cell2DIds.push_back(c);
        segment.Cell2DIds.push_back(c);
      }
      endCell0DMesh1D.Cell2DIds.push_back(c);

      mesh2D.Cell2DSetMarker(c, mesh2D.Cell2DMarker(cell2DMesh2DId));
      mesh2D.Cell2DSetState(c, true);
      mesh2D.Cell2DSetState(cell2DMesh2DId, false);

      mesh2D.Cell2DInsertUpdatedCell2D(cell2DMesh2DId, c);

      nc++;
      c++;
    }
  }
  // ***************************************************************************
  void ConformerMeshPolygon::InsertCell2DMesh2DMiddleEdgesPolygonUpdate(const Vector3d& segmentOrigin,
                                                                        const Vector3d& segmentTangent,
                                                                        ConformerMeshSegment::ConformMesh& mesh1D,
                                                                        Gedim::IMeshDAO& mesh2D,
                                                                        const list<unsigned int> cell1DMesh1DIds,
                                                                        const unsigned int& cell2DMesh2DId)
  {
    // nothing to do
    if (cell1DMesh1DIds.size() < 2)
      return;

    ConformerMeshSegment::ConformMesh::ConformMeshSegment& firstCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.front()];
    ConformerMeshSegment::ConformMesh::ConformMeshSegment& lastCell1DMesh1D = mesh1D.Segments[cell1DMesh1DIds.back()];

    const double originCurvilinearCoordinate = firstCell1DMesh1D.Points[0];
    const double endCurvilinearCoordinate = lastCell1DMesh1D.Points[1];

    ConformerMeshSegment::ConformMesh::ConformMeshPoint& originCell0DMesh1D = mesh1D.Points[originCurvilinearCoordinate];
    ConformerMeshSegment::ConformMesh::ConformMeshPoint& endCell0DMesh1D = mesh1D.Points[endCurvilinearCoordinate];

    list<unsigned int> cell2DMesh2DsToUpdate;
    mesh2D.Cell2DUpdatedCell2Ds(cell2DMesh2DId,
                                cell2DMesh2DsToUpdate);

    Output::Assert(cell2DMesh2DsToUpdate.size() == 1);

    unsigned int cell1DMesh2DToUpdate = mesh1D.Segments[cell1DMesh1DIds.front()].Edge2DIds.back();

    Output::Assert(!mesh2D.Cell1DHasUpdatedCell1Ds(cell1DMesh2DToUpdate));

    const unsigned int numberNewCell0DMesh2D = cell1DMesh1DIds.size() - 1;
    const unsigned int numberNewCell1DMesh2D = cell1DMesh1DIds.size();
    const unsigned int numberNewCell2DMesh2D = 1;

    vector<unsigned int> newCell0DMesh2Ds(numberNewCell0DMesh2D + 2);
    vector<unsigned int> newCell1DMesh2Ds(numberNewCell1DMesh2D);
    vector<unsigned int> newCell2DMesh2Ds(numberNewCell2DMesh2D);

    newCell0DMesh2Ds[0] = mesh2D.Cell1DOrigin(cell1DMesh2DToUpdate);
    newCell0DMesh2Ds[numberNewCell0DMesh2D + 1] = mesh2D.Cell1DEnd(cell1DMesh2DToUpdate);

    // Insert new Cell0Ds in Mesh2D
    unsigned int v = mesh2D.Cell0DAppend(numberNewCell0DMesh2D);

    unsigned int counter = 0;
    for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
    {
      if (counter == 0)
      {
        counter++;
        continue;
      }

      ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
      if (segment.Edge2DIds.back() != cell1DMesh2DToUpdate)
        throw runtime_error("Something goes wrong with update of " + to_string(cell1DMesh2DToUpdate) + " cell1DMesh2D in cell1DMesh1D " +to_string(cell1DMesh1DId));

      // insert only Cell1DMesh1D origin
      const double& originCurvilinearCoordinate = segment.Points[0];

      // add new Cell0Ds on Mesh2D
      const Vector3d newCell0DCoordinates = segmentOrigin + originCurvilinearCoordinate * segmentTangent;
      mesh2D.Cell0DInsertCoordinates(v, newCell0DCoordinates);
      mesh2D.Cell0DSetMarker(v, mesh2D.Cell1DMarker(cell1DMesh2DToUpdate));
      mesh2D.Cell0DSetState(v, true);

      // update mesh 1D
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];
      origin.Vertex2DIds.push_back(v);

      newCell0DMesh2Ds[counter] = v;

      counter++;
      v++;
    }

    // Insert new Cell1Ds and Cell2D in Mesh2D
    unsigned int e = mesh2D.Cell1DAppend(numberNewCell1DMesh2D);
    unsigned int c = mesh2D.Cell2DAppend(numberNewCell2DMesh2D);

    counter  = 0;
    for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
    {
      ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
      if (segment.Edge2DIds.back() != cell1DMesh2DToUpdate)
        throw runtime_error("Something goes wrong with update of " + to_string(cell1DMesh2DToUpdate) + " cell1DMesh2D in cell1DMesh1D " +to_string(cell1DMesh1DId));

      // add new Cell1Ds on Mesh2D
      mesh2D.Cell1DInsertExtremes(e,
                                  newCell0DMesh2Ds[counter],
                                  newCell0DMesh2Ds[counter + 1]);
      mesh2D.Cell1DSetState(e, true);

      // update conform mesh 1D
      const double& originCurvilinearCoordinate = segment.Points[0];
      const double& endCurvilinearCoordinate = segment.Points[1];
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& end = mesh1D.Points[endCurvilinearCoordinate];

      origin.Edge2DIds.push_back(e);
      end.Edge2DIds.push_back(e);
      segment.Edge2DIds.push_back(e);

      // update mesh2D
      mesh2D.Cell1DSetState(cell1DMesh2DToUpdate, false);
      mesh2D.Cell1DInsertUpdatedCell1D(cell1DMesh2DToUpdate, e);

      mesh2D.Cell1DSetMarker(e, mesh2D.Cell1DMarker(cell1DMesh2DToUpdate));

      mesh2D.Cell1DInitializeNeighbourCell2Ds(e,
                                              mesh2D.Cell1DNumberNeighbourCell2D(cell1DMesh2DToUpdate));

      for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(cell1DMesh2DToUpdate); n++)
      {
        if (!mesh2D.Cell1DHasNeighbourCell2D(cell1DMesh2DToUpdate, n))
          continue;

        mesh2D.Cell1DInsertNeighbourCell2D(e,
                                           n,
                                           mesh2D.Cell1DNeighbourCell2D(cell1DMesh2DToUpdate, n));
      }

      newCell1DMesh2Ds[counter] = e;
      counter++;
      e++;
    }

    // add new cell2D to mesh2D
    counter = 0;
    for (const unsigned int& cell2DMesh2DToUpdate : cell2DMesh2DsToUpdate)
    {
      vector<unsigned int> newCell2DVertices, newCell2DEdges;
      Cell2DMesh2DUpdatedPointsAndEdges(mesh2D,
                                        cell2DMesh2DToUpdate,
                                        cell1DMesh2DToUpdate,
                                        newCell0DMesh2Ds,
                                        newCell1DMesh2Ds,
                                        newCell2DVertices,
                                        newCell2DEdges);

      // cerr<< "Update Cell2D "<< cell2DMesh2DToUpdate<< " with new Cell2D "<< c<< " with "<<
      //        newCell2DVertices<< " vertices and "<<
      //        newCell2DEdges<< " edges"<<
      //        endl;

      mesh2D.Cell2DAddVertices(c, newCell2DVertices);
      mesh2D.Cell2DAddEdges(c, newCell2DEdges);

      mesh2D.Cell2DSetMarker(c, mesh2D.Cell2DMarker(cell2DMesh2DToUpdate));
      mesh2D.Cell2DSetState(c, true);
      mesh2D.Cell2DSetState(cell2DMesh2DToUpdate, false);

      mesh2D.Cell2DInsertUpdatedCell2D(cell2DMesh2DToUpdate, c);

      // update all edges neighbours
      for (unsigned int e = 0; e < mesh2D.Cell2DNumberEdges(c); e++)
      {
        const unsigned int edgeId = mesh2D.Cell2DEdge(c, e);
        for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(edgeId); n++)
        {
          if (!mesh2D.Cell1DHasNeighbourCell2D(edgeId, n))
            continue;

          if (mesh2D.Cell1DNeighbourCell2D(edgeId, n) == cell2DMesh2DToUpdate)
            mesh2D.Cell1DInsertNeighbourCell2D(edgeId,
                                               n,
                                               c);
        }
      }

      // update conform mesh 1D
      for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
      {
        ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
        const double& originCurvilinearCoordinate = segment.Points[0];
        ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];

        segment.Cell2DIds.push_back(c);
        origin.Cell2DIds.push_back(c);
      }
      mesh1D.Points[mesh1D.Segments[cell1DMesh1DIds.back()].Points[1]].Cell2DIds.push_back(c);

      newCell2DMesh2Ds[counter] = c;

      counter++;
      c++;
    }

    // cerr<< "Cell1DMesh2D "<< cell1DMesh2DToUpdate<< " updated with "<< newCell0DMesh2Ds<< " vertices and "<< newCell1DMesh2Ds<< " edges and " << newCell2DMesh2Ds<< " cells"<< endl;
  }
  // ***************************************************************************
  void ConformerMeshPolygon::InsertCell2DMesh2DMiddleEdgesPolygonCreation(const Vector3d& segmentOrigin,
                                                                          const Vector3d& segmentTangent,
                                                                          ConformerMeshSegment::ConformMesh& mesh1D,
                                                                          Gedim::IMeshDAO& mesh2D,
                                                                          const list<unsigned int> cell1DMesh1DIds,
                                                                          const unsigned int& cell2DMesh2DId)
  {
    // nothing to do
    if (cell1DMesh1DIds.size() < 2)
      return;

    list<unsigned int> cell2DMesh2DsToUpdate;
    mesh2D.Cell2DUpdatedCell2Ds(cell2DMesh2DId,
                                cell2DMesh2DsToUpdate);

    unsigned int cell1DMesh2DToUpdate = mesh1D.Segments[cell1DMesh1DIds.front()].Edge2DIds.front();

    Output::Assert(!mesh2D.Cell1DHasUpdatedCell1Ds(cell1DMesh2DToUpdate));

    const unsigned int numberNewCell0DMesh2D = cell1DMesh1DIds.size() - 1;
    const unsigned int numberNewCell1DMesh2D = cell1DMesh1DIds.size();
    const unsigned int numberNewCell2DMesh2D = cell2DMesh2DsToUpdate.size();

    // cerr<< "Insert "<< numberNewCell1DMesh2D<< " middle edges in old cell2D "<< cell2DMesh2DId<< " updated with "<< cell2DMesh2DsToUpdate<< " new cell2Ds"<< endl;

    vector<unsigned int> newCell0DMesh2Ds(numberNewCell0DMesh2D + 2);
    vector<unsigned int> newCell1DMesh2Ds(numberNewCell1DMesh2D);
    vector<unsigned int> newCell2DMesh2Ds(numberNewCell2DMesh2D);

    newCell0DMesh2Ds[0] = mesh2D.Cell1DOrigin(cell1DMesh2DToUpdate);
    newCell0DMesh2Ds[numberNewCell0DMesh2D + 1] = mesh2D.Cell1DEnd(cell1DMesh2DToUpdate);

    // Insert new Cell0Ds in Mesh2D
    unsigned int v = mesh2D.Cell0DAppend(numberNewCell0DMesh2D);

    unsigned int counter = 0;
    for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
    {
      if (counter == 0)
      {
        counter++;
        continue;
      }

      ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
      if (segment.Edge2DIds.front() != cell1DMesh2DToUpdate)
        throw runtime_error("Something goes wrong with update of " + to_string(cell1DMesh2DToUpdate) + " cell1DMesh2D in cell1DMesh1D " +to_string(cell1DMesh1DId));

      // insert only Cell1DMesh1D origin
      const double& originCurvilinearCoordinate = segment.Points[0];

      // add new Cell0Ds on Mesh2D
      const Vector3d newCell0DCoordinates = segmentOrigin + originCurvilinearCoordinate * segmentTangent;
      mesh2D.Cell0DInsertCoordinates(v, newCell0DCoordinates);
      mesh2D.Cell0DSetMarker(v, mesh2D.Cell1DMarker(cell1DMesh2DToUpdate));
      mesh2D.Cell0DSetState(v, true);

      // update mesh 1D
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];
      origin.Vertex2DIds.push_back(v);

      newCell0DMesh2Ds[counter] = v;

      counter++;
      v++;
    }

    // Insert new Cell1Ds and Cell2D in Mesh2D
    unsigned int e = mesh2D.Cell1DAppend(numberNewCell1DMesh2D);
    unsigned int c = mesh2D.Cell2DAppend(numberNewCell2DMesh2D);

    counter  = 0;
    for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
    {
      ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
      if (segment.Edge2DIds.front() != cell1DMesh2DToUpdate)
        throw runtime_error("Something goes wrong with update of " + to_string(cell1DMesh2DToUpdate) + " cell1DMesh2D in cell1DMesh1D " +to_string(cell1DMesh1DId));

      // add new Cell1Ds on Mesh2D
      mesh2D.Cell1DInsertExtremes(e,
                                  newCell0DMesh2Ds[counter],
                                  newCell0DMesh2Ds[counter + 1]);
      mesh2D.Cell1DSetState(e, true);

      // update conform mesh 1D
      const double& originCurvilinearCoordinate = segment.Points[0];
      const double& endCurvilinearCoordinate = segment.Points[1];
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];
      ConformerMeshSegment::ConformMesh::ConformMeshPoint& end = mesh1D.Points[endCurvilinearCoordinate];

      origin.Edge2DIds.push_back(e);
      end.Edge2DIds.push_back(e);
      segment.Edge2DIds.push_back(e);

      // update mesh2D
      mesh2D.Cell1DSetState(cell1DMesh2DToUpdate, false);
      mesh2D.Cell1DInsertUpdatedCell1D(cell1DMesh2DToUpdate, e);

      mesh2D.Cell1DSetMarker(e, mesh2D.Cell1DMarker(cell1DMesh2DToUpdate));

      mesh2D.Cell1DInitializeNeighbourCell2Ds(e,
                                              mesh2D.Cell1DNumberNeighbourCell2D(cell1DMesh2DToUpdate));

      for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(cell1DMesh2DToUpdate); n++)
      {
        if (!mesh2D.Cell1DHasNeighbourCell2D(cell1DMesh2DToUpdate, n))
          continue;

        mesh2D.Cell1DInsertNeighbourCell2D(e,
                                           n,
                                           mesh2D.Cell1DNeighbourCell2D(cell1DMesh2DToUpdate, n));
      }

      newCell1DMesh2Ds[counter] = e;
      counter++;
      e++;
    }

    // add new cell2D to mesh2D
    counter = 0;
    for (const unsigned int& cell2DMesh2DToUpdate : cell2DMesh2DsToUpdate)
    {
      vector<unsigned int> newCell2DVertices, newCell2DEdges;
      Cell2DMesh2DUpdatedPointsAndEdges(mesh2D,
                                        cell2DMesh2DToUpdate,
                                        cell1DMesh2DToUpdate,
                                        newCell0DMesh2Ds,
                                        newCell1DMesh2Ds,
                                        newCell2DVertices,
                                        newCell2DEdges);

      // cerr<< "Update Cell2D "<< cell2DMesh2DToUpdate<< " with new Cell2D "<< c<< " with "<<
      //        newCell2DVertices<< " vertices and "<<
      //        newCell2DEdges<< " edges"<<
      //        endl;

      mesh2D.Cell2DAddVertices(c, newCell2DVertices);
      mesh2D.Cell2DAddEdges(c, newCell2DEdges);

      mesh2D.Cell2DSetMarker(c, mesh2D.Cell2DMarker(cell2DMesh2DToUpdate));
      mesh2D.Cell2DSetState(c, true);
      mesh2D.Cell2DSetState(cell2DMesh2DToUpdate, false);

      mesh2D.Cell2DInsertUpdatedCell2D(cell2DMesh2DToUpdate, c);

      // update all edges neighbours
      for (unsigned int e = 0; e < mesh2D.Cell2DNumberEdges(c); e++)
      {
        const unsigned int edgeId = mesh2D.Cell2DEdge(c, e);
        for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(edgeId); n++)
        {
          if (!mesh2D.Cell1DHasNeighbourCell2D(edgeId, n))
            continue;

          if (mesh2D.Cell1DNeighbourCell2D(edgeId, n) == cell2DMesh2DToUpdate)
            mesh2D.Cell1DInsertNeighbourCell2D(edgeId,
                                               n,
                                               c);
        }
      }

      // update conform mesh 1D
      for (const unsigned int& cell1DMesh1DId : cell1DMesh1DIds)
      {
        ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[cell1DMesh1DId];
        const double& originCurvilinearCoordinate = segment.Points[0];
        ConformerMeshSegment::ConformMesh::ConformMeshPoint& origin = mesh1D.Points[originCurvilinearCoordinate];

        segment.Cell2DIds.push_back(c);
        origin.Cell2DIds.push_back(c);
      }
      mesh1D.Points[mesh1D.Segments[cell1DMesh1DIds.back()].Points[1]].Cell2DIds.push_back(c);

      newCell2DMesh2Ds[counter] = c;

      counter++;
      c++;
    }

    // cerr<< "Cell1DMesh2D "<< cell1DMesh2DToUpdate<< " updated with "<< newCell0DMesh2Ds<< " vertices and "<< newCell1DMesh2Ds<< " edges and " << newCell2DMesh2Ds<< " cells"<< endl;
  }
  // ***************************************************************************
  void ConformerMeshPolygon::UpdateCell2DNeighbours(ConformerMeshSegment::ConformMesh& mesh1D,
                                                    Gedim::IMeshDAO& mesh2D,
                                                    const unsigned int& cell2DMesh2DId)
  {
    for (unsigned int e = 0; e < mesh2D.Cell2DNumberEdges(cell2DMesh2DId); e++)
    {
      const unsigned int cell1DId = mesh2D.Cell2DEdge(cell2DMesh2DId, e);

      if (mesh2D.Cell1DIsActive(cell1DId))
        continue;

      for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(cell1DId); n++)
      {
        if (!mesh2D.Cell1DHasNeighbourCell2D(cell1DId, n))
          continue;

        if (mesh2D.Cell1DNeighbourCell2D(cell1DId, n) == cell2DMesh2DId)
          continue;

        const unsigned int cell2DNeigh = mesh2D.Cell1DNeighbourCell2D(cell1DId, n);

        // update neighbour cell with a new cell 2D
        list<unsigned int> newEdges;
        if (!mesh2D.Cell1DUpdatedCell1Ds(cell1DId,
                                         newEdges))
        {
          // edge is updated, update the neighbour cell2Ds
          // cerr<< "Cell neigh "<< cell2DNeigh<< " Edge "<< cell1DId<< " with "<< newEdges<< endl;

          // get new points to insert between origin and end of old edge
          vector<unsigned int> newCell1DPoints, newCell1DEdges;
          Cell2DMesh2DListCell1DToListCell0D(mesh2D,
                                             cell1DId,
                                             newEdges,
                                             newCell1DPoints,
                                             newCell1DEdges);
          // cerr<< "Old edge "<< cell1DId<< " "<<
          //        mesh2DReader.Cell1DOrigin(cell1DId)<< "x"<<
          //        mesh2DReader.Cell1DEnd(cell1DId)<< " "<<
          //        "newPoints: "<< newCell1DPoints<< " "<<
          //        "newEdges: "<< newCell1DEdges<<
          //        endl;

          vector<unsigned int> newCell2DVertices, newCell2DEdges;
          Cell2DMesh2DUpdatedPointsAndEdges(mesh2D,
                                            cell2DNeigh,
                                            cell1DId,
                                            newCell1DPoints,
                                            newCell1DEdges,
                                            newCell2DVertices,
                                            newCell2DEdges);
          // cerr<< "New Cell "<< mesh2D.NumberCell2D<< " with "<<
          //        newCell2DVertices<< " vertices and "<<
          //        newCell2DEdges<< " edges"<<
          //        endl;

          // add new cell2D to mesh2D
          unsigned int c = mesh2D.Cell2DAppend(1);

          mesh2D.Cell2DAddVertices(c, newCell2DVertices);
          mesh2D.Cell2DAddEdges(c, newCell2DEdges);

          mesh2D.Cell2DSetMarker(c, mesh2D.Cell2DMarker(cell2DNeigh));
          mesh2D.Cell2DSetState(c, true);
          mesh2D.Cell2DSetState(cell2DNeigh, false);

          mesh2D.Cell2DInsertUpdatedCell2D(cell2DNeigh, c);

          // update all edges neighbours
          for (unsigned int e = 0; e < mesh2D.Cell2DNumberEdges(c); e++)
          {
            const unsigned int edgeId = mesh2D.Cell2DEdge(c, e);
            for (unsigned int n = 0; n < mesh2D.Cell1DNumberNeighbourCell2D(edgeId); n++)
            {
              if (!mesh2D.Cell1DHasNeighbourCell2D(edgeId, n))
                continue;

              if (mesh2D.Cell1DNeighbourCell2D(edgeId, n) == cell2DNeigh)
                mesh2D.Cell1DInsertNeighbourCell2D(edgeId,
                                                   n,
                                                   c);
            }
          }

          // update conform mesh 1D
          for (map<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint>::iterator it = mesh1D.Points.begin();
               it != mesh1D.Points.end(); it++)
          {
            ConformerMeshSegment::ConformMesh::ConformMeshPoint& point = it->second;
            if (find(point.Cell2DIds.begin(), point.Cell2DIds.end(), cell2DNeigh) == point.Cell2DIds.end())
              continue;

            point.Cell2DIds.push_back(c);
          }

          for (unsigned int e1D = 0; e1D < mesh1D.Segments.size(); e1D++)
          {
            ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment = mesh1D.Segments[e1D];
            if (find(segment.Cell2DIds.begin(), segment.Cell2DIds.end(), cell2DNeigh) == segment.Cell2DIds.end())
              continue;

            segment.Cell2DIds.push_back(c);
          }
        }
      }
    }
  }
  // ***************************************************************************
  void ConformerMeshPolygon::CreateConformMesh(const Vector3d& segmentOrigin,
                                               const Vector3d& segmentEnd,
                                               ConformerMeshSegment::ConformMesh& mesh1D,
                                               Gedim::IMeshDAO& mesh2D,
                                               ConformerMeshPolygon::ConformMesh& meshConformedInformation)
  {
    const Vector3d segmentTangent = segmentEnd - segmentOrigin;

    // check starting and end of segment if inside Cell2DMesh2D
    CheckSegmentOriginAndEndIntersections(segmentOrigin,
                                          segmentEnd,
                                          mesh1D,
                                          mesh2D);

    unsigned int numVisitedCell1DMesh1D = 0;
    while (numVisitedCell1DMesh1D < mesh1D.Segments.size())
    {
      const ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment1D = mesh1D.Segments[numVisitedCell1DMesh1D];
      Output::Assert(segment1D.Cell2DIds.size() > 0);

      unsigned int cell2DId = segment1D.Cell2DIds.back();
      list<unsigned int> cell1DMesh1DIds;
      do
      {
        cell1DMesh1DIds.push_back(numVisitedCell1DMesh1D++);
      }
      while (numVisitedCell1DMesh1D < mesh1D.Segments.size() &&
             find(mesh1D.Segments[numVisitedCell1DMesh1D].Cell2DIds.begin(),
                  mesh1D.Segments[numVisitedCell1DMesh1D].Cell2DIds.end(),
                  cell2DId) !=
             mesh1D.Segments[numVisitedCell1DMesh1D].Cell2DIds.end());

      // generate cell2D mesh2D maps
      map<unsigned int, unsigned int> cell2DMesh2DVerticesMap, cell2DMesh2DEdgesMap;
      Cell2DMesh2DToMaps(mesh2D,
                         cell2DId,
                         cell2DMesh2DVerticesMap,
                         cell2DMesh2DEdgesMap);

      // update cell2d of mesh2d
      GeometryUtilities::SplitPolygonInput splitInput;
      Cell2DMesh2DToSplitInput(cell1DMesh1DIds,
                               mesh1D,
                               mesh2D,
                               cell2DId,
                               cell2DMesh2DVerticesMap,
                               cell2DMesh2DEdgesMap,
                               splitInput);

      GeometryUtilities::SplitPolygonResult splitResult = _geometryUtilities.SplitPolygon(splitInput);

      if (splitResult.Type == GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate ||
          splitResult.Type == GeometryUtilities::SplitPolygonResult::Types::PolygonCreation)
      {
        // create new cells
        SplitCell2DMesh2D(segmentOrigin,
                          segmentTangent,
                          cell1DMesh1DIds,
                          mesh1D,
                          mesh2D,
                          cell2DId,
                          cell2DMesh2DVerticesMap,
                          cell2DMesh2DEdgesMap,
                          splitResult);
      }

      // insert middle edges in new cells
      switch (splitResult.Type) {
        case GeometryUtilities::SplitPolygonResult::Types::NoAction:
        case GeometryUtilities::SplitPolygonResult::Types::PolygonCreation:
          InsertCell2DMesh2DMiddleEdgesPolygonCreation(segmentOrigin,
                                                       segmentTangent,
                                                       mesh1D,
                                                       mesh2D,
                                                       cell1DMesh1DIds,
                                                       cell2DId);
        break;
        case GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate:
          InsertCell2DMesh2DMiddleEdgesPolygonUpdate(segmentOrigin,
                                                     segmentTangent,
                                                     mesh1D,
                                                     mesh2D,
                                                     cell1DMesh1DIds,
                                                     cell2DId);
        break;
        default:
          throw runtime_error("splitResult.Type not supported");
      }

      // update neighbour cells
      UpdateCell2DNeighbours(mesh1D,
                             mesh2D,
                             cell2DId);

      // for (const auto& cell1DnewIds : mesh2DConformed.Mesh2DUpdated.UpdatedCell1Ds)
      // {
      //   cerr<< "Update Edge "<< cell1DnewIds.first<< " with "<< cell1DnewIds.second<< endl;
      // }
      //
      // for (const auto& cell2DnewIds : mesh2DConformed.Mesh2DUpdated.UpdatedCell2Ds)
      // {
      //   cerr<< "Update Cell "<< cell2DnewIds.first<< " with "<< cell2DnewIds.second<< endl;
      // }
    }

    mesh2D.Compress();
  }
  // ***************************************************************************
}
