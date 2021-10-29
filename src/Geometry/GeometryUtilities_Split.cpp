#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::SplitPolygonWithSegmentResult GeometryUtilities::SplitPolygonWithSegment(const GeometryUtilities::SplitPolygonInput& input) const
  {
    GeometryUtilities::SplitPolygonWithSegmentResult result;

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
        result.Type = SplitPolygonWithSegmentResult::Types::NoAction;
        return result;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // No split needed
          result.Type = SplitPolygonWithSegmentResult::Types::NoAction;
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
        result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonWithSegmentResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonWithSegmentResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == endIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonWithSegmentResult::NewVertex> newVertexTypes;
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
            result.NewEdges.push_back(SplitPolygonWithSegmentResult::NewEdge());
            SplitPolygonWithSegmentResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin new end
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(end,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end new origin
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(origin,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else
              throw runtime_error("This case is not possibile because origin is vertex type");
          }
        }

        for (map<unsigned int, SplitPolygonWithSegmentResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonWithSegmentResult::NewVertex());
          SplitPolygonWithSegmentResult::NewVertex& newVertex = result.NewVertices.back();
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
        result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];
        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonWithSegmentResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonWithSegmentResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == originIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonWithSegmentResult::NewVertex> newVertexTypes;
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
            result.NewEdges.push_back(SplitPolygonWithSegmentResult::NewEdge());
            SplitPolygonWithSegmentResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(end,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;
              continue;

              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(origin,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;
              continue;
            }
            else
              throw runtime_error("This case is not possibile because end is vertex type");
          }
        }

        for (map<unsigned int, SplitPolygonWithSegmentResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonWithSegmentResult::NewVertex());
          SplitPolygonWithSegmentResult::NewVertex& newVertex = result.NewVertices.back();
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
        result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
      }

      // check contigous edges
      for (unsigned int c = 0; c < input.AlignedEdges.size(); c++)
      {
        SplitPolygonInput::AlignedEdge alignedEdge = input.AlignedEdges[c];

        if (originIndex >= alignedEdge.OriginVertexIndex && endIndex <= alignedEdge.EndVertexIndex)
        {
          // Only update needed
          result.Type = SplitPolygonWithSegmentResult::Types::PolygonUpdate;
          break;
        }
      }

      if (result.Type == SplitPolygonWithSegmentResult::Types::PolygonUpdate)
      {
        // Update polygon
        result.NewPolygons.resize(1);
        SplitPolygonWithSegmentResult::NewPolygon& updatedPolygon = result.NewPolygons[0];

        for (unsigned int v = 0; v < input.NumberPolygonVertices; v++)
        {
          updatedPolygon.Vertices.push_back(v);

          if (v == originIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices);

          if (v == endIndex)
            updatedPolygon.Vertices.push_back(input.NumberPolygonVertices + 1);
        }

        unsigned int newEdgeNumber = input.NumberPolygonVertices;
        map<unsigned int, SplitPolygonWithSegmentResult::NewVertex> newVertexTypes;
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
            result.NewEdges.push_back(SplitPolygonWithSegmentResult::NewEdge());
            SplitPolygonWithSegmentResult::NewEdge& newEdge = result.NewEdges.back();
            newEdge.OriginId = origin;
            newEdge.EndId = end;
            newEdge.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdge.Cell2DNeighbours = { 0 };

            if (origin < input.NumberPolygonVertices)
            {
              // old origin
              newEdge.OldEdgeId = origin;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(end,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[end];
              newVertex.Type = (end == 4) ? SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin :
                                            SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else if (end < input.NumberPolygonVertices)
            {
              // old end
              newEdge.OldEdgeId = (end == 0) ? input.NumberPolygonVertices - 1 : end - 1;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(origin,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newVertex = newVertexTypes[origin];
              newVertex.Type = (origin == 4) ? SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin :
                                               SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;
              continue;
            }
            else
            {
              // new origin and new end
              newEdge.OldEdgeId = input.Segment.Origin.Index;
              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(origin,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newOrigin = newVertexTypes[origin];
              newOrigin.Type = (origin == 4) ? SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin :
                                               SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;

              newVertexTypes.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(end,
                                                                                                SplitPolygonWithSegmentResult::NewVertex()));

              SplitPolygonWithSegmentResult::NewVertex& newEnd = newVertexTypes[end];
              newEnd.Type = (end == 4) ? SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin :
                                         SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;
              continue;
            }
          }
        }

        for (map<unsigned int, SplitPolygonWithSegmentResult::NewVertex>::const_iterator it = newVertexTypes.begin();
             it != newVertexTypes.end(); it++)
        {
          result.NewVertices.push_back(SplitPolygonWithSegmentResult::NewVertex());
          SplitPolygonWithSegmentResult::NewVertex& newVertex = result.NewVertices.back();
          newVertex = it->second;
        }
        return result;
      }
    }

    // segment is on not contigous edges and generates two new polygons
    list<SplitPolygonWithSegmentResult::NewPolygon> polygons;

    map<unsigned int, SplitPolygonWithSegmentResult::NewVertex> newVertices;
    map<unsigned int, SplitPolygonWithSegmentResult::NewEdge> newEdges;
    vector<bool> visited(input.NumberPolygonVertices, false);
    bool allVerticesVisited = false;
    unsigned int v = 0;
    do
    {
      // starting new polygon creation
      unsigned int startingVertex = v;
      polygons.push_back(SplitPolygonWithSegmentResult::NewPolygon());
      SplitPolygonWithSegmentResult::NewPolygon& newPolygon = polygons.back();

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

            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdge = newEdges[input.NumberPolygonVertices];
            newEdge.OriginId = input.Segment.Origin.Index;
            newEdge.EndId = input.Segment.End.Index;
            newEdge.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdge.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
          }
          else
          {
            // end segment is on an edge, add one new vertex and three edges in polygon
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 1);

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices];
            newEnd.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;


            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.Segment.Origin.Index;
            newEdgeOne.EndId = input.NumberPolygonVertices;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.NumberPolygonVertices;
            newEdgeTwo.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.End.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.Segment.End.Index;
            newEdgeThree.EndId = input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdge = newEdges[input.NumberPolygonVertices];
            newEdge.OriginId = input.Segment.Origin.Index;
            newEdge.EndId = input.Segment.End.Index;
            newEdge.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdge.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
          }
          else
          {
            // origin segment is on edge, add one new vertex and three edges in polygon
            newPolygon.Vertices.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices);
            newPolygon.Edges.push_back(input.NumberPolygonVertices + 2);

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;

            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.Segment.End.Index;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;
            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices + 1,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices + 1];
            newEnd.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;


            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 3,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 4,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.NumberPolygonVertices + 1;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeFour = newEdges[input.NumberPolygonVertices + 3];
            newEdgeFour.OriginId = input.Segment.End.Index;
            newEdgeFour.EndId = input.NumberPolygonVertices + 1;
            newEdgeFour.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeFour.OldEdgeId = input.Segment.End.Index;
            newEdgeFour.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeFive = newEdges[input.NumberPolygonVertices + 4];
            newEdgeFive.OriginId = input.NumberPolygonVertices + 1;
            newEdgeFive.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeFive.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;


            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.Segment.End.Index;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1, 0 } : vector<unsigned int>{ 0, 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newOrigin = newVertices[input.NumberPolygonVertices];
            newOrigin.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin;
            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices + 1,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices + 1];
            newEnd.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;

            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 3,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 4,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.NumberPolygonVertices;
            newEdgeOne.EndId = input.NumberPolygonVertices + 1;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.Segment.Origin.Index;
            newEdgeTwo.EndId = input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.Origin.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.NumberPolygonVertices;
            newEdgeThree.EndId = (input.Segment.Origin.Index + 1) % input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeThree.OldEdgeId = input.Segment.Origin.Index;
            newEdgeThree.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeFour = newEdges[input.NumberPolygonVertices + 3];
            newEdgeFour.OriginId = input.Segment.End.Index;
            newEdgeFour.EndId = input.NumberPolygonVertices + 1;
            newEdgeFour.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeFour.OldEdgeId = input.Segment.End.Index;
            newEdgeFour.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0 } : vector<unsigned int>{ 1 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeFive = newEdges[input.NumberPolygonVertices + 4];
            newEdgeFive.OriginId = input.NumberPolygonVertices + 1;
            newEdgeFive.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeFive.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

            newVertices.insert(pair<unsigned int,SplitPolygonWithSegmentResult::NewVertex>(input.NumberPolygonVertices,
                                                                                           SplitPolygonWithSegmentResult::NewVertex()));

            SplitPolygonWithSegmentResult::NewVertex& newEnd = newVertices[input.NumberPolygonVertices];
            newEnd.Type = SplitPolygonWithSegmentResult::NewVertex::Types::SegmentEnd;

            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 1,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            newEdges.insert(pair<unsigned int, SplitPolygonWithSegmentResult::NewEdge>(input.NumberPolygonVertices + 2,
                                                                                       SplitPolygonWithSegmentResult::NewEdge()));
            SplitPolygonWithSegmentResult::NewEdge& newEdgeOne = newEdges[input.NumberPolygonVertices];
            newEdgeOne.OriginId = input.Segment.Origin.Index;
            newEdgeOne.EndId = input.NumberPolygonVertices;
            newEdgeOne.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew;
            newEdgeOne.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 0, 1 } : vector<unsigned int>{ 1, 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeTwo = newEdges[input.NumberPolygonVertices + 1];
            newEdgeTwo.OriginId = input.NumberPolygonVertices;
            newEdgeTwo.EndId = (input.Segment.End.Index + 1) % input.NumberPolygonVertices;
            newEdgeTwo.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
            newEdgeTwo.OldEdgeId = input.Segment.End.Index;
            newEdgeTwo.Cell2DNeighbours = polygons.size() == 0 ? vector<unsigned int>{ 1 } : vector<unsigned int>{ 0 };
            SplitPolygonWithSegmentResult::NewEdge& newEdgeThree = newEdges[input.NumberPolygonVertices + 2];
            newEdgeThree.OriginId = input.Segment.End.Index;
            newEdgeThree.EndId = input.NumberPolygonVertices;
            newEdgeThree.Type = SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate;
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

    result.Type = SplitPolygonWithSegmentResult::Types::PolygonCreation;
    for (map<unsigned int, SplitPolygonWithSegmentResult::NewVertex>::const_iterator it = newVertices.begin();
         it != newVertices.end(); it++)
    {
      result.NewVertices.push_back(SplitPolygonWithSegmentResult::NewVertex());
      SplitPolygonWithSegmentResult::NewVertex& newVertex = result.NewVertices.back();
      newVertex = it->second;
    }
    for (map<unsigned int, SplitPolygonWithSegmentResult::NewEdge>::const_iterator it = newEdges.begin();
         it != newEdges.end(); it++)
    {
      result.NewEdges.push_back(SplitPolygonWithSegmentResult::NewEdge());
      SplitPolygonWithSegmentResult::NewEdge& newEdge = result.NewEdges.back();
      newEdge = it->second;
    }

    result.NewPolygons.resize(polygons.size());
    copy(polygons.begin(), polygons.end(), result.NewPolygons.begin());

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::SplitPolygonWithCircleResult GeometryUtilities::SplitPolygonWithCircle(const Eigen::MatrixXd& polygonVertices,
                                                                                            const Eigen::Vector3d& circleCenter,
                                                                                            const double& circleRadius,
                                                                                            const vector<PointCirclePositionResult>& vertexPositions,
                                                                                            const IntersectionPolygonCircleResult& polygonCircleIntersections,
                                                                                            const PolygonCirclePositionTypes& polygonCirclePosition) const
  {
    SplitPolygonWithCircleResult result;

    switch (polygonCirclePosition)
    {
      case PolygonCirclePositionTypes::PolygonOutsideCircleNoIntersection:
      case PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionOnVertex:
      case PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection:
      case PolygonCirclePositionTypes::PolygonInsideCircleOneVertexIntersection:
      case PolygonCirclePositionTypes::PolygonInsideCircleIntersectionOnlyOnVertices:
      {
        result.Type = SplitPolygonWithCircleResult::Types::NoAction;
        return result;
      }
      case PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionTangentOnEdge:
      {
        result.Type = SplitPolygonWithCircleResult::Types::PolygonUpdate;
        return result;
      }
      break;
      case PolygonCirclePositionTypes::CircleInsidePolygonNoIntersection:
      {
        result.Type = SplitPolygonWithCircleResult::Types::PolygonCreation;
        return result;
      }
      break;
      case PolygonCirclePositionTypes::CircleInsidePolygonOneIntersectionTangentOnEdge:
      {
        result.Type = SplitPolygonWithCircleResult::Types::PolygonCreation;
        return result;
      }
      break;
      case PolygonCirclePositionTypes::CirclePolygonMultipleIntersections:
      {
        const unsigned int numVertices = polygonVertices.cols();
        const unsigned int numCircleIntersections = polygonCircleIntersections.Intersections.size();
        Output::Assert(numVertices > 0 && numCircleIntersections > 0);

        list<SplitPolygonWithCircleResult::NewVertex> newVertices;
        result.PolygonVerticesNewVerticesPosition.resize(numVertices);
        result.CircleIntersectionsNewVerticesPosition.resize(numCircleIntersections);

        // order vertices with circle intersections and create new vertices
        unsigned int checkIntersection = 0;
        for (unsigned int v = 0; v < numVertices; v++)
        {
          result.PolygonVerticesNewVerticesPosition[v] = newVertices.size();
          newVertices.push_back(SplitPolygonWithCircleResult::NewVertex());
          SplitPolygonWithCircleResult::NewVertex& vertex = newVertices.back();
          vertex.Type = SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex;
          vertex.Index = v;

          // add circle intersection if in the same edge
          while (checkIntersection < numCircleIntersections &&
                 polygonCircleIntersections.Intersections[checkIntersection].Index == v)
          {
            const IntersectionPolygonCircleResult::Intersection::IndexTypes& intersectionIndexType =
                polygonCircleIntersections.Intersections[checkIntersection].IndexType;

            Output::Assert(intersectionIndexType !=
                IntersectionPolygonCircleResult::Intersection::IndexTypes::Unknown);

            if (intersectionIndexType ==
                IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex)
            {
              // vertex intersection, nothing to do
              result.CircleIntersectionsNewVerticesPosition[checkIntersection] =
                  result.PolygonVerticesNewVerticesPosition[v];
              checkIntersection++;
            }
            else if (intersectionIndexType ==
                     IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge)
            {
              // edge intersection, add it in the list
              result.CircleIntersectionsNewVerticesPosition[checkIntersection] = newVertices.size();
              newVertices.push_back(SplitPolygonWithCircleResult::NewVertex());
              SplitPolygonWithCircleResult::NewVertex& vertex = newVertices.back();
              vertex.Type = SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection;
              vertex.Index = checkIntersection;
              checkIntersection++;
            }
          }
        }

        // compute new polygons number vertices
        unsigned int numNewPolygons = numCircleIntersections + 1;
        result.NewPolygons.resize(numNewPolygons);
        vector<unsigned int> numNewPolygonVertices(numNewPolygons, 0);
        for (unsigned int c = 0; c < numCircleIntersections - 1; c++)
        {
          numNewPolygonVertices[c] = result.CircleIntersectionsNewVerticesPosition[c + 1] -
                                     result.CircleIntersectionsNewVerticesPosition[c] + 1;
          result.NewPolygons[c].Vertices.resize(numNewPolygonVertices[c]);
          result.NewPolygons[c].Edges.resize(numNewPolygonVertices[c]);
        }

        numNewPolygonVertices[numCircleIntersections - 1] =
            (newVertices.size() - result.CircleIntersectionsNewVerticesPosition[numCircleIntersections - 1]) +
            result.CircleIntersectionsNewVerticesPosition[0] + 1;
        result.NewPolygons[numCircleIntersections - 1].Vertices.resize(numNewPolygonVertices[numCircleIntersections - 1]);
        result.NewPolygons[numCircleIntersections - 1].Edges.resize(numNewPolygonVertices[numCircleIntersections - 1]);

        numNewPolygonVertices[numCircleIntersections] = numCircleIntersections;
        result.NewPolygons[numCircleIntersections].Vertices.resize(numNewPolygonVertices[numCircleIntersections]);
        result.NewPolygons[numCircleIntersections].Edges.resize(numNewPolygonVertices[numCircleIntersections]);
        result.NewPolygons[numCircleIntersections].Type = SplitPolygonWithCircleResult::NewPolygon::Types::InsideCircleAndPolygon;

        // create new polygons outside the circle polygon intersection

        // then create new polygon inside the circle polygon intersection
        vector<bool> verticesVisited(numVertices, false);
        list<unsigned int> vertexInsideCircle;

        for (unsigned int v = 0; v < numVertices; v++)
        {
          if (vertexPositions[v] == PointCirclePositionResult::Outside)
            continue;

          verticesVisited[v] = true;
          vertexInsideCircle.push_back(v);
        }

        // building the result
        result.Type = SplitPolygonWithCircleResult::Types::PolygonCreation;
        result.NewVertices = vector<SplitPolygonWithCircleResult::NewVertex>(newVertices.begin(),
                                                                             newVertices.end());

        cerr<< "RESULTS"<< endl;
        cerr<< "PolygonVerticesNewVerticesPosition: "<< result.PolygonVerticesNewVerticesPosition<< endl;
        cerr<< "CircleIntersectionsNewVerticesPosition: "<< result.CircleIntersectionsNewVerticesPosition<< endl;
        cerr<< "NewVertices: "<< endl;
        for (const auto& newVertex : result.NewVertices)
        {
          switch (newVertex.Type) {
            case Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex:
              cerr<< "\t"<< "Type: "<< "PolygonVertex"<< endl;
            break;
            case Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection:
              cerr<< "\t"<< "Type: "<< "CircleIntersection"<< endl;
            break;
            default:
              cerr<< "\t"<< "Type: "<< "Unknown"<< endl;
            break;
          }

          cerr<< "\t"<< "Index: "<< newVertex.Index<< endl;
        }

        cerr<< "NewPolygons: "<< endl;
        for (const auto& newPolygon : result.NewPolygons)
        {
          switch (newPolygon.Type)
          {
            case Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle:
              cerr<< "\t"<< "Type: "<< "InsideOnlyCircle"<< endl;
            break;
            case Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon:
              cerr<< "\t"<< "Type: "<< "InsideOnlyPolygon"<< endl;
            break;
            case Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideCircleAndPolygon:
              cerr<< "\t"<< "Type: "<< "InsideCircleAndPolygon"<< endl;
            break;
            default:
               cerr<< "\t"<< "Type: "<< "Unknown"<< endl;
            break;

          }
          cerr<< "\t"<< "Vertices: "<< newPolygon.Vertices<< endl;
          cerr<< "\t"<< "Edges: "<< newPolygon.Edges<< endl;
        }

        return result;
      }
      break;
      default:
      break;
    }

    throw runtime_error("SplitPolygonWithCircle failed");
  }
  // ***************************************************************************
}
