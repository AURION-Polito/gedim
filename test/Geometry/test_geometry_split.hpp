#ifndef __TEST_GEOMETRY_SPLIT_H
#define __TEST_GEOMETRY_SPLIT_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, SplitPolygon)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // no action
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::NoAction);
        ASSERT_EQ(result.NewVertices.size(), 0);
        ASSERT_EQ(result.NewEdges.size(), 0);
        ASSERT_EQ(result.NewPolygons.size(), 0);
      }

      // update square with one point
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 2);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 5);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 5);
      }

      // update square with two points
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 6);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 6);
      }

      // update square with two points with aligned edges
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 5;
        input.AlignedEdges.resize(1);
        input.AlignedEdges[0].OriginVertexIndex = 1;
        input.AlignedEdges[0].EndVertexIndex = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 2;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 4);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 7);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 7);
      }

      // split square in two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 0);
        ASSERT_EQ(result.NewEdges.size(), 1);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
      }

      // split square in a triangle and a trapezioid
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 4);
      }

      // split square in a two trapezioids
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 5);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 4);
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 1;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 0);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 1);
        ASSERT_EQ(newEdgeIter->OriginId, 1);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 1);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 2);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 3, 5, 2 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 3, 0 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 4, 3, 0 }));
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 2);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 0);
        ASSERT_EQ(newEdgeIter->OriginId, 0);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 0);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 1);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 4, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 1, 3, 5 }));
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 2;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 1);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 2);
        ASSERT_EQ(newEdgeIter->OriginId, 2);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 2);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 0);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 1, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 0, 3, 5 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 2, 3, 1 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 4, 3, 1 }));
      }

      // split square in a a triangle and a trapezioid
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 4, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 5, 4, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 2, 4 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 1, 4, 6 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_SPLIT_H
