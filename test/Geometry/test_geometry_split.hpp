#ifndef __TEST_GEOMETRY_SPLIT_H
#define __TEST_GEOMETRY_SPLIT_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestSplitPolygonWithSegment)
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::NoAction);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 0);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 1);
        ASSERT_EQ(newEdgeIter->OriginId, 1);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 2);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 0);
        ASSERT_EQ(newEdgeIter->OriginId, 0);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 1);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 2);
        ASSERT_EQ(newEdgeIter->OriginId, 2);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::NewEdge::Types::EdgeUpdate);
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

        Gedim::GeometryUtilities::SplitPolygonWithSegmentResult result = geometryUtility.SplitPolygonWithSegment(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithSegmentResult::Types::PolygonCreation);
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

  TEST(TestGeometryUtilities, TestSplitPolygonWithCircle)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // circle outside and no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // Polygon Inside Circle with center outside no intersections
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 10.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // Polygon Inside Circle with center outside one intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(sqrt(2.0), sqrt(2.0), 0.0);
        double circleRadius = 2.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // circle inside polygon no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
        double circleRadius = 0.125;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);

        ASSERT_ANY_THROW(geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                circleCenter,
                                                                circleRadius,
                                                                vertexPositions,
                                                                polygonCircleIntersections,
                                                                polygonPosition));
      }

      // Polygon Inside Circle with center inside no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
        double circleRadius = 10.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // Polygon inside Circle intersects only vertices
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
        double circleRadius = sqrt(2) / 2;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // Circle Outside Polygon one intersection with vertex
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 2.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::NoAction);
      }

      // Circle outside Polygon tangent to edge
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, -1.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        ASSERT_ANY_THROW(geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                circleCenter,
                                                                circleRadius,
                                                                vertexPositions,
                                                                polygonCircleIntersections,
                                                                polygonPosition));
      }

      // Circle inside Polygon tangent to edge
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.125, 0.0);
        double circleRadius = 0.125;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        ASSERT_ANY_THROW(geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                circleCenter,
                                                                circleRadius,
                                                                vertexPositions,
                                                                polygonCircleIntersections,
                                                                polygonPosition));
      }

      // Circle Outside Polygon Intersects With Multiple SubPolygons
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.55, 0.0);
        double circleRadius = 0.5;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::PolygonCreation);
        ASSERT_EQ(result.PolygonVerticesNewVerticesPosition, vector<unsigned int>({ 0, 1, 4 }));
        ASSERT_EQ(result.CircleIntersectionsNewVerticesPosition, vector<unsigned int>({ 2, 3, 5 }));
        ASSERT_EQ(result.NewVertices.size(), 6);
        ASSERT_EQ(result.NewVertices[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[0].PolygonIndex, 0);
        ASSERT_EQ(result.NewVertices[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[1].PolygonIndex, 1);
        ASSERT_EQ(result.NewVertices[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[2].IntersectionIndex, 0);
        ASSERT_EQ(result.NewVertices[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[3].IntersectionIndex, 1);
        ASSERT_EQ(result.NewVertices[4].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[4].PolygonIndex, 2);
        ASSERT_EQ(result.NewVertices[5].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[5].IntersectionIndex, 2);
        ASSERT_EQ(result.NewPolygons.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle);
        ASSERT_EQ(result.NewPolygons[0].Vertices, vector<unsigned int>({ 2, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges, vector<unsigned int>({ 0, 0 }));
        ASSERT_EQ(result.NewPolygons[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[1].Vertices, vector<unsigned int>({ 3, 4, 5 }));
        ASSERT_EQ(result.NewPolygons[1].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[2].Vertices, vector<unsigned int>({ 5, 0, 1, 2 }));
        ASSERT_EQ(result.NewPolygons[2].Edges, vector<unsigned int>({ 0, 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideCircleAndPolygon);
        ASSERT_EQ(result.NewPolygons[3].Vertices, vector<unsigned int>({ 2, 3, 5 }));
        ASSERT_EQ(result.NewPolygons[3].Edges, vector<unsigned int>({ 0, 0, 0 }));
      }

      // Circle Inside Polygon Intersects With Multiple SubPolygons
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
        double circleRadius = 0.5;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::PolygonCreation);
        ASSERT_EQ(result.PolygonVerticesNewVerticesPosition, vector<unsigned int>({ 0, 2, 5 }));
        ASSERT_EQ(result.CircleIntersectionsNewVerticesPosition, vector<unsigned int>({ 1, 3, 4, 6 }));
        ASSERT_EQ(result.NewVertices.size(), 7);
        ASSERT_EQ(result.NewVertices[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[0].PolygonIndex, 0);
        ASSERT_EQ(result.NewVertices[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[1].IntersectionIndex, 0);
        ASSERT_EQ(result.NewVertices[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[2].PolygonIndex, 1);
        ASSERT_EQ(result.NewVertices[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[3].IntersectionIndex, 1);
        ASSERT_EQ(result.NewVertices[4].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[4].IntersectionIndex, 2);
        ASSERT_EQ(result.NewVertices[5].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[5].PolygonIndex, 2);
        ASSERT_EQ(result.NewVertices[6].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[6].IntersectionIndex, 3);
        ASSERT_EQ(result.NewPolygons.size(), 5);
        ASSERT_EQ(result.NewPolygons[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[0].Vertices, vector<unsigned int>({ 1, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle);
        ASSERT_EQ(result.NewPolygons[1].Vertices, vector<unsigned int>({ 3, 4 }));
        ASSERT_EQ(result.NewPolygons[1].Edges, vector<unsigned int>({ 0, 0 }));
        ASSERT_EQ(result.NewPolygons[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[2].Vertices, vector<unsigned int>({ 4, 5, 6 }));
        ASSERT_EQ(result.NewPolygons[2].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[3].Vertices, vector<unsigned int>({ 6, 0, 1 }));
        ASSERT_EQ(result.NewPolygons[3].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[4].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideCircleAndPolygon);
        ASSERT_EQ(result.NewPolygons[4].Vertices, vector<unsigned int>({ 1, 3, 4, 6 }));
        ASSERT_EQ(result.NewPolygons[4].Edges, vector<unsigned int>({ 0, 0, 0, 0 }));
      }

      // Generic intersections with quadrilateral
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 1.0, 3.0, 0.0;
        polygonVertices.col(1)<< 3.0, 3.0 - 2.0 / sqrt(3.0), 0.0;
        polygonVertices.col(2)<< 4.0 + 1.0 / 10.0, 3.0, 0.0;
        polygonVertices.col(3)<< 3.0, 4.0, 0.0;
        Eigen::Vector3d circleCenter(3.0, 3.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes polygonPosition = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                                     circleCenter,
                                                                                                                     circleRadius,
                                                                                                                     vertexPositions,
                                                                                                                     polygonCircleIntersections);
        Gedim::GeometryUtilities::SplitPolygonWithCircleResult result = geometryUtility.SplitPolygonWithCircle(polygonVertices,
                                                                                                               circleCenter,
                                                                                                               circleRadius,
                                                                                                               vertexPositions,
                                                                                                               polygonCircleIntersections,
                                                                                                               polygonPosition);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::Types::PolygonCreation);
        ASSERT_EQ(result.PolygonVerticesNewVerticesPosition, vector<unsigned int>({ 0, 2, 5, 7 }));
        ASSERT_EQ(result.CircleIntersectionsNewVerticesPosition, vector<unsigned int>({ 1, 3, 4, 6, 7, 8 }));
        ASSERT_EQ(result.NewVertices.size(), 9);
        ASSERT_EQ(result.NewVertices[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[0].PolygonIndex, 0);
        ASSERT_EQ(result.NewVertices[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[1].IntersectionIndex, 0);
        ASSERT_EQ(result.NewVertices[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[2].PolygonIndex, 1);
        ASSERT_EQ(result.NewVertices[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[3].IntersectionIndex, 1);
        ASSERT_EQ(result.NewVertices[4].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[4].IntersectionIndex, 2);
        ASSERT_EQ(result.NewVertices[5].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::PolygonVertex);
        ASSERT_EQ(result.NewVertices[5].PolygonIndex, 2);
        ASSERT_EQ(result.NewVertices[6].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[6].IntersectionIndex, 3);
        ASSERT_EQ(result.NewVertices[7].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::Both);
        ASSERT_EQ(result.NewVertices[7].PolygonIndex, 3);
        ASSERT_EQ(result.NewVertices[7].IntersectionIndex, 4);
        ASSERT_EQ(result.NewVertices[8].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewVertex::Types::CircleIntersection);
        ASSERT_EQ(result.NewVertices[8].IntersectionIndex, 5);
        ASSERT_EQ(result.NewPolygons.size(), 7);
        ASSERT_EQ(result.NewPolygons[0].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[0].Vertices, vector<unsigned int>({ 1, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[1].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle);
        ASSERT_EQ(result.NewPolygons[1].Vertices, vector<unsigned int>({ 3, 4 }));
        ASSERT_EQ(result.NewPolygons[1].Edges, vector<unsigned int>({ 0, 0 }));
        ASSERT_EQ(result.NewPolygons[2].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[2].Vertices, vector<unsigned int>({ 4, 5, 6 }));
        ASSERT_EQ(result.NewPolygons[2].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[3].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle);
        ASSERT_EQ(result.NewPolygons[3].Vertices, vector<unsigned int>({ 6, 7 }));
        ASSERT_EQ(result.NewPolygons[3].Edges, vector<unsigned int>({ 0, 0 }));
        ASSERT_EQ(result.NewPolygons[4].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyCircle);
        ASSERT_EQ(result.NewPolygons[4].Vertices, vector<unsigned int>({ 7, 8 }));
        ASSERT_EQ(result.NewPolygons[4].Edges, vector<unsigned int>({ 0, 0 }));
        ASSERT_EQ(result.NewPolygons[5].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideOnlyPolygon);
        ASSERT_EQ(result.NewPolygons[5].Vertices, vector<unsigned int>({ 8, 0, 1 }));
        ASSERT_EQ(result.NewPolygons[5].Edges, vector<unsigned int>({ 0, 0, 0 }));
        ASSERT_EQ(result.NewPolygons[6].Type, Gedim::GeometryUtilities::SplitPolygonWithCircleResult::NewPolygon::Types::InsideCircleAndPolygon);
        ASSERT_EQ(result.NewPolygons[6].Vertices, vector<unsigned int>({ 1, 3, 4, 6, 7, 8 }));
        ASSERT_EQ(result.NewPolygons[6].Edges, vector<unsigned int>({ 0, 0, 0, 0, 0, 0 }));
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
