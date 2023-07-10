#ifndef __GEOMETRYUTILITIES_H
#define __GEOMETRYUTILITIES_H

#include "Eigen/Eigen"
#include <iostream>
#include "IOUtilities.hpp"

namespace Gedim
{
  struct GeometryUtilitiesConfig final
  {
      double Tolerance = std::numeric_limits<double>::epsilon();
  };

  /// \brief The GeometryUtilities class intersects 3D segments
  class GeometryUtilities final
  {
    private:
      const GeometryUtilitiesConfig& _configuration;

    public:
      enum struct CompareTypes
      {
        Unknown = 0,
        FirstBeforeSecond = 1,
        Coincident = 2,
        SecondBeforeFirst = 3
      };

      enum struct PolygonTypes
      {
        Unknown = 0,
        Triangle = 1,
        Quadrilateral = 2,
        Generic = 3
      };

      enum struct PointSegmentPositionTypes
      {
        Unknown = 0,
        OnSegmentLineBeforeOrigin = 1,
        OnSegmentOrigin = 2,
        InsideSegment = 3,
        OnSegmentEnd = 4,
        OnSegmentLineAfterEnd = 5,
        LeftTheSegment = 6,
        RightTheSegment = 7
      };

      enum struct PointPlanePositionTypes
      {
        Unknown = 0,
        Negative = 1,
        OnPlane = 2,
        Positive = 3
      };

      enum struct PolygonCirclePositionTypes
      {
        Unknown = 0,
        PolygonOutsideCircleNoIntersection = 1,
        PolygonOutsideCircleOneIntersectionOnVertex = 2,
        PolygonOutsideCircleOneIntersectionTangentOnEdge = 3,
        CircleInsidePolygonNoIntersection = 4,
        CircleInsidePolygonOneIntersectionTangentOnEdge = 5,
        PolygonInsideCircleNoIntersection = 6,
        PolygonInsideCircleOneVertexIntersection = 7,
        PolygonInsideCircleIntersectionOnlyOnVertices = 8,
        CirclePolygonMultipleIntersections = 9
      };

      struct IntersectionPolygonCircleResult final
      {
          struct Intersection final
          {
              enum struct Types {
                Unknown = 0,
                Secant = 1,
                Tangent = 2
              };

              enum struct IndexTypes {
                Unknown = 0,
                Vertex = 1,
                Edge = 2
              };

              Types Type = Types::Unknown;
              IndexTypes IndexType = IndexTypes::Unknown;
              unsigned int Index;
              double CurvilinearCoordinate; ///< Valid only in IndexType Edge
          };

          std::vector<Intersection> Intersections = {}; ///< ordered by edge order
      };

      struct PolygonDivisionByAngleQuadrantResult final
      {
          enum struct Types
          {
            Unknown = 0,
            ExternalOrigin = 1,
            Internal = 2,
            ExternalEnd = 3
          };

          Eigen::MatrixXd Points; ///< Coordinates of generated points
          std::vector<std::vector<unsigned int>> SubPolygons; ///< Subpolygon formed
          std::vector<Types> SubPolygonTypes; ///< SubPolygon types
      };

      struct PolygonDivisionByCircleResult final
      {
          Eigen::MatrixXd Points; ///< Coordinates of generated points
          std::vector<std::vector<unsigned int>> SubTriangles; ///< Triangle formed with sub-polygons and circle Center
          std::vector<std::vector<unsigned int>> InternalTriangles; ///< Triangle formed with circle Center and new points
          std::vector<std::vector<unsigned int>> SubPolygons; ///< Subpolygon formed
      };

      struct CircleDivisionByPolygonResult final
      {
          Eigen::MatrixXd Points; ///< Coordinates of generated points
          std::vector<std::vector<unsigned int>> SubTriangles; ///< Triangle formed with sub-polygons and circle Center
          std::vector<std::vector<unsigned int>> InternalTriangles; ///< Triangle formed with circle Center and new points
          std::vector<std::vector<unsigned int>> SubPolygons; ///< Subpolygon formed
      };

      struct SplitPolygonInput final
      {
          struct AlignedEdge
          {
              unsigned int OriginVertexIndex;
              unsigned int EndVertexIndex;
          };

          struct SplitSegment
          {
              struct Vertex final
              {
                  enum struct Types {
                    Unknown = 0,
                    Vertex = 1,
                    Edge = 2
                  };

                  Types Type;
                  unsigned int Index;
              };

              Vertex Origin;
              Vertex End;
          };

          unsigned int NumberPolygonVertices;
          std::vector<AlignedEdge> AlignedEdges;
          SplitSegment Segment;
      };

      struct SplitPolygonWithSegmentResult final
      {
          enum struct Types
          {
            Unknown = 0,
            NoAction = 1,
            PolygonUpdate = 2,
            PolygonCreation = 3
          };

          struct NewVertex final
          {
              enum struct Types
              {
                Unknown = 0,
                SegmentOrigin = 1,
                SegmentEnd = 2
              };

              Types Type = Types::Unknown;
          };

          struct NewEdge final
          {
              enum struct Types
              {
                Unknown = 0,
                EdgeNew = 1,
                EdgeUpdate = 2
              };

              Types Type = Types::Unknown;
              unsigned int OldEdgeId = 0;
              unsigned int OriginId = 0;
              unsigned int EndId = 0;
              std::vector<unsigned int> Cell2DNeighbours = {};
          };

          struct NewPolygon final
          {
              std::list<unsigned int> Vertices = {};
              std::list<unsigned int> Edges = {};
          };

          Types Type = Types::Unknown;
          std::list<NewVertex> NewVertices = {};
          std::list<NewEdge> NewEdges = {};
          std::vector<NewPolygon> NewPolygons = {};
      };

      struct SplitPolygonWithCircleResult final
      {
          enum struct Types
          {
            Unknown = 0,
            NoAction = 1,
            PolygonUpdate = 2,
            PolygonCreation = 3
          };

          struct NewVertex final
          {
              enum struct Types
              {
                Unknown = 0,
                PolygonVertex = 1,
                CircleIntersection = 2,
                Both = 3
              };

              Types Type = Types::Unknown;
              unsigned int PolygonIndex; ///< Index in polygon vertices
              unsigned int IntersectionIndex; ///< Index in circle intersections
          };

          struct NewEdge final
          {
              enum struct Types
              {
                Unknown = 0,
                Segment = 1,
                Arc = 2
              };

              enum struct ArcTypes
              {
                Unknown = 0,
                InsidePolygon = 1,
                OutsidePolygon = 2
              };

              Types Type = Types::Unknown;
              ArcTypes ArcType = ArcTypes::Unknown; ///< Valid only if Type is Arc
              std::vector<unsigned int> VertexIndices; ///< Index of vertices in NewVertices
              unsigned int PolygonIndex; ///< Index of Edge in polygon intersections
          };

          struct NewPolygon final
          {
              enum struct Types
              {
                Unknown = 0,
                InsideOnlyCircle = 1,
                InsideOnlyPolygon = 2,
                InsideCircleAndPolygon = 3
              };

              Types Type = Types::Unknown;
              std::vector<unsigned int> Vertices = {};
              std::vector<unsigned int> Edges = {};
          };

          Types Type = Types::Unknown;
          std::vector<unsigned int> PolygonVerticesNewVerticesPosition = {};
          std::vector<unsigned int> CircleIntersectionsNewVerticesPosition = {};
          std::vector<NewVertex> NewVertices = {};
          std::vector<NewEdge> NewEdges = {};
          std::vector<NewPolygon> NewPolygons = {};
      };

      struct IntersectionSegmentSegmentResult final
      {
          enum struct IntersectionLineTypes
          {
            Unknown = 0,
            OnDifferentPlanes = 1,
            CoPlanarParallel = 2,
            CoPlanarIntersecting = 3
          };

          enum struct IntersectionSegmentTypes
          {
            Unknown = 0,
            NoIntersection = 1,
            SingleIntersection = 2,
            MultipleIntersections = 3
          };

          struct IntersectionPosition
          {
              PointSegmentPositionTypes Type = PointSegmentPositionTypes::Unknown;
              double CurvilinearCoordinate = 0.0;
          };

          IntersectionLineTypes IntersectionLinesType = IntersectionLineTypes::Unknown;
          IntersectionSegmentTypes IntersectionSegmentsType = IntersectionSegmentTypes::Unknown;
          /// \brief relation between first and second intersection.
          /// Values are the indeces of the SecondSegmentIntersections vector respect the FirstSegmentIntersections vector
          /// \example in MultipleIntersection case, if SecondIntersectionRelation[0] = 1,
          /// then the second intersection point SecondSegmentIntersections[1] is equal to FirstSegmentIntersections[0] point
          std::vector<unsigned int> SecondIntersectionRelation;
          /// \brief intersections of the first segment,
          /// \note if multiple intersections are found, than the origin and the end coordinate are stored
          std::vector<IntersectionPosition> FirstSegmentIntersections;
          /// \brief intersections of the second segment,
          /// \note if multiple intersections are found, than the origin and the end coordinate are stored
          std::vector<IntersectionPosition> SecondSegmentIntersections; /// intersections of the second segment
      };

      struct IntersectionSegmentCircleResult final
      {
          enum struct Types
          {
            Unknown = 0,
            NoIntersection = 1,
            TangentIntersection = 2,
            TwoIntersections = 3
          };

          struct IntersectionPosition
          {
              PointSegmentPositionTypes Type = PointSegmentPositionTypes::Unknown;
              double CurvilinearCoordinate = 0.0;
          };

          Types Type = Types::Unknown;
          /// \brief intersections of the segment,
          std::vector<IntersectionPosition> SegmentIntersections = {};
      };

      struct IntersectionSegmentPlaneResult final
      {
          enum struct Types
          {
            Unknown = 0,
            SingleIntersection = 1,
            NoIntersection = 2,
            MultipleIntersections = 3
          };

          struct Intersection
          {
              PointSegmentPositionTypes Type = PointSegmentPositionTypes::Unknown;
              double CurvilinearCoordinate = 0.0;
          };

          Types Type = Types::Unknown; ///< The intersection type
          Intersection SingleIntersection; ///< The single intersection, available only is Type is SingleIntersection
      };

      struct IntersectionPolyhedronLineResult final
      {
          enum struct Types
          {
            Unknown = 0,
            None = 1, ///< No intersection found
            OneIntersection = 2, ///< One intersection found
            TwoIntersections = 3, ///< Two intersection found
            MultipleIntersections = 4 ///< Multiple intersection found
          };

          struct PolyhedronFaceIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Intersection = 1,
                NoIntersection = 2
              };

              Types Type = Types::Unknown;
              unsigned int LineIntersectionIndex = 0; ///< Index of line intersection collection
          };

          struct PolyhedronEdgeIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Intersection = 1,
                NoIntersection = 2
              };

              Types Type = Types::Unknown;
              unsigned int LineIntersectionIndex = 0; ///< Index of line intersection collection
          };

          struct PolyhedronVertexIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Intersection = 1,
                NoIntersection = 2
              };

              Types Type = Types::Unknown;
              unsigned int LineIntersectionIndex = 0; ///< Index of line intersection collection
          };

          struct LineIntersection
          {
              enum struct Types
              {
                Unknown = 0,
                OnVertex = 1, ///< On polyhedron vertex
                OnEdge = 2, ///< On polyhedron edge
                OnFace = 3, ///< On polyhedron face
                Inside = 4, ///< Inside polyhedron
              };

              Types PolyhedronType = Types::Unknown; ///< Type of intersection
              unsigned int PolyhedronIndex = 0; ///<  Index of the intersecting element of the Polyhedron (face, edge or vertex index)
              double CurvilinearCoordinate = 0.0; ///< Curvilinear coordinate in the line
          };

          Types Type = Types::Unknown; /// The
          std::vector<LineIntersection> LineIntersections; ///< The line intersections
          std::vector<PolyhedronVertexIntersection> PolyhedronVertexIntersections = {}; ///< Polyhedron Vertex intersections, size polyhedron num vertices
          std::vector<PolyhedronEdgeIntersection> PolyhedronEdgeIntersections = {}; ///< Polyhedron Edge intersections, size polyhedron num edges
          std::vector<PolyhedronFaceIntersection> PolyhedronFaceIntersections = {}; ///< Polyhedron Face intersections, size polyhedron num faces
      };

      struct IntersectionPolyhedronsSegmentResult final
      {
          struct IntersectionPoint final
          {
              std::vector<unsigned int> Cell3DIndices = {};
          };

          struct IntersectionSegment final
          {
              std::vector<double> Points = {};
              std::vector<unsigned int> Cell3DIndices = {};
          };

          std::map<double, IntersectionPoint> Points;
          std::vector<IntersectionSegment> Segments;
      };

      struct IntersectionPolyhedronPlaneResult final
      {
          enum struct Types
          {
            Unknown = 0,
            None = 1, ///< No intersection found
            OnVertex = 2, ///< On polyhedron vertex
            OnEdge = 3, ///< On polyhedron edge
            OnFace = 4, ///< On polyhedron face
            NewPolygon = 5 ///< New polygon intersection
          };

          struct FaceIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Intersection = 1,
                NoIntersection = 2
              };

              Types Type = Types::Unknown;
          };

          struct EdgeIntersection final
          {
              IntersectionSegmentPlaneResult Intersection; ///< Intersection between edge and plane
          };

          struct VertexIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Intersection = 1,
                NoIntersection = 2
              };

              Types Type = Types::Unknown;
          };

          struct Intersection final
          {
              enum struct Types
              {
                Unknown = 0,
                Vertex = 1,
                Edge = 2
              };

              Types Type = Types::Unknown;
              unsigned int EdgeId = 0; ///<  Edge index of the Polyhedron
              unsigned int VertexId = 0; ///<  Vertex index of the Polyhedron, available only if Type is Types::Vertex
          };

          Types Type = Types::Unknown; ///< The intersection type
          unsigned int IntersectionId = 0; ///< The geometry id of the intersection, available only with Types::OnVertex, Types::OnEdge and Types::OnFace
          std::vector<VertexIntersection> VertexIntersections = {}; ///< Vertex intersections
          std::vector<EdgeIntersection> EdgeIntersections = {}; ///< Edge intersections
          std::vector<FaceIntersection> FaceIntersections = {}; ///< Face intersections
          std::vector<Intersection> Intersections = {}; ///< The resulting intersections
          Eigen::MatrixXd IntersectionCoordinates; ///< The resulting intersection coordinates
      };

      enum struct PointCirclePositionResult
      {
        Unknown = 0,
        Outside = 1,
        OnBorder = 2,
        Inside = 3
      };


      struct PointPolygonPositionResult final
      {
          enum struct Types
          {
            Unknown = 0,
            Outside = 1,
            BorderEdge = 2,
            BorderVertex = 3,
            Inside = 4
          };

          unsigned int BorderIndex = 0; ///< index of vertex/edge of border
          Types Type = Types::Unknown;
      };

      struct LinePolygonPositionResult final
      {
          enum struct Types
          {
            Unknown = 0,
            Outside = 1,
            Intersecting = 2
          };

          struct EdgeIntersection final
          {
              enum struct Types
              {
                Unknown = 0,
                OnEdgeOrigin = 1,
                InsideEdge = 2,
                OnEdgeEnd = 3,
                Parallel = 4
              };

              Types Type = Types::Unknown;
              unsigned int Index = 0;
              double CurvilinearCoordinate = 0.0;
          };

          std::vector<EdgeIntersection> EdgeIntersections = { };
          Types Type = Types::Unknown;
      };

      struct PointPolyhedronPositionResult final
      {
          enum struct Types
          {
            Unknown = 0,
            Outside = 1,
            BorderFace = 2,
            BorderEdge = 3,
            BorderVertex = 4,
            Inside = 5
          };

          unsigned int BorderIndex = 0; ///< index of vertex/edge/face of border
          Types Type = Types::Unknown;
      };

      struct Polyhedron final
      {
          Eigen::MatrixXd Vertices; ///< vertices, size 3 x numVertices
          Eigen::MatrixXi Edges; ///< edges, size 2 x numEdges
          std::vector<Eigen::MatrixXi> Faces; ///< faces vertices and edges˝, size numFaces x 2 x numFaceVertices
      };

      struct SplitPolygonWithPlaneResult final
      {
          enum struct Types
          {
            Unknown = 0,
            Split = 1,
            Positive = 2,
            Negative = 3,
            OnPlane = 4
          };

          std::vector<unsigned int> PositiveVertices; /// vertices indices of the positive sub-polygon
          std::vector<unsigned int> NegativeVertices; /// vertices indices of the negative sub-polygon
          std::vector<unsigned int> PointsOnPlane; /// vertices indices of the points on plane
          std::vector<Eigen::Vector3d> NewVertices; /// new vertices coordinates
          std::vector<unsigned int> NewVerticesEdgeIndex; /// new vertices edge indices
          Types Type = Types::Unknown; /// type of split
      };

      struct SplitPolyhedronWithPlaneResult
      {
          enum struct Types
          {
            Unknown = 0,
            Split = 1,
            None = 2
          };

          struct NewVertices
          {
              Eigen::MatrixXd Vertices; ///< all vertices contained in the new polyhedra
              std::vector<unsigned int> NewVerticesOriginalEdge; ///< For each new vertex the index of the original edge to which is located
          };

          struct NewEdges
          {
              Eigen::MatrixXi Edges; ///< all edges contained in the new polyhedra
              std::vector<int> NewEdgesOriginalEdges; ///< indices of original edges for new edges, -1 means no original edge
          };

          struct NewFaces
          {
              std::vector<Eigen::MatrixXi> Faces;
              std::vector<int> NewFacesOriginalFaces; ///< indices of original faces for new faces, -1 means no original face
          };

          struct NewPolyhedron
          {
              std::vector<unsigned int> Vertices;
              std::vector<unsigned int> Edges;
              std::vector<unsigned int> Faces;
          };

          NewVertices Vertices;
          NewEdges Edges;
          NewFaces Faces;
          NewPolyhedron PositivePolyhedron;
          NewPolyhedron NegativePolyhedron;
          Types Type = Types::Unknown; /// type of split
      };

    private:
      /// \brief Compare two values according to tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      /// \note the interval [-tolerance, tolerance] is considered 0.0
      CompareTypes CompareValues(const double& first,
                                 const double& second,
                                 const double& tolerance) const;
    public:
      GeometryUtilities(const GeometryUtilitiesConfig& configuration);
      ~GeometryUtilities();

      /// \param first the first value
      /// \param second the second value
      /// \return the relative difference between the two values according the first
      inline double RelativeDifference(const double& first,
                                       const double& second) const
      {
        return abs(second - first) / ((first == 0.0) ? 1.0 : abs(first));
      }

      /// \brief Compare two 1D values according to tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      inline CompareTypes Compare1DValues(const double& first,
                                          const double& second) const
      {
        return CompareValues(first,
                             second,
                             _configuration.Tolerance);
      }

      /// \brief Check if two 1D values are equal according to tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      inline bool Are1DValuesEqual(const double& first,
                                   const double& second) const
      {
        return CompareValues(first,
                             second,
                             _configuration.Tolerance) == CompareTypes::Coincident;
      }

      /// \param first the first value
      /// \param second the second value
      /// \return true if first is greater than second
      inline bool IsValue1DGreater(const double& first,
                                   const double& second) const
      {
        return Compare1DValues(first,
                               second) == CompareTypes::SecondBeforeFirst;
      }

      /// \param first the first value
      /// \param second the second value
      /// \return true if first is greater or equal than second
      inline bool IsValue1DGreaterOrEqual(const double& first,
                                          const double& second) const
      {
        const CompareTypes result = Compare1DValues(first,
                                                    second);

        return result == CompareTypes::SecondBeforeFirst ||
            result == CompareTypes::Coincident;
      }

      /// \param value the value
      /// \return true if value is positive
      inline bool IsValue1DPositive(const double& value) const
      {
        return Compare1DValues(0.0,
                               value) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue1DNegative(const double& value) const
      {
        return Compare1DValues(0.0,
                               value) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue1DZero(const double& value) const
      {
        return Compare1DValues(0.0,
                               value) == CompareTypes::Coincident;
      }

      /// \brief Compare two 2D values according to squared tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      inline CompareTypes Compare2DValues(const double& first,
                                          const double& second) const
      {
        return CompareValues(first,
                             second,
                             _configuration.Tolerance *
                             _configuration.Tolerance);
      }

      /// \param value the value
      /// \return true if value is positive
      inline bool IsValue2DPositive(const double& value) const
      {
        return Compare2DValues(0.0,
                               value) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue2DNegative(const double& value) const
      {
        return Compare2DValues(0.0,
                               value) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue2DZero(const double& value) const
      {
        return Compare2DValues(0.0,
                               value) == CompareTypes::Coincident;
      }

      /// \brief Compare two 3D values according to cube tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      inline CompareTypes Compare3DValues(const double& first,
                                          const double& second) const
      {
        return CompareValues(first,
                             second,
                             _configuration.Tolerance *
                             _configuration.Tolerance *
                             _configuration.Tolerance);
      }

      /// \param value the value
      /// \return true if value is positive
      inline bool IsValue3DPositive(const double& value) const
      {
        return Compare3DValues(0.0,
                               value) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue3DNegative(const double& value) const
      {
        return Compare3DValues(0.0,
                               value) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue3DZero(const double& value) const
      {
        return Compare3DValues(0.0,
                               value) == CompareTypes::Coincident;
      }

      /// \param step the distance between each coordinate
      /// \param insertExtremes if true keeps the extremes
      /// \return the equispace coordinates between [0.0, 1.0], size 1 x numCoordinates
      std::vector<double> EquispaceCoordinates(const double& step,
                                               const bool& insertExtremes) const;

      /// \param size the number of resulting coordinates
      /// \param origin the starting curvilinear coordinate
      /// \param end the ending curvilinear coordinate
      /// \param insertExtremes if true keeps the extremes
      /// \return equispaced curvilinear coordinates in the interval [origin, end]
      /// \note if size < 2 then size will be considered as 2
      std::vector<double> EquispaceCoordinates(const unsigned int& size,
                                               const double& origin,
                                               const double& end,
                                               const bool& insertExtremes) const;

      /// \brief compute the Point distance
      /// \param firstPoint the first point
      /// \param secondPoint the second point
      /// \return the distance
      inline double PointDistance(const Eigen::Vector3d& firstPoint,
                                  const Eigen::Vector3d& secondPoint) const
      {
        return (secondPoint - firstPoint).norm();
      }

      /// \brief compute the distance between a point and a list of points
      /// \param points the point collection, size 3 x numPoints
      /// \param point the point
      /// \return the collection of distances, size 1 x numPoints
      Eigen::VectorXd PointDistances(const Eigen::MatrixXd& points,
                                     const Eigen::Vector3d& point) const;

      /// \param points the point collection, size 3 x numPoints
      /// \return the distances between the points collected in matrix, size numPoints x numPoints.
      Eigen::MatrixXd PointsDistance(const Eigen::MatrixXd& points) const;

      /// \param points the point collection, size 3 x numPoints
      /// \return the maximum distance between the points.
      double PointsMaxDistance(const Eigen::MatrixXd& points) const;

      /// \param firstPoint the first point
      /// \param secondPoint the second point
      /// \return true if the points are coincident
      inline bool PointsAreCoincident(const Eigen::Vector3d& firstPoint,
                                      const Eigen::Vector3d& secondPoint) const
      {
        return IsValue1DZero(PointDistance(firstPoint, secondPoint));
      }

      /// \brief Find a point in point list
      /// \param points the point list, size 3 x numPoints
      /// \param point the point to find
      /// \return the collection of point found
      std::vector<unsigned int> FindPointInPoints(const Eigen::MatrixXd& points,
                                                  const Eigen::Vector3d& point) const;


      /// \brief Compute the distance between a point and a line
      /// \param point a point P
      /// \param lineOrigin the line origin O
      /// \param normalToLine a normal vector n to the line, in the same plane of P and the line
      /// \return the distance d
      /// \note The distance is computed as d = n^T * (P - O) / ||n||
      inline double PointLineDistance(const Eigen::Vector3d& point,
                                      const Eigen::Vector3d& lineOrigin,
                                      const Eigen::Vector3d& normalToLine) const
      {
        return abs(normalToLine.dot(point - lineOrigin)) / normalToLine.norm();
      }

      /// \param point the point
      /// \return true if the point is 2D (z == 0)
      inline bool PointIs2D(const Eigen::Vector3d& point) const
      {
        return PointsAre2D(point);
      }

      /// \param points the points to test, size 3 x numPoints
      /// \return true if the points are 2D (z == 0)
      inline bool PointsAre2D(const Eigen::MatrixXd& points) const
      {
        Gedim::Output::Assert(points.rows() == 3 && points.cols() > 0);
        return points.row(2).isZero(_configuration.Tolerance);
      }

      /// \brief compute the Point Curvilinear Coordinate of segment
      /// \param point the point
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the curvilinear coordinate computed
      inline double PointCurvilinearCoordinate(const Eigen::Vector3d& point,
                                               const Eigen::Vector3d& segmentOrigin,
                                               const Eigen::Vector3d& segmentEnd) const
      {
        const Eigen::Vector3d segmentTangent = (segmentEnd - segmentOrigin);
        return PointLineCurvilinearCoordinate(point,
                                              segmentOrigin,
                                              segmentTangent,
                                              segmentTangent.squaredNorm());
      }

      /// \brief compute the Point Curvilinear Coordinate of line
      /// \param point the point
      /// \param lineOrigin the line origin
      /// \param lineTangent the line tangent
      /// \param lineTangentSquaredLength the line tangent length squared
      /// \return the curvilinear coordinate computed
      inline double PointLineCurvilinearCoordinate(const Eigen::Vector3d& point,
                                                   const Eigen::Vector3d& lineOrigin,
                                                   const Eigen::Vector3d& lineTangent,
                                                   const double& lineTangentSquaredLength) const
      {
        return (point - lineOrigin).dot(lineTangent) / lineTangentSquaredLength;
      }

      /// \param point the point
      /// \param lineOrigin the line origin
      /// \param lineTangent the line tangent
      /// \param lineTangentSquaredLength the line tangent length squared
      /// \return true if the point belongs on line
      bool IsPointOnLine(const Eigen::Vector3d& point,
                         const Eigen::Vector3d& lineOrigin,
                         const Eigen::Vector3d& lineTangent,
                         const double& lineTangentSquaredLength) const;

      /// \brief Compute point position respect to a segment
      /// \param point the point
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return result the point position
      /// \warning left and right point positions work only in xy plane
      PointSegmentPositionTypes PointSegmentPosition(const Eigen::Vector3d& point,
                                                     const Eigen::Vector3d& segmentOrigin,
                                                     const Eigen::Vector3d& segmentEnd) const;

      /// \brief Compute point position on a segment line given the curvilinear Coordinate
      /// \param curvilinearCoordinate the curvilinear coordinate, segment is between 0.0 and 1.0
      /// \param result the point position on the line
      PointSegmentPositionTypes PointSegmentPosition(const double& curvilinearCoordinate) const;

      /// \brief Project point on a segment line
      /// \param point the point
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the projected point curvilinear coordinate
      inline double PointSegmentProjection(const Eigen::Vector3d& point,
                                           const Eigen::Vector3d& segmentOrigin,
                                           const Eigen::Vector3d& segmentEnd) const
      {
        return PointCurvilinearCoordinate(point, segmentOrigin, segmentEnd);
      }

      /// \brief Compute point position respect to a plane
      /// \param point the point
      /// \param planeNormal the plane normal
      /// \param planeOrigin the plane origin
      /// \return result the point position
      PointPlanePositionTypes PointPlanePosition(const Eigen::Vector3d& point,
                                                 const Eigen::Vector3d& planeNormal,
                                                 const Eigen::Vector3d& planeOrigin) const;

      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the segment length
      inline double SegmentLength(const Eigen::Vector3d& segmentOrigin,
                                  const Eigen::Vector3d& segmentEnd) const
      { return (segmentEnd - segmentOrigin).norm(); }

      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the segment tangent
      inline Eigen::Vector3d SegmentTangent(const Eigen::Vector3d& segmentOrigin,
                                            const Eigen::Vector3d& segmentEnd) const
      { return segmentEnd - segmentOrigin; }

      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the segment normal normalized, rotation of the normalized tangent (x,y,0) with 90° clockwise (y, -x,0)
      /// \note the segment shall be 2D
      inline Eigen::Vector3d SegmentNormal(const Eigen::Vector3d& segmentOrigin,
                                           const Eigen::Vector3d& segmentEnd) const
      {
        Gedim::Output::Assert(PointsAre2D(segmentOrigin) && PointsAre2D(segmentEnd));
        Eigen::Vector3d tangent = SegmentTangent(segmentOrigin, segmentEnd).normalized();
        return Eigen::Vector3d(tangent.y(), -tangent.x(), 0.0);
      }

      /// \brief Compute the segment slope m of line y = m * x + q
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the segment slope
      /// \note the segment shall be 2D
      inline double SegmentSlope(const Eigen::Vector3d& segmentOrigin,
                                 const Eigen::Vector3d& segmentEnd) const
      {
        Gedim::Output::Assert(!Are1DValuesEqual(segmentEnd.x(), segmentOrigin.x()));
        return (segmentEnd.y() - segmentOrigin.y()) / (segmentEnd.x() - segmentOrigin.x());
      }

      /// \brief Compute the segment intercept q of line y = m * x + q
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \return the segment intercept
      /// \note the segment shall be 2D
      inline double SegmentIntercept(const Eigen::Vector3d& segmentOrigin,
                                     const Eigen::Vector3d& segmentEnd) const
      {
        Gedim::Output::Assert(!Are1DValuesEqual(segmentEnd.x(), segmentOrigin.x()));
        return segmentOrigin.y() -
            segmentOrigin.x() * (segmentEnd.y() - segmentOrigin.y()) / (segmentEnd.x() - segmentOrigin.x());
      }

      /// \brief Check if two spheres do not intersect
      /// \param firstSphereCenter the first sphere center
      /// \param secondSphereCenter the second sphere center
      /// \param firstSphereDiameter the first sphere diameter
      /// \param secondSphereDiameter the second sphere diameter
      /// \return true if the two segments do not intersect
      /// \note if the function returns true it does not mean that the two segments intersects
      inline bool CheckNoSpheresIntersection(const Eigen::Vector3d& firstSphereCenter,
                                             const Eigen::Vector3d& secondSphereCenter,
                                             const double& firstSphereDiameter,
                                             const double& secondSphereDiameter) const
      {
        return IsValue1DGreater(2.0 * PointDistance(firstSphereCenter, secondSphereCenter),
                                firstSphereDiameter + secondSphereDiameter);
      }

      /// \param firstLineOrigin first line origin
      /// \param firstLineTangent first line tangent
      /// \param secondLineOrigin second line origin
      /// \param secondLineTangent second line tangent
      /// \return line coplanarity
      bool AreLineCoplanar(const Eigen::Vector3d& firstLineOrigin,
                           const Eigen::Vector3d& firstLineTangent,
                           const Eigen::Vector3d& secondLineOrigin,
                           const Eigen::Vector3d& secondLineTangent) const;

      /// \brief Compute the intersection between the two segments
      /// \param firstSegmentOrigin first segment origin
      /// \param firstSegmentEnd first segment end
      /// \param secondSegmentOrigin second segment origin
      /// \param secondSegmentEnd second segment end
      /// \return the resulting intersection
      /// \note no check is performed
      IntersectionSegmentSegmentResult IntersectionSegmentSegment(const Eigen::Vector3d& firstSegmentOrigin,
                                                                  const Eigen::Vector3d& firstSegmentEnd,
                                                                  const Eigen::Vector3d& secondSegmentOrigin,
                                                                  const Eigen::Vector3d& secondSegmentEnd) const;

      /// \brief Compute the intersection between a collection of segments
      /// \param segmentsVertices the segments vertices
      /// \param segmentsTangent the segments tangent
      /// \param segmentsBarycenter the segments barycenter
      /// \param segmentsLength the segments length
      /// \return for each segment the list of intersections curvilinear coordinate
      std::vector<std::list<double>> IntersectionsBetweenSegments(const std::vector<Eigen::MatrixXd>& segmentsVertices,
                                                                  const std::vector<Eigen::Vector3d>& segmentsTangent,
                                                                  const std::vector<Eigen::Vector3d>& segmentsBarycenter,
                                                                  const std::vector<double>& segmentsLength) const;

      /// \brief Compute the intersection between the a segment and a circle
      /// \param segmentOrigin first segment origin
      /// \param segmentEnd first segment end
      /// \param circleCenter circle center
      /// \param circleRadius circle radius
      /// \return the resulting intersection
      /// \note tested only in 2D
      IntersectionSegmentCircleResult IntersectionSegmentCircle(const Eigen::Vector3d& segmentOrigin,
                                                                const Eigen::Vector3d& segmentEnd,
                                                                const Eigen::Vector3d& circleCenter,
                                                                const double& circleRadius) const;

      /// \brief Intersection between a Segment, represented by origin and end and a plane
      /// represented by the normal and a point
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segement end
      /// \param planeNormal the plane normal normalized
      /// \param planeOrigin a plane point
      /// \return the resulting intersection
      IntersectionSegmentPlaneResult IntersectionSegmentPlane(const Eigen::Vector3d& segmentOrigin,
                                                              const Eigen::Vector3d& segmentEnd,
                                                              const Eigen::Vector3d& planeNormal,
                                                              const Eigen::Vector3d& planeOrigin) const;

      /// \brief Intersection between a Polyhedron and a Plane
      /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
      /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
      /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
      /// \param planeNormal the plane normal normalized
      /// \param planeOrigin the plane origin
      /// \param planeRotationMatrix the plane rotation from 3D to 2D
      /// \param planeTranslation the plane translation vector
      /// \return the intersection result
      /// \note works only with convex polyhedra
      IntersectionPolyhedronPlaneResult IntersectionPolyhedronPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                                    const Eigen::MatrixXi& polyhedronEdges,
                                                                    const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                    const Eigen::Vector3d& planeNormal,
                                                                    const Eigen::Vector3d& planeOrigin,
                                                                    const Eigen::Matrix3d& planeRotationMatrix,
                                                                    const Eigen::Vector3d& planeTranslation) const;

      GeometryUtilities::SplitPolyhedronWithPlaneResult SplitPolyhedronWithPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                                                 const Eigen::MatrixXi& polyhedronEdges,
                                                                                 const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                 const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                 const std::vector<Eigen::MatrixXd>& polyhedronFaceEdgeTangents,
                                                                                 const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                                 const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                                                                 const Eigen::Vector3d& planeNormal,
                                                                                 const Eigen::Vector3d& planeOrigin,
                                                                                 const Eigen::Matrix3d& planeRotationMatrix,
                                                                                 const Eigen::Vector3d& planeTranslation) const;

      std::vector<GeometryUtilities::Polyhedron> SplitPolyhedronWithPlaneResultToPolyhedra(const GeometryUtilities::SplitPolyhedronWithPlaneResult& result) const;

      /// \brief Intersection between a Polyhedron and a line
      /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
      /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
      /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
      /// \param lineTangent the line tangent
      /// \param lineOrigin the line origin
      /// \return the intersection result
      IntersectionPolyhedronLineResult IntersectionPolyhedronLine(const Eigen::MatrixXd& polyhedronVertices,
                                                                  const Eigen::MatrixXi& polyhedronEdges,
                                                                  const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                  const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                  const std::vector<bool>& polyhedronFaceNormalDirections,
                                                                  const Eigen::Vector3d& lineTangent,
                                                                  const Eigen::Vector3d& lineOrigin) const;

      /// \brief Intersection between a Polyhedron and a segment
      /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
      /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
      /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \param segmentTangent the segment tangent
      /// \param polyhedronLineIntersections the intersection between the polyhedron and the line of the segment
      /// \return the intersection result
      IntersectionPolyhedronLineResult IntersectionPolyhedronSegment(const Eigen::MatrixXd& polyhedronVertices,
                                                                     const Eigen::MatrixXi& polyhedronEdges,
                                                                     const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                     const Eigen::Vector3d& segmentOrigin,
                                                                     const Eigen::Vector3d& segmentEnd,
                                                                     const Eigen::Vector3d& segmentTangent,
                                                                     const IntersectionPolyhedronLineResult& polyhedronLineIntersections) const;

      /// \brief Intersection between a collectio of Polyhedrons and a segment
      /// \param polyhedrons the polyhedron collection
      /// \param polyhedronFaceNormals polyhedron face normals
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      /// \param segmentTangent the segment tangent
      /// \return the intersection result
      IntersectionPolyhedronsSegmentResult IntersectionPolyhedronsSegment(const std::vector<Polyhedron>& polyhedrons,
                                                                          const std::vector<std::vector<Eigen::Vector3d>>& polyhedronFaceNormals,
                                                                          const std::vector<std::vector<bool>>& polyhedronFaceNormalDirections,
                                                                          const Eigen::Vector3d& segmentOrigin,
                                                                          const Eigen::Vector3d& segmentEnd,
                                                                          const Eigen::Vector3d& segmentTangent) const;

      /// \brief Check if point is inside a polygon
      /// \param point the point
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param result the resulting position
      /// \warning works only in 2D
      PointPolygonPositionResult PointPolygonPosition(const Eigen::Vector3d& point,
                                                      const Eigen::MatrixXd& polygonVertices) const;

      LinePolygonPositionResult LinePolygonPosition(const Eigen::Vector3d& lineTangent,
                                                    const Eigen::Vector3d& lineOrigin,
                                                    const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Check if point is inside a polygon
      /// \param point the point
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param result the resulting position

      /// \brief Check if point is inside a polygon
      /// \param point the point
      /// \param polyhedronFaces the polyhedron faces, size numPolyhedronFaces
      /// \param polyhedronFaceVertices the polyhedron face 3D vertices, size numPolyhedronFaces
      /// \param polyhedronFaceRotatedVertices the polyhedron face 2D vertices, size numPolyhedronFaces
      /// \param polyhedronFaceNormals the polyhedron face normals
      /// \param polyhedronFaceNormalDirections the polyhedron face normal directions
      /// \param polyhedronFaceTranslations the polyhedron face translation from 2D to 3D
      /// \param polyhedronFaceRotationMatrices the polyhedron face rotation matrix from 2D to 3D
      /// \return the point position respect the polyhedron
      PointPolyhedronPositionResult PointPolyhedronPosition(const Eigen::Vector3d& point,
                                                            const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                            const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                            const std::vector<Eigen::MatrixXd>& polyhedronFaceRotatedVertices,
                                                            const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                            const std::vector<bool>& polyhedronFaceNormalDirections,
                                                            const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                            const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const;

      /// \brief Check if point is inside a circle
      /// \param point the point
      /// \param circleCenter the circle center
      /// \param circleRadius the circle radius
      /// \param result the resulting position
      /// \note tested only in 2D
      PointCirclePositionResult PointCirclePosition(const Eigen::Vector3d& point,
                                                    const Eigen::Vector3d& circleCenter,
                                                    const double& circleRadius) const;

      /// \brief Check if points are inside a circle
      /// \param points the matrix of points (size 3 x numVertices)
      /// \param circleCenter the circle center
      /// \param circleRadius the circle radius
      /// \param result the resulting positions
      /// \note tested only in 2D
      std::vector<PointCirclePositionResult> PointCirclePositions(const Eigen::MatrixXd& points,
                                                                  const Eigen::Vector3d& circleCenter,
                                                                  const double& circleRadius) const;

      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param circleCenter the circle center
      /// \param circleRadius the circle radius
      /// \param vertexPositions the polygon vertices positions respect the circle
      /// \param polygonCircleIntersections the polygon center intersections
      /// \return the Polygon Circle reciprocal position
      /// \note tested only in 2D
      PolygonCirclePositionTypes PolygonCirclePosition(const Eigen::MatrixXd& polygonVertices,
                                                       const Eigen::Vector3d& circleCenter,
                                                       const double& circleRadius,
                                                       const std::vector<PointCirclePositionResult>& vertexPositions,
                                                       const IntersectionPolygonCircleResult& polygonCircleIntersections) const;

      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param circleCenter the circle center
      /// \param circleRadius the circle radius
      /// \return the Polygon Circle reciprocal intersections
      /// \note tested only in 2D
      IntersectionPolygonCircleResult IntersectionPolygonCircle(const Eigen::MatrixXd& polygonVertices,
                                                                const Eigen::Vector3d& circleCenter,
                                                                const double& circleRadius) const;

      /// \brief Convex Polygon simple Triangulation from the first vertex
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \return the sub-division triangulation, size 1 x 3 * numTriangles
      /// \note works only for convex polygon
      std::vector<unsigned int> PolygonTriangulationByFirstVertex(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Convex Polygon simple Triangulation from an internal point
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \param point internal polygon point
      /// \return the sub-division triangulation, size 1 x 3 * numPolygonVertices,
      /// \note the internal point index is numPolygonVertices
      std::vector<unsigned int> PolygonTriangulationByInternalPoint(const Eigen::MatrixXd& polygonVertices,
                                                                    const Eigen::Vector3d& internalPoint) const;

      /// \brief Convex Polygon sub division by angle quadrant which intersects a polygon in a curved edge
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \param circleCenter the circle center from which the curved edge derives
      /// \param circleRadius the radius of the circle from which the curved edge derives
      /// \param curvedEdgeIndex curved edge index, from 0 to numPolygonVertices
      /// \return the sub-division polygons result
      PolygonDivisionByAngleQuadrantResult PolygonOutsideCircleDivisionByAngleQuadrant(const Eigen::MatrixXd& polygonVertices,
                                                                                       const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                       const Eigen::Vector3d& circleCenter,
                                                                                       const double& circleRadius,
                                                                                       const unsigned int& curvedEdgeIndex) const;

      /// \brief Convex Polygon sub division by angle quadrant which intersects a polygon in a curved edge
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \param circleCenter the circle center from which the curved edge derives
      /// \param circleRadius the radius of the circle from which the curved edge derives
      /// \param curvedEdgeIndex curved edge index, from 0 to numPolygonVertices
      /// \return the sub-division polygons result
      PolygonDivisionByAngleQuadrantResult PolygonInsideCircleDivisionByAngleQuadrant(const Eigen::MatrixXd& polygonVertices,
                                                                                      const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                      const Eigen::Vector3d& circleCenter,
                                                                                      const double& circleRadius,
                                                                                      const unsigned int& curvedEdgeIndex) const;

      /// \brief Convex Polygon sub division from a circle which intersects a polygon in a curved edge
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \param circleCenter the circle center from which the curved edge derives
      /// \param circleRadius the radius of the circle from which the curved edge derives
      /// \param curvedEdgeIndex curved edge index, from 0 to numPolygonVertices
      /// \return the sub-division polygons result
      /// \note the polygon should be inside the angle quadrant formed by the curved edge
      /// \note otherwise use PolygonDivisionByAngleQuadrant function to split the polygon
      PolygonDivisionByCircleResult PolygonDivisionByCircle(const Eigen::MatrixXd& polygonVertices,
                                                            const Eigen::MatrixXd& polygonEdgeTangents,
                                                            const Eigen::Vector3d& circleCenter,
                                                            const double& circleRadius,
                                                            const unsigned int& curvedEdgeIndex) const;

      /// \brief Circle division from Convex Polygon sub division which intersects a polygon in a curved edge
      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \param circleCenter the circle center from which the curved edge derives
      /// \param circleRadius the radius of the circle from which the curved edge derives
      /// \param curvedEdgeIndex curved edge index, from 0 to numPolygonVertices
      /// \return the sub-division circle result
      /// \note the polygon should be inside the angle quadrant formed by the curved edge
      CircleDivisionByPolygonResult CircleDivisionByPolygon(const Eigen::MatrixXd& polygonVertices,
                                                            const Eigen::MatrixXd& polygonEdgeTangents,
                                                            const Eigen::Vector3d& circleCenter,
                                                            const double& circleRadius,
                                                            const unsigned int& curvedEdgeIndex) const;

      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \return the polygon area
      /// \note the polygon shall be 2D
      double PolygonArea(const Eigen::MatrixXd& polygonVertices) const;

      /// \param polygonVertices the polygon vertices, size 3 x numPolygonVertices
      /// \return the polygon area
      double PolygonArea3D(const Eigen::MatrixXd& polygonVertices) const;


      /// \param polygonCentroid the centroid
      /// \param polygonTriangulationPoints the internal polygon sub-triangulation
      /// \return the polygon intertia tensor \int_E (x-x_C)^2
      Eigen::Matrix3d PolygonInertia(const Eigen::Vector3d& polygonCentroid,
                                     const std::vector<Eigen::Matrix3d>& polygonTriangulationPoints) const;

      /// \brief Split a polygon with n vertices numbered from 0 to n unclockwise given a segment contained inside
      /// \param input the input data
      /// \param result the resulting split
      /// \note only indices are threated in this function, no space points
      SplitPolygonWithSegmentResult SplitPolygonWithSegment(const SplitPolygonInput& input) const;

      /// \brief Split a polygon with n vertices numbered from 0 to n unclockwise given a cirle
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param circleCenter the circle center
      /// \param circleRadius the circle radius
      /// \param vertexPositions the polygon vertices positions respect the circle
      /// \param polygonCircleIntersections the polygon center intersections
      /// \param polygonCirclePosition the polygon position respect the circle
      /// \note tested only in 2D
      /// \return the split result
      /// \note only indices are threated in this function, no space points
      SplitPolygonWithCircleResult SplitPolygonWithCircle(const Eigen::MatrixXd& polygonVertices,
                                                          const Eigen::Vector3d& circleCenter,
                                                          const double& circleRadius,
                                                          const std::vector<PointCirclePositionResult>& vertexPositions,
                                                          const IntersectionPolygonCircleResult& polygonCircleIntersections,
                                                          const PolygonCirclePositionTypes& polygonCirclePosition) const;

      /// \brief Build the subpolygon coordinates from split result
      /// \param splitResult the split result
      /// \param subPolygonIndex the subpolygon index, from 0 to SplitPolygonWithCircleResult::NewPolygons.size()
      /// \param polygonVertices the original polygon vertices
      /// \param polygonCircleIntersections the polygon circle intersection
      /// \return the resulting subpolygon coordinates
      Eigen::MatrixXd SplitPolygonWithCircleBuildSubPolygon(const SplitPolygonWithCircleResult& splitResult,
                                                            const unsigned int& subPolygonIndex,
                                                            const Eigen::MatrixXd& polygonVertices,
                                                            const Gedim::GeometryUtilities::IntersectionPolygonCircleResult& polygonCircleIntersections) const;

      /// \brief Split 3d Polygon With Plane
      /// \param polygonVertices the 3D polygon vertices
      /// \param polygonEdgeTangents the 3D polygon edge tangents
      /// \param planeNormal the plane normal
      /// \param planeOrigin the plane origin
      /// \param polygonTranslation the polygon translation vector for rotation
      /// \param polygonRotationMatrix the polygon rotation matrix from 2D to 3D
      /// \return the splitted polygons
      SplitPolygonWithPlaneResult SplitPolygonWithPlane(const Eigen::MatrixXd& polygonVertices,
                                                        const Eigen::MatrixXd& polygonEdgeTangents,
                                                        const Eigen::Vector3d& planeNormal,
                                                        const Eigen::Vector3d& planeOrigin,
                                                        const Eigen::Vector3d& polygonTranslation,
                                                        const Eigen::Matrix3d& polygonRotationMatrix) const;

      /// \brief Compute the Polygon tridimensional normalized Normal
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \return the resulting normalized normal
      Eigen::Vector3d PolygonNormal(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the Polygon edge lengths
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \return the resulting edge lengths, size 1 x numVertices
      Eigen::VectorXd PolygonEdgeLengths(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the Polygon edge tangents
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \return the resulting edge tangents, size 3 x numVertices
      Eigen::MatrixXd PolygonEdgeTangents(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the Polygon edge normals outgoing the polygon
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \return the resulting edge normals outgoing the polygon, size 3 x numVertices
      Eigen::MatrixXd PolygonEdgeNormals(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the simplex barycenter as a mean of all vertices
      /// \param vertices the matrix of vertices of the simplex (size 3 x numVertices)
      inline Eigen::Vector3d SimplexBarycenter(const Eigen::MatrixXd& vertices) const
      {
        Gedim::Output::Assert(vertices.rows() == 3);
        return vertices.rowwise().mean();
      }

      /// \brief Compute the segment barycenter as a mean of all vertices
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segment end
      inline Eigen::Vector3d SegmentBarycenter(const Eigen::Vector3d& segmentOrigin,
                                               const Eigen::Vector3d& segmentEnd) const
      {
        return SimplexBarycenter((Eigen::MatrixXd(3, 2)<< segmentOrigin, segmentEnd).finished());
      }

      /// \brief Compute the Polygon barycenter as a mean of all vertices
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      inline Eigen::Vector3d PolygonBarycenter(const Eigen::MatrixXd& polygonVertices) const
      {
        Gedim::Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);
        return SimplexBarycenter(polygonVertices);
      }

      /// \brief Compute the polyhedron barycenter as a mean of all vertices
      /// \param polyhedronVertices the matrix of vertices of the polyhedron (size 3 x numVertices)
      inline Eigen::Vector3d PolyhedronBarycenter(const Eigen::MatrixXd& polyhedronVertices) const
      {
        Gedim::Output::Assert(polyhedronVertices.rows() == 3 && polyhedronVertices.cols() > 2);
        return SimplexBarycenter(polyhedronVertices);
      }

      /// \brief Compute the Polygon centroid as described in https://en.wikipedia.org/wiki/Centroid
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param polygonArea the area of the polygon
      /// \note the polygon shall be 2D
      Eigen::Vector3d PolygonCentroid(const Eigen::MatrixXd& polygonVertices,
                                      const double& polygonArea) const;

      /// \brief Compute the Polygon centroid using polygon sub-division
      /// \param subPolygonCentroids the centroid of each subPolygon (size 3 x numSubPolygons)
      /// \param subPolygonAreas the areas of each subPolygon, size 1 x numSubPolygons
      /// \param polygonArea the total area of the polygon
      inline Eigen::Vector3d PolygonCentroid(const Eigen::MatrixXd& subPolygonCentroids,
                                             const Eigen::VectorXd& subPolygonAreas,
                                             const double& polygonArea) const
      {
        Gedim::Output::Assert(subPolygonCentroids.rows() == 3 &&
                              subPolygonCentroids.cols() > 0 &&
                              subPolygonCentroids.cols() == subPolygonAreas.size());

        return subPolygonCentroids * subPolygonAreas / polygonArea;
      }

      /// \brief Polygon Area By Integral on edges
      /// \param polygonVertices the polygon vertices, size 3 x numVertices
      /// \param edgeLengths the edge lengths, size numEdges
      /// \param edgeTangents the edge tangents, size 3 x numEdges
      /// \param edgeNormals the edge outgoint normals, size 3 x numEdges
      /// \return the polygon area
      /// \note the area is computed as integral_edges x dot n_x with gauss formula on edges of order 1
      double PolygonAreaByIntegral(const Eigen::MatrixXd& polygonVertices,
                                   const Eigen::VectorXd& edgeLengths,
                                   const Eigen::MatrixXd& edgeTangents,
                                   const Eigen::MatrixXd& edgeNormals) const;

      /// \brief Polygon Area By Integral on edges
      /// \param polygonVertices the polygon vertices, size 3 x numVertices
      /// \param edgeLengths the edge lengths, size numEdges
      /// \param edgeTangents the edge tangents, size 3 x numEdges
      /// \param edgeNormals the edge outgoint normals, size 3 x numEdges
      /// \param polygonArea the polygon area
      /// \return the polygon centroid
      /// \note the area is computed as integral_edges (x^2, y^2) dot n with gauss formula on edges of order 2
      Eigen::Vector3d PolygonCentroidByIntegral(const Eigen::MatrixXd& polygonVertices,
                                                const Eigen::VectorXd& edgeLengths,
                                                const Eigen::MatrixXd& edgeTangents,
                                                const Eigen::MatrixXd& edgeNormals,
                                                const double& polygonArea) const;

      /// \param polygonVertices the polygon vertices, size 3 x numVertices
      /// \param polygonCentroid the polygon centroid
      /// \param polygonEdgeNormals the polygon edge normals outgoing the polygon, size 3 x numEdges
      /// \return the distance between the centroid and all the polygon edges, size 1 x numEdges
      Eigen::VectorXd PolygonCentroidEdgesDistance(const Eigen::MatrixXd& polygonVertices,
                                                   const Eigen::Vector3d& polygonCentroid,
                                                   const Eigen::MatrixXd& polygonEdgeNormals) const;

      /// \param polygonCentroidEdgesDistance the polygon centroid edges distance, size 1 x numEdges
      /// \return the polygon in radius, as the minimum distance between the polygon centroid and the edges
      inline double PolygonInRadius(const Eigen::VectorXd& polygonCentroidEdgesDistance) const
      { return polygonCentroidEdgesDistance.minCoeff(); }

      /// \param polygonDiameter the polygon diameter
      /// \param polygonInRadius the polygon in radius
      /// \return the polygon aspect ratio, defined as the ratio bewteen the in and out diameter
      inline double PolygonAspectRatio(const double& polygonDiameter,
                                       const double& polygonInRadius) const
      {
        return polygonDiameter / (2.0 * polygonInRadius);
      };

      /// \brief Compute the Polygon diameter defined as the maximum distance between the vertices
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      inline double PolygonDiameter(const Eigen::MatrixXd& polygonVertices) const
      {
        return PointsMaxDistance(polygonVertices);
      }

      /// \brief Compute the translation vector of a tridimensional Polygon
      /// \param polygonVertices the vertices of the polygon unclockwise (size 3 x numVertices)
      /// \return the resulting translation vector t which corresponds to the first vertex of the polygon
      /// \note to rotate some point P from 2D to 3D use Q * P + t
      /// \note to rotate some point P from 3D to 2D use Q^T * (P - t)
      inline Eigen::Vector3d PolygonTranslation(const Eigen::MatrixXd& polygonVertices) const
      {
        return polygonVertices.col(0);
      }

      /// \brief Compute the rotation matrix and translation vector of a tridimensional Polygon
      /// \param polygonVertices the vertices of the polygon unclockwise (size 3 x numVertices)
      /// \param polygonNormal the normalized normal of the plane which contains the polygon
      /// \param polygonTranslation the translation vector t
      /// \return the resulting rotation matrix Q which rotates 2D points to 3D points
      /// \note to rotate some point P from 2D to 3D use Q * P + t
      /// \note to rotate some point P from 3D to 2D use Q^T * (P - t)
      Eigen::Matrix3d PolygonRotationMatrix(const Eigen::MatrixXd& polygonVertices,
                                            const Eigen::Vector3d& polygonNormal,
                                            const Eigen::Vector3d& polygonTranslation) const;

      /// \brief Check if Polygon is Convex
      /// \param polygonVertices the polygon vertices, size 3 x numVertices
      /// \param convexHull the convex hull vertices unclockwise
      /// \return true if polygon is convex, false otherwise
      /// \note works only in 2D-plane
      bool PolygonIsConvex(const Eigen::MatrixXd& polygonVertices,
                           const Eigen::MatrixXd& convexHull) const;

      PolygonTypes PolygonType(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the rotation matrix of a plane from 2D to 3D
      /// \param planeNormal the normalized normal of the plane
      /// \return the resulting rotation matrix Q which rotates 2D points to 3D points
      /// \note to rotate some point P from 2D to 3D use Q * P + t
      /// \note to rotate some point P from 3D to 2D use Q^T * (P - t)
      Eigen::Matrix3d PlaneRotationMatrix(const Eigen::Vector3d& planeNormal) const;

      /// \brief Compute the translation vector of a plane from 2D to 3D
      /// \param planeNormal the normalized normal of the plane
      /// \param planeOrigin the 3D plane origin
      /// \return the resulting translation vector t which translates 2D points to 3D points
      /// \note to rotate some point P from 2D to 3D use Q * P + t
      /// \note to rotate some point P from 3D to 2D use Q^T * (P - t)
      inline Eigen::Vector3d PlaneTranslation(const Eigen::Vector3d& planeOrigin) const
      { return planeOrigin; }

      /// \brief Rotate Points P using rotation matrix Q and translation t: Q * P + t
      /// \param points the points (size 3 x numPoints)
      /// \param rotationMatrix the rotation matrix, size 3x3
      /// \param translation the translation vector, size 1x3
      /// \param rotatedPoints the resulting rotated points (size 3 x numPoints) rP = Q * P + t
      inline Eigen::MatrixXd RotatePoints(const Eigen::MatrixXd& points,
                                          const Eigen::Matrix3d& rotationMatrix,
                                          const Eigen::Vector3d& translation = Eigen::Vector3d::Zero()) const
      {
        Gedim::Output::Assert(points.rows() == 3 && points.cols() > 0);
        return (rotationMatrix * points).colwise() + translation;
      }

      /// \brief Rotate Points P From 2D To 3D using rotation matrix Q and translation t: Q * P + t
      /// \param points the points (size 3 x numPoints)
      /// \param rotationMatrix the rotation matrix from 2D to 3D
      /// \param translation the translation vector
      /// \param rotatedPoints the resulting rotated points (size 3 x numPoints) rP = Q * P + t
      inline Eigen::MatrixXd RotatePointsFrom2DTo3D(const Eigen::MatrixXd& points,
                                                    const Eigen::Matrix3d& rotationMatrix,
                                                    const Eigen::Vector3d& translation = Eigen::Vector3d::Zero()) const
      {
        Gedim::Output::Assert(points.rows() == 3 && points.cols() > 0 && PointsAre2D(points));
        return (rotationMatrix * points).colwise() + translation;
      }
      /// \brief Rotate Points P From 3D To 2D using rotation matrix Q and translation t: Q * (P - t)
      /// \param points the points (size 3 x numPoints)
      /// \param rotationMatrix the rotation matrix from 3D to 2D
      /// \param translation the translation vector
      /// \param rotatedPoints the resulting rotated points (size 3 x numPoints) rP = Q * (P - t)
      inline Eigen::MatrixXd RotatePointsFrom3DTo2D(const Eigen::MatrixXd& points,
                                                    const Eigen::Matrix3d& rotationMatrix,
                                                    const Eigen::Vector3d& translation = Eigen::Vector3d::Zero()) const
      {
        Gedim::Output::Assert(points.rows() == 3 && points.cols() > 0);
        Eigen::MatrixXd rotatedPoints = rotationMatrix * (points.colwise() - translation);
        rotatedPoints.row(2).setZero();
        return rotatedPoints;
      }

      /// \brief Compute the Convex Hull of 2D points
      /// \param points the points, size 3 x numPoints
      /// \param includeCollinear include the collinear points, default true
      /// \return the convex hull indices unclockwise, size numConvexHullPoints, numConvexHullPoints <= numPoints
      /// \note works in 2D, use the Graham_scan algorithm https://en.wikipedia.org/wiki/Graham_scan
      std::vector<unsigned int> ConvexHull(const Eigen::MatrixXd& points,
                                           const bool& includeCollinear = true) const;

      /// \brief Check if a set of points are aligned to a line identified by a segment
      /// \param segmentOrigin segment origin of the line
      /// \param segmentEnd segment end of the line
      /// \param points the points, size 3 x numPoints
      /// \return true if the i-th point is aligned, size 1 x numPoints
      inline std::vector<bool> PointsAreAligned(const Eigen::Vector3d& segmentOrigin,
                                                const Eigen::Vector3d& segmentEnd,
                                                const Eigen::MatrixXd& points) const
      {
        return PointsAreOnLine(points,
                               segmentOrigin,
                               SegmentTangent(segmentOrigin, segmentEnd));
      }

      /// \brief Check if a set of points are on a line
      /// \param points the points, size 3 x numPoints
      /// \param lineTangent the line tangent
      /// \param lineOrigin the line origin
      /// \return true if the i-th point is aligned, size 1 x numPoints
      std::vector<bool> PointsAreOnLine(const Eigen::MatrixXd& points,
                                        const Eigen::Vector3d& lineOrigin,
                                        const Eigen::Vector3d& lineTangent) const;

      /// \brief Check if a point is aligned to a line identified by a segment
      /// \param segmentOrigin segment origin of the line
      /// \param segmentEnd segment end of the line
      /// \param point the point
      /// \return true if the point is aligned
      inline bool PointIsAligned(const Eigen::Vector3d& segmentOrigin,
                                 const Eigen::Vector3d& segmentEnd,
                                 const Eigen::Vector3d& point) const
      { return PointsAreAligned(segmentOrigin, segmentEnd, point)[0]; }

      /// \brief Check if a point is on a line
      /// \param point the point
      /// \param lineTangent the line tangent
      /// \param lineOrigin the line origin
      /// \return true if the point is aligned
      inline bool PointIsOnLine(const Eigen::Vector3d& point,
                                const Eigen::Vector3d& lineOrigin,
                                const Eigen::Vector3d& lineTangent) const
      { return PointsAreOnLine(point, lineOrigin, lineTangent)[0]; }

      /// \brief Extract the circumscribed unaligned points (minimum 2) in a set of points
      /// \param points the points, size 3 x numPoints
      /// \return the unaligned points indices unclockwise, size numUnalignedPoints, 2 <= numUnalignedPoints <= numPoints
      std::vector<unsigned int> UnalignedPoints(const Eigen::MatrixXd& points) const;

      /// \param points the points, size 3 x numPoints
      /// \param filter indices unclockwise, size numFilterPoints, numFilterPoints <= numPoints
      /// \return the points coordinates filtered, size 3 x numFilterPoints
      Eigen::MatrixXd ExtractPoints(const Eigen::MatrixXd& points,
                                    const std::vector<unsigned int>& filter) const;

      /// \param points the points, size 3 x numPoints
      /// \param pointsTriangulation the polygon sub-division triangulation, size 1 x 3 * numTriangles
      /// \return the triangles coordinates, size 1 x numTriangles
      std::vector<Eigen::Matrix3d> ExtractTriangulationPoints(const Eigen::MatrixXd& points,
                                                              const std::vector<unsigned int>& pointsTriangulation) const;

      /// \param points the points, size 3 x numPoints
      /// \param externalPoint the external point coordinates
      /// \param pointsTriangulation the polygon sub-division triangulation, size 1 x 3 * numTriangles
      /// \return the triangles coordinates, size 1 x numTriangles
      std::vector<Eigen::Matrix3d> ExtractTriangulationPointsByInternalPoint(const Eigen::MatrixXd& points,
                                                                             const Eigen::Vector3d& internalPoint,
                                                                             const std::vector<unsigned int>& pointsTriangulation) const;

      /// \brief Create 2D Ellipse approximation with 2D polygon
      /// \param axisMajorLength the ellipse axis major length
      /// \param axisMinorLength the ellipse axis minor length
      /// \param resolution the number of points on each ellipse quadrant
      /// \return the polygon which approximate the ellipse
      /// \note the ellipse is centered in the origin and parallel to xy-axis
      Eigen::MatrixXd CreateEllipse(const double& axisMajorLength,
                                    const double& axisMinorLength,
                                    const unsigned int& resolution) const;

      /// \brief Create a triangle with points
      Eigen::MatrixXd CreateTriangle(const Eigen::Vector3d& p1,
                                     const Eigen::Vector3d& p2,
                                     const Eigen::Vector3d& p3) const;
      /// \brief Create a parallelogram with origin and dimension
      /// \param origin the origin
      /// \param lengthVector the length vector
      /// \param widthVector the width vector
      Eigen::MatrixXd CreateParallelogram(const Eigen::Vector3d& origin,
                                          const Eigen::Vector3d& lengthVector,
                                          const Eigen::Vector3d& widthVector) const;
      /// \brief Create a rectangle with origin and dimensions parallel to axis
      /// \param origin the origin
      /// \param base the base length
      /// \param height the height length
      inline Eigen::MatrixXd CreateRectangle(const Eigen::Vector3d& origin,
                                             const double& base,
                                             const double& height) const
      {
        return CreateParallelogram(origin,
                                   Eigen::Vector3d(base, 0.0, 0.0),
                                   Eigen::Vector3d(0.0, height, 0.0));
      }

      /// \brief Create a square with origin and dimensions parallel to axis
      /// \param origin the origin
      /// \param edgeLength the edge length
      inline Eigen::MatrixXd CreateSquare(const Eigen::Vector3d& origin,
                                          const double& edgeLength) const
      {
        return CreateRectangle(origin,
                               edgeLength,
                               edgeLength);
      }

      /// \brief Create a Tetrahedron with origin and dimension
      /// \param origin the origin
      /// \param lengthVector the length vector
      /// \param heightVector the heigth vector
      /// \param widthVector the width vector
      /// \return the tetrahedron created
      Polyhedron CreateTetrahedronWithOrigin(const Eigen::Vector3d& origin,
                                             const Eigen::Vector3d& lengthVector,
                                             const Eigen::Vector3d& heightVector,
                                             const Eigen::Vector3d& widthVector) const;

      /// \brief Create a Tetrahedron with the four vertices
      /// \param v1 the first vertex
      /// \param v2 the second vertex
      /// \param v3 the third vertex
      /// \param v4 the fourth vertex
      /// \return the tetrahedron created
      Polyhedron CreateTetrahedronWithVertices(const Eigen::Vector3d& v1,
                                               const Eigen::Vector3d& v2,
                                               const Eigen::Vector3d& v3,
                                               const Eigen::Vector3d& v4) const;

      /// \brief Create a Parallelepiped with origin and dimension
      /// \param origin the origin
      /// \param lengthVector the length vector
      /// \param heightVector the heigth vector
      /// \param widthVector the width vector
      /// \return the parallelepiped created
      Polyhedron CreateParallelepipedWithOrigin(const Eigen::Vector3d& origin,
                                                const Eigen::Vector3d& lengthVector,
                                                const Eigen::Vector3d& heightVector,
                                                const Eigen::Vector3d& widthVector) const;

      /// \brief Create a Cube with origin aligned to axis
      /// \param origin the origin
      /// \param edgeLength the edge length
      /// \return the cube created
      inline Polyhedron CreateCubeWithOrigin(const Eigen::Vector3d& origin,
                                             const double& edgeLength) const
      {
        return CreateParallelepipedWithOrigin(origin,
                                              Eigen::Vector3d(edgeLength, 0.0, 0.0),
                                              Eigen::Vector3d(0.0, 0.0, edgeLength),
                                              Eigen::Vector3d(0.0, edgeLength, 0.0));
      }

      /// \brief Compute the Polyhedron Volume
      /// \param polyhedronRotatedFaceTriangulationPoints polyhedron face triangulation points 2D, size numPolyhedronFaces x numTrianglesPerFace
      /// \param polyhedronFaceNormals polyhedron face normals, size numPolyhedronFaces
      /// \param polyhedronFaceNormalDirections polyhedron face normal directions, size numPolyhedronFaces
      /// \param polyhedronFaceTranslations polyhedron face translation vector from 2D to 3D
      /// \param polyhedronFaceRotationMatrices polyhedron face rotation matrix from 2D to 3D
      /// \return the polyhedron volume
      /// \note use the divergence theorem, with F = 1/3 (x, y, z), see https://en.wikipedia.org/wiki/Divergence_theorem
      double PolyhedronVolume(const std::vector<std::vector<Eigen::Matrix3d>>& polyhedronFaceRotatedTriangulationPoints,
                              const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                              const std::vector<bool>& polyhedronFaceNormalDirections,
                              const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                              const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const;

      /// \brief Compute the Polyhedron centroid
      /// \param polyhedronRotatedFaceTriangulationPoints polyhedron face triangulation points 2D, size numPolyhedronFaces x numTrianglesPerFace
      /// \param polyhedronFaceNormals polyhedron face normals, size numPolyhedronFaces
      /// \param polyhedronFaceNormalDirections polyhedron face normal directions, size numPolyhedronFaces
      /// \param polyhedronFaceTranslations polyhedron face translation vector from 2D to 3D
      /// \param polyhedronFaceRotationMatrices polyhedron face rotation matrix from 2D to 3D
      /// \param polyhedronVolume the polyhedron volume
      /// \return the polyhedron centroid
      /// \note use the divergence theorem, with F_x = 1/2 (x^2, 0, 0), F_y = 1/2 (0, y^2, 0), F_z = 1/2 (0, 0, z^2), see https://en.wikipedia.org/wiki/Divergence_theorem
      Eigen::Vector3d PolyhedronCentroid(const std::vector<std::vector<Eigen::Matrix3d>>& polyhedronFaceRotatedTriangulationPoints,
                                         const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                         const std::vector<bool>& polyhedronFaceNormalDirections,
                                         const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                         const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                         const double& polyhedronVolume) const;

      /// \brief Compute the Polyhedron diameter defined as the maximum distance between the vertices
      /// \param polyhedronVertices the matrix of vertices of the polyhedron (size 3 x numVertices)
      inline double PolyhedronDiameter(const Eigen::MatrixXd& polyhedronVertices) const
      {
        return PointsMaxDistance(polyhedronVertices);
      }

      /// \brief Compute Polyhedron Faces Vertices
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \return for each edge the tangent, size 3xnumEdges
      Eigen::MatrixXd PolyhedronEdgeTangents(const Eigen::MatrixXd& polyhedronVertices,
                                             const Eigen::MatrixXi& polyhedronEdges) const;

      /// \brief Compute Polyhedron Faces Vertices
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \param polyhedronFaces the polyhedron faces
      /// \return for each face the vertices, size 1xnumFaces
      std::vector<Eigen::MatrixXd> PolyhedronFaceVertices(const Eigen::MatrixXd& polyhedronVertices,
                                                          const std::vector<Eigen::MatrixXi>& polyhedronFaces) const;

      /// \brief Compute Polyhedron Faces Edge Direction
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \param polyhedronFaces the polyhedron faces
      /// \return for each face the edge direction compare to polyhedron edge directions, size 1xnumFaces
      std::vector<std::vector<bool>> PolyhedronFaceEdgeDirections(const Eigen::MatrixXd& polyhedronVertices,
                                                                  const Eigen::MatrixXi& polyhedronEdges,
                                                                  const std::vector<Eigen::MatrixXi>& polyhedronFaces) const;

      /// \brief Compute Polyhedron Faces Edge Tangents
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \param polyhedronFaces the polyhedron faces
      /// \param polyhedronEdgeTangents for each polyhedron edge the tangent, size 3xnumEdges
      /// \return for each face the edge tangents, size 1xnumFaces
      std::vector<Eigen::MatrixXd> PolyhedronFaceEdgeTangents(const Eigen::MatrixXd& polyhedronVertices,
                                                              const Eigen::MatrixXi& polyhedronEdges,
                                                              const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                              const std::vector<std::vector<bool>>& polyhedronFaceEdgeDirections,
                                                              const Eigen::MatrixXd& polyhedronEdgeTangents) const;

      /// \brief Compute Polyhedron Faces Rotation matrix
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \return for each polyhedron face the rotation matrix from 2D to 3D
      std::vector<Eigen::Matrix3d> PolyhedronFaceRotationMatrices(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                  const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                  const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations) const;

      /// \brief Compute Polyhedron Faces translation vectors
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \return for each polyhedron face the translation vector
      std::vector<Eigen::Vector3d> PolyhedronFaceTranslations(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices) const;

      /// \brief Compute Polyhedron Faces Rotated Vertices 2D
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronFaceTranslations the polyhedron face translations from 2D to 3D
      /// \param polyhedronFaceRotationMatrices the polyhedron face rotation matrix from 2D to 3D
      /// \return for each face the 2D vertices, size 1xnumFaces
      std::vector<Eigen::MatrixXd> PolyhedronFaceRotatedVertices(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                 const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                 const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const;


      /// \brief Compute Polyhedron Faces Normals
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \return for each polyhedron face the normal
      std::vector<Eigen::Vector3d> PolyhedronFaceNormals(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices) const;

      /// \brief Compute Polyhedron Faces barycenters
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \return for each polyhedron face the barycenter
      std::vector<Eigen::Vector3d> PolyhedronFaceBarycenter(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices) const;

      /// \brief Check if Polyhedron is Convex
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \param polyhedronFaceInternalPoints the polyhedron face internal points
      /// \param polyhedronFaceNormals the normal of each face
      /// \param polyhedronFaceNormalDirections the normal outgoing direction
      /// \param pointInsidePolyhedron a point inside polyhedron
      /// \return true if polyhedron is convex, false otherwise
      /// \warning still not working
      bool PolyhedronIsConvex(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                              const std::vector<Eigen::MatrixXd>& polyhedronFaceRotatedVertices,
                              const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                              const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                              const std::vector<bool>& polyhedronFaceNormalDirections,
                              const Eigen::Vector3d& polyhedronInternalPoint) const;

      /// \brief Compute Polyhedron Face Normal Directions
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \param pointInsidePolyhedron a point inside polyhedron
      /// \param polyhedronFaceNormals the normal of each face
      /// \return true if the face has normal outgoing
      std::vector<bool> PolyhedronFaceNormalDirections(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                       const Eigen::Vector3d& pointInsidePolyhedron,
                                                       const std::vector<Eigen::Vector3d>& polyhedronFaceNormals) const;
      /// \brief Polyhedron Face Triangulations of each face
      /// \param polyhedronFaces the polyhedron faces
      /// \param localFaceTriangulations the local faces triangulations indices, size 1xnumFaces x (3xnumTriangles)
      /// \return for each face the triangulation indices by first vertex, size 1xnumFaces x (3xnumTriangles)
      std::vector<std::vector<unsigned int>> PolyhedronFaceTriangulations(const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                          const std::vector<std::vector<unsigned int>>& localFaceTriangulations) const;
      /// \brief Polyhedron Face Triangulations by first vertex of each face
      /// \param polyhedronFaces the polyhedron faces
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \return for each face the triangulation indices by first vertex, size 1xnumFaces x (3xnumTriangles)
      std::vector<std::vector<unsigned int>> PolyhedronFaceTriangulationsByFirstVertex(const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                       const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices) const;

      std::vector<std::vector<Eigen::Matrix3d>> PolyhedronFaceTriangulationPointsByFirstVertex(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                               const std::vector<std::vector<unsigned int>>& polyhedronFaceTriangulations) const;

      std::vector<std::vector<Eigen::Matrix3d>> PolyhedronFaceTriangulationPointsByInternalPoint(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                                 const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                                                 const std::vector<std::vector<unsigned int>>& polyhedronFaceTriangulations) const;

      /// \brief Polyhedron Face Triangulations by internal point of each face
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronFaces the polyhedron faces
      /// \param polyhedronFaceVertices the polyhedron faces vertices
      /// \param polyhedronFaceInternalPoints the polyhedron face internal points
      /// \return for each face the triangulation indices by first vertex, size 1xnumFaces x (3xnumTriangles)
      /// \note the internal point index is polyhedronVertices.size()
      std::vector<std::vector<unsigned int>> PolyhedronFaceTriangulationsByInternalPoint(const Eigen::MatrixXd& polyhedronVertices,
                                                                                         const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                         const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                         const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints) const;

      /// \brief Polyhedron Tetrahedrons By Face Triangulations
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronFaces the polyhedron faces
      /// \param faceTriangulations the triangulation on face vertices
      /// \param polyhedronInternalPoint a polyhedron internal point
      /// \return the polyhedron tetrahedrons indices, size 1x(4*numTetrahedrons)
      /// \note the polyhedron internal point index is polyhedronVertices.size() + f
      std::vector<unsigned int> PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                           const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                           const std::vector<std::vector<unsigned int>>& polyhedronFaceTriangulations,
                                                                           const Eigen::Vector3d& polyhedronInternalPoint) const;

      /// \brief Polyhedron Tetrahedrons By Face Triangulations with face internal points
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronFaces the polyhedron faces
      /// \param faceTriangulations the triangulation on face vertices by internal points
      /// \param polyhedronFaceInternalPoints the polyhedron face internal points
      /// \param polyhedronInternalPoint a polyhedron internal point
      /// \return the polyhedron tetrahedrons indices, size 1x(4*numTetrahedrons)
      /// \note the polyhedron face internal points are polyhedronVertices.size() + f
      /// \note the polyhedron internal point index is polyhedronVertices.size() + polyhedronFaceInternalPoints.size()
      std::vector<unsigned int> PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                           const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                           const std::vector<std::vector<unsigned int>>& polyhedronFaceTriangulations,
                                                                           const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                           const Eigen::Vector3d& polyhedronInternalPoint) const;

      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronInternalPoint a polyhedron internal point
      /// \param pointTetrahedrons the polyhedron sub-division tetrahedrons, size 1 x 4 * numTetra
      /// \return the tetrahedrons coordinates, size 1 x numTetra
      std::vector<Eigen::MatrixXd> ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                            const Eigen::Vector3d& polyhedronInternalPoint,
                                                            const std::vector<unsigned int>& pointTetrahedrons) const;
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronInternalPoint a polyhedron internal point
      /// \param polyhedronFaceInternalPoints the polyhedron face internal points
      /// \param pointTetrahedrons the polyhedron sub-division tetrahedrons, size 1 x 4 * numTetra
      /// \return the tetrahedrons coordinates, size 1 x numTetra
      std::vector<Eigen::MatrixXd> ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                            const Eigen::Vector3d& polyhedronInternalPoint,
                                                            const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                            const std::vector<unsigned int>& pointTetrahedrons) const;

      /// \brief Get Polyhedron Coordinate System
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \return the four vertices indices forming a coordinate system for the polyhedron, size 1x4
      std::vector<unsigned int> PolyhedronCoordinateSystem(const Eigen::MatrixXd& polyhedronVertices,
                                                           const Eigen::MatrixXi& polyhedronEdges);

      /// \brief Export Polyhedron To VTU
      /// \param polyhedronVertices the polyhedron vertices
      /// \param polyhedronEdges the polyhedron edges
      /// \param polyhedronFaces the polyhedron faces
      /// \param exportFolder the folder in which to export
      void ExportPolyhedronToVTU(const Eigen::MatrixXd& polyhedronVertices,
                                 const Eigen::MatrixXi& polyhedronEdges,
                                 const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                 const std::string& exportFolder) const;
  };
}

#endif // __GEOMETRYUTILITIES_H
