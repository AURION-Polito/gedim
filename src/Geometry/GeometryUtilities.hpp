#ifndef __GEOMETRYUTILITIES_H
#define __GEOMETRYUTILITIES_H

#include "Eigen"
#include <iostream>

using namespace std;

namespace Gedim
{
  struct GeometryUtilitiesConfig final
  {
      double Tolerance = numeric_limits<double>::epsilon();
  };

  /// \brief The GeometryUtilities class intersects 3D segments
  class GeometryUtilities final
  {
    private:
      const GeometryUtilitiesConfig& _configuration;

    public:
      enum struct CompareTypes {
        Unknown = 0,
        FirstBeforeSecond = 1,
        Coincident = 2,
        SecondBeforeFirst = 3
      };

      enum struct PointSegmentPositionTypes {
        Unknown = 0,
        OnSegmentLineBeforeOrigin = 1,
        OnSegmentOrigin = 2,
        InsideSegment = 3,
        OnSegmentEnd = 4,
        OnSegmentLineAfterEnd = 5,
        LeftTheSegment = 6,
        RightTheSegment = 7
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
          vector<AlignedEdge> AlignedEdges;
          SplitSegment Segment;
      };

      struct SplitPolygonResult final
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
              vector<unsigned int> Cell2DNeighbours = {};
          };

          struct NewPolygon final
          {
              list<unsigned int> Vertices = {};
              list<unsigned int> Edges = {};
          };

          Types Type = Types::Unknown;
          list<NewVertex> NewVertices = {};
          list<NewEdge> NewEdges = {};
          vector<NewPolygon> NewPolygons = {};
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
              double CurvilinearCoordinate = -1.0;
          };

          IntersectionLineTypes IntersectionLinesType = IntersectionLineTypes::Unknown;
          IntersectionSegmentTypes IntersectionSegmentsType = IntersectionSegmentTypes::Unknown;
          /// \brief relation between first and second intersection.
          /// Values are the indeces of the SecondSegmentIntersections vector respect the FirstSegmentIntersections vector
          /// \example in MultipleIntersection case, if SecondIntersectionRelation[0] = 1,
          /// then the second intersection point SecondSegmentIntersections[1] is equal to FirstSegmentIntersections[0] point
          vector<unsigned int> SecondIntersectionRelation;
          /// \brief intersections of the first segment,
          /// \note if multiple intersections are found, than the origin and the end coordinate are stored
          vector<IntersectionPosition> FirstSegmentIntersections;
          /// \brief intersections of the second segment,
          /// \note if multiple intersections are found, than the origin and the end coordinate are stored
          vector<IntersectionPosition> SecondSegmentIntersections; /// intersections of the second segment
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
              double CurvilinearCoordinate = -1.0;
          };

          Types Type = Types::Unknown; ///< The intersection type
          Intersection SingleIntersection; ///< The single intersection, available only is Type is SingleIntersection
      };

      struct PointPolygonPositionResult final
      {
          enum struct PositionTypes
          {
            Unknown = 0,
            Outside = 1,
            BorderEdge = 2,
            BorderVertex = 3,
            Inside = 4
          };

          unsigned int BorderIndex = 0; ///< index of vertex/edge of border
          PositionTypes PositionType = PositionTypes::Unknown;
      };

    private:
      /// \brief Compare two values according to tolerance
      /// \param first the first value
      /// \param second the second value
      /// \return the result
      CompareTypes CompareValues(const double& first,
                                 const double& second,
                                 const double& tolerance) const;
    public:

      GeometryUtilities(const GeometryUtilitiesConfig& configuration);
      ~GeometryUtilities();

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

      /// \param value the value
      /// \return true if value is positive
      inline bool IsValue1DPositive(const double& value) const
      {
        return Compare1DValues(value,
                               0.0) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue1DNegative(const double& value) const
      {
        return Compare1DValues(value,
                               0.0) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue1DZero(const double& value) const
      {
        return Compare1DValues(value,
                               0.0) == CompareTypes::Coincident;
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
        return Compare2DValues(value,
                               0.0) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue2DNegative(const double& value) const
      {
        return Compare2DValues(value,
                               0.0) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue2DZero(const double& value) const
      {
        return Compare2DValues(value,
                               0.0) == CompareTypes::Coincident;
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
        return Compare3DValues(value,
                               0.0) == CompareTypes::SecondBeforeFirst;
      }

      /// \param value the value
      /// \return true if value is negative
      inline bool IsValue3DNegative(const double& value) const
      {
        return Compare3DValues(value,
                               0.0) == CompareTypes::FirstBeforeSecond;
      }

      /// \param value the value
      /// \return true if value is zero
      inline bool IsValue3DZero(const double& value) const
      {
        return Compare3DValues(value,
                               0.0) == CompareTypes::Coincident;
      }

      /// \brief compute the Point distance
      /// \param firstPoint the first point
      /// \param secondPoint the second point
      /// \return the distance
      inline double PointDistance(const Eigen::Vector3d& firstPoint,
                                  const Eigen::Vector3d& secondPoint) const
      {
        return (secondPoint - firstPoint).norm();
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
        return (point - segmentOrigin).dot(segmentEnd - segmentOrigin) / (segmentEnd - segmentOrigin).squaredNorm();
      }

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

      /// \brief Intersection between a Segment, represented by origin and end and a plane
      /// represented by the normal and a point
      /// \param segmentOrigin the segment origin
      /// \param segmentEnd the segement end
      /// \param planeNormal the plane normal
      /// \param planeOrigin a plane point
      /// \return the resulting intersection
      IntersectionSegmentPlaneResult IntersectionSegmentPlane(const Eigen::Vector3d& segmentOrigin,
                                                              const Eigen::Vector3d& segmentEnd,
                                                              const Eigen::Vector3d& planeNormal,
                                                              const Eigen::Vector3d& planeOrigin) const;

      /// \brief Check if point is inside a polygon
      /// \param point the point
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param result the resulting position
      PointPolygonPositionResult PointPolygonPosition(const Eigen::Vector3d& point,
                                                      const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Split a polygon with n vertices numbered from 0 to n unclockwise given a segment contained inside
      /// \param input the input data
      /// \param result the resulting split
      /// \note only indices are threated in this function, no space points
      SplitPolygonResult SplitPolygon(const SplitPolygonInput& input) const;

      /// \brief Compute the Polygon tridimensional normalized Normal
      /// \param polygonVertices the matrix of vertices of the polygon (size 3 x numVertices)
      /// \param normal the resulting normalized normal
      Eigen::Vector3d PolygonNormal(const Eigen::MatrixXd& polygonVertices) const;

      /// \brief Compute the rotation matrix and translation vector of a tridimensional Polygon
      /// \param polygonVertices the vertices of the polygon (size 3 x numVertices)
      /// \param normal the normalized normal of the plane which contains the polygon
      /// \param rotationMatrix the resulting rotation matrix Q which rotates 2D points to 3D points
      /// \param translation the resulting translation vector t which corresponds to the first vertex of the polygon
      /// \note to rotate some point P from 2D to 3D use Q * P + t
      /// \note to rotate some point P from 3D to 2D use Q^T * (P - t)
      void PolygonRotation(const Eigen::MatrixXd& polygonVertices,
                           const Eigen::Vector3d& normal,
                           Eigen::Matrix3d& rotationMatrix,
                           Eigen::Vector3d& translation) const;

      /// \brief Rotate Points P From 2D To 3D using rotation matrix Q and translation t: Q * P + t
      /// \param points the points (size 3 x numPoints)
      /// \param rotationMatrix the rotation matrix from 2D to 3D
      /// \param translation the translation vector
      /// \param rotatedPoints the resulting rotated points (size 3 x numPoints) rP = Q * P + t
      inline Eigen::MatrixXd RotatePointsFrom2DTo3D(const Eigen::MatrixXd& points,
                                                    const Eigen::Matrix3d& rotationMatrix,
                                                    const Eigen::Vector3d& translation) const
      {
        return (rotationMatrix * points).colwise() + translation;
      }
      /// \brief Rotate Points P From 3D To 2D using rotation matrix Q and translation t: Q * (P - t)
      /// \param points the points (size 3 x numPoints)
      /// \param rotationMatrix the rotation matrix from 3D to 2D
      /// \param translation the translation vector
      /// \param rotatedPoints the resulting rotated points (size 3 x numPoints) rP = Q * (P - t)
      inline Eigen::MatrixXd RotatePointsFrom3DTo2D(const Eigen::MatrixXd& points,
                                                    const Eigen::Matrix3d& rotationMatrix,
                                                    const Eigen::Vector3d& translation) const
      {
        Eigen::MatrixXd rotatedPoints = rotationMatrix * (points.colwise() - translation);
        rotatedPoints.row(2).setZero();
        return rotatedPoints;
      }

      ///
      /// \brief Compute the Convex Hull of 2D points
      /// \param points the points, size 3 x numPoints
      /// \return the convex hull indices, size numConvexHullPoints, numConvexHullPoints <= numPoints
      /// \note works in 2D, use the Gift wrapping algorithm (see https://en.wikipedia.org/wiki/Gift_wrapping_algorithm)
      vector<unsigned int> ConvexHull(const Eigen::MatrixXd& points);
  };
}

#endif // __GEOMETRYUTILITIES_H
