#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"
#include <unordered_set>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // **************************************************************
  GeometryUtilities::GeometryUtilities(const GeometryUtilitiesConfig& configuration) :
    _configuration(configuration)
  {
  }
  GeometryUtilities::~GeometryUtilities()
  {
  }
  // ***************************************************************************
  vector<double> GeometryUtilities::EquispaceCoordinates(const double& step,
                                                         const bool& insertExtremes) const
  {
    Output::Assert(IsValue1DPositive(step) &&
                   Compare1DValues(step, 1.0) != CompareTypes::SecondBeforeFirst);

    VectorXd generated = VectorXd::LinSpaced(static_cast<unsigned int>(1.0 / step + 0.5) + 1, 0.0, 1.0);
    vector<double> coordinates;
    if (insertExtremes)
      coordinates.resize(generated.size());
    else
      coordinates.resize(generated.size() - 2);

    for (unsigned int c = 0; c < coordinates.size(); c++)
      coordinates[c] = insertExtremes ? generated[c] : generated[c + 1];

    return coordinates;
  }
  // ***************************************************************************
  std::vector<double> GeometryUtilities::EquispaceCoordinates(const unsigned int& numLinSpace,
                                                              const double& origin,
                                                              const double& end,
                                                              const bool& insertExtremes) const
  {
    VectorXd generated = VectorXd::LinSpaced(numLinSpace, origin, end);
    vector<double> coordinates;
    if (insertExtremes)
      coordinates.resize(generated.size());
    else
      coordinates.resize(generated.size() - 2);

    for (unsigned int c = 0; c < coordinates.size(); c++)
      coordinates[c] = insertExtremes ? generated[c] : generated[c + 1];

    return coordinates;
  }
  // ***************************************************************************
  GeometryUtilities::CompareTypes GeometryUtilities::CompareValues(const double& first,
                                                                   const double& second,
                                                                   const double& tolerance) const
  {
    double relativeValue = (first == 0.0 ||
                            second == 0.0) ? 1.0 :
                                             abs(first);
    double difference = second - first;

    if (abs(difference) <= tolerance * relativeValue)
      return CompareTypes::Coincident;
    else if (difference < -tolerance * relativeValue)
      return CompareTypes::SecondBeforeFirst;
    else
      return CompareTypes::FirstBeforeSecond;
  }
  // ***************************************************************************
  Matrix3d GeometryUtilities::PlaneRotationMatrix(const Eigen::Vector3d& planeNormal) const
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
  vector<unsigned int> GeometryUtilities::ConvexHull(const Eigen::MatrixXd& points) const
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

    Output::Assert(points.rows() == 3 && points.cols() > 2 && PointsAre2D(points));

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
      unsigned int numPointTested = 0;

      while (numPointTested < numPoints) // to avoid aligned points
      {
        for (unsigned int j = 0; j < numPoints; j++)
        {
          if (endpoint == pointOnHull)
          {
            endpoint = j;
            continue;
          }

          const  GeometryUtilities::PointSegmentPositionTypes pointPosition =  PointSegmentPosition(points.col(j),
                                                                                                    points.col(pointOnHull),
                                                                                                    points.col(endpoint));

          if (pointPosition == GeometryUtilities::PointSegmentPositionTypes::RightTheSegment)
          {
            endpoint = j;
            numPointTested = 0;
            break;
          }

          numPointTested++;
        }
      }
      pointOnHull = endpoint;
    }
    while (pointOnHull != convexHull.front());

    return vector<unsigned int>(convexHull.begin(), convexHull.end());
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::ExtractPoints(const Eigen::MatrixXd& points,
                                            const vector<unsigned int>& filter) const
  {
    Eigen::MatrixXd extraction(3, filter.size());
    for (unsigned int c = 0; c < filter.size(); c++)
      extraction.col(c) = points.col(filter[c]);

    return extraction;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::ExtractTriangulationPoints(const Eigen::MatrixXd& points,
                                                                 const vector<unsigned int>& pointsTriangulation) const
  {
    const unsigned int numTriangles = pointsTriangulation.size() / 3;
    vector<Matrix3d> triangulations(numTriangles);

    for (unsigned int t = 0; t < numTriangles; t++)
    {
      Eigen::Matrix3d& triangleVertices = triangulations[t];
      triangleVertices.col(0)<< points.col(pointsTriangulation[3 * t]);
      triangleVertices.col(1)<< points.col(pointsTriangulation[3 * t + 1]);
      triangleVertices.col(2)<< points.col(pointsTriangulation[3 * t + 2]);
    }

    return triangulations;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::ExtractTriangulationPointsByInternalPoint(const Eigen::MatrixXd& points,
                                                                                const Eigen::Vector3d& internalPoint,
                                                                                const vector<unsigned int>& pointsTriangulation) const
  {
    const unsigned int numTriangles = pointsTriangulation.size() / 3;
    vector<Matrix3d> triangulations(numTriangles);

    for (unsigned int t = 0; t < numTriangles; t++)
    {
      Eigen::Matrix3d& triangleVertices = triangulations[t];
      triangleVertices.col(0)<< internalPoint;
      triangleVertices.col(1)<< points.col(pointsTriangulation[3 * t + 1]);
      triangleVertices.col(2)<< points.col(pointsTriangulation[3 * t + 2]);
    }

    return triangulations;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::CreateTriangle(const Eigen::Vector3d& p1,
                                             const Eigen::Vector3d& p2,
                                             const Eigen::Vector3d& p3) const
  {
    MatrixXd vertices(3, 3);
    vertices.col(0)<< p1;
    vertices.col(1)<< p2;
    vertices.col(2)<< p3;

    return vertices;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::CreateParallelogram(const Eigen::Vector3d& origin,
                                                  const Eigen::Vector3d& lengthVector,
                                                  const Eigen::Vector3d& widthVector) const
  {
    MatrixXd vertices(3, 4);
    vertices.col(0)<< origin;
    vertices.col(1)<< origin + lengthVector;
    vertices.col(2)<< origin + lengthVector + widthVector;
    vertices.col(3)<< origin + widthVector;
    return vertices;
  }
  // ***************************************************************************
}
