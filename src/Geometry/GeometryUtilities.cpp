#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::GeometryUtilities(const GeometryUtilitiesConfig& configuration) :
    _configuration(configuration)
  {
  }
  GeometryUtilities::~GeometryUtilities()
  {
  }

  Eigen::VectorXd GeometryUtilities::PointDistances(const Eigen::MatrixXd& points,
                                                   const Eigen::Vector3d& point) const
  {
    return (points.colwise() - point).colwise().norm();
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::FindPointInPoints(const Eigen::MatrixXd& points,
                                                            const Eigen::Vector3d& point) const
  {
    VectorXd pointDistances = PointDistances(points,
                                             point);
    list<unsigned int> indices;
    for (unsigned int p = 0; p < pointDistances.size(); p++)
    {
      if (IsValue1DZero(pointDistances[p]))
        indices.push_back(p);
    }

    return vector<unsigned int>(indices.begin(), indices.end());
  }
  // ***************************************************************************
  GeometryUtilities::CompareTypes GeometryUtilities::CompareValues(const double& first,
                                                                   const double& second,
                                                                   const double& tolerance) const
  {
    double maxFirstSecondOne = std::max( { 1.0, abs(first) , abs(second) } ) ;
    double difference = second - first;

    if (abs(difference) <= tolerance * maxFirstSecondOne)
      return CompareTypes::Coincident;
    else if (difference < -tolerance * maxFirstSecondOne)
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
  Vector3d GeometryUtilities::PlaneTranslation(const Eigen::Vector3d& planeNormal,
                                               const Eigen::Vector3d& planeOrigin,
                                               const Eigen::Matrix3d& planeRotationMatrix) const
  {
    Vector3d planeOrigin2D = RotatePointsFrom3DTo2D(planeOrigin,
                                                    planeRotationMatrix.transpose());

    return (planeOrigin - RotatePointsFrom2DTo3D(planeOrigin2D,
                                                 planeRotationMatrix));
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

    Output::Assert(points.rows() == 3 && points.cols() > 0 && PointsAre2D(points));
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
      for (unsigned int j = 0; j < numPoints; j++)
      {
        if ((endpoint == pointOnHull) ||
            PointSegmentPosition(points.col(j),
                                 points.col(pointOnHull),
                                 points.col(endpoint)) == GeometryUtilities::PointSegmentPositionTypes::RightTheSegment)
          endpoint = j;
      }
      pointOnHull = endpoint;
    }
    while (pointOnHull != convexHull.front());

    return vector<unsigned int>(convexHull.begin(), convexHull.end());
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::UnalignedPoints(const Eigen::MatrixXd& points) const
  {
    Output::Assert(points.rows() == 3 && points.cols() > 0);
    const unsigned int numPoints = points.cols();

    list<unsigned int> unalignedPoints;

    unalignedPoints.push_back(0);
    for (unsigned int p = 0; p < numPoints; p++)
    {

    }

    return vector<unsigned int>(unalignedPoints.begin(), unalignedPoints.end());
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::ExtractConvexHull(const Eigen::MatrixXd& points,
                                                const vector<unsigned int>& convexHull) const
  {
    Eigen::MatrixXd extraction(3, convexHull.size());
    for (unsigned int c = 0; c < convexHull.size(); c++)
      extraction.col(c) = points.col(convexHull[c]);

    return extraction;
  }
  // ***************************************************************************
}
