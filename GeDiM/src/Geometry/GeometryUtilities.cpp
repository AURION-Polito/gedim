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

    return insertExtremes ?
          EquispaceCoordinates(static_cast<unsigned int>(1.0 / step + 0.5) + 1, 0.0, 1.0, true) :
          EquispaceCoordinates(static_cast<unsigned int>(1.0 / step + 0.5) - 1, 0.0, 1.0, false);
  }
  // ***************************************************************************
  std::vector<double> GeometryUtilities::EquispaceCoordinates(const unsigned int& size,
                                                              const double& origin,
                                                              const double& end,
                                                              const bool& insertExtremes) const
  {
    if (size == 0)
      return { };
    else if (size == 1 && insertExtremes)
      throw invalid_argument("size is not valid with false insertExtremes");

    const VectorXd generated = insertExtremes ? VectorXd::LinSpaced(size,
                                                                    origin,
                                                                    end) :
                                                VectorXd::LinSpaced(size + 2,
                                                                    origin,
                                                                    end);

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
    const double max_tolerance = std::max(abs(tolerance),
                                          std::numeric_limits<double>::epsilon());

    double relativeValue = (first == 0.0 ||
                            second == 0.0) ? 1.0 :
                                             abs(first);
    double difference = second - first;

    if (abs(difference) <= max_tolerance * relativeValue)
      return CompareTypes::Coincident;
    else if (difference < -max_tolerance * relativeValue)
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
  vector<unsigned int> GeometryUtilities::ConvexHull(const Eigen::MatrixXd& points,
                                                     const bool& includeCollinear) const
  {
    Output::Assert(points.rows() == 3 && PointsAre2D(points));

    list<unsigned int> convexHull;
    const unsigned int numPoints = points.cols();

    if (numPoints == 1)
      return vector<unsigned int> { 0 };

    struct pt
    {
        double x;
        double y;
        unsigned int id;
    };

    vector<pt> structPoints(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
      structPoints[p] = { points.col(p).x(), points.col(p).y(), p };

    pt p0 = { std::numeric_limits<double>::max(),
              std::numeric_limits<double>::max(),
              0
            };
    for (const auto& p : structPoints)
    {
      if (IsValue1DGreater(p0.y, p.y))
        p0 = p;

      if (IsValue1DZero(abs(p0.y - p.y)) &&
          IsValue1DGreater(p0.x, p.x))
        p0 = p;
    }

    sort(structPoints.begin(),
         structPoints.end(),
         [&p0, this](const pt& a, const pt& b)
    {
      const double orientation = p0.x * (a.y - b.y) +
                                 a.x * (b.y - p0.y) +
                                 b.x * (p0.y - a.y); // < 0 clockwise, > 0 counter-clockwise

      if (IsValue1DZero(orientation))
        return IsValue1DGreater((p0.x-b.x) * (p0.x-b.x) +
                                (p0.y-b.y) * (p0.y-b.y),
                                (p0.x-a.x) * (p0.x-a.x) +
                                (p0.y-a.y) * (p0.y-a.y));
      return IsValue1DNegative(orientation);
    });

    if (includeCollinear)
    {
      int i = (int)structPoints.size() - 1;

      while (i >= 0 &&
             IsValue1DZero(p0.x * (structPoints[i].y - structPoints.back().y) +
                           structPoints[i].x * (structPoints.back().y - p0.y) +
                           structPoints.back().x * (p0.y - structPoints[i].y)))
        i--;
      reverse(structPoints.begin() + i + 1,
              structPoints.end());
    }

    vector<pt> st;
    for (int i = 0; i < (int)structPoints.size(); i++)
    {
      while (st.size() > 1)
      {
        double orientation = st[st.size()-2].x * (st.back().y - structPoints[i].y) +
                             st.back().x * (structPoints[i].y - st[st.size()-2].y) +
            structPoints[i].x*(st[st.size()-2].y - st.back().y); // < 0 clockwise, > 0 counter-clockwise

        if (IsValue1DNegative(orientation) || (includeCollinear && IsValue1DZero(orientation)))
          break;

        st.pop_back();
      }

      st.push_back(structPoints[i]);
    }

    structPoints = st;

    vector<unsigned int> result(structPoints.size());
    for (unsigned int p = 0; p < structPoints.size(); p++)
      result[result.size() - 1 - p] = structPoints[p].id;

    return result;
  }
  //  // ***************************************************************************
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
  MatrixXd GeometryUtilities::CreateEllipse(const double& axisMajorLength,
                                            const double& axisMinorLength,
                                            const unsigned int& resolution) const
  {
    const vector<double> ellipseXPoints = EquispaceCoordinates(resolution,
                                                               0.0,
                                                               axisMajorLength,
                                                               false);
    const unsigned int numInternalPoints = ellipseXPoints.size();
    vector<double> ellipseYPoints(numInternalPoints, 0.0);

    for (unsigned int p = 0; p < numInternalPoints; p++)
    {
      ellipseYPoints[p] =  axisMinorLength * sqrt(1.0 -
                                                  ellipseXPoints.at(p) *
                                                  ellipseXPoints.at(p) /
                                                  (axisMajorLength *
                                                   axisMajorLength));
    }

    Eigen::MatrixXd vertices(3, 4 + 4 * numInternalPoints);

    vertices.col(0)<< Vector3d(axisMajorLength, 0.0, 0.0);
    vertices.col(numInternalPoints + 1)<< Vector3d(0.0, axisMinorLength, 0.0);
    vertices.col(2 * (numInternalPoints + 1))<< Vector3d(-axisMajorLength, 0.0, 0.0);
    vertices.col(3 * (numInternalPoints + 1))<< Vector3d(0.0, -axisMinorLength, 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(v + 1)<< Vector3d(ellipseXPoints[numInternalPoints - 1 - v], ellipseYPoints[numInternalPoints - 1 - v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(numInternalPoints + 1 + v + 1)<< Vector3d(-ellipseXPoints[v], ellipseYPoints[v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(2 * (numInternalPoints + 1) + v + 1)<< Vector3d(-ellipseXPoints[numInternalPoints - 1 - v], -ellipseYPoints[numInternalPoints - 1 - v], 0.0);

    for (unsigned int v = 0; v < numInternalPoints; v++)
      vertices.col(3 * (numInternalPoints + 1) + v + 1)<< Vector3d(ellipseXPoints[v], -ellipseYPoints[v], 0.0);

    return vertices;
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
