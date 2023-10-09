#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"
#include "CommonUtilities.hpp"
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
  std::vector<double> GeometryUtilities::RandomCoordinates(const unsigned int size,
                                                           const bool insertExtremes,
                                                           const unsigned int seed) const
  {
    if (size == 0)
      return { };
    else if (size == 1 && insertExtremes)
      throw invalid_argument("size is not valid with false insertExtremes");

    const unsigned int max_value = 1.0 / (3.0 * _configuration.Tolerance);

    Gedim::Output::Assert(size <= max_value);

    const unsigned int num_random_points = insertExtremes ? (size - 2) :
                                                            size;
    std::vector<unsigned int> random_points = num_random_points == 0 ? std::vector<unsigned int>() :
                                                                       Gedim::Utilities::RandomArrayNoRepetition(num_random_points,
                                                                                                                 max_value,
                                                                                                                 seed);
    std::sort(random_points.begin(),
              random_points.end());

    list<double> coordinates;

    if (insertExtremes)
      coordinates.push_back(0.0);

    for (const unsigned int random_point : random_points)
    {
      if (insertExtremes &&
          (random_point == 0 ||
           random_point == max_value))
        continue;

      coordinates.push_back(static_cast<double>(random_point) /
                            static_cast<double>(max_value));
    }

    if (insertExtremes)
      coordinates.push_back(1.0);

    if (coordinates.size() < size)
    {
      const double middlePoint = 0.5 * (*coordinates.begin() +
                                        *std::next(coordinates.begin()));
      coordinates.insert(std::next(coordinates.begin()),
                         middlePoint);
    }

    if (coordinates.size() < size)
    {
      const double middlePoint = 0.5 * (*coordinates.rbegin() +
                                        *std::next(coordinates.rbegin()));
      coordinates.insert(std::prev(coordinates.end()),
                         middlePoint);
    }

    return std::vector<double>(coordinates.begin(),
                               coordinates.end());
  }
  // ***************************************************************************
  GeometryUtilities::CompareTypes GeometryUtilities::CompareValues(const double& first,
                                                                   const double& second,
                                                                   const double& tolerance) const
  {
    const double max_tolerance = std::max(abs(tolerance),
                                          std::numeric_limits<double>::epsilon());

    double relativeValue = (abs(first) <= max_tolerance ||
                            abs(second) <= max_tolerance) ? 1.0 :
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
  vector<unsigned int> GeometryUtilities::OLD_ConvexHull(const Eigen::MatrixXd& points,
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

    std::cout<< "OLD P0 "<< p0.id<< std::endl;

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

    std::cout<< "OLD B P: ";
    for (auto p : structPoints)
      std::cout<< " "<< p.id;
    std::cout<< std::endl;

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

    std::cout<< "OLD P: ";
    for (auto p : structPoints)
      std::cout<< " "<< p.id;
    std::cout<< std::endl;

    vector<pt> st;
    for (int i = 0; i < (int)structPoints.size(); i++)
    {
      std::cout<< "OLD Point "<< structPoints[i].id<< " - ";
      std::cout<< "st: ";
      for (unsigned int k = 0; k < st.size(); k++)
        std::cout<< " "<< st[k].id;
      std::cout<< std::endl;

      while (st.size() > 1)
      {
        double orientation = st[st.size()-2].x * (st.back().y - structPoints[i].y) +
                             st.back().x * (structPoints[i].y - st[st.size()-2].y) +
            structPoints[i].x*(st[st.size()-2].y - st.back().y); // < 0 clockwise, > 0 counter-clockwise

        if (IsValue1DNegative(orientation) || (includeCollinear && IsValue1DZero(orientation)))
          break;

        std::cout.precision(4);
        std::cout<< std::scientific<< "\t OLD rm: "<< st[st.size() - 1].id<< " value "<< orientation<< std::endl;

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
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::ConvexHull(const Eigen::MatrixXd& points,
                                                     const bool& includeCollinear) const
  {
    Output::Assert(points.rows() == 3 && PointsAre2D(points));

    list<unsigned int> convexHull;
    const unsigned int numPoints = points.cols();

    if (numPoints == 1)
      return vector<unsigned int> { 0 };

    struct pt final
    {
        double x;
        double y;
        unsigned int id;
        inline Eigen::Vector3d Point() const { return Eigen::Vector3d(x, y, 0.0); }
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

      if (Are1DValuesEqual(p0.y, p.y) &&
          IsValue1DGreater(p0.x, p.x))
        p0 = p;
    }

    std::cout<< "P0 "<< p0.id<< std::endl;

    sort(structPoints.begin(),
         structPoints.end(),
         [&p0, this](const pt& a, const pt& b)
    {
      const double norm_b_a = (a.x - b.x) * (a.x - b.x) +
                              (a.y - b.y) * (a.y - b.y);
      const double norm_b_p0 = (p0.x - b.x) * (p0.x - b.x) +
                               (p0.y - b.y) * (p0.y - b.y);

      const double orientation = p0.x * (a.y - b.y) +
                                 a.x * (b.y - p0.y) +
                                 b.x * (p0.y - a.y); // < 0 clockwise, > 0 counter-clockwise

      const double orientation_polar = PolarAngle(a.Point(),
                                                  b.Point(),
                                                  p0.Point(),
                                                  norm_b_a,
                                                  norm_b_p0);

      if(Compare1DValues(0.0,
                         orientation) !=
         Compare1DValues(0.0,
                         orientation_polar))
      {
        std::cout<< "HERE 1 "<< std::endl;
        std::cout.precision(16);
        std::cout<< std::scientific<< "\t 0"<< " value "<< orientation<< " polar "<< orientation_polar<< std::endl;
        std::cout<< std::scientific<< "\t 0"<< " norm_b_a "<< norm_b_a<< " norm_b_p0 "<< norm_b_p0<< std::endl;
        std::cout<< scientific<< a.Point().transpose()<< "\n"<< b.Point().transpose()<< "\n"<< p0.Point().transpose()<< std::endl;
      }

      const double norm_a_p0 = (p0.x - a.x) * (p0.x - a.x) +
                               (p0.y - a.y) * (p0.y - a.y);

      if (IsValue1DZero(orientation))
        return IsValue1DGreater(norm_b_p0,
                                norm_a_p0);
      return IsValue1DNegative(orientation);
    });

    std::cout<< "B P: ";
    for (auto p : structPoints)
      std::cout<< " "<< p.id;
    std::cout<< std::endl;

    if (includeCollinear)
    {
      int i = (int)structPoints.size() - 1;

      double norm_b_a = (structPoints[i].x - structPoints.back().x) * (structPoints[i].x - structPoints.back().x) +
                        (structPoints[i].y - structPoints.back().y) * (structPoints[i].y - structPoints.back().y);
      double norm_b_p0 = (p0.x - structPoints.back().x) * (p0.x - structPoints.back().x) +
                         (p0.y - structPoints.back().y) * (p0.y - structPoints.back().y);
      double orientation = p0.x * (structPoints[i].y - structPoints.back().y) +
                           structPoints[i].x * (structPoints.back().y - p0.y) +
                           structPoints.back().x * (p0.y - structPoints[i].y);
      double orientation_polar = PolarAngle(structPoints[i].Point(),
                                            structPoints.back().Point(),
                                            p0.Point(),
                                            norm_b_a,
                                            norm_b_p0);

      if(Compare1DValues(0.0,
                         orientation) !=
         Compare1DValues(0.0,
                         orientation_polar))
      {
        std::cout<< "HERE 2 "<< std::endl;
        std::cout.precision(16);
        std::cout<< std::scientific<< "\t 0"<< " value "<< orientation<< " polar "<< orientation_polar<< std::endl;
        std::cout<< std::scientific<< "\t 0"<< " norm_b_a "<< norm_b_a<< " norm_b_p0 "<< norm_b_p0<< std::endl;
      }

      while (i >= 0 &&
             IsValue1DZero(p0.x * (structPoints[i].y - structPoints.back().y) +
                           structPoints[i].x * (structPoints.back().y - p0.y) +
                           structPoints.back().x * (p0.y - structPoints[i].y)))
      {
        i--;

        norm_b_a = (structPoints[i].x - structPoints.back().x) * (structPoints[i].x - structPoints.back().x) +
                   (structPoints[i].y - structPoints.back().y) * (structPoints[i].y - structPoints.back().y);
        norm_b_p0 = (p0.x - structPoints.back().x) * (p0.x - structPoints.back().x) +
                    (p0.y - structPoints.back().y) * (p0.y - structPoints.back().y);
        orientation = p0.x * (structPoints[i].y - structPoints.back().y) +
                      structPoints[i].x * (structPoints.back().y - p0.y) +
                      structPoints.back().x * (p0.y - structPoints[i].y);
        orientation_polar = PolarAngle(structPoints[i].Point(),
                                       structPoints.back().Point(),
                                       p0.Point(),
                                       norm_b_a,
                                       norm_b_p0);

        if(Compare1DValues(0.0,
                           orientation) !=
           Compare1DValues(0.0,
                           orientation_polar))
        {
          std::cout<< "HERE 3 "<< std::endl;
          std::cout.precision(16);
          std::cout<< std::scientific<< "\t 0"<< " value "<< orientation<< " polar "<< orientation_polar<< std::endl;
          std::cout<< std::scientific<< "\t 0"<< " norm_b_a "<< norm_b_a<< " norm_b_p0 "<< norm_b_p0<< std::endl;
        }
      }

      reverse(structPoints.begin() + i + 1,
              structPoints.end());
    }

    std::cout<< "P: ";
    for (auto p : structPoints)
      std::cout<< " "<< p.id;
    std::cout<< std::endl;

    vector<pt> st;
    for (unsigned int i = 0; i < structPoints.size(); i++)
    {
      std::cout<< "Point "<< structPoints[i].id<< " - ";
      std::cout<< "st: ";
      for (unsigned int k = 0; k < st.size(); k++)
        std::cout<< " "<< st[k].id;
      std::cout<< std::endl;

      while (st.size() > 1)
      {
        const double norm_b_a = (st.back().x - structPoints[i].x) * (st.back().x - structPoints[i].x) +
                                (st.back().y - structPoints[i].y) * (st.back().y - structPoints[i].y);
        const double norm_b_c = (st[st.size() - 2].x - structPoints[i].x) * (st[st.size() - 2].x - structPoints[i].x) +
            (st[st.size() - 2].y - structPoints[i].y) * (st[st.size() - 2].y - structPoints[i].y);
        double orientation = st[st.size()-2].x * (st.back().y - structPoints[i].y) +
                             st.back().x * (structPoints[i].y - st[st.size()-2].y) +
            structPoints[i].x*(st[st.size()-2].y - st.back().y); // < 0 clockwise, > 0 counter-clockwise


        const double orientation_polar = PolarAngle(st.back().Point(),
                                                    structPoints[i].Point(),
                                                    st[st.size() - 2].Point(),
            norm_b_a,
            norm_b_c);

        if(Compare1DValues(0.0,
                           orientation) !=
           Compare1DValues(0.0,
                           orientation_polar))
        {
          std::cout<< "HERE 4 "<< std::endl;
          std::cout.precision(16);
          std::cout<< std::scientific<< "\t 0"<< " value "<< orientation<< " polar "<< orientation_polar<< std::endl;
          std::cout<< std::scientific<< "\t 0"<< " norm_b_a "<< norm_b_a<< " norm_b_p0 "<< norm_b_c<< std::endl;
        }

        if (IsValue1DNegative(orientation) ||
            (includeCollinear && IsValue1DZero(orientation)))
          break;

        std::cout.precision(4);
        std::cout<< std::scientific<< "\t rm: "<< st[st.size() - 1].id<< " value "<< orientation<< " polar "<< orientation_polar<< std::endl;

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
