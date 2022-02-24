#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  Eigen::Vector3d GeometryUtilities::PolygonNormal(const MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    Vector3d normal;

    normal.setZero();
    const unsigned int& numVertices = polygonVertices.cols();

    for (unsigned int i = 0; i < numVertices; i++)
    {
      Vector3d edge = polygonVertices.col((i + 1) % numVertices) - polygonVertices.col(i);
      Vector3d edgePrevious = polygonVertices.col((i - 1) % numVertices) - polygonVertices.col(i);
      normal.noalias() += edge.cross(edgePrevious);
    }

    return normal.normalized();
  }
  // ***************************************************************************
  VectorXd GeometryUtilities::PolygonEdgeLengths(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    VectorXd edgeLengths(numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeLengths[v] = SegmentLength(polygonVertices.col(v),
                                     polygonVertices.col((v + 1) % numVertices));
    }

    return edgeLengths;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolygonEdgeTangents(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    MatrixXd edgeTangents(3, numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeTangents.col(v) = SegmentTangent(polygonVertices.col(v),
                                           polygonVertices.col((v + 1) % numVertices));
    }

    return edgeTangents;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolygonEdgeNormals(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    MatrixXd edgeNormals(3, numVertices);
    for (unsigned int v = 0; v < numVertices; v++)
    {
      edgeNormals.col(v) = SegmentNormal(polygonVertices.col(v),
                                         polygonVertices.col((v + 1) % numVertices));
    }

    return edgeNormals;
  }
  // ***************************************************************************
  Vector3d GeometryUtilities::PolygonCentroid(const Eigen::MatrixXd& polygonVertices,
                                              const double& polygonArea) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    Eigen::Vector3d centroid;
    centroid.setZero();

    const unsigned int& numVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d& vertex = polygonVertices.col(v);
      const Eigen::Vector3d& nextVertex = polygonVertices.col((v + 1) % numVertices);
      centroid.x() += (vertex.x() + nextVertex.x()) *
                      (vertex.x() * nextVertex.y() - nextVertex.x() * vertex.y());

      centroid.y() += (vertex.y() + nextVertex.y()) *
                      (vertex.x() * nextVertex.y() - nextVertex.x() * vertex.y());
    }

    centroid /= 6.0 * polygonArea;

    return centroid;
  }
  // ***************************************************************************
  double GeometryUtilities::PolygonDiameter(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();

    double diameter = 0.0;

    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d& vertexOne = polygonVertices.col(v);
      for (unsigned int w = v + 1; w < numVertices; w++)
      {
        const Eigen::Vector3d& vertexTwo = polygonVertices.col(w);
        double distance = PointDistance(vertexOne, vertexTwo);
        if (Compare1DValues(diameter, distance) == CompareTypes::FirstBeforeSecond)
          diameter = distance;
      }
    }

    return diameter;
  }
  // ***************************************************************************
  Matrix3d GeometryUtilities::PolygonRotationMatrix(const Eigen::MatrixXd& polygonVertices,
                                                    const Eigen::Vector3d& polygonNormal,
                                                    const Eigen::Vector3d& polygonTranslation) const
  {
    Output::Assert(Compare1DValues(polygonNormal.norm(), 1.0) == CompareTypes::Coincident);

    const unsigned int& numVertices = polygonVertices.cols();
    MatrixXd Z, W;
    Z.setZero(3, numVertices);
    W.setZero(3, numVertices);

    const Vector3d V1mV0 = polygonVertices.col(1) - polygonTranslation;
    double normVectorOne = V1mV0.norm();
    Z.col(0)<< V1mV0;
    W.col(0)<< normVectorOne, 0.0, 0.0;
    for (unsigned int i = 2; i < numVertices; i++)
    {
      const Vector3d VimV0 = polygonVertices.col(i) - polygonTranslation;
      Z.col(i - 1)<< VimV0;

      const double normVectorI = VimV0.norm();
      const double cosTheta = VimV0.dot(V1mV0) / (normVectorOne * normVectorI);

      if (Compare1DValues(cosTheta, 1.0) == CompareTypes::Coincident)
        W.col(i - 1)<< normVectorI, 0.0, 0.0;
      else if (Compare1DValues(cosTheta, -1.0) == CompareTypes::Coincident)
        W.col(i - 1)<< -normVectorI, 0.0, 0.0;
      else if (Compare1DValues(cosTheta, 0.0) == CompareTypes::Coincident)
        W.col(i - 1)<< 0.0, normVectorI, 0.0;
      else
        W.col(i - 1) << normVectorI * cosTheta, normVectorI * sqrt(1.0 - cosTheta*cosTheta), 0;
    }
    Z.col(numVertices - 1) = polygonNormal;
    W.col(numVertices - 1)<< 0.0, 0.0, 1.0;
    MatrixXd H = W * Z.transpose();
    JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);

    return svd.matrixV() * (svd.matrixU()).transpose();
  }
  // ***************************************************************************
  bool GeometryUtilities::PolygonIsConvex(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
    for (unsigned int v = 0; v < numVertices; v++)
    {
      const Eigen::Vector3d edgeOrigin = polygonVertices.col(v);
      const Eigen::Vector3d edgeEnd = polygonVertices.col((v + 1) % numVertices);
      const Eigen::Vector3d vertexNextEdge = polygonVertices.col((v + 2) % numVertices);

      if (PointSegmentPosition(vertexNextEdge,
                               edgeOrigin,
                               edgeEnd) == PointSegmentPositionTypes::RightTheSegment)
        return false;
    }

    return true;
  }
  // ***************************************************************************
  GeometryUtilities::PolygonTypes GeometryUtilities::PolygonType(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
    if (numVertices == 3)
      return PolygonTypes::Triangle;
    else if (numVertices == 4)
      return PolygonTypes::Quadrilateral;
    else
      return PolygonTypes::Generic;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolygonTriangulationByFirstVertex(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    list<unsigned int> triangleList;

    const unsigned int numPolygonVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      const unsigned int nextVertex = (v + 1) % numPolygonVertices;
      const unsigned int nextNextVertex = (v + 2) % numPolygonVertices;

      if (nextNextVertex == 0)
        break;

      triangleList.push_back(0);
      triangleList.push_back(nextVertex);
      triangleList.push_back(nextNextVertex);
    }

    Output::Assert(triangleList.size() % 3 == 0);

    return vector<unsigned int>(triangleList.begin(), triangleList.end());
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolygonTriangulationByInternalPoint(const Eigen::MatrixXd& polygonVertices,
                                                                              const Eigen::Vector3d& point) const
  {
    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);

    const unsigned int numPolygonVertices = polygonVertices.cols();
    vector<unsigned int> triangles(3 * numPolygonVertices);

    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      triangles[3 * v] = numPolygonVertices;
      triangles[3 * v + 1] = v;
      triangles[3 * v + 2] = (v + 1) % numPolygonVertices;
    }

    Output::Assert(triangles.size() % 3 == 0);

    return triangles;
  }
  // ***************************************************************************
  GeometryUtilities::PolygonDivisionByAngleQuadrantResult GeometryUtilities::PolygonOutsideCircleDivisionByAngleQuadrant(const Eigen::MatrixXd& polygonVertices,
                                                                                                                         const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                                                         const Eigen::Vector3d& circleCenter,
                                                                                                                         const double& circleRadius,
                                                                                                                         const unsigned int& curvedEdgeIndex) const
  {
    PolygonDivisionByAngleQuadrantResult result;

    list<Vector3d> newCoordinates;

    // insert polygon and circle center in output
    const unsigned int& numPolygonVertices = polygonVertices.cols();

    for (unsigned int p = 0; p < numPolygonVertices; p++)
      newCoordinates.push_back(polygonVertices.col(p));

    const unsigned int curvedEdgeOriginIndex = curvedEdgeIndex;
    const unsigned int curvedEdgeEndIndex = (curvedEdgeIndex + 1) % numPolygonVertices;
    const Vector3d& curvedEdgeOrigin = polygonVertices.col(curvedEdgeOriginIndex);
    const Vector3d& curvedEdgeEnd = polygonVertices.col(curvedEdgeEndIndex);

    // check if the polygon is contained in the circular arc
    int curvedEdgeOriginEdgeIntersectionIndex = -1;
    int curvedEdgeEndEdgeIntersectionIndex = -1;
    Eigen::Vector3d curvedEdgeOriginEdgeIntersection;
    Eigen::Vector3d curvedEdgeEndEdgeIntersection;

    // Intersect with curved edge origin
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      if (e == curvedEdgeIndex)
        continue;

      const unsigned int edgeOriginIndex = e;
      const unsigned int edgeEndIndex = (e + 1) % numPolygonVertices;

      if (edgeEndIndex == curvedEdgeIndex)
        continue;

      const Vector3d& edgeOrigin = polygonVertices.col(edgeOriginIndex);
      const Vector3d& edgeEnd = polygonVertices.col(edgeEndIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(e);

      IntersectionSegmentSegmentResult resultOrigin = IntersectionSegmentSegment(circleCenter,
                                                                                 curvedEdgeOrigin,
                                                                                 edgeOrigin,
                                                                                 edgeEnd);
      if (resultOrigin.IntersectionLinesType !=
          IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting)
        continue;

      if (resultOrigin.IntersectionSegmentsType !=
          IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection)
        continue;

      if (resultOrigin.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::OnSegmentOrigin &&
          resultOrigin.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::InsideSegment)
        continue;

      curvedEdgeOriginEdgeIntersectionIndex = e;
      curvedEdgeOriginEdgeIntersection = edgeOrigin +
                                         resultOrigin.SecondSegmentIntersections[0].CurvilinearCoordinate *
                                         edgeTangent;
    }

    // Intersect with curved edge end
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      if (e == curvedEdgeIndex)
        continue;

      const unsigned int edgeOriginIndex = e;
      const unsigned int edgeEndIndex = (e + 1) % numPolygonVertices;

      if (edgeOriginIndex == curvedEdgeIndex)
        continue;

      const Vector3d& edgeOrigin = polygonVertices.col(edgeOriginIndex);
      const Vector3d& edgeEnd = polygonVertices.col(edgeEndIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(e);

      IntersectionSegmentSegmentResult resultEnd = IntersectionSegmentSegment(circleCenter,
                                                                              curvedEdgeEnd,
                                                                              edgeOrigin,
                                                                              edgeEnd);
      if (resultEnd.IntersectionLinesType !=
          IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting)
        continue;

      if (resultEnd.IntersectionSegmentsType !=
          IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection)
        continue;

      if (resultEnd.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::OnSegmentEnd &&
          resultEnd.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::InsideSegment)
        continue;

      curvedEdgeEndEdgeIntersectionIndex = e;
      curvedEdgeEndEdgeIntersection = edgeOrigin +
                                      resultEnd.SecondSegmentIntersections[0].CurvilinearCoordinate *
                                      edgeTangent;
    }

    unsigned int numSubPolygons = 1;
    int curvedEdgeOriginVertexIntersectionIndex = -1;
    int curvedEdgeEndVertexIntersectionIndex = -1;

    if (curvedEdgeOriginEdgeIntersectionIndex != -1)
    {
      curvedEdgeOriginVertexIntersectionIndex = newCoordinates.size();
      newCoordinates.push_back(curvedEdgeOriginEdgeIntersection);
      numSubPolygons++;
    }

    if (curvedEdgeEndEdgeIntersectionIndex != -1)
    {
      curvedEdgeEndVertexIntersectionIndex = newCoordinates.size();
      newCoordinates.push_back(curvedEdgeEndEdgeIntersection);
      numSubPolygons++;
    }

    vector<list<unsigned int>> subPolygons(numSubPolygons);
    result.SubPolygonTypes.resize(numSubPolygons);
    unsigned int subPolygonCounter = 0;

    // create origin subpolygon
    if (curvedEdgeOriginVertexIntersectionIndex != -1)
    {
      result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::ExternalOrigin;

      list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
      subPolygon.push_back(curvedEdgeOriginIndex);
      subPolygon.push_back(curvedEdgeOriginVertexIntersectionIndex);
      unsigned int vertexIndex = (curvedEdgeOriginEdgeIntersectionIndex + 1) % numPolygonVertices;
			while (vertexIndex != curvedEdgeOriginIndex)
      {
        subPolygon.push_back(vertexIndex);
        vertexIndex = (vertexIndex + 1) % numPolygonVertices;
      }
      subPolygonCounter++;
    }

    // create end subpolygon
    if (curvedEdgeEndVertexIntersectionIndex != -1)
    {
      result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::ExternalEnd;

      list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
      subPolygon.push_back(curvedEdgeEndIndex);
      unsigned int vertexIndex = (curvedEdgeEndIndex + 1) % numPolygonVertices;
			while (vertexIndex != curvedEdgeEndEdgeIntersectionIndex)
      {
        subPolygon.push_back(vertexIndex);
        vertexIndex = (vertexIndex + 1) % numPolygonVertices;
      }
			subPolygon.push_back(curvedEdgeEndEdgeIntersectionIndex);
      subPolygon.push_back(curvedEdgeEndVertexIntersectionIndex);

      subPolygonCounter++;
    }

    // create central subPolygon
    result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::Internal;

    list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
    subPolygon.push_back(curvedEdgeOriginIndex);
    subPolygon.push_back(curvedEdgeEndIndex);
    if (curvedEdgeEndVertexIntersectionIndex != -1)
      subPolygon.push_back(curvedEdgeEndVertexIntersectionIndex);
    const unsigned int middleVerticesStart = (curvedEdgeEndVertexIntersectionIndex != -1) ?  (curvedEdgeEndEdgeIntersectionIndex + 1) % numPolygonVertices :
                                                                                             (curvedEdgeEndIndex + 1) % numPolygonVertices;
    const unsigned int middleVerticesEnd = (curvedEdgeOriginVertexIntersectionIndex != -1) ?  (curvedEdgeOriginEdgeIntersectionIndex + 1) % numPolygonVertices :
                                                                                              curvedEdgeOriginIndex;

    unsigned int middleVertexIndex = middleVerticesStart;
    while (middleVertexIndex != middleVerticesEnd)
    {
      subPolygon.push_back(middleVertexIndex);
      middleVertexIndex = (middleVertexIndex + 1) % numPolygonVertices;
    }

    if (curvedEdgeOriginVertexIntersectionIndex != -1)
      subPolygon.push_back(curvedEdgeOriginVertexIntersectionIndex);

    // convert output
    result.Points.setZero(3, newCoordinates.size());
    unsigned int counter = 0;
    for (const Vector3d& point : newCoordinates)
      result.Points.col(counter++)<< point;

    result.SubPolygons.resize(subPolygons.size());
    for (unsigned int s = 0; s < subPolygons.size(); s++)
      result.SubPolygons[s] = vector<unsigned int>(subPolygons[s].begin(),
                                                   subPolygons[s].end());

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::PolygonDivisionByAngleQuadrantResult GeometryUtilities::PolygonInsideCircleDivisionByAngleQuadrant(const Eigen::MatrixXd& polygonVertices,
                                                                                                                        const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                                                        const Eigen::Vector3d& circleCenter,
                                                                                                                        const double& circleRadius,
                                                                                                                        const unsigned int& curvedEdgeIndex) const
  {
    PolygonDivisionByAngleQuadrantResult result;

    list<Vector3d> newCoordinates;

    // insert polygon and circle center in output
    const unsigned int& numPolygonVertices = polygonVertices.cols();

    for (unsigned int p = 0; p < numPolygonVertices; p++)
      newCoordinates.push_back(polygonVertices.col(p));

    const unsigned int curvedEdgeOriginIndex = curvedEdgeIndex;
    const unsigned int curvedEdgeEndIndex = (curvedEdgeIndex + 1) % numPolygonVertices;
    const Vector3d& curvedEdgeOrigin = polygonVertices.col(curvedEdgeOriginIndex);
    const Vector3d& curvedEdgeEnd = polygonVertices.col(curvedEdgeEndIndex);

    // check if the polygon is contained in the circular arc
    int curvedEdgeOriginEdgeIntersectionIndex = -1;
    int curvedEdgeEndEdgeIntersectionIndex = -1;
    Eigen::Vector3d curvedEdgeOriginEdgeIntersection;
    Eigen::Vector3d curvedEdgeEndEdgeIntersection;

    // Intersect with curved edge origin
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      if (e == curvedEdgeIndex)
        continue;

      const unsigned int edgeOriginIndex = e;
      const unsigned int edgeEndIndex = (e + 1) % numPolygonVertices;

			if (edgeEndIndex == curvedEdgeIndex)
        continue;

      const Vector3d& edgeOrigin = polygonVertices.col(edgeOriginIndex);
      const Vector3d& edgeEnd = polygonVertices.col(edgeEndIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(e);

      IntersectionSegmentSegmentResult resultOrigin = IntersectionSegmentSegment(circleCenter,
                                                                                 curvedEdgeOrigin,
                                                                                 edgeOrigin,
                                                                                 edgeEnd);
      if (resultOrigin.IntersectionLinesType !=
          IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting)
        continue;

      if (resultOrigin.IntersectionSegmentsType !=
          IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
        continue;

      if (resultOrigin.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::OnSegmentOrigin &&
          resultOrigin.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::InsideSegment)
        continue;

      curvedEdgeOriginEdgeIntersectionIndex = e;
      curvedEdgeOriginEdgeIntersection = edgeOrigin +
                                         resultOrigin.SecondSegmentIntersections[0].CurvilinearCoordinate *
                                         edgeTangent;
    }

    // Intersect with curved edge end
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      if (e == curvedEdgeIndex)
        continue;

      const unsigned int edgeOriginIndex = e;
      const unsigned int edgeEndIndex = (e + 1) % numPolygonVertices;

      if (edgeOriginIndex == curvedEdgeIndex)
        continue;

      const Vector3d& edgeOrigin = polygonVertices.col(edgeOriginIndex);
      const Vector3d& edgeEnd = polygonVertices.col(edgeEndIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(e);

      IntersectionSegmentSegmentResult resultEnd = IntersectionSegmentSegment(circleCenter,
                                                                              curvedEdgeEnd,
                                                                              edgeOrigin,
                                                                              edgeEnd);
      if (resultEnd.IntersectionLinesType !=
          IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting)
        continue;

      if (resultEnd.IntersectionSegmentsType !=
          IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
        continue;

      if (resultEnd.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::OnSegmentEnd &&
          resultEnd.SecondSegmentIntersections[0].Type !=
          PointSegmentPositionTypes::InsideSegment)
        continue;

      curvedEdgeEndEdgeIntersectionIndex = e;
      curvedEdgeEndEdgeIntersection = edgeOrigin +
                                      resultEnd.SecondSegmentIntersections[0].CurvilinearCoordinate *
                                      edgeTangent;
    }

    unsigned int numSubPolygons = 1;
    int curvedEdgeOriginVertexIntersectionIndex = -1;
    int curvedEdgeEndVertexIntersectionIndex = -1;

    if (curvedEdgeOriginEdgeIntersectionIndex != -1)
    {
      curvedEdgeOriginVertexIntersectionIndex = newCoordinates.size();
      newCoordinates.push_back(curvedEdgeOriginEdgeIntersection);
      numSubPolygons++;
    }

    if (curvedEdgeEndEdgeIntersectionIndex != -1)
    {
      curvedEdgeEndVertexIntersectionIndex = newCoordinates.size();
      newCoordinates.push_back(curvedEdgeEndEdgeIntersection);
      numSubPolygons++;
    }

    vector<list<unsigned int>> subPolygons(numSubPolygons);
    result.SubPolygonTypes.resize(numSubPolygons);
    unsigned int subPolygonCounter = 0;

    // create origin subpolygon
    if (curvedEdgeOriginVertexIntersectionIndex != -1)
    {
      result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::ExternalOrigin;

      list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
      subPolygon.push_back(curvedEdgeOriginIndex);
      subPolygon.push_back(curvedEdgeOriginVertexIntersectionIndex);
      unsigned int vertexIndex = (curvedEdgeOriginEdgeIntersectionIndex + 1) % numPolygonVertices;
			while (vertexIndex != curvedEdgeOriginIndex)
      {
        subPolygon.push_back(vertexIndex);
        vertexIndex = (vertexIndex + 1) % numPolygonVertices;
      }
      subPolygonCounter++;
    }

    // create end subpolygon
    if (curvedEdgeEndVertexIntersectionIndex != -1)
    {
      result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::ExternalEnd;

      list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
      subPolygon.push_back(curvedEdgeEndIndex);
      unsigned int vertexIndex = (curvedEdgeEndIndex + 1) % numPolygonVertices;
			while (vertexIndex != curvedEdgeEndEdgeIntersectionIndex)
      {
        subPolygon.push_back(vertexIndex);
        vertexIndex = (vertexIndex + 1) % numPolygonVertices;
      }
			subPolygon.push_back(curvedEdgeEndEdgeIntersectionIndex);
      subPolygon.push_back(curvedEdgeEndVertexIntersectionIndex);

      subPolygonCounter++;
    }

    // create central subPolygon
    result.SubPolygonTypes[subPolygonCounter] = PolygonDivisionByAngleQuadrantResult::Types::Internal;

    list<unsigned int>& subPolygon = subPolygons[subPolygonCounter];
    subPolygon.push_back(curvedEdgeOriginIndex);
    subPolygon.push_back(curvedEdgeEndIndex);
    if (curvedEdgeEndVertexIntersectionIndex != -1)
      subPolygon.push_back(curvedEdgeEndVertexIntersectionIndex);
    const unsigned int middleVerticesStart = (curvedEdgeEndVertexIntersectionIndex != -1) ?  (curvedEdgeEndEdgeIntersectionIndex + 1) % numPolygonVertices :
                                                                                             (curvedEdgeEndIndex + 1) % numPolygonVertices;
    const unsigned int middleVerticesEnd = (curvedEdgeOriginVertexIntersectionIndex != -1) ?  (curvedEdgeOriginEdgeIntersectionIndex + 1) % numPolygonVertices :
                                                                                              curvedEdgeOriginIndex;

    unsigned int middleVertexIndex = middleVerticesStart;
    while (middleVertexIndex != middleVerticesEnd)
    {
      subPolygon.push_back(middleVertexIndex);
      middleVertexIndex = (middleVertexIndex + 1) % numPolygonVertices;
    }

    if (curvedEdgeOriginVertexIntersectionIndex != -1)
      subPolygon.push_back(curvedEdgeOriginVertexIntersectionIndex);

    // convert output
    result.Points.setZero(3, newCoordinates.size());
    unsigned int counter = 0;
    for (const Vector3d& point : newCoordinates)
      result.Points.col(counter++)<< point;

    result.SubPolygons.resize(subPolygons.size());
    for (unsigned int s = 0; s < subPolygons.size(); s++)
      result.SubPolygons[s] = vector<unsigned int>(subPolygons[s].begin(),
                                                   subPolygons[s].end());

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::PolygonDivisionByCircleResult GeometryUtilities::PolygonDivisionByCircle(const Eigen::MatrixXd& polygonVertices,
                                                                                              const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                              const Eigen::Vector3d& circleCenter,
                                                                                              const double& circleRadius,
                                                                                              const unsigned int& curvedEdgeIndex) const
  {
    PolygonDivisionByCircleResult result;

    list<Vector3d> newCoordinates;

    // insert polygon and circle center in output
    const unsigned int& numPolygonVertices = polygonVertices.cols();
    const unsigned int cirlceCenterIndex = numPolygonVertices;

    for (unsigned int p = 0; p < numPolygonVertices; p++)
      newCoordinates.push_back(polygonVertices.col(p));
    newCoordinates.push_back(circleCenter);

    const unsigned int curvedEdgeOriginIndex = curvedEdgeIndex;
    const unsigned int curvedEdgeEndIndex = (curvedEdgeIndex + 1) % numPolygonVertices;
    const Vector3d& curvedEdgeOrigin = polygonVertices.col(curvedEdgeOriginIndex);
    const Vector3d& curvedEdgeEnd = polygonVertices.col(curvedEdgeEndIndex);

    // intersects all the lines from circleCenter to each vertex not belonging to curved edge
    unsigned int edgeNumber = curvedEdgeEndIndex;
    unsigned int lastCheck = (curvedEdgeIndex == 0) ?  numPolygonVertices - 1 : curvedEdgeIndex - 1;
    vector<int> newVerticesIndicesPerEdge(numPolygonVertices, -1);
    newVerticesIndicesPerEdge[curvedEdgeIndex] = curvedEdgeIndex;

    while (edgeNumber != lastCheck)
    {
      const unsigned int originEdgeIndex = edgeNumber;
      const unsigned int endEdgeIndex = (edgeNumber + 1) % numPolygonVertices;

      const Eigen::Vector3d& edgeOrigin = polygonVertices.col(originEdgeIndex);
      const Eigen::Vector3d& edgeEnd = polygonVertices.col(endEdgeIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(edgeNumber);

      // check if the edge is aligned to the angle quadrant
      if (PointIsAligned(edgeOrigin,
                         edgeEnd,
                         circleCenter))
      {
        edgeNumber++;
        if (edgeNumber == numPolygonVertices)
          edgeNumber = 0;

        continue;
      }

      // check if the next edge is aligned to the angle quadrant
      const unsigned int nextEdgeEndIndex = (endEdgeIndex + 1) % numPolygonVertices;

      if (PointIsAligned(edgeEnd,
                         polygonVertices.col(nextEdgeEndIndex),
                         circleCenter))
      {
        edgeNumber++;
        if (edgeNumber == numPolygonVertices)
          edgeNumber = 0;

        continue;
      }

      // intersecting the line with the circle
      const Eigen::Vector3d endEdgeCenterTangent = SegmentTangent(edgeEnd,
                                                                  circleCenter);

      const IntersectionSegmentCircleResult intersection = IntersectionSegmentCircle(edgeEnd,
                                                                                     circleCenter,
                                                                                     circleCenter,
                                                                                     circleRadius);

      Output::Assert(intersection.Type == IntersectionSegmentCircleResult::Types::TwoIntersections);

      double intersectionCoordinate = -1.0;
      for (unsigned int i = 0; i < 2; i++)
      {
        if (intersection.SegmentIntersections[i].Type != PointSegmentPositionTypes::InsideSegment)
          continue;

        intersectionCoordinate = intersection.SegmentIntersections[i].CurvilinearCoordinate;
        break;
      }

      Output::Assert(IsValue1DPositive(intersectionCoordinate));

      const Eigen::Vector3d intersectionPoint = edgeEnd + intersectionCoordinate * endEdgeCenterTangent;
      newVerticesIndicesPerEdge[edgeNumber] = newCoordinates.size();
      newCoordinates.push_back(intersectionPoint);

      edgeNumber++;
      if (edgeNumber == numPolygonVertices)
        edgeNumber = 0;
    }

    // convert newCoordinates
    result.Points.setZero(3, newCoordinates.size());
    unsigned int counter = 0;
    for (const Vector3d& point : newCoordinates)
      result.Points.col(counter++)<< point;

    // create sub-polygons
    list<list<unsigned int>> subPolygonsIndices;
    subPolygonsIndices.push_back(list<unsigned int>());

    edgeNumber = curvedEdgeEndIndex;
    while (edgeNumber != curvedEdgeIndex)
    {
      list<unsigned int>& subPolygonIndices = subPolygonsIndices.back();

      const unsigned int originEdgeIndex = edgeNumber;
      const unsigned int endEdgeIndex = (edgeNumber + 1) % numPolygonVertices;

      subPolygonIndices.push_back(originEdgeIndex);

      if (endEdgeIndex == curvedEdgeIndex)
      {
        subPolygonIndices.push_back(endEdgeIndex);
        break;
      }

      if (newVerticesIndicesPerEdge[edgeNumber] != -1)
      {
        subPolygonIndices.push_back(endEdgeIndex);
        subPolygonIndices.push_back(newVerticesIndicesPerEdge[edgeNumber]);
        subPolygonsIndices.push_back(list<unsigned int> { (unsigned int)newVerticesIndicesPerEdge[edgeNumber] });
      }

      edgeNumber++;
      if (edgeNumber == numPolygonVertices)
        edgeNumber = 0;
    }

    // convert sub-polygons
    result.SubPolygons.resize(subPolygonsIndices.size());
    unsigned int subPolygonsCounter = 0;
    for (const list<unsigned int>& subPolygonIndices : subPolygonsIndices)
      result.SubPolygons[subPolygonsCounter++] = vector<unsigned int>(subPolygonIndices.begin(), subPolygonIndices.end());

    // create sub triangles
    result.SubTriangles.resize(result.SubPolygons.size(),
                               vector<unsigned int> { cirlceCenterIndex,
                                                      (unsigned int)newCoordinates.size(),
                                                      (unsigned int)newCoordinates.size() });

    // if more then one sub-polygon
    if (result.SubPolygons.size() > 1)
    {
      // first sub triangle
      const vector<unsigned int>& subPolygon = result.SubPolygons[0];
      vector<unsigned int>& subTriangle = result.SubTriangles[0];

      if (subPolygon.size() == 3)
      {
        subTriangle[1] = subPolygon[0];
        subTriangle[2] = subPolygon[1];
      }
      else
      {
        subTriangle[1] = subPolygon[1];
        subTriangle[2] = subPolygon[subPolygon.size() - 2];
      }

      // internal sub triangles
      for (unsigned int sp = 1; sp < result.SubPolygons.size() - 1; sp++)
      {
        const vector<unsigned int>& subPolygon = result.SubPolygons[sp];
        vector<unsigned int>& subTriangle = result.SubTriangles[sp];

        if (subPolygon.size() == 3 &&
            sp == 0 &&
            result.SubPolygons.size() > 1) // first subPolygon
        {
          subTriangle[1] = subPolygon[0];
          subTriangle[2] = subPolygon[1];
        }
        else if (subPolygon.size() == 3) // last subPolygon
        {
          subTriangle[1] = subPolygon[1];
          subTriangle[2] = subPolygon[2];
        }
        else
        {
          subTriangle[1] = subPolygon[1];
          subTriangle[2] = subPolygon[subPolygon.size() - 2];
        }
      }

      // last sub triangle
      const vector<unsigned int>& lastSubPolygon = result.SubPolygons[result.SubPolygons.size() - 1];
      vector<unsigned int>& lastSubTriangle = result.SubTriangles[result.SubPolygons.size() - 1];
      if (lastSubPolygon.size() == 3)
      {
        lastSubTriangle[1] = lastSubPolygon[1];
        lastSubTriangle[2] = lastSubPolygon[2];
      }
      else
      {
        lastSubTriangle[1] = lastSubPolygon[1];
        lastSubTriangle[2] = lastSubPolygon[lastSubPolygon.size() - 2];
      }
    }
    else
    {
      // just one sub-triangle
      const vector<unsigned int>& subPolygon = result.SubPolygons[0];
      vector<unsigned int>& subTriangle = result.SubTriangles[0];

      if (subPolygon.size() == 3)
      {
        // check which point index is aligned
        if (PointIsAligned(circleCenter,
                           curvedEdgeEnd,
                           polygonVertices.col(subPolygon[1])))
        {
          subTriangle[1] = subPolygon[1];
          subTriangle[2] = subPolygon[2];
        }
        else
        {
          subTriangle[1] = subPolygon[0];
          subTriangle[2] = subPolygon[1];
        }
      }
      else
      {
        subTriangle[1] = subPolygon[1];
        subTriangle[2] = subPolygon[subPolygon.size() - 2];
      }
    }

    // Internal triangles
    result.InternalTriangles.resize(result.SubPolygons.size(),
                                    vector<unsigned int> { cirlceCenterIndex,
                                                           (unsigned int)newCoordinates.size(),
                                                           (unsigned int)newCoordinates.size() });

    for (unsigned int sp = 0; sp < result.SubPolygons.size(); sp++)
    {
      const vector<unsigned int>& subPolygon = result.SubPolygons[sp];
      vector<unsigned int>& internalTriangle = result.InternalTriangles[sp];

      internalTriangle[1] = subPolygon[0];
      internalTriangle[2] = subPolygon[subPolygon.size() - 1];
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::CircleDivisionByPolygonResult GeometryUtilities::CircleDivisionByPolygon(const Eigen::MatrixXd& polygonVertices,
                                                                                              const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                              const Eigen::Vector3d& circleCenter,
                                                                                              const double& circleRadius,
                                                                                              const unsigned int& curvedEdgeIndex) const
  {
    CircleDivisionByPolygonResult result;

    list<Vector3d> newCoordinates;

    // insert polygon and circle center in output
    const unsigned int& numPolygonVertices = polygonVertices.cols();
    const unsigned int cirlceCenterIndex = numPolygonVertices;

    for (unsigned int p = 0; p < numPolygonVertices; p++)
      newCoordinates.push_back(polygonVertices.col(p));
    newCoordinates.push_back(circleCenter);

    const unsigned int curvedEdgeOriginIndex = curvedEdgeIndex;
    const unsigned int curvedEdgeEndIndex = (curvedEdgeIndex + 1) % numPolygonVertices;
    const Vector3d& curvedEdgeOrigin = polygonVertices.col(curvedEdgeOriginIndex);
    const Vector3d& curvedEdgeEnd = polygonVertices.col(curvedEdgeEndIndex);

    // intersects all the lines from circleCenter to each vertex not belonging to curved edge
    unsigned int edgeNumber = curvedEdgeEndIndex;
    unsigned int triangleVertexNumber = curvedEdgeEndIndex;
    const unsigned int lastCheck = (curvedEdgeIndex == 0) ?  numPolygonVertices - 1 : curvedEdgeIndex - 1;

    list<vector<unsigned int>> subTriangles;
    list<vector<unsigned int>> internalTriangles;
    list<list<unsigned int>> subPolygons;

    while (edgeNumber != lastCheck)
    {
      const unsigned int originEdgeIndex = edgeNumber;
      const unsigned int endEdgeIndex = (edgeNumber + 1) % numPolygonVertices;

      const Eigen::Vector3d& edgeOrigin = polygonVertices.col(originEdgeIndex);
      const Eigen::Vector3d& edgeEnd = polygonVertices.col(endEdgeIndex);
      const Eigen::Vector3d& edgeTangent = polygonEdgeTangents.col(edgeNumber);

      // check if the edge is aligned to the angle quadrant
      if (PointIsAligned(edgeOrigin,
                         edgeEnd,
                         circleCenter))
      {
        edgeNumber++;
        if (edgeNumber == numPolygonVertices)
          edgeNumber = 0;

        continue;
      }

      // check if the next edge is aligned to the angle quadrant
      const unsigned int nextEdgeEndIndex = (endEdgeIndex + 1) % numPolygonVertices;

      if (PointIsAligned(edgeEnd,
                         polygonVertices.col(nextEdgeEndIndex),
                         circleCenter))
        break;

      // intersecting the line with the circle
      const Eigen::Vector3d endEdgeCenterTangent = SegmentTangent(edgeEnd,
                                                                  circleCenter);

      const IntersectionSegmentCircleResult intersection = IntersectionSegmentCircle(edgeEnd,
                                                                                     circleCenter,
                                                                                     circleCenter,
                                                                                     circleRadius);

      Output::Assert(intersection.Type == IntersectionSegmentCircleResult::Types::TwoIntersections);

      double intersectionCoordinate = 0.0;
      for (unsigned int i = 0; i < 2; i++)
      {
        if (intersection.SegmentIntersections[i].Type != PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
          continue;

        intersectionCoordinate = intersection.SegmentIntersections[i].CurvilinearCoordinate;
        break;
      }

      Output::Assert(IsValue1DNegative(intersectionCoordinate));

      // create new coordinate
      const Eigen::Vector3d intersectionPoint = edgeEnd + intersectionCoordinate * endEdgeCenterTangent;
      newCoordinates.push_back(intersectionPoint);

      // create new sub-triangle
      subTriangles.push_back(vector<unsigned int>({ cirlceCenterIndex, 0, 0 }));
      internalTriangles.push_back(vector<unsigned int>({ cirlceCenterIndex, 0, 0 }));
      subPolygons.push_back(list<unsigned int>());

      vector<unsigned int>& subTriangle = subTriangles.back();
      vector<unsigned int>& internalTriangle = internalTriangles.back();
      list<unsigned int>& subPolygon = subPolygons.back();

      subTriangle[1] = newCoordinates.size() - 1;
      subTriangle[2] = triangleVertexNumber;

      internalTriangle[1] = endEdgeIndex;
      internalTriangle[2] = originEdgeIndex;

      subPolygon.push_back(triangleVertexNumber);
      if (originEdgeIndex != triangleVertexNumber)
        subPolygon.push_back(originEdgeIndex);
      subPolygon.push_back(endEdgeIndex);
      subPolygon.push_back(newCoordinates.size() - 1);

      triangleVertexNumber = newCoordinates.size() - 1;

      edgeNumber++;
      if (edgeNumber == numPolygonVertices)
        edgeNumber = 0;
    }

    // insert last sub-triangle
    subTriangles.push_back(vector<unsigned int>({ cirlceCenterIndex, 0, 0 }));
    internalTriangles.push_back(vector<unsigned int>({ cirlceCenterIndex, 0, 0 }));
    subPolygons.push_back(list<unsigned int>());

    vector<unsigned int>& subTriangle = subTriangles.back();
    vector<unsigned int>& internalTriangle = internalTriangles.back();
    list<unsigned int>& subPolygon = subPolygons.back();

    subTriangle[1] = curvedEdgeOriginIndex;
    subTriangle[2] = triangleVertexNumber;

    internalTriangle[1] = (edgeNumber + 1) % numPolygonVertices;
    internalTriangle[2] = edgeNumber;

    subPolygon.push_back(triangleVertexNumber);
    subPolygon.push_back(edgeNumber);
    if ((edgeNumber + 1) % numPolygonVertices != curvedEdgeOriginIndex)
      subPolygon.push_back((edgeNumber + 1) % numPolygonVertices);
    subPolygon.push_back(curvedEdgeOriginIndex);

    // convert newCoordinates
    result.Points.setZero(3, newCoordinates.size());
    unsigned int counter = 0;
    for (const Vector3d& point : newCoordinates)
      result.Points.col(counter++)<< point;

    // convert sub-triangles
    result.SubTriangles = vector<vector<unsigned int>>(subTriangles.begin(),
                                                       subTriangles.end());

    // convert sub-triangles
    result.InternalTriangles = vector<vector<unsigned int>>(internalTriangles.begin(),
                                                            internalTriangles.end());

    // convert sub-polygons
    result.SubPolygons.resize(subPolygons.size());
    unsigned int sp = 0;
    for (list<unsigned int> subPolygon : subPolygons)
      result.SubPolygons[sp++] = vector<unsigned int>(subPolygon.begin(), subPolygon.end());

    return result;
  }
  // ***************************************************************************
  double GeometryUtilities::PolygonArea(const Eigen::MatrixXd& polygonVertices) const
  {
    Output::Assert(PointsAre2D(polygonVertices) && polygonVertices.cols() > 2);

    const unsigned int& numVertices = polygonVertices.cols();
    double area = 0.0;
    Eigen::VectorXd xPoints(numVertices + 1);
    Eigen::VectorXd yPoints(numVertices + 1);
    xPoints<< polygonVertices.row(0).transpose(), polygonVertices(0, 0);
    yPoints<< polygonVertices.row(1).transpose(), polygonVertices(1, 0);

    return 0.5 * (xPoints.segment(0, numVertices).dot(yPoints.segment(1, numVertices)) -
                  xPoints.segment(1, numVertices).dot(yPoints.segment(0, numVertices)));
  }
  // ***************************************************************************
  GeometryUtilities::PolygonCirclePositionTypes GeometryUtilities::PolygonCirclePosition(const Eigen::MatrixXd& polygonVertices,
                                                                                         const Eigen::Vector3d& circleCenter,
                                                                                         const double& circleRadius,
                                                                                         const vector<PointCirclePositionResult>& vertexPositions,
                                                                                         const IntersectionPolygonCircleResult& polygonCircleIntersections) const
  {
    // check circle center position respect the polygon
    GeometryUtilities::PointPolygonPositionResult centerPosition = PointPolygonPosition(circleCenter,
                                                                                        polygonVertices);
    Output::Assert(centerPosition.Type != Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Unknown);

    const unsigned int& numVertices = polygonVertices.cols();
    const unsigned int& numIntersections = polygonCircleIntersections.Intersections.size();

    // compute polygon vertices position respect the circle
    bool oneVertexOutsideCircle = false;
    for (unsigned int v = 0; v < numVertices; v++)
    {
      if (vertexPositions[v] == PointCirclePositionResult::Outside)
        oneVertexOutsideCircle = true;
    }

    switch (centerPosition.Type)
    {
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside:
      {
        // the circle is outside the polygon
        if (oneVertexOutsideCircle &&
            numIntersections == 0) // vertices are far from the circle, no intersection found
          return PolygonCirclePositionTypes::PolygonOutsideCircleNoIntersection;
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 0) // all vertices are inside the circle and there are no intersections then the polygon is inside the circle
          return PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection;
        else if (oneVertexOutsideCircle &&
                 numIntersections == 1) // at least one vertex is outside the circle and only one intersection found, the circle is outside and touches the polygon in one vertex
        {
          if (polygonCircleIntersections.Intersections[0].IndexType == IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex)
            return PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionOnVertex;
          else if (polygonCircleIntersections.Intersections[0].IndexType == IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge &&
                   polygonCircleIntersections.Intersections[0].Type == IntersectionPolygonCircleResult::Intersection::Types::Tangent)
            return PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionTangentOnEdge;
        }
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 1) // all vertices are inside the circle and one intersection found
          return PolygonCirclePositionTypes::PolygonInsideCircleOneVertexIntersection;
        else if (numIntersections > 1)
          return PolygonCirclePositionTypes::CirclePolygonMultipleIntersections;
      }
      break;
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex:
      case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside:
      {
        // the circle center is inside the polygon
        if (oneVertexOutsideCircle &&
            numIntersections == 0) // if at least one vertex is ouside the circle and there are no intersections than the circle is inside the polygon
          return PolygonCirclePositionTypes::CircleInsidePolygonNoIntersection;
        else if (!oneVertexOutsideCircle &&
                 numIntersections == 0) // if at least a vertex is not outside the circle and there are no intersections then the polygon is inside the circle
          return PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection;
        else if (numIntersections == 1) // only one intersection found, the circle touch the polygon tangent in one edge
          return PolygonCirclePositionTypes::CircleInsidePolygonOneIntersectionTangentOnEdge;
        else if (numIntersections > 1)
        {
          // check if the circle intersects only the vertices (a sort of circumscribed circle)
          bool onlyVerticesIntersections = true;
          for (unsigned int i = 0; i < numIntersections; i++)
          {
            if (polygonCircleIntersections.Intersections[i].IndexType !=
                IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex)
            {
              onlyVerticesIntersections = false;
              break;
            }
          }

          if (onlyVerticesIntersections)
            return PolygonCirclePositionTypes::PolygonInsideCircleIntersectionOnlyOnVertices;

          return PolygonCirclePositionTypes::CirclePolygonMultipleIntersections;
        }
      }
      break;
      default:
      break;
    }

    throw runtime_error("PolygonCirclePosition failed");
  }
  // ***************************************************************************
}
