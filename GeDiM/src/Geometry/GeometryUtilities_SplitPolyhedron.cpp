#include "GeometryUtilities.hpp"
#include <string>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::SplitPolygonWithPlaneResult GeometryUtilities::SplitPolygonWithPlane(const Eigen::MatrixXd& polygonVertices,
                                                                                          const Eigen::MatrixXd& polygonEdgeTangents,
                                                                                          const Eigen::Vector3d& planeNormal,
                                                                                          const Eigen::Vector3d& planeOrigin,
                                                                                          const Eigen::Vector3d& polygonTranslation,
                                                                                          const Eigen::Matrix3d& polygonRotationMatrix) const
  {
    SplitPolygonWithPlaneResult result;
    
    bool positiveUsed = false;
    bool negativeUsed = false;
    
    std::list<unsigned int> positiveVerticesIndex;
    std::list<unsigned int> negativeVerticesIndex;
    std::list<unsigned int> pointsOnPlaneIndex;
    std::list<Eigen::Vector3d> newVertices;
    std::list<unsigned int> newVerticesEdgeIndex;
    
    // note: for convex polygons only
    for (unsigned int i = 0; i < polygonVertices.cols(); i++)
    {
      // Place the vertex in the front or back polygon (or both)
      const unsigned int& vertexIndex = i;
      PointPlanePositionTypes vertexPlanePosition = PointPlanePosition(PointPlaneDistance(polygonVertices.col(vertexIndex),
                                                                                          planeNormal,
                                                                                          planeOrigin));
      switch (vertexPlanePosition)
      {
        case PointPlanePositionTypes::Positive:
          positiveVerticesIndex.push_back(vertexIndex);
          positiveUsed = true;
          break;
          
        case PointPlanePositionTypes::OnPlane:
          positiveVerticesIndex.push_back(vertexIndex);
          negativeVerticesIndex.push_back(vertexIndex);
          pointsOnPlaneIndex.push_back(vertexIndex);
          break;
          
        case PointPlanePositionTypes::Negative:
          negativeVerticesIndex.push_back(vertexIndex);
          negativeUsed = true;
          break;

        default:
          throw runtime_error("Unsupported plane position");
      }
      
      // Add a vertex to both polygons where edges intersect the plane
      unsigned int nextVertexIndex = (i + 1) % polygonVertices.cols();
      
      PointPlanePositionTypes nextVertexPlanePosition = PointPlanePosition(PointPlaneDistance(polygonVertices.col(nextVertexIndex),
                                                                                              planeNormal,
                                                                                              planeOrigin));
      if ((vertexPlanePosition == PointPlanePositionTypes::Positive &&
           nextVertexPlanePosition == PointPlanePositionTypes::Negative) ||
          (vertexPlanePosition == PointPlanePositionTypes::Negative &&
           nextVertexPlanePosition == PointPlanePositionTypes::Positive))
      {
        IntersectionSegmentPlaneResult edgePlaneIntersection = IntersectionSegmentPlane(polygonVertices.col(vertexIndex),
                                                                                        polygonVertices.col(nextVertexIndex),
                                                                                        planeNormal,
                                                                                        planeOrigin);
        Output::Assert(edgePlaneIntersection.Type == IntersectionSegmentPlaneResult::Types::SingleIntersection &&
                       edgePlaneIntersection.SingleIntersection.Type == PointSegmentPositionTypes::InsideSegment);
        
        const unsigned int newVertexIndex = polygonVertices.cols() +
                                            newVertices.size();
        
        newVertices.push_back(polygonVertices.col(vertexIndex) +
                              edgePlaneIntersection.SingleIntersection.CurvilinearCoordinate *
                              polygonEdgeTangents.col(vertexIndex));
        newVerticesEdgeIndex.push_back(vertexIndex);
        pointsOnPlaneIndex.push_back(newVertexIndex);
        positiveVerticesIndex.push_back(newVertexIndex);
        negativeVerticesIndex.push_back(newVertexIndex);
      }
    }
    
    vector<unsigned int> positiveVertices(positiveVerticesIndex.begin(), positiveVerticesIndex.end());
    vector<unsigned int> negativeVertices(negativeVerticesIndex.begin(), negativeVerticesIndex.end());
    result.PointsOnPlane = vector<unsigned int>(pointsOnPlaneIndex.begin(), pointsOnPlaneIndex.end());
    result.NewVertices = vector<Eigen::Vector3d>(newVertices.begin(), newVertices.end());
    result.NewVerticesEdgeIndex = vector<unsigned int>(newVerticesEdgeIndex.begin(), newVerticesEdgeIndex.end());
    result.Type = (positiveUsed) ? (negativeUsed ? SplitPolygonWithPlaneResult::Types::Split :
                                                   SplitPolygonWithPlaneResult::Types::Positive) :
                                   negativeUsed ? SplitPolygonWithPlaneResult::Types::Negative :
                                                  SplitPolygonWithPlaneResult::Types::OnPlane;

    // sort points with correct order
    Eigen::MatrixXd globalVertices(3, polygonVertices.cols() + result.NewVertices.size());
    globalVertices.block(0, 0, 3, polygonVertices.cols()) = polygonVertices;
    for (unsigned int v = 0; v < result.NewVertices.size(); v++)
      globalVertices.col(polygonVertices.cols() + v)<< result.NewVertices[v];

    if (result.Type == SplitPolygonWithPlaneResult::Types::Split ||
        result.Type == SplitPolygonWithPlaneResult::Types::Positive)
    {
      const Eigen::MatrixXd positive3DVertices = ExtractPoints(globalVertices,
                                                               positiveVertices);
      const Eigen::MatrixXd positive2DVertices = RotatePointsFrom3DTo2D(positive3DVertices,
                                                                        polygonRotationMatrix.transpose(),
                                                                        polygonTranslation);
      const vector<unsigned int> convexHull = ConvexHull(positive2DVertices);
      Output::Assert(convexHull.size() == positiveVertices.size());

      result.PositiveVertices.resize(convexHull.size());
      for (unsigned int c = 0; c < convexHull.size(); c++)
        result.PositiveVertices[c] = positiveVertices[convexHull[c]];
    }

    if (result.Type == SplitPolygonWithPlaneResult::Types::Split ||
        result.Type == SplitPolygonWithPlaneResult::Types::Negative)
    {
      const Eigen::MatrixXd negative3DVertices = ExtractPoints(globalVertices,
                                                               negativeVertices);
      const Eigen::MatrixXd negative2DVertices = RotatePointsFrom3DTo2D(negative3DVertices,
                                                                        polygonRotationMatrix.transpose(),
                                                                        polygonTranslation);

      const vector<unsigned int> convexHull = ConvexHull(negative2DVertices);
      Output::Assert(convexHull.size() == negativeVertices.size());

      result.NegativeVertices.resize(convexHull.size());
      for (unsigned int c = 0; c < convexHull.size(); c++)
        result.NegativeVertices[c] = negativeVertices[convexHull[c]];
    }

    return result;
  }
  // ***************************************************************************
  GeometryUtilities::SplitPolyhedronWithPlaneResult GeometryUtilities::SplitPolyhedronWithPlane(const Eigen::MatrixXd& polyhedronVertices,
                                                                                                const Eigen::MatrixXi& polyhedronEdges,
                                                                                                const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                                const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                                const vector<Eigen::MatrixXd>& polyhedronFaceEdgeTangents,
                                                                                                const vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                                                const vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                                                                                const Eigen::Vector3d& planeNormal,
                                                                                                const Eigen::Vector3d& planeOrigin,
                                                                                                const Eigen::Matrix3d& planeRotationMatrix,
                                                                                                const Eigen::Vector3d& planeTranslation) const
  {
    std::list<Eigen::Vector3d> newVertices;
    std::unordered_map<unsigned int, unsigned int> newVerticesByEdgeIndex;
    std::unordered_map<unsigned int, unsigned int> newVerticesEdgeIndex;
    std::set<unsigned int> pointsOnPlaneIndices;
    std::set<unsigned int> positivePolyhedronVerticesIndices;
    std::set<unsigned int> negativePolyhedronVerticesIndices;
    std::list<Eigen::MatrixXi> positivePolyhedronFaces;
    std::list<Eigen::MatrixXi> negativePolyhedronFaces;
    std::list<int> positivePolyhedronOriginalFacesIndices;
    std::list<int> negativePolyhedronOriginalFacesIndices;
    bool positiveUsed = false, negativeUsed = false;
    
    const unsigned int numVertices = polyhedronVertices.cols();

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const unsigned int numFaceVertices = polyhedronFaces[f].cols();
      const Eigen::MatrixXd& faceVertices = polyhedronFaceVertices[f];
      Eigen::MatrixXd faceEdgeTangents = polyhedronFaceEdgeTangents[f];
      SplitPolygonWithPlaneResult splitFaceByPlane = SplitPolygonWithPlane(faceVertices,
                                                                           faceEdgeTangents,
                                                                           planeNormal,
                                                                           planeOrigin,
                                                                           polyhedronFaceTranslations[f],
                                                                           polyhedronFaceRotationMatrices[f]);
      
      switch (splitFaceByPlane.Type)
      {
        case SplitPolygonWithPlaneResult::Types::Positive:
        {
          const unsigned int& positiveFaceVertices = splitFaceByPlane.PositiveVertices.size();
          
          positivePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, positiveFaceVertices));
          positivePolyhedronOriginalFacesIndices.push_back(f);

          Eigen::MatrixXi& positiveFace = positivePolyhedronFaces.back();
          for (unsigned int v = 0; v < positiveFaceVertices; v++)
          {
            const unsigned int polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]);

            if (positivePolyhedronVerticesIndices.find(polyhedronVertexIndex) == positivePolyhedronVerticesIndices.end())
              positivePolyhedronVerticesIndices.insert(polyhedronVertexIndex);

            positiveFace(0, v) = polyhedronVertexIndex;
          }
          
          positiveUsed = true;
        }
          break;
          
        case SplitPolygonWithPlaneResult::Types::OnPlane:
        case SplitPolygonWithPlaneResult::Types::Split:
        {
          vector<unsigned int> newVerticesIndices(splitFaceByPlane.NewVertices.size());
          for (unsigned int nv = 0; nv < splitFaceByPlane.NewVertices.size(); nv++)
          {
            const unsigned int polyhedronEdgeIndex = polyhedronFaces[f](1, splitFaceByPlane.NewVerticesEdgeIndex[nv]);
            if (newVerticesByEdgeIndex.find(polyhedronEdgeIndex) == newVerticesByEdgeIndex.end())
            {
              // Add new vertex on the list
              const unsigned int polyhedronNewVertexIndex = polyhedronVertices.cols() + newVertices.size();
              newVertices.push_back(splitFaceByPlane.NewVertices[nv]);
              newVerticesByEdgeIndex.insert(make_pair(polyhedronEdgeIndex, polyhedronNewVertexIndex));
              newVerticesEdgeIndex.insert(make_pair(polyhedronNewVertexIndex, polyhedronEdgeIndex));
            }

            newVerticesIndices[nv] = newVerticesByEdgeIndex.at(polyhedronEdgeIndex);
          }

          for (unsigned int p = 0; p < splitFaceByPlane.PointsOnPlane.size(); p++)
          {
            unsigned int polyhedronVertexIndex = 0;

            if (splitFaceByPlane.PointsOnPlane[p] < numFaceVertices)
              polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.PointsOnPlane[p]);
            else
              polyhedronVertexIndex = newVerticesIndices[splitFaceByPlane.PointsOnPlane[p] - numFaceVertices];

            if (pointsOnPlaneIndices.find(polyhedronVertexIndex) == pointsOnPlaneIndices.end())
              pointsOnPlaneIndices.insert(polyhedronVertexIndex);
          }

          const unsigned int& positiveFaceVertices = splitFaceByPlane.PositiveVertices.size();
          positivePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, positiveFaceVertices));
          positivePolyhedronOriginalFacesIndices.push_back(f);

          Eigen::MatrixXi& positiveFace = positivePolyhedronFaces.back();
          for (unsigned int v = 0; v < positiveFaceVertices; v++)
          {
            unsigned int polyhedronVertexIndex = 0;

            if (splitFaceByPlane.PositiveVertices[v] < numFaceVertices)
              polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.PositiveVertices[v]);
            else
              polyhedronVertexIndex = newVerticesIndices[splitFaceByPlane.PositiveVertices[v] - numFaceVertices];

            if (positivePolyhedronVerticesIndices.find(polyhedronVertexIndex) == positivePolyhedronVerticesIndices.end())
              positivePolyhedronVerticesIndices.insert(polyhedronVertexIndex);
            positiveFace(0, v) = polyhedronVertexIndex;
          }

          const unsigned int& negativeFaceVertices = splitFaceByPlane.NegativeVertices.size();
          negativePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, negativeFaceVertices));
          negativePolyhedronOriginalFacesIndices.push_back(f);

          Eigen::MatrixXi& negativeFace = negativePolyhedronFaces.back();
          for (unsigned int v = 0; v < negativeFaceVertices; v++)
          {
            unsigned int polyhedronVertexIndex = 0;

            if (splitFaceByPlane.NegativeVertices[v] < numFaceVertices)
              polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]);
            else
              polyhedronVertexIndex = newVerticesIndices[splitFaceByPlane.NegativeVertices[v] - numFaceVertices];

            if (negativePolyhedronVerticesIndices.find(polyhedronVertexIndex) == negativePolyhedronVerticesIndices.end())
              negativePolyhedronVerticesIndices.insert(polyhedronVertexIndex);
            negativeFace(0, v) = polyhedronVertexIndex;
          }

          if (splitFaceByPlane.Type == SplitPolygonWithPlaneResult::Types::Split)
          {
            positiveUsed = true;
            negativeUsed = true;
          }
        }
          break;

        case SplitPolygonWithPlaneResult::Types::Negative:
        {
          const unsigned int& negativeFaceVertices = splitFaceByPlane.NegativeVertices.size();

          negativePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, negativeFaceVertices));
          negativePolyhedronOriginalFacesIndices.push_back(f);

          Eigen::MatrixXi& negativeFace = negativePolyhedronFaces.back();
          for (unsigned int v = 0; v < negativeFaceVertices; v++)
          {
            const unsigned int polyhedronVertexIndex = polyhedronFaces[f](0, splitFaceByPlane.NegativeVertices[v]);

            if (negativePolyhedronVerticesIndices.find(polyhedronVertexIndex) == negativePolyhedronVerticesIndices.end())
              negativePolyhedronVerticesIndices.insert(polyhedronVertexIndex);

            negativeFace(0, v) = polyhedronVertexIndex;
          }

          negativeUsed = true;
        }
          break;
        default:
          throw runtime_error("Unsupported split polygon with plane type");
      }
    }

    SplitPolyhedronWithPlaneResult result;

    // no split is necessary, return original polyhedron
    if (!positiveUsed || !negativeUsed)
    {
      result.Type = SplitPolyhedronWithPlaneResult::Types::None;

      return result;
    }

    // Split is done, create the last face to fill in the holes inside the new polyhedra
    result.Type = SplitPolyhedronWithPlaneResult::Types::Split;

    // Create the last face to fill in the holes inside the new polyhedra
    vector<Eigen::Vector3d> newVerticesVector = vector<Eigen::Vector3d>(newVertices.begin(),
                                                                        newVertices.end());
    vector<unsigned int> pointsOnPlaneIndicesVector = vector<unsigned int>(pointsOnPlaneIndices.begin(),
                                                                           pointsOnPlaneIndices.end());
    Eigen::MatrixXd pointsOnPlane3D(3, pointsOnPlaneIndicesVector.size());

    for (unsigned int p = 0; p < pointsOnPlaneIndicesVector.size(); p++)
    {
      const unsigned int pointOnPlaneIndex = pointsOnPlaneIndicesVector[p];
      if (pointOnPlaneIndex < numVertices)
        pointsOnPlane3D.col(p)<< polyhedronVertices.col(pointOnPlaneIndex);
      else
        pointsOnPlane3D.col(p)<< newVerticesVector[pointOnPlaneIndex - numVertices];
    }

    Eigen::MatrixXd pointsOnPlane2D = RotatePointsFrom3DTo2D(pointsOnPlane3D,
                                                             planeRotationMatrix.transpose(),
                                                             planeTranslation);

    vector<unsigned int> convexHull = ConvexHull(pointsOnPlane2D);

    positivePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, convexHull.size()));
    positivePolyhedronOriginalFacesIndices.push_back(-1);
    Eigen::MatrixXi& positiveNewFace = positivePolyhedronFaces.back();
    for (unsigned int v = 0; v < convexHull.size(); v++)
    {
      const unsigned int polyhedronVertexIndex = pointsOnPlaneIndicesVector[convexHull[v]];
      positiveNewFace(0, v) = polyhedronVertexIndex;
    }

    negativePolyhedronFaces.push_back(Eigen::MatrixXi::Zero(2, convexHull.size()));
    negativePolyhedronOriginalFacesIndices.push_back(-1);
    Eigen::MatrixXi& negativeNewFace = negativePolyhedronFaces.back();
    for (unsigned int v = 0; v < convexHull.size(); v++)
    {
      const unsigned int polyhedronVertexIndex = pointsOnPlaneIndicesVector[convexHull[v]];
      negativeNewFace(0, v) = polyhedronVertexIndex;
    }

    // Craete new edges
    unordered_map<string, unsigned int> originalEdges;
    for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      const string key = to_string(polyhedronEdges(0, e)) + "-" +
                         to_string(polyhedronEdges(1, e));

      originalEdges.insert(make_pair(key, e));
    }

    unordered_map<string, unsigned int> newEdges;
    unordered_map<string, vector<unsigned int>> newEdgesOriginEnd;
    unordered_map<unsigned int, int> newEdgesOriginalEdges;
    set<unsigned int> positiveEdges;
    set<unsigned int> negativeEdges;

    // create edges of the new face first
    for (unsigned int v = 0; v < positiveNewFace.cols(); v++)
    {
      const unsigned int edgeOrigin = positiveNewFace(0, v);
      const unsigned int edgeEnd = positiveNewFace(0, (v + 1) % positiveNewFace.cols());
      const string edgeFrom = to_string(edgeOrigin) + "-" + to_string(edgeEnd);
      const string edgeTo = to_string(edgeEnd) + "-" + to_string(edgeOrigin);

      unsigned int newEdgeIndex = newEdges.size();

      if (originalEdges.find(edgeFrom) != originalEdges.end())
      {
        newEdges.insert(make_pair(edgeFrom, newEdges.size()));
        newEdgesOriginEnd.insert(make_pair(edgeFrom, vector<unsigned int>({ edgeOrigin, edgeEnd })));
        newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, originalEdges.at(edgeFrom)));
      }
      else if (originalEdges.find(edgeTo) != originalEdges.end())
      {
        newEdges.insert(make_pair(edgeTo, newEdges.size()));
        newEdgesOriginEnd.insert(make_pair(edgeTo, vector<unsigned int>({ edgeEnd, edgeOrigin })));
        newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, originalEdges.at(edgeTo)));
      }
      else
      {
        newEdges.insert(make_pair(edgeFrom, newEdges.size()));
        newEdgesOriginEnd.insert(make_pair(edgeFrom, vector<unsigned int>({ edgeOrigin, edgeEnd })));
        newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, -1));
      }
    }

    // crate all the other edges
    for (unsigned int p = 0; p < 2; p++)
    {
      list<Eigen::MatrixXi>& polyhedronFaces = (p == 0) ? positivePolyhedronFaces :
                                                          negativePolyhedronFaces;
      set<unsigned int>& polyhedronEdges = (p == 0) ? positiveEdges :
                                                      negativeEdges;

      for (Eigen::MatrixXi& newFace : polyhedronFaces)
      {
        for (unsigned int v = 0; v < newFace.cols(); v++)
        {
          const unsigned int edgeOrigin = newFace(0, v);
          const unsigned int edgeEnd = newFace(0, (v + 1) % newFace.cols());
          const string edgeFrom = to_string(edgeOrigin) + "-" + to_string(edgeEnd);
          const string edgeTo = to_string(edgeEnd) + "-" + to_string(edgeOrigin);

          unsigned int newEdgeIndex = 0;

          if (newEdges.find(edgeFrom) != newEdges.end())
            newEdgeIndex = newEdges.at(edgeFrom);
          else if (newEdges.find(edgeTo) != newEdges.end())
            newEdgeIndex = newEdges.at(edgeTo);
          else if (originalEdges.find(edgeFrom) != originalEdges.end())
          {
            newEdges.insert(make_pair(edgeFrom, newEdges.size()));
            newEdgesOriginEnd.insert(make_pair(edgeFrom, vector<unsigned int>({ edgeOrigin, edgeEnd })));
            newEdgeIndex = newEdges.at(edgeFrom);
            newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, originalEdges.at(edgeFrom)));
          }
          else if (originalEdges.find(edgeTo) != originalEdges.end())
          {
            newEdges.insert(make_pair(edgeTo, newEdges.size()));
            newEdgesOriginEnd.insert(make_pair(edgeTo, vector<unsigned int>({ edgeEnd, edgeOrigin })));
            newEdgeIndex = newEdges.at(edgeTo);
            newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, originalEdges.at(edgeTo)));
          }
          else
          {
            newEdges.insert(make_pair(edgeFrom, newEdges.size()));
            newEdgesOriginEnd.insert(make_pair(edgeFrom, vector<unsigned int>({ edgeOrigin, edgeEnd })));
            newEdgeIndex = newEdges.at(edgeFrom);

            if (newVerticesEdgeIndex.find(edgeOrigin) != newVerticesEdgeIndex.end() &&
                newVerticesEdgeIndex.find(edgeEnd) == newVerticesEdgeIndex.end())
              newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, newVerticesEdgeIndex.at(edgeOrigin)));
            else if (newVerticesEdgeIndex.find(edgeOrigin) == newVerticesEdgeIndex.end() &&
                     newVerticesEdgeIndex.find(edgeEnd) != newVerticesEdgeIndex.end())
              newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, newVerticesEdgeIndex.at(edgeEnd)));
            else
              newEdgesOriginalEdges.insert(make_pair(newEdgeIndex, -1));
          }

          polyhedronEdges.insert(newEdgeIndex);
          newFace(1, v) = newEdgeIndex;
        }
      }
    }

    // Create new polyhedra vertices
    result.Vertices.Vertices.setZero(3, numVertices + newVerticesVector.size());
    result.Vertices.NewVerticesOriginalEdge.resize(newVerticesVector.size());

    result.Vertices.Vertices.block(0, 0, 3, numVertices) = polyhedronVertices;
    for (unsigned int v = 0; v < newVerticesVector.size(); v++)
      result.Vertices.Vertices.col(numVertices + v)<< newVerticesVector[v];

    for (unordered_map<unsigned int, unsigned int>::const_iterator it = newVerticesByEdgeIndex.begin();
         it != newVerticesByEdgeIndex.end();
         it++)
      result.Vertices.NewVerticesOriginalEdge[it->second - numVertices] = it->first;

    result.PositivePolyhedron.Vertices.reserve(positivePolyhedronVerticesIndices.size());
    for (const unsigned int& v : positivePolyhedronVerticesIndices)
      result.PositivePolyhedron.Vertices.push_back(v);

    result.NegativePolyhedron.Vertices.reserve(negativePolyhedronVerticesIndices.size());
    for (const unsigned int& v : negativePolyhedronVerticesIndices)
      result.NegativePolyhedron.Vertices.push_back(v);

    // Create new polyhedra edges
    result.Edges.Edges.resize(2, newEdges.size());
    result.Edges.NewEdgesOriginalEdges.resize(newEdges.size(), -1);

    {
      for (const auto& it : newEdges)
      {
        const string& newEdgeKey = it.first;
        const unsigned int& newEdgeIndex = it.second;
        const vector<unsigned int>& newEdgeOriginEnd = newEdgesOriginEnd.at(newEdgeKey);
        const int& originalEdge = newEdgesOriginalEdges.at(newEdgeIndex);

        result.Edges.Edges(0, newEdgeIndex) = newEdgeOriginEnd[0];
        result.Edges.Edges(1, newEdgeIndex) = newEdgeOriginEnd[1];
        result.Edges.NewEdgesOriginalEdges[newEdgeIndex] = originalEdge;
      }
    }

    result.PositivePolyhedron.Edges.reserve(positiveEdges.size());
    for (const unsigned int& e : positiveEdges)
      result.PositivePolyhedron.Edges.push_back(e);

    result.NegativePolyhedron.Edges.reserve(negativeEdges.size());
    for (const unsigned int& e : negativeEdges)
      result.NegativePolyhedron.Edges.push_back(e);

    // Create new polyhedra faces
    vector<int> positivePolyhedronOriginalFaces = vector<int>(positivePolyhedronOriginalFacesIndices.begin(),
                                                              positivePolyhedronOriginalFacesIndices.end());
    vector<int> negativePolyhedronOriginalFaces = vector<int>(negativePolyhedronOriginalFacesIndices.begin(),
                                                              negativePolyhedronOriginalFacesIndices.end());

    result.Faces.Faces.resize(positivePolyhedronFaces.size() + negativePolyhedronFaces.size() - 1);
    result.Faces.NewFacesOriginalFaces.resize(positivePolyhedronFaces.size() + negativePolyhedronFaces.size() - 1, -1);
    result.PositivePolyhedron.Faces.resize(positivePolyhedronFaces.size());
    result.NegativePolyhedron.Faces.resize(negativePolyhedronFaces.size());

    {
      unsigned int f = 0;
      unsigned int nf = 0;
      unsigned int pf = 0;
      for (const Eigen::MatrixXi& newFace : positivePolyhedronFaces)
      {
        if (positivePolyhedronOriginalFaces[pf] < 0)
        {
          result.Faces.Faces[positivePolyhedronFaces.size() + negativePolyhedronFaces.size() - 2] = newFace;
          result.PositivePolyhedron.Faces[pf] = positivePolyhedronFaces.size() + negativePolyhedronFaces.size() - 2;
          pf++;
          continue;
        }

        result.Faces.Faces[f] = newFace;
        result.Faces.NewFacesOriginalFaces[f] = positivePolyhedronOriginalFaces[pf];
        result.PositivePolyhedron.Faces[pf] = f;
        pf++;
        f++;
      }
      for (const Eigen::MatrixXi& newFace : negativePolyhedronFaces)
      {
        if (negativePolyhedronOriginalFaces[nf] < 0)
        {
          result.NegativePolyhedron.Faces[nf] = positivePolyhedronFaces.size() + negativePolyhedronFaces.size() - 2;
          nf++;
          continue;
        }

        result.Faces.Faces[f] = newFace;
        result.Faces.NewFacesOriginalFaces[f] = negativePolyhedronOriginalFaces[nf];
        result.NegativePolyhedron.Faces[nf] = f;
        nf++;
        f++;
      }
    }

    //    cerr<< "RESULT"<< endl;
    //    cerr<< "**> result.Vertices.Vertices:\n"<< result.Vertices.Vertices<< endl;
    //    cerr<< "**> result.Vertices.NewVerticesOriginalEdge:\n"<< result.Vertices.NewVerticesOriginalEdge<< endl;
    //    cerr<< "**> result.Edges.Edges:\n"<< result.Edges.Edges<< endl;
    //    cerr<< "**> result.Edges.NewEdgesOriginalEdges:\n"<< result.Edges.NewEdgesOriginalEdges<< endl;
    //    cerr<< "**> result.Faces.Faces:\n"<< result.Faces.Faces<< endl;
    //    cerr<< "**> result.Faces.NewFacesOriginalFaces:\n"<< result.Faces.NewFacesOriginalFaces<< endl;
    //    cerr<< "**> result.PositivePolyhedron.Vertices:\n"<< result.PositivePolyhedron.Vertices<< endl;
    //    cerr<< "**> result.PositivePolyhedron.Edges:\n"<< result.PositivePolyhedron.Edges<< endl;
    //    cerr<< "**> result.PositivePolyhedron.Faces:\n"<< result.PositivePolyhedron.Faces<< endl;
    //    cerr<< "**> result.NegativePolyhedron.Vertices:\n"<< result.NegativePolyhedron.Vertices<< endl;
    //    cerr<< "**> result.NegativePolyhedron.Edges:\n"<< result.NegativePolyhedron.Edges<< endl;
    //    cerr<< "**> result.NegativePolyhedron.Faces:\n"<< result.NegativePolyhedron.Faces<< endl;

    return result;
  }
  // ***************************************************************************
  vector<GeometryUtilities::Polyhedron> GeometryUtilities::SplitPolyhedronWithPlaneResultToPolyhedra(const SplitPolyhedronWithPlaneResult& result) const
  {
    Output::Assert(result.Type == SplitPolyhedronWithPlaneResult::Types::Split);

    vector<GeometryUtilities::Polyhedron> polyhedra(2);

    for (unsigned int p = 0; p < 2; p++)
    {
      const SplitPolyhedronWithPlaneResult::NewPolyhedron& polyhedron = (p == 0) ? result.PositivePolyhedron :
                                                                                   result.NegativePolyhedron;

      unordered_map<unsigned int, unsigned int> mapVertices;
      polyhedra[p].Vertices = ExtractPoints(result.Vertices.Vertices,
                                            polyhedron.Vertices);
      for (unsigned int v = 0; v < polyhedron.Vertices.size(); v++)
        mapVertices.insert(make_pair(polyhedron.Vertices[v], v));

      unordered_map<unsigned int, unsigned int> mapEdges;
      polyhedra[p].Edges.resize(2, polyhedron.Edges.size());
      for (unsigned int e = 0; e < polyhedron.Edges.size(); e++)
        mapEdges.insert(make_pair(polyhedron.Edges[e], e));

      for (unsigned int e = 0; e < polyhedron.Edges.size(); e++)
      {
        const unsigned int edgeOrigin = result.Edges.Edges(0, polyhedron.Edges[e]);
        const unsigned int edgeEnd = result.Edges.Edges(1, polyhedron.Edges[e]);
        polyhedra[p].Edges(0, e) = mapVertices.at(edgeOrigin);
        polyhedra[p].Edges(1, e) = mapVertices.at(edgeEnd);
      }

      polyhedra[p].Faces.resize(polyhedron.Faces.size());
      for (unsigned int f = 0; f < polyhedron.Faces.size(); f++)
      {
        const Eigen::MatrixXi& face = result.Faces.Faces[polyhedron.Faces[f]];
        polyhedra[p].Faces[f].resize(2, face.cols());

        for (unsigned int v = 0; v < face.cols(); v++)
        {
          polyhedra[p].Faces[f](0, v) = mapVertices.at(face(0, v));
          polyhedra[p].Faces[f](1, v) = mapEdges.at(face(1, v));
        }
      }
    }

    return polyhedra;
  }
  // ***************************************************************************
}
