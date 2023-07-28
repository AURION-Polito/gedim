#include "MeshUtilities.hpp"

#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"
#include <numeric>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::FillMesh2D(const MatrixXd& cell0Ds,
                                 const MatrixXi& cell1Ds,
                                 const vector<MatrixXi>& cell2Ds,
                                 IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    // Create Cell0Ds
    Output::Assert(cell0Ds.rows() == 3);
    const unsigned int& numCell0Ds = cell0Ds.cols();
    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell0DsInsertCoordinates(cell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
      mesh.Cell0DSetState(v, true);

    // Create Cell1Ds
    Output::Assert(cell1Ds.rows() == 2);
    unsigned int numCell1Ds = cell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell1DsInsertExtremes(cell1Ds);
    for (int e = 0; e < cell1Ds.cols(); e++)
      mesh.Cell1DSetState(e, true);

    // Create Cell2Ds
    const unsigned int& numCell2Ds = cell2Ds.size();
    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const MatrixXi& polygon = cell2Ds[f];
      Output::Assert(polygon.rows() == 2);
      const unsigned int& numVertices = polygon.cols();

      mesh.Cell2DInitializeVertices(f, numVertices);
      mesh.Cell2DInitializeEdges(f, numVertices);

      for (unsigned int v = 0; v < numVertices; v++)
        mesh.Cell2DInsertVertex(f, v, polygon(0, v));
      for (unsigned int e = 0; e < numVertices; e++)
        mesh.Cell2DInsertEdge(f, e, polygon(1, e));

      mesh.Cell2DSetState(f, true);
    }
  }
  // ***************************************************************************
  MeshUtilities::ComputeMesh2DCell1DsResult MeshUtilities::ComputeMesh2DCell1Ds(const Eigen::MatrixXd& cell0Ds,
                                                                                const std::vector<Eigen::VectorXi>& cell2Ds) const
  {
    ComputeMesh2DCell1DsResult result;

    const unsigned int numVertices = cell0Ds.cols();
    const unsigned int numCell2Ds = cell2Ds.size();

    Eigen::SparseMatrix<unsigned int> edges;
    edges.resize(numVertices,
                 numVertices);

    std::list<Eigen::Triplet<unsigned int>> triplets;
    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];
        triplets.push_back(Eigen::Triplet<unsigned int>(origin, end, 1));
        triplets.push_back(Eigen::Triplet<unsigned int>(end, origin, 1));
      }
    }

    edges.setFromTriplets(triplets.begin(), triplets.end());
    edges.makeCompressed();

    unsigned int numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          it.valueRef() = 1 + numEdges++;
      }
    }

    result.Cell1Ds.resize(2, numEdges);

    numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          result.Cell1Ds.col(numEdges++)<< it.row(), it.col();
      }
    }

    result.Cell2Ds.resize(numCell2Ds);

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      Eigen::MatrixXi& cell2D = result.Cell2Ds.at(c);
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      cell2D.resize(2, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];

        cell2D(0, v) = origin;
        cell2D(1, v) = (origin < end) ? edges.coeff(origin, end) - 1 : edges.coeff(end, origin) - 1;
      }
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh2D(const CheckMesh2DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 2);

    // check Cell0D are 2D
    if (configuration.Cell0D_CheckCoordinates2D)
      Output::Assert(geometryUtilities.PointsAre2D(convexMesh.Cell0DsCoordinates()));

    // check Cell0D duplications
    if (configuration.Cell0D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell0DTotalNumber(); p1++)
      {
        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell0DTotalNumber(); p2++)
        {
          Output::Assert(!geometryUtilities.PointsAreCoincident(convexMesh.Cell0DCoordinates(p1),
                                                                convexMesh.Cell0DCoordinates(p2)));
        }
      }
    }

    if (configuration.Cell1D_CheckDuplications)
    {
      for (unsigned int e1 = 0; e1 < convexMesh.Cell1DTotalNumber(); e1++)
      {
        Output::Assert(convexMesh.Cell1DExists(convexMesh.Cell1DOrigin(e1),
                                               convexMesh.Cell1DEnd(e1)));
        Output::Assert(!convexMesh.Cell1DExists(convexMesh.Cell1DEnd(e1),
                                                convexMesh.Cell1DOrigin(e1)));

        for (unsigned int e2 = e1 + 1; e2 < convexMesh.Cell1DTotalNumber(); e2++)
        {
          Output::Assert(!(convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DEnd(e2)));
          Output::Assert(!(convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DEnd(e2)));
        }
      }
    }

    if (configuration.Cell1D_CheckNeighbours)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(convexMesh.Cell1DNumberNeighbourCell2D(e) == 2);

        if (convexMesh.Cell1DHasNeighbourCell2D(e, 0))
        {
          const unsigned int cell2DRight = convexMesh.Cell1DNeighbourCell2D(e, 0);
          const vector<unsigned int> cell2DEdges = convexMesh.Cell2DEdges(cell2DRight);

          // check edge orientation
          vector<unsigned int>::const_iterator it = std::find(cell2DEdges.begin(), cell2DEdges.end(), e);
          Output::Assert(it != cell2DEdges.end());

          const unsigned int cell2DEdgeIndex = std::distance(cell2DEdges.begin(), it);
          const unsigned int edgeOrigin = convexMesh.Cell2DVertex(cell2DRight,
                                                                  (cell2DEdgeIndex + 1) % cell2DEdges.size());
          const unsigned int edgeEnd = convexMesh.Cell2DVertex(cell2DRight,
                                                               cell2DEdgeIndex);

          Output::Assert(convexMesh.Cell1DExists(edgeOrigin,
                                                 edgeEnd) &&
                         convexMesh.Cell1DByExtremes(edgeOrigin,
                                                     edgeEnd) == e);
        }

        if (convexMesh.Cell1DHasNeighbourCell2D(e, 1))
        {
          const unsigned int cell2DLeft = convexMesh.Cell1DNeighbourCell2D(e, 1);
          const vector<unsigned int> cell2DEdges = convexMesh.Cell2DEdges(cell2DLeft);

          // check edge orientation
          vector<unsigned int>::const_iterator it = std::find(cell2DEdges.begin(), cell2DEdges.end(), e);
          Output::Assert(it != cell2DEdges.end());

          const unsigned int cell2DEdgeIndex = std::distance(cell2DEdges.begin(), it);
          const unsigned int edgeOrigin = convexMesh.Cell2DVertex(cell2DLeft,
                                                                  cell2DEdgeIndex);
          const unsigned int edgeEnd = convexMesh.Cell2DVertex(cell2DLeft,
                                                               (cell2DEdgeIndex + 1) % cell2DEdges.size());

          Output::Assert(convexMesh.Cell1DExists(edgeOrigin,
                                                 edgeEnd) &&
                         convexMesh.Cell1DByExtremes(edgeOrigin,
                                                     edgeEnd) == e);
        }
      }
    }

    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(geometryUtilities.IsValue1DPositive(
                         geometryUtilities.SegmentLength(convexMesh.Cell1DOriginCoordinates(e),
                                                         convexMesh.Cell1DEndCoordinates(e))));
      }
    }

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        for (unsigned int v = 0; v < convexMesh.Cell2DNumberVertices(p); v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % convexMesh.Cell2DNumberVertices(p));
          Output::Assert(convexMesh.Cell1DExists(eO, eE) || convexMesh.Cell1DExists(eE, eO));
          const unsigned int edgeFromVertices = convexMesh.Cell1DExists(eO, eE) ? convexMesh.Cell1DByExtremes(eO, eE) :
                                                                                  convexMesh.Cell1DByExtremes(eE, eO);
          Output::Assert(convexMesh.Cell2DEdge(p, v) == edgeFromVertices);
        }
      }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell2DTotalNumber(); p1++)
      {
        vector<unsigned int> cell2D1Vertices = convexMesh.Cell2DVertices(p1);
        sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
        vector<unsigned int> cell2D1Edges = convexMesh.Cell2DEdges(p1);
        sort(cell2D1Edges.begin(), cell2D1Edges.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell2DTotalNumber(); p2++)
        {
          vector<unsigned int> cell2D2Vertices = convexMesh.Cell2DVertices(p2);
          sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
          vector<unsigned int> cell2D2Edges = convexMesh.Cell2DEdges(p2);
          sort(cell2D2Edges.begin(), cell2D2Edges.end());

          Output::Assert(cell2D1Vertices.size() != cell2D2Vertices.size() || !equal(cell2D1Vertices.begin(),
                                                                                    cell2D1Vertices.end(),
                                                                                    cell2D2Vertices.begin()));
          Output::Assert(cell2D1Edges.size() != cell2D2Edges.size() || !equal(cell2D1Edges.begin(),
                                                                              cell2D1Edges.end(),
                                                                              cell2D2Edges.begin()));
        }
      }
    }

    if (configuration.Cell2D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const Eigen::MatrixXd cell2DVertices = convexMesh.Cell2DVerticesCoordinates(p);
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);
        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(convexCell2DUnalignedVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(convexCell2DUnalignedVertices,
                                                                                   convexHull);

        Output::Assert(geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices,
                                                         convexHullVertices));
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        Output::Assert(geometryUtilities.IsValue2DPositive(
                         geometryUtilities.PolygonArea(convexMesh.Cell2DVerticesCoordinates(p))));
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::Mesh2DFromPolygon(const Eigen::MatrixXd& polygonVertices,
                                        const vector<unsigned int> vertexMarkers,
                                        const vector<unsigned int> edgeMarkers,
                                        IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);
    const unsigned int numPolygonVertices = polygonVertices.cols();
    Output::Assert(vertexMarkers.size() == numPolygonVertices);
    Output::Assert(edgeMarkers.size() == numPolygonVertices);

    // Create Cell0Ds
    const unsigned int& numCell0Ds = numPolygonVertices;
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polygonVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numPolygonVertices;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                e,
                                (e + 1) % numPolygonVertices);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, edgeMarkers[e]);
    }

    // Create Cell2Ds
    const unsigned int& numCell2Ds = 1;
    mesh.Cell2DsInitialize(numCell2Ds);

    mesh.Cell2DInitializeVertices(0, numPolygonVertices);
    mesh.Cell2DInitializeEdges(0, numPolygonVertices);

    for (unsigned int v = 0; v < numPolygonVertices; v++)
      mesh.Cell2DInsertVertex(0, v, v);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
      mesh.Cell2DInsertEdge(0, e, e);

    mesh.Cell2DSetState(0, true);
    mesh.Cell2DSetMarker(0, 0);

    // Create Cell1D neighbours
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInitializeNeighbourCell2Ds(e, 2);
      mesh.Cell1DInsertNeighbourCell2D(e, 1, 0);
    }
  }
  // ***************************************************************************
  vector<unsigned int> MeshUtilities::MeshCell2DRoots(const IMeshDAO& mesh) const
  {
    vector<unsigned int> rootCell2Ds(mesh.Cell2DTotalNumber());

    for (unsigned int cc = 0; cc < mesh.Cell2DTotalNumber(); cc++)
    {
      unsigned int rootCell = cc;
      while (mesh.Cell2DHasOriginalCell2D(rootCell))
        rootCell = mesh.Cell2DOriginalCell2D(rootCell);

      rootCell2Ds[cc] = rootCell;
    }

    return rootCell2Ds;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(convexMesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell2DTotalNumber(); c++)
    {
      const unsigned int& domainCell2DIndex = c;

      // Extract original cell2D geometric information
      vector<Eigen::Matrix3d> convexCell2DTriangulationPoints;
      double convexCell2DArea;
      Eigen::Vector3d convexCell2DCentroid;

      const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // compute original cell2D triangulation
      const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
      const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                            convexCell2DUnalignedVerticesFilter);

      const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
      vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
      for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
        convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

      convexCell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                     convexCell2DTriangulation);

      const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();
      unsigned int cell2DTriangulationSize = numConvexCell2DTriangulation;

      // compute original cell2D area and centroids
      Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
      Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
      for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
      {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
      }

      convexCell2DArea = convexCell2DTriangulationAreas.sum();
      convexCell2DCentroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                               convexCell2DTriangulationAreas,
                                                               convexCell2DArea);

      result.Cell2DsVertices[c] = convexCell2DVertices;

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints.size(); cct++)
        result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cct];

      result.Cell2DsAreas[c] = convexCell2DArea;
      result.Cell2DsCentroids[c] = convexCell2DCentroid;
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      result.Cell2DsEdgeDirections[c].resize(convexMesh.Cell2DNumberEdges(domainCell2DIndex));
      for (unsigned int e = 0; e < convexMesh.Cell2DNumberEdges(domainCell2DIndex); e++)
      {
        const unsigned int origin = convexMesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = convexMesh.Cell2DVertex(domainCell2DIndex,
                                                         (e + 1) % convexMesh.Cell2DNumberEdges(domainCell2DIndex));

        result.Cell2DsEdgeDirections[c][e] = convexMesh.Cell1DExists(origin,
                                                                     end);
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const IMeshDAO& convexMesh,
                                                                            const vector<vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      const unsigned int& domainCell2DIndex = c;
      const vector<unsigned int>& domainConvexCell2DIndices = meshCell2DToConvexCell2DIndices[domainCell2DIndex];
      const unsigned int& numConvexCells = domainConvexCell2DIndices.size();

      // Get domain cell2D geometry information
      map<unsigned int, unsigned int> cell2DVerticesToPosition;
      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(domainCell2DIndex); v++)
        cell2DVerticesToPosition.insert(pair<unsigned int, unsigned int>(mesh.Cell2DVertex(domainCell2DIndex,
                                                                                           v),
                                                                         v));

      // Extract original cell2D geometric information
      unsigned int cell2DTriangulationSize = 0;
      vector<vector<Eigen::Matrix3d>> convexCell2DTriangulationPoints(numConvexCells);
      Eigen::VectorXd convexCell2DAreas(numConvexCells);
      Eigen::MatrixXd convexCell2DCentroids(3, numConvexCells);

      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        const unsigned int& domainConvexCell2DIndex = domainConvexCell2DIndices[cc];
        const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainConvexCell2DIndex);

        // compute original cell2D triangulation
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);

        const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
        vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
        for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
          convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

        convexCell2DTriangulationPoints[cc] = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                           convexCell2DTriangulation);

        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints[cc].size();
        cell2DTriangulationSize += numConvexCell2DTriangulation;

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cc][cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cc][cct],
                                                                                          convexCell2DTriangulationAreas[cct]);
        }

        convexCell2DAreas[cc] = convexCell2DTriangulationAreas.sum();
        convexCell2DCentroids.col(cc) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                          convexCell2DTriangulationAreas,
                                                                          convexCell2DAreas[cc]);
      }

      result.Cell2DsVertices[c] = mesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints[cc].size(); cct++)
          result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cc][cct];
      }

      result.Cell2DsAreas[c] = convexCell2DAreas.sum();
      result.Cell2DsCentroids[c] = geometryUtilities.PolygonCentroid(convexCell2DCentroids,
                                                                     convexCell2DAreas,
                                                                     result.Cell2DsAreas[c]);
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      result.Cell2DsEdgeDirections[c].resize(mesh.Cell2DNumberEdges(domainCell2DIndex));
      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(domainCell2DIndex); e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(domainCell2DIndex,
                                                   (e + 1) % mesh.Cell2DNumberEdges(domainCell2DIndex));

        result.Cell2DsEdgeDirections[c][e] = mesh.Cell1DExists(origin,
                                                               end);
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell1DCell2DNeighbours(IMeshDAO& mesh) const
  {
    // Initialize cell1D neighbours
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
      mesh.Cell1DInitializeNeighbourCell2Ds(c1D, 2);

    // Compute Cell1D neighbours starting from cell2Ds
    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
      const unsigned int numCell2DEdges = mesh.Cell2DNumberEdges(c2D);
      for (unsigned int e = 0; e < numCell2DEdges; e++)
      {
        const unsigned int cell1D = mesh.Cell2DEdge(c2D, e);
        const unsigned int edgeOrigin =  mesh.Cell2DVertex(c2D, e);
        const unsigned int edgeEnd =  mesh.Cell2DVertex(c2D, (e + 1) % numCell2DEdges);

        if (mesh.Cell1DExists(edgeOrigin,
                              edgeEnd)) // left cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           1,
                                           c2D);
        }
        else // right cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           0,
                                           c2D);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateRectangleMesh(const Eigen::Vector3d& rectangleOrigin,
                                          const Eigen::Vector3d& rectangleBaseTangent,
                                          const Eigen::Vector3d& rectangleHeightTangent,
                                          const vector<double>& baseMeshCurvilinearCoordinates,
                                          const vector<double>& heightMeshCurvilinearCoordinates,
                                          IMeshDAO& mesh) const
  {
    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    const unsigned int numCell0Ds = numBasePoints * numHeightPoints;
    const unsigned int numCell1Ds = numHeightPoints * (numBasePoints - 1) + numBasePoints * (numHeightPoints - 1);
    const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

    mesh.InitializeDimension(2);

    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell2DsInitialize(numCell2Ds);

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const Eigen::Vector3d coordinate = rectangleOrigin +
                                           baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                           heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;
        const unsigned int marker = 1 * (b == 0 && h == 0) +
                                    2 * (b == (numBasePoints - 1) && h == 0) +
                                    4 * (b == 0 && h == (numHeightPoints - 1)) +
                                    3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                    5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                    7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                    8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                    6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

        mesh.Cell0DSetState(cell0DIndex, true);
        mesh.Cell0DInsertCoordinates(cell0DIndex,
                                     coordinate);

        mesh.Cell0DSetMarker(cell0DIndex, marker);
        cell0DIndex++;
      }
    }

    // create cell1Ds
    unsigned int cell1DIndex = 0;

    // create horizontal cell1Ds
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + 1;

        const unsigned int marker = 5 * (h == 0) +
                                    7 * (h == (numHeightPoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create vertical cell1Ds
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

        const unsigned int marker = 8 * (b == 0) +
                                    6 * (b == (numBasePoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create cell2Ds
    unsigned int cell2DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DHorizontalIndex = b + h * (numBasePoints - 1);
        const unsigned int cell1DVerticalIndex = cell0DIndex + numHeightPoints * (numBasePoints - 1);

        vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                cell0DIndex + 1,
                                                cell0DIndex + numBasePoints + 1,
                                                cell0DIndex + numBasePoints };
        vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                             cell1DVerticalIndex + 1,
                                             cell1DHorizontalIndex + (numBasePoints - 1),
                                             cell1DVerticalIndex
                                           };

        mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
        mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

        mesh.Cell2DSetState(cell2DIndex, true);

        mesh.Cell2DSetMarker(cell2DIndex, 0);

        cell2DIndex++;
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateRectanglePlusHangingNodesMesh(const Eigen::Vector3d& rectangleOrigin,
                                                          const Eigen::Vector3d& rectangleBaseTangent,
                                                          const Eigen::Vector3d& rectangleHeightTangent,
                                                          const vector<double>& baseMeshCurvilinearCoordinates,
                                                          const vector<double>& heightMeshCurvilinearCoordinates,
                                                          const vector<unsigned int>& numberOfAddedVerticesForEachRectangle,
                                                          const GeometryUtilities& geometryUtilities,
                                                          IMeshDAO& mesh) const
  {
    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    unsigned int totalAddedHangNodes = std::accumulate(numberOfAddedVerticesForEachRectangle.begin(), numberOfAddedVerticesForEachRectangle.end(), 0);

    if ( totalAddedHangNodes == 0)
    {
      MeshUtilities::CreateRectangleMesh(rectangleOrigin,
                                         rectangleBaseTangent,
                                         rectangleHeightTangent,
                                         baseMeshCurvilinearCoordinates,
                                         heightMeshCurvilinearCoordinates,
                                         mesh);
    }
    else
    {

      const unsigned int numCell0Ds = numBasePoints * numHeightPoints
                                      + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                      + ceil(numBasePoints * 0.5) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[1]
                                      + (numHeightPoints / 2) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2]
                                      + (numBasePoints / 2) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[3];

      const unsigned int numCell1Ds = ceil((numHeightPoints * 0.5)) *
                                      (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1)
          + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
          * ( numberOfAddedVerticesForEachRectangle[1] + 1)
          + (numHeightPoints / 2) * (numBasePoints - 1) *
          (numberOfAddedVerticesForEachRectangle[2] + 1)
          + (numBasePoints / 2) * (numHeightPoints - 1) *
          (numberOfAddedVerticesForEachRectangle[3] + 1);

      const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

      mesh.InitializeDimension(2);

      mesh.Cell0DsInitialize(numCell0Ds);
      mesh.Cell1DsInitialize(numCell1Ds);
      mesh.Cell2DsInitialize(numCell2Ds);

      // create cell0Ds (# numBasePoints * numHeightPoints)
      unsigned int cell0DIndex = 0;
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b++)
        {
          const Eigen::Vector3d coordinate = rectangleOrigin +
                                             baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                             heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;
          const unsigned int marker = 1 * (b == 0 && h == 0) +
                                      2 * (b == (numBasePoints - 1) && h == 0) +
                                      4 * (b == 0 && h == (numHeightPoints - 1)) +
                                      3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                      5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                      7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                      8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                      6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

          mesh.Cell0DSetState(cell0DIndex, true);
          mesh.Cell0DInsertCoordinates(cell0DIndex,
                                       coordinate);

          mesh.Cell0DSetMarker(cell0DIndex, marker);
          cell0DIndex++;
        }
      }

      // create hanging nodes onto the base of rectangle
      // (# numberOfAddedNodesToFirstNEdges * (numBasePoints - 1) * ceil((numHeightPoints * 0.5)))
      for (unsigned int h = 0; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < (numBasePoints - 1); b++)
        {

          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[0],
              baseMeshCurvilinearCoordinates[b],
              baseMeshCurvilinearCoordinates[b+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               curvilinearPoints[s] * rectangleBaseTangent +
                                               heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;



            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 0; h < (numHeightPoints - 1); h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b = b + 2)
        {

          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[1],
              heightMeshCurvilinearCoordinates[h],
              heightMeshCurvilinearCoordinates[h+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                               curvilinearPoints[s] * rectangleHeightTangent;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 1; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < (numBasePoints - 1); b++)
        {
          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[2],
              baseMeshCurvilinearCoordinates[b],
              baseMeshCurvilinearCoordinates[b+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               curvilinearPoints[s] * rectangleBaseTangent +
                                               heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;



            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 0; h < (numHeightPoints-1); h++)
      {
        for (unsigned int b = 1; b < numBasePoints; b = b + 2)
        {
          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[3],
              heightMeshCurvilinearCoordinates[h],
              heightMeshCurvilinearCoordinates[h+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                               curvilinearPoints[s] * rectangleHeightTangent;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      // create cell1Ds
      unsigned int cell1DIndex = 0;

      // create horizontal cell1Ds on b % 2 == 0
      unsigned int cell0DIndexBaseHangsEvenB = numBasePoints * numHeightPoints - 1;
      for (unsigned int h = 0; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {
          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsEvenB + 1;

            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsEvenB++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + 1;

          const unsigned int marker = 5 * (h == 0) +
                                      7 * (h == (numHeightPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create horizontal cell1Ds on b % 2 == 1
      unsigned int cell0DIndexBaseHangsOddB = numBasePoints * numHeightPoints - 1
                                              + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                              + ceil(numBasePoints * 0.5) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[1];

      for (unsigned int h = 1; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsOddB + 1;

            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsOddB++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + 1;

          const unsigned int marker = 5 * (h == 0) +
                                      7 * (h == (numHeightPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create vertical cell1Ds on h % 2 == 0
      unsigned int cell0DIndexBaseHangsEvenH = numBasePoints * numHeightPoints - 1
                                               + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b = b + 2)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsEvenH + 1;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsEvenH++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

          const unsigned int marker = 8 * (b == 0) +
                                      6 * (b == (numBasePoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create vertical cell1Ds on h % 2 == 1
      unsigned int cell0DIndexBaseHangsOddH = numBasePoints * numHeightPoints - 1
                                              + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                              + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                              * numberOfAddedVerticesForEachRectangle[1]
                                              + (numHeightPoints / 2) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2];

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int b = 1; b < numBasePoints; b = b + 2)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsOddH + 1;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsOddH++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

          const unsigned int marker = 8 * (b == 0) +
                                      6 * (b == (numBasePoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create cell2Ds
      unsigned int cell2DIndex = 0;

      cell0DIndexBaseHangsEvenB = numBasePoints * numHeightPoints - 1;

      cell0DIndexBaseHangsEvenH = cell0DIndexBaseHangsEvenB
                                  + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];

      cell0DIndexBaseHangsOddB = cell0DIndexBaseHangsEvenH
                                 + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                 * (numberOfAddedVerticesForEachRectangle[1]);

      cell0DIndexBaseHangsOddH = cell0DIndexBaseHangsOddB
                                 + (numHeightPoints / 2) * (numBasePoints - 1) *
                                 (numberOfAddedVerticesForEachRectangle[2]);

      unsigned int cell1DHorizontalEvenB = 0;
      unsigned int cell1DHorizontalOddB = ceil((numHeightPoints * 0.5)) *
                                          (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1);
      unsigned int cell1DVerticalEvenH = cell1DHorizontalOddB
                                         + (numHeightPoints / 2) * (numBasePoints - 1) *
                                         (numberOfAddedVerticesForEachRectangle[2] + 1);
      unsigned int cell1DVerticalOddH = cell1DVerticalEvenH
                                        + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                        * ( numberOfAddedVerticesForEachRectangle[1] + 1);

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        if ( h != 0)
        {
          if (h % 2 == 0)
          {
            cell0DIndexBaseHangsEvenB -= (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];
            cell1DHorizontalEvenB -= (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1);
          }
          else
          {
            cell0DIndexBaseHangsOddB -= (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2];
            cell1DHorizontalOddB -= (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[2] + 1);
          }
        }

        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {

          if ( b != 0)
          {
            if (b % 2 == 0)
            {
              cell0DIndexBaseHangsEvenH -= numberOfAddedVerticesForEachRectangle[1];
              cell1DVerticalEvenH -= numberOfAddedVerticesForEachRectangle[1] + 1;
            }
            else
            {
              cell0DIndexBaseHangsOddH -= numberOfAddedVerticesForEachRectangle[3];
              cell1DVerticalOddH -= numberOfAddedVerticesForEachRectangle[3] + 1;
            }
          }

          const unsigned int cell0DIndex = b + h * numBasePoints;

          vector<unsigned int> cell2DVertices(4 + totalAddedHangNodes);
          vector<unsigned int> cell2DEdges(4 + totalAddedHangNodes);

          unsigned int countV = 0;
          unsigned int countE = 0;
          cell2DVertices[countV] = cell0DIndex;
          countV++;

          //unsigned int numberOfAddedNodesTolastEdges;

          if (h % 2 == 0)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenB + 1;
              cell0DIndexBaseHangsEvenB++;
              cell2DEdges[countE++] = cell1DHorizontalEvenB;
              cell1DHorizontalEvenB++;
            }
            cell2DEdges[countE++] = cell1DHorizontalEvenB;
            cell1DHorizontalEvenB++;
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddB + 1;
              cell0DIndexBaseHangsOddB++;
              cell2DEdges[countE++] = cell1DHorizontalOddB;
              cell1DHorizontalOddB++;
            }

            cell2DEdges[countE++] = cell1DHorizontalOddB;
            cell1DHorizontalOddB++;
          }

          cell2DVertices[countV] = cell0DIndex + 1;
          countV++;

          if(b % 2 == 0)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddH + 1;
              cell0DIndexBaseHangsOddH++;
              cell2DEdges[countE++] = cell1DVerticalOddH;
              cell1DVerticalOddH++;
            }
            cell2DEdges[countE++] = cell1DVerticalOddH;
            cell1DVerticalOddH++;
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenH + 1;
              cell0DIndexBaseHangsEvenH++;
              cell2DEdges[countE++] = cell1DVerticalEvenH;
              cell1DVerticalEvenH++;
            }
            cell2DEdges[countE++] = cell1DVerticalEvenH;
            cell1DVerticalEvenH++;
          }

          cell2DVertices[countV] = cell0DIndex + numBasePoints + 1;
          countV++;

          if (h % 2 == 1)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenB + (numberOfAddedVerticesForEachRectangle[0] - s);
              cell2DEdges[countE++] = cell1DHorizontalEvenB + (numberOfAddedVerticesForEachRectangle[0] - s);
            }

            cell2DEdges[countE++] = cell1DHorizontalEvenB;
            cell1DHorizontalEvenB += numberOfAddedVerticesForEachRectangle[0] + 1;
            cell0DIndexBaseHangsEvenB += numberOfAddedVerticesForEachRectangle[0];
          }
          else
          {

            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddB + (numberOfAddedVerticesForEachRectangle[2] - s);
              cell2DEdges[countE++] = cell1DHorizontalOddB + (numberOfAddedVerticesForEachRectangle[2] - s);
            }

            cell2DEdges[countE++] = cell1DHorizontalOddB;
            cell1DHorizontalOddB += numberOfAddedVerticesForEachRectangle[2] + 1;
            cell0DIndexBaseHangsOddB += numberOfAddedVerticesForEachRectangle[2];
          }

          cell2DVertices[countV] = cell0DIndex + numBasePoints;
          countV++;

          if(b % 2 == 1)
          {

            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddH + (numberOfAddedVerticesForEachRectangle[3] - s);
              cell2DEdges[countE++] = cell1DVerticalOddH + (numberOfAddedVerticesForEachRectangle[3] - s);
            }

            cell2DEdges[countE++] = cell1DVerticalOddH;
            cell1DVerticalOddH += numberOfAddedVerticesForEachRectangle[3] + 1;
            cell0DIndexBaseHangsOddH += numberOfAddedVerticesForEachRectangle[3];
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenH + (numberOfAddedVerticesForEachRectangle[1] - s);
              cell2DEdges[countE++] = cell1DVerticalEvenH + (numberOfAddedVerticesForEachRectangle[1] - s);
            }

            cell2DEdges[countE++] = cell1DVerticalEvenH;
            cell1DVerticalEvenH += numberOfAddedVerticesForEachRectangle[1] + 1;
            cell0DIndexBaseHangsEvenH += numberOfAddedVerticesForEachRectangle[1];
          }

          mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
          mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

          mesh.Cell2DSetState(cell2DIndex, true);

          mesh.Cell2DSetMarker(cell2DIndex, 0);

          cell2DIndex++;
        }
      }
    }
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::SplitCell2D(const unsigned int& cell2DIndex,
                                                       const std::vector<Eigen::MatrixXi> subCell2Ds,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numSubCells = subCell2Ds.size();
    unsigned int newCell2DsStartingIndex = mesh.Cell2DAppend(numSubCells);

    vector<unsigned int> newCell2DsIndex(numSubCells);

    mesh.Cell2DSetState(cell2DIndex, false);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell2DsIndex[c] = newCell2DsStartingIndex + c;

      const unsigned int& newCell2DIndex = newCell2DsIndex[c];
      mesh.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                     subCell2Ds[c]);

      mesh.Cell2DSetMarker(newCell2DIndex, mesh.Cell2DMarker(cell2DIndex));
      mesh.Cell2DSetState(newCell2DIndex, true);

      mesh.Cell2DInsertUpdatedCell2D(cell2DIndex, newCell2DIndex);

      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(newCell2DIndex); e++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(newCell2DIndex, e);

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
            continue;

          if (mesh.Cell1DNeighbourCell2D(cell1DIndex, n) == cell2DIndex)
            mesh.Cell1DInsertNeighbourCell2D(cell1DIndex,
                                             n,
                                             newCell2DIndex);
        }
      }
    }

    return newCell2DsIndex;
  }
  // ***************************************************************************
  unsigned int MeshUtilities::FindTrianglesCommonVertex(const std::vector<unsigned int>& trianglesIndex,
                                                        const IMeshDAO& mesh) const
  {
    Output::Assert(trianglesIndex.size() > 1);

    std::vector<unsigned int> firstTriangleVertices = mesh.Cell2DVertices(trianglesIndex[0]);
    Output::Assert(firstTriangleVertices.size() == 3);

    std::vector<unsigned int> secondTriangleVertices = mesh.Cell2DVertices(trianglesIndex[1]);
    Output::Assert(secondTriangleVertices.size() == 3);

    std::vector<unsigned int> intersection;

    std::sort(firstTriangleVertices.begin(),
              firstTriangleVertices.end());
    std::sort(secondTriangleVertices.begin(),
              secondTriangleVertices.end());

    std::set_intersection(firstTriangleVertices.begin(),
                          firstTriangleVertices.end(),
                          secondTriangleVertices.begin(),
                          secondTriangleVertices.end(),
                          back_inserter(intersection));

    Output::Assert(intersection.size() == 1);

    // check other triangles
    for (unsigned int t = 2; t < trianglesIndex.size(); t++)
    {
      std::vector<unsigned int> triangleVertices = mesh.Cell2DVertices(trianglesIndex[t]);
      Output::Assert(triangleVertices.size() == 3);
      Output::Assert(find(triangleVertices.begin(),
                          triangleVertices.end(),
                          intersection[0]) != triangleVertices.end());
    }

    return intersection[0];
  }
  // ***************************************************************************
  void MeshUtilities::CreateConcaveMeshFromTriangularMesh(const std::vector<std::vector<unsigned int> >& trianglesIndicsToAgglomerate,
                                                          IMeshDAO& mesh) const
  {

  }
  // ***************************************************************************
  void MeshUtilities::CreateTriangularMesh(const Eigen::MatrixXd& polygonVertices,
                                           const double& maxTriangleArea,
                                           IMeshDAO& mesh,
                                           const string& options) const
  {
    TriangleInterface triangleInterface;

    triangleInterface.CreateMesh(polygonVertices,
                                 maxTriangleArea,
                                 mesh,
                                 options);
  }
  // ***************************************************************************
  void MeshUtilities::ChangePolygonMeshMarkers(const Eigen::MatrixXd& polygonVertices,
                                               const vector<unsigned int>& cell0DMarkers,
                                               const vector<unsigned int>& cell1DMarkers,
                                               IMeshDAO& mesh) const
  {
    Output::Assert(mesh.Dimension() == 2);

    const unsigned int numVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      if (mesh.Cell0DMarker(v) == 0)
        continue;

      const unsigned int newMarker = (mesh.Cell0DMarker(v) <= numVertices) ?
                                       cell0DMarkers.at(mesh.Cell0DMarker(v) - 1) :
                                       cell1DMarkers.at(mesh.Cell0DMarker(v) - numVertices - 1);

      mesh.Cell0DSetMarker(v,
                           newMarker);
    }

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      if (mesh.Cell1DMarker(e) == 0)
        continue;

      const unsigned int newMarker = cell1DMarkers.at(mesh.Cell1DMarker(e) - numVertices - 1);

      mesh.Cell1DSetMarker(e,
                           newMarker);
    }
  }
  // ***************************************************************************
}
