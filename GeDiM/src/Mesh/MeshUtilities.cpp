#include "MeshUtilities.hpp"

#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"

using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  MeshUtilities::MeshUtilities()
  {
  }
  MeshUtilities::~MeshUtilities()
  {
  }
  // ***************************************************************************
  void MeshUtilities::ExtractActiveMesh(IMeshDAO& mesh,
                                        ExtractActiveMeshData& extractionData) const
  {
    // remove inactive Cell0Ds
    unsigned int numNewCell0Ds = 0;
    list<unsigned int> cell0DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
      if (!mesh.Cell0DIsActive(c))
      {
        cell0DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell0DToOldCell0D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell0Ds, c));
      extractionData.OldCell0DToNewCell0D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell0Ds));
      numNewCell0Ds++;
    }

    unsigned int removedCell0Ds = 0;
    for (const unsigned int& c : cell0DIdToRemove)
    {
      mesh.Cell0DRemove(c - removedCell0Ds);
      removedCell0Ds++;
    }

    // remove inactive Cell1D
    unsigned int numNewCell1Ds = 0;
    list<unsigned int> cell1DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
      if (!mesh.Cell1DIsActive(c))
      {
        cell1DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell1DToOldCell1D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell1Ds, c));
      extractionData.OldCell1DToNewCell1D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell1Ds));
      numNewCell1Ds++;
    }

    unsigned int removedCell1Ds = 0;
    for (const unsigned int& c : cell1DIdToRemove)
    {
      mesh.Cell1DRemove(c - removedCell1Ds);
      removedCell1Ds++;
    }

    // remove inactive Cell2Ds
    unsigned int numNewCell2Ds = 0;
    list<unsigned int> cell2DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
      {
        cell2DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell2DToOldCell2D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell2Ds, c));
      extractionData.OldCell2DToNewCell2D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell2Ds));
      numNewCell2Ds++;
    }

    unsigned int removedCell2Ds = 0;
    for (const unsigned int& c : cell2DIdToRemove)
    {
      mesh.Cell2DRemove(c - removedCell2Ds);
      removedCell2Ds++;
    }

    // remove inactive Cell3Ds
    unsigned int numNewCell3Ds = 0;
    list<unsigned int> cell3DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
      if (!mesh.Cell3DIsActive(c))
      {
        cell3DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell3DToOldCell3D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell3Ds, c));
      extractionData.OldCell3DToNewCell3D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell3Ds));
      numNewCell3Ds++;
    }

    unsigned int removedCell3Ds = 0;
    for (const unsigned int& c : cell3DIdToRemove)
    {
      mesh.Cell3DRemove(c - removedCell3Ds);
      removedCell3Ds++;
    }

    mesh.Compress();
  }
  // ***************************************************************************
  void MeshUtilities::FillMesh1D(const GeometryUtilities& geometryUtilities,
                                 const Vector3d& segmentOrigin,
                                 const Vector3d& segmentTangent,
                                 const vector<double>& coordinates,
                                 IMeshDAO& mesh) const
  {
    if (coordinates.size() == 0)
      return;

    mesh.InitializeDimension(1);

    const unsigned int& numCell0Ds = coordinates.size();
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
    {
      mesh.Cell0DSetId(c, c);
      mesh.Cell0DSetState(c, true);
      mesh.Cell0DInsertCoordinates(c, segmentOrigin + coordinates[c] * segmentTangent);
    }

    const unsigned int numCell1Ds = numCell0Ds - 1;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      mesh.Cell1DSetId(e, e);
      mesh.Cell1DInsertExtremes(e,
                                e,
                                e + 1);
      mesh.Cell1DSetState(e, true);
    }
  }
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
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetId(v, v);
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, cell0Ds.col(v));
    }

    // Create Cell1Ds
    Output::Assert(cell1Ds.rows() == 2);
    unsigned int numCell1Ds = cell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    for (int e = 0; e < cell1Ds.cols(); e++)
    {
      mesh.Cell1DSetId(e, e);
      mesh.Cell1DInsertExtremes(e,
                                cell1Ds(0, e),
                                cell1Ds(1, e));
      mesh.Cell1DSetState(e, true);
    }

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

      mesh.Cell2DSetId(f, f);
      mesh.Cell2DSetState(f, true);
    }
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh2D(const CheckMesh2DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 2);

    // check Cell0D are 2D
    if (configuration.Cell0D_CheckCoordinates2D)
      Output::Assert(geometryUtilities.PointsAre2D(convexMesh.Cell0DCoordinates()));

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

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        for (unsigned int v = 0; v < convexMesh.Cell2DNumberVertices(p); v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % convexMesh.Cell2DNumberVertices(p));
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
      mesh.Cell0DSetId(v, v);
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polygonVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numPolygonVertices;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DSetId(e, e);
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

    mesh.Cell2DSetId(0, 0);
    mesh.Cell2DSetState(0, true);
    mesh.Cell2DSetMarker(0, 0);

    // Create Cell1D neighbours
    for (int e = 0; e < numPolygonVertices; e++)
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
  MeshUtilities::MeshGeometricData MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                          const IMeshDAO& convexMesh) const
  {
    MeshGeometricData result;

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
  MeshUtilities::MeshGeometricData MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                          const IMeshDAO& mesh,
                                                                          const IMeshDAO& convexMesh,
                                                                          const vector<vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const
  {
    MeshGeometricData result;

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
    for (int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
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

        mesh.Cell0DSetId(cell0DIndex, cell0DIndex);
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

        mesh.Cell1DSetId(cell1DIndex, cell1DIndex);
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

        mesh.Cell1DSetId(cell1DIndex, cell1DIndex);
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

        mesh.Cell2DSetId(cell2DIndex, cell2DIndex);
        mesh.Cell2DSetState(cell2DIndex, true);

        mesh.Cell2DSetMarker(cell2DIndex, 0);

        cell2DIndex++;
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateTriangularMesh(const Eigen::MatrixXd& polygonVertices,
                                           const double& maxTriangleArea,
                                           IMeshDAO& mesh) const
  {
    TriangleInterface triangleInterface;

    triangleInterface.CreateMesh(polygonVertices,
                                 maxTriangleArea,
                                 mesh);
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
  void MeshUtilities::ExportMeshToVTU(const IMeshDAO& mesh,
                                      const string& exportFolder,
                                      const string& fileName) const
  {
    // Export Cell0Ds
    {
      Gedim::VTKUtilities vtpUtilities;
      for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
      {
        vector<double> id(1, mesh.Cell0DId(g));
        vector<double> marker(1, mesh.Cell0DMarker(g));

        vtpUtilities.AddPoint(mesh.Cell0DCoordinates(g),
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Marker",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(marker.size()),
                                  marker.data()
                                }
                              });
      }

      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell0Ds.vtu");
    }

    // Export Cell1Ds
    {
      Gedim::VTKUtilities vtpUtilities;
      for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
      {
        vector<double> id(1, mesh.Cell1DId(g));
        vector<double> marker(1, mesh.Cell1DMarker(g));

        vtpUtilities.AddSegment(mesh.Cell1DCoordinates(g),
                                {
                                  {
                                    "Id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(id.size()),
                                    id.data()
                                  },
                                  {
                                    "Marker",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(marker.size()),
                                    marker.data()
                                  }
                                });
      }

      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell1Ds.vtu");
    }

    // Export Cell2Ds
    {
      Gedim::VTKUtilities vtpUtilities;
      for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
      {
        vector<double> id(1, mesh.Cell2DId(g));
        vector<double> marker(1, mesh.Cell2DMarker(g));

        vtpUtilities.AddPolygon(mesh.Cell2DVerticesCoordinates(g),
                                {
                                  {
                                    "Id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(id.size()),
                                    id.data()
                                  },
                                  {
                                    "Marker",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(marker.size()),
                                    marker.data()
                                  }
                                });
      }

      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell2Ds.vtu");
    }
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron MeshUtilities::MeshCell3DToPolyhedron(const IMeshDAO& mesh,
                                                                      const unsigned int& cell3DIndex) const
  {
    GeometryUtilities::Polyhedron polyhedron;

    unordered_map<unsigned int, unsigned int> cell0DIndexToVertexIndex;
    unordered_map<unsigned int, unsigned int> cell1DIndexToEdgeIndex;

    polyhedron.Vertices = mesh.Cell3DVerticesCoordinates(cell3DIndex);
    polyhedron.Edges.resize(2, mesh.Cell3DNumberEdges(cell3DIndex));
    polyhedron.Faces.resize(mesh.Cell3DNumberFaces(cell3DIndex));

    for (unsigned int v = 0; v < polyhedron.Vertices.cols(); v++)
    {
      cell0DIndexToVertexIndex.insert(make_pair(mesh.Cell3DVertex(cell3DIndex,
                                                                  v),
                                                v));
    }

    for (unsigned int e = 0; e < polyhedron.Edges.cols(); e++)
    {
      const unsigned int cell1DIndex = mesh.Cell3DEdge(cell3DIndex,
                                                       e);
      cell1DIndexToEdgeIndex.insert(make_pair(cell1DIndex,
                                              e));

      polyhedron.Edges(0, e) = cell0DIndexToVertexIndex.at(mesh.Cell1DOrigin(cell1DIndex));
      polyhedron.Edges(1, e) = cell0DIndexToVertexIndex.at(mesh.Cell1DEnd(cell1DIndex));
    }

    for (unsigned int f = 0; f < polyhedron.Faces.size(); f++)
    {
      const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                       f);
      const unsigned int numFaceVertices = mesh.Cell2DNumberVertices(cell2DIndex);

      polyhedron.Faces[f].resize(2, numFaceVertices);
      for (unsigned int v = 0; v < numFaceVertices; v++)
      {
        polyhedron.Faces[f](0, v) = cell0DIndexToVertexIndex.at(mesh.Cell2DVertex(cell2DIndex,
                                                                                  v));
        polyhedron.Faces[f](1, v) = cell1DIndexToEdgeIndex.at(mesh.Cell2DEdge(cell2DIndex ,
                                                                              v));
      }
    }

    return polyhedron;
  }
  // ***************************************************************************
}
