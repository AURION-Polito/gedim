#include "MeshUtilities.hpp"
#include <iostream>
#include <fstream>

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
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cct],
                                                                                        convexCell2DTriangulationAreas[cct]);
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
}
