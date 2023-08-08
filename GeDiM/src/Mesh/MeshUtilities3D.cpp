#include "MeshUtilities.hpp"

#include "TetgenInterface.hpp"
#include "VTKUtilities.hpp"
#include "MapTetrahedron.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::CheckMesh3D(const CheckMesh3DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 3);

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
        Output::Assert(convexMesh.Cell1DByExtremes(convexMesh.Cell1DOrigin(e1),
                                                   convexMesh.Cell1DEnd(e1)) ==
                       e1);
        Output::Assert(convexMesh.Cell1DByExtremes(convexMesh.Cell1DEnd(e1),
                                                   convexMesh.Cell1DOrigin(e1)) ==
                       convexMesh.Cell1DTotalNumber());

        for (unsigned int e2 = e1 + 1; e2 < convexMesh.Cell1DTotalNumber(); e2++)
        {
          Output::Assert(!(convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DEnd(e2)));
          Output::Assert(!(convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DEnd(e2)));
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
        const unsigned int cell2DNumEdges = convexMesh.Cell2DNumberEdges(p);
        for (unsigned int v = 0; v < cell2DNumEdges; v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % cell2DNumEdges);

          const unsigned int edgeFromVerticesOE = convexMesh.Cell2DFindEdgeByExtremes(p,
                                                                                      eO,
                                                                                      eE);
          const unsigned int edgeFromVerticesEO = convexMesh.Cell2DFindEdgeByExtremes(p,
                                                                                      eE,
                                                                                      eO);

          Output::Assert((edgeFromVerticesOE < cell2DNumEdges && edgeFromVerticesOE == v) ||
                         (edgeFromVerticesEO < cell2DNumEdges && edgeFromVerticesEO == v));
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
        const Eigen::MatrixXd cell2DVertices3D = convexMesh.Cell2DVerticesCoordinates(p);
        const Eigen::Vector3d cell2DNormal = geometryUtilities.PolygonNormal(cell2DVertices3D);
        const Eigen::Vector3d cell2DTranslation = geometryUtilities.PolygonTranslation(cell2DVertices3D);
        const Eigen::Matrix3d cell2DRotationMatrix = geometryUtilities.PolygonRotationMatrix(cell2DVertices3D,
                                                                                             cell2DNormal,
                                                                                             cell2DTranslation);
        const Eigen::MatrixXd cell2DVertices2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D,
                                                                                          cell2DRotationMatrix.transpose(),
                                                                                          cell2DTranslation);
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices2D);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices2D,
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
        const Eigen::MatrixXd cell2DVertices3D = convexMesh.Cell2DVerticesCoordinates(p);
        const Eigen::Vector3d cell2DNormal = geometryUtilities.PolygonNormal(cell2DVertices3D);
        const Eigen::Vector3d cell2DTranslation = geometryUtilities.PolygonTranslation(cell2DVertices3D);
        const Eigen::Matrix3d cell2DRotationMatrix = geometryUtilities.PolygonRotationMatrix(cell2DVertices3D,
                                                                                             cell2DNormal,
                                                                                             cell2DTranslation);
        const Eigen::MatrixXd cell2DVertices2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D,
                                                                                          cell2DRotationMatrix.transpose(),
                                                                                          cell2DTranslation);

        Output::Assert(geometryUtilities.IsValue2DPositive(
                         geometryUtilities.PolygonArea(cell2DVertices2D)));
      }
    }

    if (configuration.Cell3D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell3DTotalNumber(); p1++)
      {
        vector<unsigned int> cell3D1Vertices = convexMesh.Cell3DVertices(p1);
        sort(cell3D1Vertices.begin(), cell3D1Vertices.end());
        vector<unsigned int> cell3D1Edges = convexMesh.Cell3DEdges(p1);
        sort(cell3D1Edges.begin(), cell3D1Edges.end());
        vector<unsigned int> cell3D1Faces = convexMesh.Cell3DFaces(p1);
        sort(cell3D1Faces.begin(), cell3D1Faces.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell3DTotalNumber(); p2++)
        {
          vector<unsigned int> cell3D2Vertices = convexMesh.Cell3DVertices(p2);
          sort(cell3D2Vertices.begin(), cell3D2Vertices.end());
          vector<unsigned int> cell3D2Edges = convexMesh.Cell3DEdges(p2);
          sort(cell3D2Edges.begin(), cell3D2Edges.end());
          vector<unsigned int> cell3D2Faces = convexMesh.Cell3DFaces(p2);
          sort(cell3D2Faces.begin(), cell3D2Faces.end());

          Output::Assert(cell3D1Vertices.size() != cell3D2Vertices.size() || !equal(cell3D1Vertices.begin(),
                                                                                    cell3D1Vertices.end(),
                                                                                    cell3D2Vertices.begin()));
          Output::Assert(cell3D1Edges.size() != cell3D2Edges.size() || !equal(cell3D1Edges.begin(),
                                                                              cell3D1Edges.end(),
                                                                              cell3D2Edges.begin()));
          Output::Assert(cell3D1Faces.size() != cell3D2Faces.size() || !equal(cell3D1Faces.begin(),
                                                                              cell3D1Faces.end(),
                                                                              cell3D2Faces.begin()));
        }
      }
    }


    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);
      }
    }

    if (configuration.Cell3D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);

        const GeometryUtilities::PointPolyhedronPositionResult polyhedronBarycenterPosition = geometryUtilities.PointPolyhedronPosition(polyhedronBarycenter,
                                                                                                                                        polyhedron.Faces,
                                                                                                                                        polyhedronFace3DVertices,
                                                                                                                                        polyhedronFace2DVertices,
                                                                                                                                        polyhedronFaceNormals,
                                                                                                                                        polyhedronFaceNormalDirections,
                                                                                                                                        polyhedronFaceTranslations,
                                                                                                                                        polyhedronFaceRotationMatrices);

        Output::Assert(polyhedronBarycenterPosition.Type ==
                       GeometryUtilities::PointPolyhedronPositionResult::Types::Inside);

        Output::Assert(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                            polyhedronFace2DVertices,
                                                            polyhedronFaceBarycenters,
                                                            polyhedronFaceNormals,
                                                            polyhedronFaceNormalDirections,
                                                            polyhedronBarycenter));
      }
    }

    if (configuration.Cell3D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                      polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices,
                                                                                                                                                         polyhedronFaceTriangulations);

        Output::Assert(geometryUtilities.IsValue3DPositive(
                         geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                            polyhedronFaceNormals,
                                                            polyhedronFaceNormalDirections,
                                                            polyhedronFaceTranslations,
                                                            polyhedronFaceRotationMatrices)));
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::Mesh3DFromPolyhedron(const Eigen::MatrixXd& polyhedronVertices,
                                           const Eigen::MatrixXi& polyhedronEdges,
                                           const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                           const vector<unsigned int> vertexMarkers,
                                           const vector<unsigned int> edgeMarkers,
                                           const vector<unsigned int> faceMarkers,
                                           IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(3);

    Output::Assert(polyhedronVertices.rows() == 3 && polyhedronVertices.cols() > 3);
    const unsigned int numVertices = polyhedronVertices.cols();
    const unsigned int numEdges = polyhedronEdges.cols();
    const unsigned int numFaces = polyhedronFaces.size();

    Output::Assert(vertexMarkers.size() == numVertices);
    Output::Assert(edgeMarkers.size() == numEdges);
    Output::Assert(faceMarkers.size() == numFaces);

    // Create Cell0Ds
    const unsigned int& numCell0Ds = numVertices;
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polyhedronVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numEdges;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                polyhedronEdges(0, e),
                                polyhedronEdges(1, e));
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, edgeMarkers[e]);
    }


    // Create Cell2Ds
    const unsigned int& numCell2Ds = numFaces;
    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const unsigned int numCell2DVertices = polyhedronFaces.at(f).cols();
      mesh.Cell2DInitializeVertices(f, numCell2DVertices);
      mesh.Cell2DInitializeEdges(f, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
        mesh.Cell2DInsertVertex(f, v, polyhedronFaces.at(f)(0, v));
      for (unsigned int e = 0; e < numCell2DVertices; e++)
        mesh.Cell2DInsertEdge(f, e, polyhedronFaces.at(f)(1, e));

      mesh.Cell2DSetState(f, true);
      mesh.Cell2DSetMarker(f, faceMarkers[f]);
    }

    // Create Cell3Ds
    const unsigned int& numCell3Ds = 1;
    mesh.Cell3DsInitialize(numCell3Ds);

    mesh.Cell3DInitializeVertices(0, numVertices);
    mesh.Cell3DInitializeEdges(0, numEdges);
    mesh.Cell3DInitializeFaces(0, numFaces);

    for (unsigned int v = 0; v < numVertices; v++)
      mesh.Cell3DInsertVertex(0, v, v);
    for (unsigned int e = 0; e < numEdges; e++)
      mesh.Cell3DInsertEdge(0, e, e);
    for (unsigned int f = 0; f < numFaces; f++)
      mesh.Cell3DInsertFace(0, f, f);

    mesh.Cell3DSetState(0, true);
    mesh.Cell3DSetMarker(0, 0);
  }
  // ***************************************************************************
  void MeshUtilities::SetMeshMarkersOnPlane(const GeometryUtilities& geometryUtilities,
                                            const Eigen::Vector3d& planeNormal,
                                            const Eigen::Vector3d& planeOrigin,
                                            const unsigned int& marker,
                                            IMeshDAO& mesh) const
  {
    // set cell0Ds markers
    std::vector<bool> vertices_on_plane(mesh.Cell0DTotalNumber(),
                                        false);
    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      if (geometryUtilities.IsPointOnPlane(mesh.Cell0DCoordinates(v),
                                           planeNormal,
                                           planeOrigin))
      {
        vertices_on_plane[v] = true;
        mesh.Cell0DSetMarker(v,
                             marker);
      }
    }

    // set cell1Ds markers
    for (unsigned int s = 0; s < mesh.Cell1DTotalNumber(); s++)
    {
      const Eigen::VectorXi extremes = mesh.Cell1DExtremes(s);

      if (vertices_on_plane[extremes[0]] &&
          vertices_on_plane[extremes[1]])
      {
        mesh.Cell1DSetMarker(s,
                             marker);
      }
    }

    // set cell2Ds markers
    for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
    {
      const vector<unsigned int> extremes = mesh.Cell2DVertices(p);

      bool isOnPlane = true;
      for (unsigned int v = 0; v < extremes.size(); v++)
      {
        if (!vertices_on_plane[extremes[v]])
        {
          isOnPlane = false;
          break;
        }
      }

      if (isOnPlane)
      {
        mesh.Cell2DSetMarker(p,
                             marker);
      }
    }
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData3D MeshUtilities::FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData3D result;

    result.Cell3DsVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdges.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(convexMesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  convexMesh.Cell3DTotalNumber(); c++)
    {
      const GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                              c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsEdgeTangents[c] = geometryUtilities.PolyhedronEdgeTangents(result.Cell3DsVertices[c],
                                                                               result.Cell3DsEdges[c]);
      result.Cell3DsEdgeDirections[c].resize(convexMesh.Cell3DNumberEdges(c));
      for (unsigned int e = 0; e < convexMesh.Cell3DNumberEdges(c); e++)
      {
        const unsigned int meshOrigin = convexMesh.Cell3DVertex(c, polyhedron.Edges(0, e));
        const unsigned int meshEnd = convexMesh.Cell3DVertex(c, polyhedron.Edges(1, e));

        result.Cell3DsEdgeDirections[c][e] = (convexMesh.Cell3DFindEdgeByExtremes(c,
                                                                                  meshOrigin,
                                                                                  meshEnd) == e);
      }

      result.Cell3DsFaces3DVertices[c] = geometryUtilities.PolyhedronFaceVertices(result.Cell3DsVertices[c],
                                                                                  result.Cell3DsFaces[c]);
      result.Cell3DsFacesTranslations[c] = geometryUtilities.PolyhedronFaceTranslations(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesNormals[c] = geometryUtilities.PolyhedronFaceNormals(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesRotationMatrices[c] = geometryUtilities.PolyhedronFaceRotationMatrices(result.Cell3DsFaces3DVertices[c],
                                                                                                result.Cell3DsFacesNormals[c],
                                                                                                result.Cell3DsFacesTranslations[c]);



      /// TODO: qui secondo me questo non è corretto, occore utilizzare le informazioni della mesh
      const unsigned int numFaces = result.Cell3DsFaces[c].size();
      result.Cell3DsFacesEdgeDirections[c].resize(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
      {
        const unsigned int numFaceEdges = polyhedron.Faces[f].cols();
        const unsigned int cell2DIndex = convexMesh.Cell3DFace(c, f);

        result.Cell3DsFacesEdgeDirections[c][f].resize(numFaceEdges);
        for (unsigned int e = 0; e < numFaceEdges; e++)
        {
          const unsigned int faceEdgeOrigin = polyhedron.Faces[f](0, e);
          const unsigned int faceEdgeEnd = polyhedron.Faces[f](0, (e + 1) % numFaceEdges);

          const unsigned int meshOrigin = convexMesh.Cell3DVertex(c, faceEdgeOrigin);
          const unsigned int meshEnd = convexMesh.Cell3DVertex(c, faceEdgeEnd);

          result.Cell3DsFacesEdgeDirections[c][f][e] = (convexMesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                                                            meshOrigin,
                                                                                            meshEnd) == e);
        }
      }

      result.Cell3DsDiameters[c] = geometryUtilities.PolyhedronDiameter(result.Cell3DsVertices[c]);

      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    result.Cell3DsFaces3DVertices[c]);

      result.Cell3DsFaces2DVertices[c] = geometryUtilities.PolyhedronFaceRotatedVertices(result.Cell3DsFaces3DVertices[c],
                                                                                         result.Cell3DsFacesTranslations[c],
                                                                                         result.Cell3DsFacesRotationMatrices[c]);

      result.Cell3DsFaces2DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces2DVertices[c],
                                                                                                          polyhedronFaceTriangulations);


      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdge2DTangents[c].resize(numFaces);

      for(unsigned int f = 0; f < numFaces; f++)
      {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d>& convexCell2DTriangulationPoints = result.Cell3DsFaces2DTriangulations[c][f];
        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
        }

        result.Cell3DsFacesAreas[c][f] = convexCell2DTriangulationAreas.sum();
        result.Cell3DsFaces2DCentroids[c][f] = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                                 convexCell2DTriangulationAreas,
                                                                                 result.Cell3DsFacesAreas[c][f]);
        result.Cell3DsFacesDiameters[c][f] = geometryUtilities.PolygonDiameter(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdgeLengths[c][f] = geometryUtilities.PolygonEdgeLengths(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DNormals[c][f] = geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);
      }

      result.Cell3DsFacesNormalDirections[c] = geometryUtilities.PolyhedronFaceNormalDirections(result.Cell3DsFaces3DVertices[c],
                                                                                                geometryUtilities.PolyhedronBarycenter(result.Cell3DsVertices[c]),
                                                                                                result.Cell3DsFacesNormals[c]);

      result.Cell3DsVolumes[c] = geometryUtilities.PolyhedronVolume(result.Cell3DsFaces2DTriangulations[c],
                                                                    result.Cell3DsFacesNormals[c],
                                                                    result.Cell3DsFacesNormalDirections[c],
                                                                    result.Cell3DsFacesTranslations[c],
                                                                    result.Cell3DsFacesRotationMatrices[c]);

      result.Cell3DsCentroids[c] = geometryUtilities.PolyhedronCentroid(result.Cell3DsFaces2DTriangulations[c],
                                                                        result.Cell3DsFacesNormals[c],
                                                                        result.Cell3DsFacesNormalDirections[c],
                                                                        result.Cell3DsFacesTranslations[c],
                                                                        result.Cell3DsFacesRotationMatrices[c],
                                                                        result.Cell3DsVolumes[c]);

      vector<unsigned int> polyhedronTetrahedrons = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(result.Cell3DsVertices[c],
                                                                                                                 result.Cell3DsFaces[c],
                                                                                                                 polyhedronFaceTriangulations,
                                                                                                                 result.Cell3DsCentroids[c]);

      result.Cell3DsTetrahedronPoints[c] = geometryUtilities.ExtractTetrahedronPoints(result.Cell3DsVertices[c],
                                                                                      result.Cell3DsCentroids[c],
                                                                                      polyhedronTetrahedrons);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData3D MeshUtilities::FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const IMeshDAO& convexMesh,
                                                                            const std::vector<std::vector<unsigned int>>& meshCell3DToConvexCell3DIndices) const
  {
    MeshGeometricData3D result;

    result.Cell3DsVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdges.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(mesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  mesh.Cell3DTotalNumber(); c++)
    {
      const vector<unsigned int>& convexCell3DIndices = meshCell3DToConvexCell3DIndices[c];
      const unsigned int& numConvexCell3Ds = convexCell3DIndices.size();

      const GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh,
                                                                              c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsEdgeTangents[c] = geometryUtilities.PolyhedronEdgeTangents(result.Cell3DsVertices[c],
                                                                               result.Cell3DsEdges[c]);
      result.Cell3DsEdgeDirections[c].resize(mesh.Cell3DNumberEdges(c));
      for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(c); e++)
      {
        const unsigned int meshOrigin = mesh.Cell3DVertex(c, polyhedron.Edges(0, e));
        const unsigned int meshEnd = mesh.Cell3DVertex(c, polyhedron.Edges(1, e));

        result.Cell3DsEdgeDirections[c][e] = (mesh.Cell3DFindEdgeByExtremes(c,
                                                                            meshOrigin,
                                                                            meshEnd) == e);
      }

      result.Cell3DsFaces3DVertices[c] = geometryUtilities.PolyhedronFaceVertices(result.Cell3DsVertices[c],
                                                                                  result.Cell3DsFaces[c]);
      result.Cell3DsFacesTranslations[c] = geometryUtilities.PolyhedronFaceTranslations(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesNormals[c] = geometryUtilities.PolyhedronFaceNormals(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesRotationMatrices[c] = geometryUtilities.PolyhedronFaceRotationMatrices(result.Cell3DsFaces3DVertices[c],
                                                                                                result.Cell3DsFacesNormals[c],
                                                                                                result.Cell3DsFacesTranslations[c]);

      const unsigned int numFaces = result.Cell3DsFaces[c].size();

      result.Cell3DsFacesEdgeDirections[c].resize(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
      {
        const unsigned int numFaceEdges = polyhedron.Faces[f].cols();
        const unsigned int cell2DIndex = convexMesh.Cell3DFace(c, f);

        result.Cell3DsFacesEdgeDirections[c][f].resize(numFaceEdges);
        for (unsigned int e = 0; e < numFaceEdges; e++)
        {
          const unsigned int faceEdgeOrigin = polyhedron.Faces[f](0, e);
          const unsigned int faceEdgeEnd = polyhedron.Faces[f](0, (e + 1) % numFaceEdges);

          const unsigned int meshOrigin = mesh.Cell3DVertex(c, faceEdgeOrigin);
          const unsigned int meshEnd = mesh.Cell3DVertex(c, faceEdgeEnd);

          result.Cell3DsFacesEdgeDirections[c][f][e] = (mesh.Cell2DFindEdgeByExtremes(cell2DIndex,
                                                                                      meshOrigin,
                                                                                      meshEnd) == e);
        }
      }

      result.Cell3DsDiameters[c] = geometryUtilities.PolyhedronDiameter(result.Cell3DsVertices[c]);

      result.Cell3DsFaces2DVertices[c] = geometryUtilities.PolyhedronFaceRotatedVertices(result.Cell3DsFaces3DVertices[c],
                                                                                         result.Cell3DsFacesTranslations[c],
                                                                                         result.Cell3DsFacesRotationMatrices[c]);

      // fix orientation for concave cells
      for (unsigned int f = 0; f < numFaces; f++)
      {
        if (geometryUtilities.IsValue1DNegative(result.Cell3DsFacesRotationMatrices[c][f].determinant()))
          result.Cell3DsFaces2DVertices[c][f].block(0, 1, 3, result.Cell3DsFaces2DVertices[c][f].cols() - 1).rowwise().reverseInPlace();
      }

      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(polyhedron.Faces,
                                                                                                                                    result.Cell3DsFaces2DVertices[c]);

      result.Cell3DsFaces2DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces2DVertices[c],
                                                                                                          polyhedronFaceTriangulations);

      std::vector<Eigen::Vector3d> faceInternalPoints(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
        faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(result.Cell3DsFaces2DTriangulations[c][f][0]);


      result.Cell3DsFacesNormalDirections[c] = geometryUtilities.PolyhedronFaceNormalDirections(result.Cell3DsFaces3DVertices[c],
                                                                                                faceInternalPoints,
                                                                                                result.Cell3DsFaces2DVertices[c],
                                                                                                result.Cell3DsFacesNormals[c],
                                                                                                result.Cell3DsFacesTranslations[c],
                                                                                                result.Cell3DsFacesRotationMatrices[c]);
      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdge2DTangents[c].resize(numFaces);

      for(unsigned int f = 0; f < numFaces; f++)
      {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d>& convexCell2DTriangulationPoints = result.Cell3DsFaces2DTriangulations[c][f];
        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
        }

        result.Cell3DsFacesAreas[c][f] = convexCell2DTriangulationAreas.sum();
        result.Cell3DsFaces2DCentroids[c][f] = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                                 convexCell2DTriangulationAreas,
                                                                                 result.Cell3DsFacesAreas[c][f]);
        result.Cell3DsFacesDiameters[c][f] = geometryUtilities.PolygonDiameter(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdgeLengths[c][f] = geometryUtilities.PolygonEdgeLengths(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DNormals[c][f] = geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);
      }

      unsigned int cell3DTetrahedronsSize = 0;
      std::vector<std::vector<Eigen::MatrixXd>> convexCell3DTetrahedronsPoints(numConvexCell3Ds);
      Eigen::VectorXd convexCell3DsVolume(numConvexCell3Ds);
      Eigen::MatrixXd convexCell3DsCentroid(3, numConvexCell3Ds);

      for (unsigned int cc = 0; cc < numConvexCell3Ds; cc++)
      {
        const unsigned int convexCell3DIndex = convexCell3DIndices[cc];

        const GeometryUtilities::Polyhedron convexCell3DPolyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                                            convexCell3DIndex);

        const std::vector<Eigen::MatrixXd> convexCell3DFaces3DVertices = geometryUtilities.PolyhedronFaceVertices(convexCell3DPolyhedron.Vertices,
                                                                                                                  convexCell3DPolyhedron.Faces);
        const std::vector<Eigen::Vector3d> convexCell3DFacesTranslation = geometryUtilities.PolyhedronFaceTranslations(convexCell3DFaces3DVertices);
        const std::vector<Eigen::Vector3d> convexCell3DFacesNormal = geometryUtilities.PolyhedronFaceNormals(convexCell3DFaces3DVertices);
        const std::vector<Eigen::Matrix3d> convexCell3DFacesRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(convexCell3DFaces3DVertices,
                                                                                                                                convexCell3DFacesNormal,
                                                                                                                                convexCell3DFacesTranslation);

        const std::vector<Eigen::MatrixXd> convexCell3DFaces2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(convexCell3DFaces3DVertices,
                                                                                                                         convexCell3DFacesTranslation,
                                                                                                                         convexCell3DFacesRotationMatrices);

        const std::vector<std::vector<unsigned int>> convexCell3DFacesTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(convexCell3DPolyhedron.Faces,
                                                                                                                                                   convexCell3DFaces3DVertices);
        const std::vector<std::vector<Eigen::Matrix3d>> convexCell3DFaces2DTriangulations = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(convexCell3DFaces2DVertices,
                                                                                                                                                       convexCell3DFacesTriangulations);

        const std::vector<bool> convexCell3DFacesNormalDirection = geometryUtilities.PolyhedronFaceNormalDirections(convexCell3DFaces3DVertices,
                                                                                                                    geometryUtilities.PolyhedronBarycenter(convexCell3DPolyhedron.Vertices),
                                                                                                                    convexCell3DFacesNormal);

        convexCell3DsVolume[cc] = geometryUtilities.PolyhedronVolume(convexCell3DFaces2DTriangulations,
                                                                     convexCell3DFacesNormal,
                                                                     convexCell3DFacesNormalDirection,
                                                                     convexCell3DFacesTranslation,
                                                                     convexCell3DFacesRotationMatrices);

        convexCell3DsCentroid.col(cc) = geometryUtilities.PolyhedronCentroid(convexCell3DFaces2DTriangulations,
                                                                             convexCell3DFacesNormal,
                                                                             convexCell3DFacesNormalDirection,
                                                                             convexCell3DFacesTranslation,
                                                                             convexCell3DFacesRotationMatrices,
                                                                             convexCell3DsVolume[cc]);

        const vector<unsigned int> convexCell3DTetrahedrons = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(convexCell3DPolyhedron.Vertices,
                                                                                                                           convexCell3DPolyhedron.Faces,
                                                                                                                           convexCell3DFacesTriangulations,
                                                                                                                           convexCell3DsCentroid.col(cc));

        convexCell3DTetrahedronsPoints.at(cc) = geometryUtilities.ExtractTetrahedronPoints(convexCell3DPolyhedron.Vertices,
                                                                                           convexCell3DsCentroid.col(cc),
                                                                                           convexCell3DTetrahedrons);
        cell3DTetrahedronsSize += convexCell3DTetrahedronsPoints.at(cc).size();
      }

      result.Cell3DsVolumes[c] = convexCell3DsVolume.sum();
      result.Cell3DsCentroids[c] = convexCell3DsCentroid * convexCell3DsVolume / result.Cell3DsVolumes[c];

      result.Cell3DsTetrahedronPoints[c].resize(cell3DTetrahedronsSize);
      unsigned int tetraCounter = 0;
      for (unsigned int cc = 0; cc < numConvexCell3Ds; cc++)
      {
        for (unsigned int cct = 0; cct < convexCell3DTetrahedronsPoints[cc].size(); cct++)
          result.Cell3DsTetrahedronPoints[c][tetraCounter++] =
              convexCell3DTetrahedronsPoints[cc][cct];
      }
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell2DCell3DNeighbours(IMeshDAO& mesh) const
  {
    // Compute Cell2D neighbours starting from cell3Ds
    std::vector<std::list<unsigned int>> cell2DsNeighbours(mesh.Cell2DTotalNumber());
    for (unsigned int c3D = 0; c3D < mesh.Cell3DTotalNumber(); c3D++)
    {
      const unsigned int numCell3DFaces = mesh.Cell3DNumberFaces(c3D);
      for (unsigned int f = 0; f < numCell3DFaces; f++)
      {
        const unsigned int cell2D = mesh.Cell3DFace(c3D, f);
        cell2DsNeighbours[cell2D].push_back(c3D);
      }
    }

    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
      mesh.Cell2DInitializeNeighbourCell3Ds(c2D, cell2DsNeighbours[c2D].size());

      unsigned int n = 0;
      for (const auto& cell3DIndex : cell2DsNeighbours[c2D])
        mesh.Cell2DInsertNeighbourCell3D(c2D,
                                         n++,
                                         cell3DIndex);
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
  MeshUtilities::VTPPolyhedron MeshUtilities::MeshCell3DToVTPPolyhedron(const IMeshDAO& mesh,
                                                                        const unsigned int& cell3DIndex) const
  {
    VTPPolyhedron vtpPolyhedron;

    unordered_map<unsigned int, unsigned int> cell0DIndexToVertexIndex;

    vtpPolyhedron.Vertices = mesh.Cell3DVerticesCoordinates(cell3DIndex);
    vtpPolyhedron.PolyhedronFaces.resize(mesh.Cell3DNumberFaces(cell3DIndex));

    for (unsigned int v = 0; v < vtpPolyhedron.Vertices.cols(); v++)
    {
      cell0DIndexToVertexIndex.insert(make_pair(mesh.Cell3DVertex(cell3DIndex,
                                                                  v),
                                                v));
    }

    for (unsigned int f = 0; f < vtpPolyhedron.PolyhedronFaces.size(); f++)
    {
      const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                       f);
      const unsigned int numFaceVertices = mesh.Cell2DNumberVertices(cell2DIndex);

      vtpPolyhedron.PolyhedronFaces[f].resize(numFaceVertices);
      for (unsigned int v = 0; v < numFaceVertices; v++)
      {
        vtpPolyhedron.PolyhedronFaces[f][v] = cell0DIndexToVertexIndex.at(mesh.Cell2DVertex(cell2DIndex,
                                                                                            v));
      }
    }

    return vtpPolyhedron;
  }
  // ***************************************************************************
}
