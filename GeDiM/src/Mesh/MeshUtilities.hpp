#ifndef __MeshUtilities_H
#define __MeshUtilities_H

#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"

namespace Gedim
{
  /// \brief MeshUtilities
  /// \copyright See top level LICENSE file for details.
  class MeshUtilities final
  {
    public:
      struct CheckMesh2DConfiguration final
      {
          bool Cell0D_CheckCoordinates2D = true;
          bool Cell0D_CheckDuplications = true;
          bool Cell1D_CheckDuplications = true;
          bool Cell1D_CheckNeighbours = true;
          bool Cell1D_CheckMeasure = true;
          bool Cell2D_CheckEdges = true;
          bool Cell2D_CheckDuplications = true;
          bool Cell2D_CheckConvexity = true;
          bool Cell2D_CheckMeasure = true;
      };

      struct CheckMesh3DConfiguration final
      {
          bool Cell0D_CheckDuplications = true;
          bool Cell1D_CheckDuplications = true;
          bool Cell1D_CheckMeasure = true;
          bool Cell2D_CheckEdges = true;
          bool Cell2D_CheckDuplications = true;
          bool Cell2D_CheckConvexity = true;
          bool Cell2D_CheckMeasure = true;
          bool Cell3D_CheckDuplications = true;
          bool Cell3D_CheckConvexity = true;
          bool Cell3D_CheckMeasure = true;
      };

      struct ExtractActiveMeshData final
      {
          std::unordered_map<unsigned int, unsigned int> OldCell0DToNewCell0D; ///< each pair is {old Cell0D index, new Cell0D index}
          std::unordered_map<unsigned int, unsigned int> OldCell1DToNewCell1D; ///< each pair is {old Cell1D index, new Cell1D index}
          std::unordered_map<unsigned int, unsigned int> OldCell2DToNewCell2D; ///< each pair is {old Cell2D index, new Cell2D index}
          std::unordered_map<unsigned int, unsigned int> OldCell3DToNewCell3D; ///< each pair is {old Cell3D index, new Cell3D index}
          std::unordered_map<unsigned int, unsigned int> NewCell0DToOldCell0D; ///< each pair is {new Cell0D index, old Cell0D index}
          std::unordered_map<unsigned int, unsigned int> NewCell1DToOldCell1D; ///< each pair is {new Cell1D index, old Cell1D index}
          std::unordered_map<unsigned int, unsigned int> NewCell2DToOldCell2D; ///< each pair is {new Cell2D index, old Cell2D index}
          std::unordered_map<unsigned int, unsigned int> NewCell3DToOldCell3D; ///< each pair is {new Cell3D index, old Cell3D index}
      };

      struct ComputeMesh2DCell1DsResult final
      {
          Eigen::MatrixXi Cell1Ds; /// Cell1Ds vertices, size 2 x Cell1DTotalNumber()
          std::vector<Eigen::MatrixXi> Cell2Ds; ///< Cell2Ds vertices and edges, size Cell2DTotalNumber()x2xCell2DNumberVertices()
      };

      struct MeshGeometricData1D final
      {
          std::vector<Eigen::MatrixXd> Cell1DsVertices; ///< cell1D vertices coordinates
          std::vector<Eigen::Vector3d> Cell1DsTangents; ///< cell1D tangents
          std::vector<double> Cell1DsLengths; ///< cell1D lengths
          std::vector<double> Cell1DsSquaredLengths; ///< cell1D squared lengths
          std::vector<Eigen::Vector3d> Cell1DsCentroids; ///< cell1D centroids
      };

      struct MeshGeometricData2D final
      {
          std::vector<Eigen::MatrixXd> Cell2DsVertices; ///< cell2D vertices coordinates
          std::vector<std::vector<Eigen::Matrix3d>> Cell2DsTriangulations; ///< cell2D triangulations
          std::vector<double> Cell2DsAreas; ///< cell2D areas
          std::vector<Eigen::Vector3d> Cell2DsCentroids; ///< cell2D centroids
          std::vector<double> Cell2DsDiameters; ///< cell2D diameters
          std::vector<std::vector<bool>> Cell2DsEdgeDirections; ///< cell2D edge directions
          std::vector<Eigen::VectorXd> Cell2DsEdgeLengths; ///< cell2D edge lengths
          std::vector<Eigen::MatrixXd> Cell2DsEdgeTangents; ///< cell2D edge tangents
          std::vector<Eigen::MatrixXd> Cell2DsEdgeNormals; ///< cell2D edge normals
      };

      struct MeshGeometricData3D final
      {
          std::vector<Eigen::MatrixXd> Cell3DsVertices;
          std::vector<Eigen::MatrixXi> Cell3DsEdges;
          std::vector<std::vector<Eigen::MatrixXi>> Cell3DsFaces;
          std::vector<double> Cell3DsVolumes;
          std::vector<double> Cell3DsDiameters;
          std::vector<Eigen::Vector3d> Cell3DsCentroids;
          std::vector<Eigen::MatrixXd> Cell3DsEdgeTangents;
          std::vector<std::vector<bool>> Cell3DsEdgeDirections;
          std::vector<std::vector<Eigen::MatrixXd>> Cell3DsTetrahedronPoints;
          std::vector<std::vector<Eigen::Vector3d>> Cell3DsFacesTranslations;
          std::vector<std::vector<Eigen::Matrix3d>> Cell3DsFacesRotationMatrices;
          std::vector<std::vector<Eigen::Vector3d>> Cell3DsFacesNormals;
          std::vector<std::vector<bool>> Cell3DsFacesNormalDirections;
          std::vector<std::vector<std::vector<bool>>> Cell3DsFacesEdgeDirections;

          std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFaces3DVertices; ///< faces vertices 3D coordinates
          std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFaces2DVertices; ///< faces vertices 2D coordinates
          std::vector<std::vector<std::vector<Eigen::Matrix3d>>> Cell3DsFaces2DTriangulations; ///< faces triangulations 2D
          std::vector<std::vector<double>> Cell3DsFacesAreas; ///< faces areas
          std::vector<std::vector<Eigen::Vector3d>> Cell3DsFaces2DCentroids; ///< faces centroids
          std::vector<std::vector<double>> Cell3DsFacesDiameters; ///< faces diameters
          std::vector<std::vector<Eigen::VectorXd>> Cell3DsFacesEdgeLengths; ///< faces edge lengths
          std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdge2DTangents; ///< faces edge tangents
          std::vector<std::vector<Eigen::MatrixXd>> Cell3DsFacesEdge2DNormals; ///< faces edge normals
      };

      struct VTPPolyhedron final
      {
          Eigen::MatrixXd Vertices; /// size 3xnumVertices
          std::vector<std::vector<unsigned int>> PolyhedronFaces; /// size numFaces x numFaceVertices
      };

      struct AgglomerateTrianglesResult final
      {
          std::vector<unsigned int> RemovedEdges;
          std::vector<unsigned int> VerticesIndex;
          std::vector<unsigned int> EdgesIndex;
      };

      struct AgglomerateMeshFromTriangularMeshResult final
      {
          struct ConcaveCell2D
          {
              unsigned int Cell2DIndex;
              std::vector<unsigned int> ConvexCell2DsIndex;
          };

          std::vector<ConcaveCell2D> ConcaveCell2Ds;
          std::vector<unsigned int> RemovedCell1Ds;
          std::vector<unsigned int> RemovedCell2Ds;
      };

      struct AgglomerationInformation final
      {
          std::vector<unsigned int> OriginalCell0DToAgglomeratedCell0Ds;
          std::vector<unsigned int> OriginalCell1DToAgglomeratedCell1Ds;
          std::vector<unsigned int> OriginalCell2DToAgglomeratedCell2Ds;
          std::vector<unsigned int> AgglomeratedCell0DToOriginalCell0Ds;
          std::vector<std::vector<unsigned int>> AgglomeratedCell1DToOriginalCell1Ds;
          std::vector<std::vector<unsigned int>> AgglomeratedCell2DToOriginalCell2Ds;
      };

    public:
      MeshUtilities() { };
      ~MeshUtilities() { };

      /// \brief Extract Active Cells from mesh
      /// \note the resulting mesh has no inactive elements
      void ExtractActiveMesh(IMeshDAO& mesh,
                             ExtractActiveMeshData& extractionData) const;

      /// \brief Fill Mesh 1D From segment Coordinates
      /// \param segmentOrigin the segment origin
      /// \param segmentTangent the segment tangent vector
      /// \param coordinates relative coordinates between [0.0, 1.0]
      /// \param mesh the resulting mesh
      void FillMesh1D(const GeometryUtilities& geometryUtilities,
                      const Eigen::Vector3d& segmentOrigin,
                      const Eigen::Vector3d& segmentTangent,
                      const std::vector<double>& coordinates,
                      IMeshDAO& mesh) const;

      /// \brief Fill a Mesh 2D with vertices, edges and polygons
      /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
      /// \param cell1Ds the origin and end as Eigen MatrixXd of cell1Ds, size 2xCell1DTotalNumber()
      /// \param cell2Ds the vertices and edges indices of the cell2Ds ordered counterclockwise, size Cell2DTotalNumber()x2xCell2DNumberVertices()
      void FillMesh2D(const Eigen::MatrixXd& cell0Ds,
                      const Eigen::MatrixXi& cell1Ds,
                      const std::vector<Eigen::MatrixXi>& cell2Ds,
                      IMeshDAO& mesh) const;

      /// \brief Compute edges in a Mesh 2D with vertices and polygons
      /// \param cell0Ds the coordinates as Eigen MatrixXd of cell0Ds, size 3xCell0DTotalNumber()
      /// \param cell2Ds the vertices indices of the cell2Ds ordered counterclockwise, size Cell2DTotalNumber()xCell2DNumberVertices()
      /// \return the Cell1Ds data
      ComputeMesh2DCell1DsResult ComputeMesh2DCell1Ds(const Eigen::MatrixXd& cell0Ds,
                                                      const std::vector<Eigen::VectorXi>& cell2Ds) const;

      /// \brief Check Mesh2D correctness
      /// \param geometryUtilities the geometry utilities
      /// \param convexMesh a convex 2D mesh
      void CheckMesh2D(const CheckMesh2DConfiguration& configuration,
                       const GeometryUtilities& geometryUtilities,
                       const IMeshDAO& convexMesh) const;

      /// \brief Check Mesh3D correctness
      /// \param geometryUtilities the geometry utilities
      /// \param convexMesh a convex 3D mesh
      void CheckMesh3D(const CheckMesh3DConfiguration& configuration,
                       const GeometryUtilities& geometryUtilities,
                       const IMeshDAO& convexMesh) const;

      /// \brief Create a Mesh 1D with a segment
      /// \param segmentVertices the segment coordinates, size 3x2
      /// \param vertexMarkers mesh markers of vertices, size 1xNumPolygonVertices()
      void Mesh1DFromSegment(const GeometryUtilities& geometryUtilities,
                             const Eigen::MatrixXd& segmentVertices,
                             const std::vector<unsigned int> vertexMarkers,
                             IMeshDAO& mesh) const;

      /// \brief Create a Mesh 2D with a polygon
      /// \param polygonVertices the polygon coordinates, size 3xNumPolygonVertices()
      /// \param vertexMarkers mesh markers of vertices, size 1xNumPolygonVertices()
      /// \param edgeMarkers mesh markers of edges, size 1xNumPolygonVertices()
      void Mesh2DFromPolygon(const Eigen::MatrixXd& polygonVertices,
                             const std::vector<unsigned int> vertexMarkers,
                             const std::vector<unsigned int> edgeMarkers,
                             IMeshDAO& mesh) const;

      /// \brief Create a Mesh 3D with a polyhedron
      /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
      /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
      /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
      /// \param vertexMarkers mesh markers of vertices, size 1xnumVertices
      /// \param edgeMarkers mesh markers of edges, size 1xnumEdges
      /// \param faceMarkers mesh markers of faces, size 1xnumFaces
      void Mesh3DFromPolyhedron(const Eigen::MatrixXd& polyhedronVertices,
                                const Eigen::MatrixXi& polyhedronEdges,
                                const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                const std::vector<unsigned int> vertexMarkers,
                                const std::vector<unsigned int> edgeMarkers,
                                const std::vector<unsigned int> faceMarkers,
                                IMeshDAO& mesh) const;

      /// \brief Set the marker on all the mesh 3D elements lyin on the plane
      /// \param geometryUtilities the geometry utilities
      /// \param planeNormal the plane normal
      /// \param planeOrigin the plane origin
      /// \param marker the marker
      /// \param mesh the mesh
      void SetMeshMarkersOnPlane(const GeometryUtilities& geometryUtilities,
                                 const Eigen::Vector3d& planeNormal,
                                 const Eigen::Vector3d& planeOrigin,
                                 const unsigned int& marker,
                                 IMeshDAO& mesh) const;

      /// \brief Extract the mesh Cell2D Roots
      /// \param mesh the mesh
      /// \return the root cell for each cell2D, size 1xCell2DTotalNumber()
      std::vector<unsigned int> MeshCell2DRoots(const IMeshDAO& mesh) const;

      /// \brief Fill Mesh1D Geometric Data given a mesh with convex mesh cells
      /// \param convexMesh the convex mesh
      /// \return the MeshGeometricData computed
      MeshGeometricData1D FillMesh1DGeometricData(const GeometryUtilities& geometryUtilities,
                                                  const IMeshDAO& convexMesh) const;

      /// \brief Fill Mesh2D Geometric Data given a mesh with convex mesh cells
      /// \param convexMesh the convex mesh
      /// \return the MeshGeometricData computed
      MeshGeometricData2D FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                  const IMeshDAO& convexMesh) const;

      /// \brief Fill Mesh2D Geometric Data starting given a mesh with non convex mesh cells and its convex sub-mesh cells
      /// \param mesh the mesh
      /// \param convexMesh the convex mesh cells of mesh
      /// \param meshCell2DToConvexCell2DIndices the collection of convex cell2Ds for each mesh cell2D
      /// \return the MeshGeometricData computed
      MeshGeometricData2D FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                  const IMeshDAO& mesh,
                                                  const IMeshDAO& convexMesh,
                                                  const std::vector<std::vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const;

      /// \brief Fill Mesh3D Geometric Data given a mesh with convex mesh cells
      /// \param convexMesh the convex mesh
      /// \return the MeshGeometricData computed
      MeshGeometricData3D FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                  const IMeshDAO& convexMesh) const;

      /// \brief Fill Mesh3D Geometric Data starting given a mesh with non convex mesh cells and its convex sub-mesh cells
      /// \param mesh the mesh
      /// \param convexMesh the convex mesh
      /// \param meshCell2DToConvexCell2DIndices the collection of convex cell2Ds for each mesh cell2D
      /// \param meshCell3DToConvexCell3DIndices the collection of convex cell3Ds for each mesh cell3D
      /// \return the MeshGeometricData computed
      MeshGeometricData3D FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                  const IMeshDAO& mesh,
                                                  const IMeshDAO& convexMesh,
                                                  const std::vector<std::vector<unsigned int>>& meshCell2DToConvexCell2DIndices,
                                                  const std::vector<std::vector<unsigned int>>& meshCell2DToConvexCell3DIndices,
                                                  const std::vector<std::vector<unsigned int>>& meshCell3DToConvexCell3DIndices) const;

      /// \brief Compute Cell1D Cell2DNeighbours with given mesh data
      /// \param mesh the resulting mesh
      void ComputeCell1DCell2DNeighbours(IMeshDAO& mesh) const;

      /// \brief Compute Cell2D Cell3DNeighbours with given mesh data
      /// \param mesh the resulting mesh
      void ComputeCell2DCell3DNeighbours(IMeshDAO& mesh) const;

      /// \brief Crete rectange Mesh on rectangle base x height
      /// \param rectangleOrigin the rectangle origin point
      /// \param rectangleBaseTangent the rectangle base tangent vector
      /// \param rectangleHeightTangent the rectangle height tangent vector
      /// \param baseMeshCurvilinearCoordinates the base mesh 1D curvilinear coordinates
      /// \param heightMeshCurvilinearCoordinates the height mesh 1D curvilinear coordinates
      /// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 * numVertices } for cell1Ds
      void CreateRectangleMesh(const Eigen::Vector3d& rectangleOrigin,
                               const Eigen::Vector3d& rectangleBaseTangent,
                               const Eigen::Vector3d& rectangleHeightTangent,
                               const std::vector<double>& baseMeshCurvilinearCoordinates,
                               const std::vector<double>& heightMeshCurvilinearCoordinates,
                               IMeshDAO& mesh) const;

      void CreateParallelepipedMesh(const Eigen::Vector3d& rectangleOrigin,
                                    const Eigen::Vector3d& rectangleLengthTangent,
                                    const Eigen::Vector3d& rectangleHeightTangent,
                                    const Eigen::Vector3d& rectangleWidthTangent,
                                    const std::vector<double>& lengthMeshCurvilinearCoordinates,
                                    const std::vector<double>& heightMeshCurvilinearCoordinates,
                                    const std::vector<double>& widthMeshCurvilinearCoordinates,
                                    IMeshDAO& mesh) const;

      /// \brief Create triangular mesh on 2D polygon
      /// \param polygonVertices the 2D polygon vertices, size 3xnumVertices
      /// \param maxTriangleArea the maximum triangular area
      /// \param options mesh options, see https://www.cs.cmu.edu/~quake/triangle.switch.html
      /// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 * numVertices } for cell1Ds
      /// \note use triangle library
      void CreateTriangularMesh(const Eigen::MatrixXd& polygonVertices,
                                const double& maxTriangleArea,
                                IMeshDAO& mesh,
                                const std::string& options = "-QDzpqnea") const;

      /// \brief Create tetrahedral mesh on 3D polyhedron
      /// \param polyhedronVertices the polyhedron vertices, size 3 x numVertices
      /// \param polyhedronEdges the polyhedron edges, size 2 x numEdges
      /// \param polyhedronFaces the polyhedron face vertices and edges, size numFaces x 2 x numVertices
      /// \param maxTetrahedronVolume the maximum tetrahedron area
      /// \param options mesh options, see https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#cmd-q
      /// \note markers on border are set as { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 * numVertices } for cell1Ds
      /// \note use tetgen library
      void CreateTetrahedralMesh(const Eigen::MatrixXd& polyhedronVertices,
                                 const Eigen::MatrixXi& polyhedronEdges,
                                 const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                 const double& maxTetrahedronVolume,
                                 IMeshDAO& mesh,
                                 const std::string& options = "Qpqfezna") const;

      /// \brief Import 3D mesh from OVM file
      void ImportOpenVolumeMesh(const std::string& ovmFilePath,
                                IMeshDAO& mesh,
                                std::vector<std::vector<bool>>& meshCell3DsFacesOrientation) const;

      /// \brief Export 3D mesh to OVM file
      void ExportMeshToOpenVolume(const IMeshDAO& mesh,
                                  const std::vector<std::vector<bool>>& meshCell3DsFacesOrientation,
                                  const std::string& ovmFilePath) const;

      /// \brief Change Polygon Mesh Markers from { 1, 2, 3, 4, ..., numVertices } for cell0Ds and { 5, 6, 7, 8, ..., 2 * numVertices } for cell1Ds to cell0DMarkers and cell1DMarkers
      /// \param polygonVertices the 2D polygon vertices, size 3xnumVertices
      /// \param cell0DMarkers the new cell0D markers, size 1xnumPolygonVertices
      /// \param cell1DMarkers the new cell1D markers, size 1xnumPolygonVertices
      /// \param mesh the mesh
      void ChangePolygonMeshMarkers(const Eigen::MatrixXd& polygonVertices,
                                    const std::vector<unsigned int>& cell0DMarkers,
                                    const std::vector<unsigned int>& cell1DMarkers,
                                    IMeshDAO& mesh) const;

      /// \brief Export Mesh To VTU
      /// \param mesh the mesh
      /// \param exportFolder the folder in which the mesh is exported
      void ExportMeshToVTU(const IMeshDAO& mesh,
                           const std::string& exportFolder,
                           const std::string& fileName) const;

      /// \brief Export Cell2D To VTU
      /// \param mesh the mesh
      /// \param cell2DIndex the cell2D index
      /// \param cell2DVertices the cell2D vertices
      /// \param cell2DTriangulations the cell2D triangulation
      /// \param cell2DArea the cell2D area
      /// \param cell2DCentroid the cell2D centroid
      /// \param exportFolder the folder in which to export
      void ExportCell2DToVTU(const IMeshDAO& mesh,
                             const unsigned int& cell2DIndex,
                             const Eigen::MatrixXd& cell2DVertices,
                             const std::vector<Eigen::Matrix3d>& cell2DTriangulations,
                             const double& cell2DArea,
                             const Eigen::Vector3d& cell2DCentroid,
                             const std::string& exportFolder) const;

      /// \brief Convert a mesh cell3D to a geometric polydheron
      /// \param mesh a mesh
      /// \param cell3DIndex the cell3D index
      /// \return polyhedron from mesh 3D cell
      GeometryUtilities::Polyhedron MeshCell3DToPolyhedron(const IMeshDAO& mesh,
                                                           const unsigned int& cell3DIndex) const;

      /// \brief Convert a mesh cell3D to a VTP polydheron
      /// \param mesh a mesh
      /// \param cell3DIndex the cell3D index
      /// \return VTP polyhedron from mesh 3D cell
      MeshUtilities::VTPPolyhedron MeshCell3DToVTPPolyhedron(const IMeshDAO& mesh,
                                                             const unsigned int& cell3DIndex) const;

      void CreateRectanglePlusHangingNodesMesh(const Eigen::Vector3d& rectangleOrigin,
                                               const Eigen::Vector3d& rectangleBaseTangent,
                                               const Eigen::Vector3d& rectangleHeightTangent,
                                               const std::vector<double>& baseMeshCurvilinearCoordinates,
                                               const std::vector<double>& heightMeshCurvilinearCoordinates,
                                               const std::vector<unsigned int>& numberOfAddedVerticesForEachRectangle,
                                               const GeometryUtilities& geometryUtilities,
                                               IMeshDAO& mesh) const;

      /// \brief Split cell2D into subcells
      /// \param cell1DIndex the index of Cell1D from 0 to Cell1DTotalNumber()
      /// \param subCell1Ds the list of sub-cells 1D mesh vertices indices, size 2 x numSubCells)
      /// \param mesh the mesh to update
      /// \return the list of new cell1Ds indices, from 0 to Cell1DTotalNumber()
      std::vector<unsigned int> SplitCell1D(const unsigned int& cell1DIndex,
                                            const Eigen::MatrixXi subCell1Ds,
                                            IMeshDAO& mesh) const;

      /// \brief Split cell2D into subcells
      /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
      /// \param subCell2Ds the list of sub-cells 2D mesh vertices and edges indices, size numSubCells x (2 x numVertices)
      /// \param mesh the mesh to update
      /// \return the list of new cell2Ds indices, from 0 to Cell2DTotalNumber()
      std::vector<unsigned int> SplitCell2D(const unsigned int& cell2DIndex,
                                            const std::vector<Eigen::MatrixXi> subCell2Ds,
                                            IMeshDAO& mesh) const;

      void CreateRandomlyDeformedQuadrilaterals(const GeometryUtilities& geometryUtilities,
                                                const Eigen::Vector3d& rectangleOrigin,
                                                const Eigen::Vector3d& rectangleBaseTangent,
                                                const Eigen::Vector3d& rectangleHeightTangent,
                                                const unsigned int& numQuadrilateralsBaseTangent,
                                                const unsigned int& numQuadrilateralsHeightTangent,
                                                const double& maxDeformingPercentageBase,
                                                const double& maxDeformingPercentageHeight,
                                                IMeshDAO& mesh) const;

      void CreateDistortedQuadrilaterals(const GeometryUtilities& geometryUtilities,
                                         const Eigen::Vector3d& rectangleOrigin,
                                         const Eigen::Vector3d& rectangleBaseTangent,
                                         const Eigen::Vector3d& rectangleHeightTangent,
                                         const unsigned int& numQuadrilateralsBaseTangent,
                                         const unsigned int& numQuadrilateralsHeightTangent,
                                         IMeshDAO& mesh) const;
      /// \brief Given a set of Cell2Ds find the common Cell0Ds
      /// \param cell2DsIndex the cell2Ds index
      /// \param mesh the mesh
      /// \return the Cell0D indices
      std::vector<unsigned int> FindCell2DsCommonVertices(const std::vector<unsigned int>& cell2DsIndex,
                                                          const IMeshDAO& mesh) const;

      /// \brief Given a set of Cell2Ds find the common Cell1Ds
      /// \param cell2DsIndex the cell2Ds index
      /// \param mesh the mesh
      /// \return the Cell1D indices
      std::vector<unsigned int> FindCell2DsCommonEdges(const std::vector<unsigned int>& cell2DsIndex,
                                                       const IMeshDAO& mesh) const;

      /// \brief Agglomerate Triangles with one vertex in common
      /// \param trianglesIndexToAgglomerate the cell2Ds triangular index in the mesh
      /// \param triangularMesh the triangular mesh
      /// \return the agglomearted polygon indices
      /// \note the triangular index shall be done counterclockwise
      AgglomerateTrianglesResult AgglomerateTriangles(const std::vector<unsigned int>& trianglesIndexToAgglomerate,
                                                      IMeshDAO& triangularMesh) const;

      AgglomerateMeshFromTriangularMeshResult AgglomerateMeshFromTriangularMesh(const std::vector<std::vector<unsigned int>>& trianglesIndicesToAgglomerate,
                                                                                IMeshDAO& triangularMesh) const;

      /// \brief Import Agglomeration mesh Information From file Csv
      /// \param geometryUtilities the geometry utilities
      /// \param originalMesh the original mesh
      /// \param agglomeratedMesh the agglomerated mesh
      /// \param fileName the csv file name
      /// \param separator the csv file separator
      /// \param originalCell0DToAgglomeratedCell0Ds original Cell0Ds to agglomerated Cell0Ds
      /// \param originalCell1DToAgglomeratedCell1Ds original Cell1Ds to agglomerated Cell1Ds
      /// \param originalCell2DToAgglomeratedCell2Ds original Cell2Ds to agglomerated Cell2Ds
      /// \param agglomeratedCell0DToOriginalCell0Ds agglomerated Cell0Ds to original Cell0Ds
      /// \param agglomeratedCell1DToOriginalCell1Ds agglomerated Cell1Ds to original Cell1Ds
      /// \param agglomeratedCell2DToOriginalCell2Ds agglomerated Cell2Ds to original Cell2Ds
      AgglomerationInformation ImportAgglomerationInformationFromCsv(const Gedim::GeometryUtilities geometryUtilities,
                                                                     const Gedim::IMeshDAO& originalMesh,
                                                                     const Gedim::IMeshDAO& agglomeratedMesh,
                                                                     const std::string& fileName,
                                                                     const char& separator) const;

      /// \brief Export mesh to csv file
      void ExportMeshToCsv(const IMeshDAO& mesh,
                           const char& separator,
                           const std::string& exportFolderPath) const;

      /// \brief Export 2D concave mesh to csv file
      void ExportConcaveMesh2DToCsv(const IMeshDAO& mesh,
                                    const std::vector<std::vector<unsigned int>>& convexCell2DsIndex,
                                    const char& separator,
                                    const std::string& exportFolderPath) const;
  };

}

#endif // __MeshUtilities_H
