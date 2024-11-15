#include "MeshUtilities.hpp"

#include "TetgenInterface.hpp"
#include "VTKUtilities.hpp"
#include "MapTetrahedron.hpp"
#include <numeric>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::FillMesh3D(const Eigen::MatrixXd& cell0Ds,
                                 const Eigen::MatrixXi& cell1Ds,
                                 const std::vector<Eigen::MatrixXi>& cell2Ds,
                                 const std::vector<Mesh3DPolyhedron>& cell3Ds,
                                 IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(3);

    // Create Cell0Ds
    if (cell0Ds.size() > 0)
    {
      Output::Assert(cell0Ds.rows() == 3);
      const unsigned int numCell0Ds = cell0Ds.cols();
      mesh.Cell0DsInitialize(numCell0Ds);
      mesh.Cell0DsInsertCoordinates(cell0Ds);
      for (unsigned int v = 0; v < numCell0Ds; v++)
        mesh.Cell0DSetState(v, true);
    }

    // Create Cell1Ds
    if (cell1Ds.size() > 0)
    {
      Output::Assert(cell1Ds.rows() == 2);
      unsigned int numCell1Ds = cell1Ds.cols();
      mesh.Cell1DsInitialize(numCell1Ds);
      mesh.Cell1DsInsertExtremes(cell1Ds);
      for (int e = 0; e < cell1Ds.cols(); e++)
        mesh.Cell1DSetState(e, true);
    }

    // Create Cell2Ds
    if (cell2Ds.size() > 0)
    {
      const unsigned int& numCell2Ds = cell2Ds.size();
      mesh.Cell2DsInitialize(numCell2Ds);
      std::vector<unsigned int> cell2DsNumVertices(numCell2Ds), cell2DsNumEdges(numCell2Ds);

      for (unsigned int f = 0; f < numCell2Ds; f++)
      {
        const unsigned int numVertices = cell2Ds[f].cols();
        cell2DsNumVertices[f] = numVertices;
        cell2DsNumEdges[f] = numVertices;
      }

      mesh.Cell2DsInitializeVertices(cell2DsNumVertices);
      mesh.Cell2DsInitializeEdges(cell2DsNumEdges);

      for (unsigned int f = 0; f < numCell2Ds; f++)
      {
        const MatrixXi& polygon = cell2Ds[f];
        Output::Assert(polygon.rows() == 2);
        const unsigned int& numVertices = polygon.cols();

        for (unsigned int v = 0; v < numVertices; v++)
          mesh.Cell2DInsertVertex(f, v, polygon(0, v));
        for (unsigned int e = 0; e < numVertices; e++)
          mesh.Cell2DInsertEdge(f, e, polygon(1, e));

        mesh.Cell2DSetState(f, true);
      }
    }

    // Create Cell3Ds
    if (cell3Ds.size() > 0)
    {
      const unsigned int& numCell3Ds = cell3Ds.size();
      mesh.Cell3DsInitialize(numCell3Ds);
      std::vector<unsigned int> cell3DsNumVertices(numCell3Ds);
      std::vector<unsigned int> cell3DsNumEdges(numCell3Ds);
      std::vector<unsigned int> cell3DsNumFaces(numCell3Ds);

      for (unsigned int c = 0; c < numCell3Ds; c++)
      {
        cell3DsNumVertices[c] = cell3Ds[c].VerticesIndex.size();
        cell3DsNumEdges[c] = cell3Ds[c].EdgesIndex.size();
        cell3DsNumFaces[c] = cell3Ds[c].FacesIndex.size();
      }

      mesh.Cell3DsInitializeVertices(cell3DsNumVertices);
      mesh.Cell3DsInitializeEdges(cell3DsNumEdges);
      mesh.Cell3DsInitializeFaces(cell3DsNumFaces);

      for (unsigned int c = 0; c < numCell3Ds; c++)
      {
        const Mesh3DPolyhedron& polyhedron = cell3Ds[c];

        for (unsigned int v = 0; v < polyhedron.VerticesIndex.size(); v++)
          mesh.Cell3DInsertVertex(c, v, polyhedron.VerticesIndex.at(v));
        for (unsigned int e = 0; e < polyhedron.EdgesIndex.size(); e++)
          mesh.Cell3DInsertEdge(c, e, polyhedron.EdgesIndex.at(e));
        for (unsigned int f = 0; f < polyhedron.FacesIndex.size(); f++)
          mesh.Cell3DInsertFace(c, f, polyhedron.FacesIndex.at(f));

        mesh.Cell3DSetState(c, true);
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh3D(const CheckMesh3DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& mesh) const
  {
    Output::Assert(mesh.Dimension() == 3);

    // check Cell0D duplications
    if (configuration.Cell0D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < mesh.Cell0DTotalNumber(); p1++)
      {
        if (!mesh.Cell0DIsActive(p1))
          continue;

        for (unsigned int p2 = p1 + 1; p2 < mesh.Cell0DTotalNumber(); p2++)
        {
          if (!mesh.Cell0DIsActive(p2))
            continue;

          if (geometryUtilities.PointsAreCoincident(mesh.Cell0DCoordinates(p1),
                                                    mesh.Cell0DCoordinates(p2)))
          {
            throw std::runtime_error("Cell0D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are coincident");
          }
        }
      }
    }

    if (configuration.Cell1D_CheckDuplications)
    {
      for (unsigned int e1 = 0; e1 < mesh.Cell1DTotalNumber(); e1++)
      {
        if (!mesh.Cell1DIsActive(e1))
          continue;

        if (mesh.Cell1DByExtremes(mesh.Cell1DOrigin(e1),
                                  mesh.Cell1DEnd(e1)) !=
            e1)
        {
          throw std::runtime_error("Cell1D " +
                                   std::to_string(e1) +
                                   " has wrong extremes");
        }

        if (mesh.Cell1DByExtremes(mesh.Cell1DEnd(e1),
                                  mesh.Cell1DOrigin(e1)) !=
            mesh.Cell1DTotalNumber())
        {
          throw std::runtime_error("Cell1D " +
                                   std::to_string(e1) +
                                   " has wrong extremes");
        }

        for (unsigned int e2 = e1 + 1; e2 < mesh.Cell1DTotalNumber(); e2++)
        {
          if (!mesh.Cell1DIsActive(e2))
            continue;

          if ((mesh.Cell1DOrigin(e1) == mesh.Cell1DOrigin(e2) &&
               mesh.Cell1DEnd(e1) == mesh.Cell1DEnd(e2)))
          {
            throw std::runtime_error("Cell1D " +
                                     std::to_string(e1) +
                                     " and " +
                                     std::to_string(e2) +
                                     " are duplicated");
          }
          if ((mesh.Cell1DEnd(e1) == mesh.Cell1DOrigin(e2) &&
               mesh.Cell1DOrigin(e1) == mesh.Cell1DEnd(e2)))
          {
            throw std::runtime_error("Cell1D " +
                                     std::to_string(e1) +
                                     " and " +
                                     std::to_string(e2) +
                                     " are duplicated");
          }
        }
      }
    }

    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
      {
        if (!mesh.Cell1DIsActive(e))
          continue;

        if (geometryUtilities.IsValueZero(
              geometryUtilities.SegmentLength(mesh.Cell1DOriginCoordinates(e),
                                              mesh.Cell1DEndCoordinates(e)),
              geometryUtilities.Tolerance1D()))
        {
          throw std::runtime_error("Cell1D " +
                                   std::to_string(e) +
                                   " is zero");
        }
      }
    }

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
      {
        if (!mesh.Cell2DIsActive(p))
          continue;

        const unsigned int cell2DNumEdges = mesh.Cell2DNumberEdges(p);

        for (unsigned int v = 0; v < cell2DNumEdges; v++)
        {
          const unsigned int eO = mesh.Cell2DVertex(p, v);
          const unsigned int eE = mesh.Cell2DVertex(p, (v + 1) % cell2DNumEdges);

          const unsigned int edgeFromVerticesOE = mesh.Cell2DFindEdgeByExtremes(p,
                                                                                eO,
                                                                                eE);
          const unsigned int edgeFromVerticesEO = mesh.Cell2DFindEdgeByExtremes(p,
                                                                                eE,
                                                                                eO);

          if (!((edgeFromVerticesOE < cell2DNumEdges && edgeFromVerticesOE == v) ||
                (edgeFromVerticesEO < cell2DNumEdges && edgeFromVerticesEO == v)))
          {
            throw std::runtime_error("Cell2D " +
                                     std::to_string(p) +
                                     " has wrong edge " +
                                     std::to_string(v));
          }
        }
      }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < mesh.Cell2DTotalNumber(); p1++)
      {
        if (!mesh.Cell2DIsActive(p1))
          continue;

        vector<unsigned int> cell2D1Vertices = mesh.Cell2DVertices(p1);
        sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
        vector<unsigned int> cell2D1Edges = mesh.Cell2DEdges(p1);
        sort(cell2D1Edges.begin(), cell2D1Edges.end());

        for (unsigned int p2 = p1 + 1; p2 < mesh.Cell2DTotalNumber(); p2++)
        {
          if (!mesh.Cell2DIsActive(p2))
            continue;

          vector<unsigned int> cell2D2Vertices = mesh.Cell2DVertices(p2);
          sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
          vector<unsigned int> cell2D2Edges = mesh.Cell2DEdges(p2);
          sort(cell2D2Edges.begin(), cell2D2Edges.end());

          if (!(cell2D1Vertices.size() != cell2D2Vertices.size() ||
                !equal(cell2D1Vertices.begin(),
                       cell2D1Vertices.end(),
                       cell2D2Vertices.begin())))
          {
            throw std::runtime_error("Cell2D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are duplicated");
          }

          if (!(cell2D1Edges.size() != cell2D2Edges.size() ||
                !equal(cell2D1Edges.begin(),
                       cell2D1Edges.end(),
                       cell2D2Edges.begin())))
          {
            throw std::runtime_error("Cell2D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are duplicated");
          }
        }
      }
    }

    if (configuration.Cell2D_CheckConvexity)
    {
      for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
      {
        if (!mesh.Cell2DIsActive(p))
          continue;

        const Eigen::MatrixXd cell2DVertices3D = mesh.Cell2DVerticesCoordinates(p);
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

        if (!geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices,
                                               convexHullVertices))
        {
          throw std::runtime_error("Cell2D " +
                                   std::to_string(p) +
                                   " is not convex ");
        }

        if (geometryUtilities.PolygonOrientation(convexHull) ==
            Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
        {
          throw std::runtime_error("Cell2D " +
                                   std::to_string(p) +
                                   " is not CounterClockwise");
        }
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
      {
        if (!mesh.Cell2DIsActive(p))
          continue;

        const Eigen::MatrixXd cell2DVertices3D = mesh.Cell2DVerticesCoordinates(p);
        const Eigen::Vector3d cell2DNormal = geometryUtilities.PolygonNormal(cell2DVertices3D);
        const Eigen::Vector3d cell2DTranslation = geometryUtilities.PolygonTranslation(cell2DVertices3D);
        const Eigen::Matrix3d cell2DRotationMatrix = geometryUtilities.PolygonRotationMatrix(cell2DVertices3D,
                                                                                             cell2DNormal,
                                                                                             cell2DTranslation);
        const Eigen::MatrixXd cell2DVertices2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D,
                                                                                          cell2DRotationMatrix.transpose(),
                                                                                          cell2DTranslation);

        if (geometryUtilities.IsValueZero(
              geometryUtilities.PolygonArea(cell2DVertices2D),
              geometryUtilities.Tolerance2D()))
        {
          throw std::runtime_error("Cell2D " +
                                   std::to_string(p) +
                                   " is zero");
        }
      }
    }

    if (configuration.Cell3D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < mesh.Cell3DTotalNumber(); p1++)
      {
        if (!mesh.Cell3DIsActive(p1))
          continue;

        vector<unsigned int> cell3D1Vertices = mesh.Cell3DVertices(p1);
        sort(cell3D1Vertices.begin(), cell3D1Vertices.end());
        vector<unsigned int> cell3D1Edges = mesh.Cell3DEdges(p1);
        sort(cell3D1Edges.begin(), cell3D1Edges.end());
        vector<unsigned int> cell3D1Faces = mesh.Cell3DFaces(p1);
        sort(cell3D1Faces.begin(), cell3D1Faces.end());

        for (unsigned int p2 = p1 + 1; p2 < mesh.Cell3DTotalNumber(); p2++)
        {
          if (!mesh.Cell3DIsActive(p2))
            continue;

          vector<unsigned int> cell3D2Vertices = mesh.Cell3DVertices(p2);
          sort(cell3D2Vertices.begin(), cell3D2Vertices.end());
          vector<unsigned int> cell3D2Edges = mesh.Cell3DEdges(p2);
          sort(cell3D2Edges.begin(), cell3D2Edges.end());
          vector<unsigned int> cell3D2Faces = mesh.Cell3DFaces(p2);
          sort(cell3D2Faces.begin(), cell3D2Faces.end());

          if (!(cell3D1Vertices.size() != cell3D2Vertices.size() ||
                !equal(cell3D1Vertices.begin(),
                       cell3D1Vertices.end(),
                       cell3D2Vertices.begin())))
          {
            throw std::runtime_error("Cell3D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are duplicated");
          }

          if (!(cell3D1Edges.size() != cell3D2Edges.size() ||
                !equal(cell3D1Edges.begin(),
                       cell3D1Edges.end(),
                       cell3D2Edges.begin())))
          {
            throw std::runtime_error("Cell3D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are duplicated");
          }

          if (!(cell3D1Faces.size() != cell3D2Faces.size() ||
                !equal(cell3D1Faces.begin(),
                       cell3D1Faces.end(),
                       cell3D2Faces.begin())))
          {
            throw std::runtime_error("Cell3D " +
                                     std::to_string(p1) +
                                     " and " +
                                     std::to_string(p2) +
                                     " are duplicated");
          }
        }
      }
    }

    if (configuration.Cell3D_CheckConvexity)
    {
      for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
      {
        if (!mesh.Cell3DIsActive(p))
          continue;

        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh,
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

        if (!geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                  polyhedronFace2DVertices,
                                                  polyhedronFaceBarycenters,
                                                  polyhedronFaceNormals,
                                                  polyhedronFaceNormalDirections,
                                                  polyhedronBarycenter))
        {
          throw std::runtime_error("Cell3D " +
                                   std::to_string(p) +
                                   " is not convex");
        }
      }
    }

    if (configuration.Cell3D_CheckMeasure)
    {
      for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
      {
        if (!mesh.Cell3DIsActive(p))
          continue;

        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh,
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

        if (geometryUtilities.IsValueZero(
              geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                   polyhedronFaceNormals,
                                                                   polyhedronFaceNormalDirections,
                                                                   polyhedronFaceTranslations,
                                                                   polyhedronFaceRotationMatrices),
              geometryUtilities.Tolerance3D()))
        {
          throw std::runtime_error("Cell3D " +
                                   std::to_string(p) +
                                   " is zero");
        }
      }
    }
  }
  // ***************************************************************************
  MeshUtilities::ComputeMesh3DAlignedCell1DsResult MeshUtilities::ComputeMesh3DAlignedCell1Ds(const std::vector<std::vector<std::vector<unsigned int>>>& cell3DsAlignedEdgesVertices,
                                                                                              const std::vector<std::vector<std::vector<unsigned int>>>& cell3DsAlignedEdgesEdges,
                                                                                              const IMeshDAO& mesh) const
  {
    ComputeMesh3DAlignedCell1DsResult result;

    struct AlignedCell1D final
    {
        unsigned int Index;
        std::vector<unsigned int> Cell0Ds;
        std::vector<unsigned int> Cell1Ds;
        std::list<unsigned int> Cell3Ds;
    };

    unsigned int alignedCell1DIndex = 0;
    std::vector<std::unordered_map<unsigned int, AlignedCell1D>> alignedCell1DsByOrigin(mesh.Cell0DTotalNumber());
    std::vector<std::list<std::pair<unsigned int, unsigned int>>> cell0DsAlignedCell1DsIndex(mesh.Cell0DTotalNumber());
    std::vector<std::list<std::pair<unsigned int, unsigned int>>> cell1DsAlignedCell1DsIndex(mesh.Cell1DTotalNumber());
    result.Cell3DsAlignedCell1DsIndex.resize(mesh.Cell3DTotalNumber());

    for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
    {
      if (!mesh.Cell3DIsActive(cell3DIndex))
        continue;

      const std::vector<std::vector<unsigned int>>& cell3DAlignedEdgesVertices = cell3DsAlignedEdgesVertices.at(cell3DIndex);
      const std::vector<std::vector<unsigned int>>& cell3DAlignedEdgesEdges = cell3DsAlignedEdgesEdges.at(cell3DIndex);
      result.Cell3DsAlignedCell1DsIndex[cell3DIndex].resize(2, cell3DAlignedEdgesVertices.size());

      for (unsigned int ae = 0; ae < cell3DAlignedEdgesVertices.size(); ae++)
      {
        const std::vector<unsigned int>& cell3DAlignedEdgeVertices = cell3DAlignedEdgesVertices.at(ae);
        const std::vector<unsigned int>& cell3DAlignedEdgeEdges = cell3DAlignedEdgesEdges.at(ae);

        const unsigned int alignedEdgeNumVertices = cell3DAlignedEdgeVertices.size();
        const bool direction = mesh.Cell3DVertex(cell3DIndex,
                                                 cell3DAlignedEdgeVertices.at(0)) <
                               mesh.Cell3DVertex(cell3DIndex,
                                                 cell3DAlignedEdgeVertices.at(alignedEdgeNumVertices - 1));

        const unsigned int alignedEdgeOrigin = direction ?
                                                 cell3DAlignedEdgeVertices.at(0) :
                                                 cell3DAlignedEdgeVertices.at(alignedEdgeNumVertices - 1);
        const unsigned int alignedEdgeEnd = direction ?
                                              cell3DAlignedEdgeVertices.at(alignedEdgeNumVertices - 1) :
                                              cell3DAlignedEdgeVertices.at(0);
        const unsigned int cell0DOrigin = mesh.Cell3DVertex(cell3DIndex,
                                                            alignedEdgeOrigin);
        const unsigned int cell0DEnd = mesh.Cell3DVertex(cell3DIndex,
                                                         alignedEdgeEnd);

        const auto alignedEdge = alignedCell1DsByOrigin[cell0DOrigin].find(cell0DEnd);
        if (alignedEdge ==
            alignedCell1DsByOrigin[cell0DOrigin].end())
        {
          std::vector<unsigned int> alignedCell1D_Cell0Ds(alignedEdgeNumVertices);
          std::vector<unsigned int> alignedCell1D_Cell1Ds(alignedEdgeNumVertices - 1);

          for (unsigned int ae_v = 0; ae_v < alignedEdgeNumVertices - 1; ae_v++)
          {
            const unsigned int direction_vertex_index = direction ? ae_v :
                                                                    (alignedEdgeNumVertices - 1) - ae_v;
            const unsigned int direction_edge_index = direction ? ae_v :
                                                                  (alignedEdgeNumVertices - 2) - ae_v;

            alignedCell1D_Cell0Ds[ae_v] =  mesh.Cell3DVertex(cell3DIndex,
                                                             cell3DAlignedEdgeVertices.at(direction_vertex_index));
            alignedCell1D_Cell1Ds[ae_v] =  mesh.Cell3DEdge(cell3DIndex,
                                                           cell3DAlignedEdgeEdges.at(direction_edge_index));

            cell0DsAlignedCell1DsIndex[alignedCell1D_Cell0Ds[ae_v]].push_back(std::make_pair(alignedCell1DIndex, ae_v));
            cell1DsAlignedCell1DsIndex[alignedCell1D_Cell1Ds[ae_v]].push_back(std::make_pair(alignedCell1DIndex, ae_v));
          }
          alignedCell1D_Cell0Ds[alignedEdgeNumVertices - 1] = cell0DEnd;
          cell0DsAlignedCell1DsIndex[cell0DEnd].push_back(std::make_pair(alignedCell1DIndex, alignedEdgeNumVertices - 1));

          result.Cell3DsAlignedCell1DsIndex[cell3DIndex](0, ae) = alignedCell1DIndex;
          result.Cell3DsAlignedCell1DsIndex[cell3DIndex](1, ae) = 0;
          alignedCell1DsByOrigin[cell0DOrigin].insert(std::pair<unsigned int,
                                                      AlignedCell1D>(
                                                        cell0DEnd,
                                                        {
                                                          alignedCell1DIndex++,
                                                          alignedCell1D_Cell0Ds,
                                                          alignedCell1D_Cell1Ds,
                                                          { cell3DIndex }
                                                        }));
        }
        else
        {
          result.Cell3DsAlignedCell1DsIndex[cell3DIndex](0, ae) = alignedEdge->second.Index;
          result.Cell3DsAlignedCell1DsIndex[cell3DIndex](1, ae) = alignedEdge->second.Cell3Ds.size();
          alignedEdge->second.Cell3Ds.push_back(cell3DIndex);
        }
      }
    }

    result.Cell0DsAlignedCell1DsIndex.resize(mesh.Cell0DTotalNumber());
    for (unsigned int c0D = 0; c0D < mesh.Cell0DTotalNumber(); c0D++)
    {
      result.Cell0DsAlignedCell1DsIndex[c0D].resize(2, cell0DsAlignedCell1DsIndex[c0D].size());
      unsigned int col = 0;
      for (const auto& pair : cell0DsAlignedCell1DsIndex[c0D])
        result.Cell0DsAlignedCell1DsIndex[c0D].col(col++)<< pair.first, pair.second;
    }

    result.Cell1DsAlignedCell1DsIndex.resize(mesh.Cell1DTotalNumber());
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
    {
      result.Cell1DsAlignedCell1DsIndex[c1D].resize(2, cell1DsAlignedCell1DsIndex[c1D].size());
      unsigned int col = 0;
      for (const auto& pair : cell1DsAlignedCell1DsIndex[c1D])
        result.Cell1DsAlignedCell1DsIndex[c1D].col(col++)<< pair.first, pair.second;
    }

    result.AlignedCell1Ds.resize(2, alignedCell1DIndex);
    result.AlignedCell1Ds_SubCell0Ds.resize(alignedCell1DIndex);
    result.AlignedCell1Ds_SubCell1Ds.resize(alignedCell1DIndex);
    result.AlignedCell1Ds_Cell3Ds.resize(alignedCell1DIndex);

    for (unsigned int cell0DIndex = 0; cell0DIndex < mesh.Cell0DTotalNumber(); cell0DIndex++)
    {
      for (const auto& alignedEdge : alignedCell1DsByOrigin[cell0DIndex])
      {
        result.AlignedCell1Ds_SubCell0Ds[alignedEdge.second.Index] = alignedEdge.second.Cell0Ds;
        result.AlignedCell1Ds_SubCell1Ds[alignedEdge.second.Index] = alignedEdge.second.Cell1Ds;
        result.AlignedCell1Ds.col(alignedEdge.second.Index)<< alignedEdge.second.Cell0Ds.at(0), alignedEdge.second.Cell0Ds.at(alignedEdge.second.Cell0Ds.size() - 1);
        result.AlignedCell1Ds_Cell3Ds[alignedEdge.second.Index] = std::vector<unsigned int>(alignedEdge.second.Cell3Ds.begin(),
                                                                                            alignedEdge.second.Cell3Ds.end());
      }
    }

    //    std::cout.precision(16);
    //    std::cout<< result.AlignedCell1Ds<< std::endl;
    //    std::cout<< result.AlignedCell1Ds_SubCell0Ds<< std::endl;
    //    std::cout<< result.AlignedCell1Ds_SubCell1Ds<< std::endl;
    //    std::cout<< result.AlignedCell1Ds_Cell3Ds<< std::endl;
    //    std::cout<< result.Cell0DsAlignedCell1DsIndex<< std::endl;
    //    std::cout<< result.Cell1DsAlignedCell1DsIndex<< std::endl;
    //    std::cout<< result.Cell3DsAlignedCell1DsIndex<< std::endl;

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::CheckMeshGeometricData3D(const CheckMeshGeometricData3DConfiguration& configuration,
                                               const GeometryUtilities& geometryUtilities,
                                               const IMeshDAO& mesh,
                                               const MeshGeometricData3D& geometricData) const
  {
    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
          const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex, f);
          for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e++)
            Output::Assert(geometryUtilities.IsValuePositive(geometricData.Cell3DsFacesEdgeLengths[cell3DIndex][f][e],
                                                             geometryUtilities.Tolerance1D()));
        }
      }
    }

    if (configuration.Cell1D_CheckNormals)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
          const double area = geometryUtilities.PolygonAreaByBoundaryIntegral(geometricData.Cell3DsFaces2DVertices[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdgeLengths[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdge2DTangents[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdge2DNormals[cell3DIndex][f]);

          if (!geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                area,
                                                geometryUtilities.Tolerance2D()))
          {
            std::cout.precision(16);
            std::cout<< "Cell1D_CheckNormals cell3DIndex "<< cell3DIndex<< " cell2DIndex "<< mesh.Cell3DFace(cell3DIndex, f)<< std::endl;
            std::cout<< scientific<< "Areas: "<< area<< " "<< geometricData.Cell3DsFacesAreas[cell3DIndex][f]<< std::endl;
            std::cout<< scientific<< "Difference: "<< (area - geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
            std::cout<< scientific<< "Tolerance: "<< geometryUtilities.Tolerance2D() * std::max(area, geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
          }

          Output::Assert(geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                          area,
                                                          geometryUtilities.Tolerance2D()));
        }
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
          Output::Assert(geometryUtilities.IsValuePositive(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                           geometryUtilities.Tolerance2D()));
        }
      }
    }

    if (configuration.Cell2D_CheckTriangles)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
          const double area = geometryUtilities.PolygonAreaByInternalIntegral(geometricData.Cell3DsFaces2DTriangulations[cell3DIndex][f]);

          if (!geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                area,
                                                geometryUtilities.Tolerance2D()))
          {
            std::cout.precision(16);
            std::cout<< "Cell2D_CheckTriangles cell3DIndex "<< cell3DIndex<< " cell2DIndex "<< mesh.Cell3DFace(cell3DIndex, f)<< std::endl;
            std::cout<< scientific<< "Areas: "<< area<< " "<< geometricData.Cell3DsFacesAreas[cell3DIndex][f]<< std::endl;
            std::cout<< scientific<< "Difference: "<< (area - geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
            std::cout<< scientific<< "Tolerance: "<< geometryUtilities.Tolerance2D() * std::max(area, geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
          }

          Output::Assert(geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                          area,
                                                          geometryUtilities.Tolerance2D()));
        }
      }
    }

    if (configuration.Cell2D_CheckNormals)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
        {
          const double area = geometryUtilities.PolygonAreaByBoundaryIntegral(geometricData.Cell3DsFaces2DVertices[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdgeLengths[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdge2DTangents[cell3DIndex][f],
                                                                              geometricData.Cell3DsFacesEdge2DNormals[cell3DIndex][f]);

          if (!geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                area,
                                                geometryUtilities.Tolerance2D()))
          {
            std::cout.precision(16);
            std::cout<< "Cell2D_CheckNormals cell3DIndex "<< cell3DIndex<< " cell2DIndex "<< mesh.Cell3DFace(cell3DIndex, f)<< std::endl;
            std::cout<< scientific<< "Areas: "<< area<< " "<< geometricData.Cell3DsFacesAreas[cell3DIndex][f]<< std::endl;
            std::cout<< scientific<< "Difference: "<< (area - geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
            std::cout<< scientific<< "Tolerance: "<< geometryUtilities.Tolerance2D() * std::max(area, geometricData.Cell3DsFacesAreas[cell3DIndex][f])<< std::endl;
          }

          Output::Assert(geometryUtilities.AreValuesEqual(geometricData.Cell3DsFacesAreas[cell3DIndex][f],
                                                          area,
                                                          geometryUtilities.Tolerance2D()));
        }


        const double volume = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(geometricData.Cell3DsFaces2DTriangulations[cell3DIndex],
                                                                                   geometricData.Cell3DsFacesNormals[cell3DIndex],
                                                                                   geometricData.Cell3DsFacesNormalDirections[cell3DIndex],
                                                                                   geometricData.Cell3DsFacesTranslations[cell3DIndex],
                                                                                   geometricData.Cell3DsFacesRotationMatrices[cell3DIndex]);

        if (!geometryUtilities.AreValuesEqual(geometricData.Cell3DsVolumes[cell3DIndex],
                                              volume,
                                              geometryUtilities.Tolerance3D()))
        {
          std::cout.precision(16);
          std::cout<< "Cell2D_CheckNormals cell3DIndex "<< cell3DIndex<< std::endl;
          std::cout<< scientific<< "Volumes: "<< volume<< " "<< geometricData.Cell3DsVolumes[cell3DIndex]<< std::endl;
          std::cout<< scientific<< "Difference: "<< (volume - geometricData.Cell3DsVolumes[cell3DIndex])<< std::endl;
          std::cout<< scientific<< "Tolerance: "<< geometryUtilities.Tolerance3D() * std::max(volume, geometricData.Cell3DsVolumes[cell3DIndex])<< std::endl;
        }

        Output::Assert(geometryUtilities.AreValuesEqual(geometricData.Cell3DsVolumes[cell3DIndex],
                                                        volume,
                                                        geometryUtilities.Tolerance3D()));
      }
    }

    if (configuration.Cell3D_CheckMeasure)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        Output::Assert(geometryUtilities.IsValuePositive(geometricData.Cell3DsVolumes[cell3DIndex],
                                                         geometryUtilities.Tolerance3D()));
      }
    }

    if (configuration.Cell3D_CheckTetrahedra)
    {
      for (unsigned int cell3DIndex = 0; cell3DIndex < mesh.Cell3DTotalNumber(); cell3DIndex++)
      {
        if (!mesh.Cell3DIsActive(cell3DIndex))
          continue;

        const double volume = geometryUtilities.PolyhedronVolumeByInternalIntegral(geometricData.Cell3DsTetrahedronPoints[cell3DIndex]);

        if (!geometryUtilities.AreValuesEqual(geometricData.Cell3DsVolumes[cell3DIndex],
                                              volume,
                                              geometryUtilities.Tolerance3D()))
        {
          std::cout.precision(16);
          std::cout<< "Cell3D_CheckTetrahedra cell3DIndex "<< cell3DIndex<< std::endl;
          std::cout<< scientific<< "Volume: "<< volume<< " "<< geometricData.Cell3DsVolumes[cell3DIndex]<< std::endl;
          std::cout<< scientific<< "Difference: "<< (volume - geometricData.Cell3DsVolumes[cell3DIndex])<< std::endl;
          std::cout<< scientific<< "Tolerance: "<< geometryUtilities.Tolerance3D() * std::max(volume, geometricData.Cell3DsVolumes[cell3DIndex])<< std::endl;
        }

        Output::Assert(geometryUtilities.AreValuesEqual(geometricData.Cell3DsVolumes[cell3DIndex],
                                                        volume,
                                                        geometryUtilities.Tolerance3D()));
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
  void MeshUtilities::SetMeshMarkersByFaceNormal(const GeometryUtilities& geometryUtilities,
                                                 const Eigen::Vector3d& normal,
                                                 const std::vector<Eigen::Vector3d>& cell2Ds_normal,
                                                 const unsigned int& marker,
                                                 IMeshDAO& mesh) const
  {
    const double normal_norm = normal.norm();

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); ++f)
    {
      unsigned int num_cell3D_neighs = 0;

      for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(f); ++n)
      {
        if (!mesh.Cell2DHasNeighbourCell3D(f, n))
          continue;

        num_cell3D_neighs++;
      }

      if (num_cell3D_neighs != 1)
        continue;

      const auto& cell2D_normal = cell2Ds_normal.at(f);

      if (!geometryUtilities.AreValuesEqual(normal_norm * cell2D_normal.norm(),
                                            std::abs(cell2D_normal.dot(normal)),
                                            geometryUtilities.Tolerance1DSquared()))
        continue;

      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(f); ++v)
      {
        mesh.Cell0DSetMarker(mesh.Cell2DVertex(f, v),
                             marker);
      }

      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(f); ++e)
      {
        mesh.Cell1DSetMarker(mesh.Cell2DEdge(f, e),
                             marker);
      }

      mesh.Cell2DSetMarker(f,
                           marker);
    }
  }
  // ***************************************************************************
  void MeshUtilities::SetMeshMarkersOnPolygon(const GeometryUtilities& geometryUtilities,
                                              const Eigen::Vector3d& polygon_plane_normal,
                                              const Eigen::Vector3d& polygon_plane_origin,
                                              const MatrixXd& polygon_vertices_2D,
                                              const Vector3d& polygon_translation,
                                              const Matrix3d& polygon_rotation_matrix,
                                              const unsigned int& marker,
                                              IMeshDAO& mesh) const
  {
    // set cell0Ds markers
    std::vector<bool> vertices_on_polygon(mesh.Cell0DTotalNumber(),
                                          false);

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      const Eigen::Vector3d vertex = mesh.Cell0DCoordinates(v);

      if (!geometryUtilities.IsPointOnPlane(vertex,
                                            polygon_plane_normal,
                                            polygon_plane_origin))
        continue;

      const Eigen::Vector3d vertex_2D = geometryUtilities.RotatePointsFrom3DTo2D(vertex,
                                                                                 polygon_rotation_matrix.transpose(),
                                                                                 polygon_translation);

      if (!geometryUtilities.IsPointInsidePolygon(vertex_2D,
                                                  polygon_vertices_2D))
        continue;


      vertices_on_polygon[v] = true;
      mesh.Cell0DSetMarker(v,
                           marker);
    }

    // set cell1Ds markers
    for (unsigned int s = 0; s < mesh.Cell1DTotalNumber(); s++)
    {
      const Eigen::VectorXi extremes = mesh.Cell1DExtremes(s);

      if (vertices_on_polygon[extremes[0]] &&
          vertices_on_polygon[extremes[1]])
      {
        mesh.Cell1DSetMarker(s,
                             marker);
      }
    }

    // set cell2Ds markers
    for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
    {
      const vector<unsigned int> extremes = mesh.Cell2DVertices(p);

      bool is_on_polygon = true;
      for (unsigned int v = 0; v < extremes.size(); v++)
      {
        if (!vertices_on_polygon[extremes[v]])
        {
          is_on_polygon = false;
          break;
        }
      }

      if (!is_on_polygon)
        continue;

      mesh.Cell2DSetMarker(p,
                           marker);
    }
  }
  // ***************************************************************************
  void MeshUtilities::SetMeshMarkersOnPolygon(const GeometryUtilities& geometryUtilities,
                                              const Eigen::Vector3d& polygon_plane_normal,
                                              const Eigen::Vector3d& polygon_plane_origin,
                                              const Eigen::MatrixXd& polygon_vertices_2D,
                                              const Eigen::Vector3d& polygon_translation,
                                              const Eigen::Matrix3d& polygon_rotation_matrix,
                                              const std::vector<Eigen::Vector3d>& cell1Ds_centroid,
                                              const std::vector<Eigen::Vector3d>& cell2Ds_centroid,
                                              const unsigned int& marker,
                                              IMeshDAO& mesh) const
  {
    // set cell0Ds markers
    std::vector<bool> vertices_on_polygon(mesh.Cell0DTotalNumber(),
                                          false);

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      const Eigen::Vector3d vertex = mesh.Cell0DCoordinates(v);

      if (!geometryUtilities.IsPointOnPlane(vertex,
                                            polygon_plane_normal,
                                            polygon_plane_origin))
        continue;

      const Eigen::Vector3d vertex_2D = geometryUtilities.RotatePointsFrom3DTo2D(vertex,
                                                                                 polygon_rotation_matrix.transpose(),
                                                                                 polygon_translation);

      if (!geometryUtilities.IsPointInsidePolygon_RayCasting(vertex_2D,
                                                             polygon_vertices_2D))
        continue;


      vertices_on_polygon[v] = true;
      mesh.Cell0DSetMarker(v,
                           marker);
    }

    // set cell1Ds markers
    for (unsigned int s = 0; s < mesh.Cell1DTotalNumber(); s++)
    {
      const Eigen::VectorXi extremes = mesh.Cell1DExtremes(s);

      if (!vertices_on_polygon[extremes[0]] ||
          !vertices_on_polygon[extremes[1]])
        continue;

      const Eigen::Vector3d& cell1D_centroid = cell1Ds_centroid.at(s);

      const Eigen::Vector3d cell1D_centroid_2D = geometryUtilities.RotatePointsFrom3DTo2D(cell1D_centroid,
                                                                                          polygon_rotation_matrix.transpose(),
                                                                                          polygon_translation);


      if (!geometryUtilities.IsPointInsidePolygon_RayCasting(cell1D_centroid_2D,
                                                             polygon_vertices_2D))
        continue;

      mesh.Cell1DSetMarker(s,
                           marker);
    }

    // set cell2Ds markers
    for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
    {
      const vector<unsigned int> extremes = mesh.Cell2DVertices(p);

      bool is_on_polygon = true;
      for (unsigned int v = 0; v < extremes.size(); v++)
      {
        if (!vertices_on_polygon[extremes[v]])
        {
          is_on_polygon = false;
          break;
        }
      }

      if (!is_on_polygon)
        continue;

      const Eigen::Vector3d& cell2D_centroid = cell2Ds_centroid.at(p);

      const Eigen::Vector3d cell2D_centroid_2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2D_centroid,
                                                                                          polygon_rotation_matrix.transpose(),
                                                                                          polygon_translation);


      if (!geometryUtilities.IsPointInsidePolygon_RayCasting(cell2D_centroid_2D,
                                                             polygon_vertices_2D))
        continue;

      mesh.Cell2DSetMarker(p,
                           marker);
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
    result.Cell3DsBoundingBox.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgesCentroid.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeLengths.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalGlobalDirection.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangentsGlobalDirection.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces3DTriangulations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges3DCentroid.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges2DCentroid.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge3DTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(convexMesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  convexMesh.Cell3DTotalNumber(); c++)
    {
      if (!convexMesh.Cell3DIsActive(c))
        continue;

      const GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                              c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsBoundingBox[c] = geometryUtilities.PointsBoundingBox(result.Cell3DsVertices[c]);

      result.Cell3DsEdgesCentroid[c] = geometryUtilities.PolyhedronEdgesCentroid(result.Cell3DsVertices[c],
                                                                                 result.Cell3DsEdges[c]);
      result.Cell3DsEdgeLengths[c] = geometryUtilities.PolyhedronEdgesLength(result.Cell3DsVertices[c],
                                                                             result.Cell3DsEdges[c]);

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

      for (unsigned int f = 0; f < numFaces; f++)
      {
        const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(result.Cell3DsFaces2DVertices[c][f],
                                                                                 false);


        if (geometryUtilities.PolygonOrientation(faceConvexHull) ==
            Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
        {
          const auto cell2DIndex = convexMesh.Cell3DFace(c,
                                                         f);
          std::cout<< "Cell3D "<< c<< " face "<< cell2DIndex<< std::endl;
          std::cout<< "Cell3D vertices "<< convexMesh.Cell3DVertices(c)<< std::endl;
          std::cout<< "Cell2D vertices "<< convexMesh.Cell2DVertices(cell2DIndex)<< std::endl;

          std::cout<< "Cell2D hull vertices";
          for (const auto faceConvexHullIndex : faceConvexHull)
            std::cout<< " "<< convexMesh.Cell2DVertex(cell2DIndex,
                                                      faceConvexHullIndex);
          std::cout<< std::endl;
        }

        Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                       Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
      }

      result.Cell3DsFaces3DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces3DVertices[c],
                                                                                                          polyhedronFaceTriangulations);
      result.Cell3DsFaces2DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces2DVertices[c],
                                                                                                          polyhedronFaceTriangulations);


      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdges3DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdges2DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdge3DTangents[c].resize(numFaces);
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
        result.Cell3DsFacesEdges3DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdges2DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge3DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DNormals[c][f] = geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);
      }

      result.Cell3DsFacesNormalDirections[c] = geometryUtilities.PolyhedronFaceNormalDirections(result.Cell3DsFaces3DVertices[c],
                                                                                                geometryUtilities.PolyhedronBarycenter(result.Cell3DsVertices[c]),
                                                                                                result.Cell3DsFacesNormals[c]);


      result.Cell3DsFacesTangents[c] = geometryUtilities.PolyhedronFaceTangents(result.Cell3DsFaces3DVertices[c],
                                                                                result.Cell3DsFacesNormals[c],
                                                                                result.Cell3DsFacesNormalDirections[c]);

      result.Cell3DsFacesNormalGlobalDirection[c].resize(numFaces);
      result.Cell3DsFacesTangentsGlobalDirection[c].resize(numFaces);

      for (unsigned int f = 0; f < numFaces; ++f)
      {
        const unsigned int cell2D_index = convexMesh.Cell3DFace(c, f);

        Gedim::Output::Assert(convexMesh.Cell2DNumberNeighbourCell3D(cell2D_index) > 0);
        unsigned int first_cell3D_neigh_position = 0;
        for (unsigned int f_n = 0; f_n < convexMesh.Cell2DNumberNeighbourCell3D(cell2D_index); ++f_n)
        {
          if (!convexMesh.Cell2DHasNeighbourCell3D(cell2D_index, f_n))
            continue;

          first_cell3D_neigh_position = f_n;
          break;
        }

        result.Cell3DsFacesNormalGlobalDirection[c][f] =
            convexMesh.Cell2DNeighbourCell3D(cell2D_index,
                                             first_cell3D_neigh_position) == c;

        result.Cell3DsFacesTangentsGlobalDirection[c][f][0] =
            result.Cell3DsFacesEdgeDirections[c][f][0];
        result.Cell3DsFacesTangentsGlobalDirection[c][f][1] =
            result.Cell3DsFacesTangentsGlobalDirection[c][f][0] ==
            result.Cell3DsFacesNormalGlobalDirection[c][f];
      }

      result.Cell3DsVolumes[c] = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(result.Cell3DsFaces2DTriangulations[c],
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
    result.Cell3DsBoundingBox.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgesCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeLengths.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalGlobalDirection.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangentsGlobalDirection.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces3DTriangulations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges3DCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges2DCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge3DTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(mesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  mesh.Cell3DTotalNumber(); c++)
    {
      if (!mesh.Cell3DIsActive(c))
        continue;

      const vector<unsigned int>& convexCell3DIndices = meshCell3DToConvexCell3DIndices[c];
      const unsigned int& numConvexCell3Ds = convexCell3DIndices.size();

      const GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh,
                                                                              c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsBoundingBox[c] = geometryUtilities.PointsBoundingBox(result.Cell3DsVertices[c]);

      result.Cell3DsEdgesCentroid[c] = geometryUtilities.PolyhedronEdgesCentroid(result.Cell3DsVertices[c],
                                                                                 result.Cell3DsEdges[c]);
      result.Cell3DsEdgeLengths[c] = geometryUtilities.PolyhedronEdgesLength(result.Cell3DsVertices[c],
                                                                             result.Cell3DsEdges[c]);

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
        const unsigned int cell2DIndex = mesh.Cell3DFace(c, f);

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
      std::vector<bool> face2DCCWOrientation(numFaces);
      std::vector<Eigen::MatrixXd> face2DVerticesCCW(numFaces);
      std::vector<std::vector<unsigned int>> face2DsCCW(numFaces);

      for (unsigned int f = 0; f < numFaces; f++)
      {
        const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(result.Cell3DsFaces2DVertices[c][f],
                                                                                 false);

        switch (geometryUtilities.PolygonOrientation(faceConvexHull))
        {
          case Gedim::GeometryUtilities::PolygonOrientations::Clockwise:
          {
            face2DCCWOrientation[f] = false;
            face2DsCCW[f] = geometryUtilities.ChangePolygonOrientation(result.Cell3DsFaces2DVertices[c][f].cols());
            face2DVerticesCCW[f] = geometryUtilities.ExtractPoints(result.Cell3DsFaces2DVertices[c][f],
                                                                   face2DsCCW[f]);
          }
            break;
          case Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise:
          {
            face2DCCWOrientation[f] = true;
            face2DVerticesCCW[f] = result.Cell3DsFaces2DVertices[c][f];

            face2DsCCW[f].resize(result.Cell3DsFaces2DVertices[c][f].cols());
            std::iota(face2DsCCW[f].begin(), face2DsCCW[f].end(), 0);
          }
            break;
          default:
            throw std::runtime_error("Not managed face orientation case");
        }
      }

      for (unsigned int f = 0; f < numFaces; f++)
      {
        const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(face2DVerticesCCW[f],
                                                                                 false);

        Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                       Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
      }

      {
        const vector<vector<unsigned int>> polyhedronFaceTriangulationsCCW = geometryUtilities.PolyhedronFaceTriangulationsByEarClipping(numFaces,
                                                                                                                                         face2DVerticesCCW);

        vector<vector<unsigned int>> polyhedronFaceTriangulations = polyhedronFaceTriangulationsCCW;
        for (unsigned int f = 0; f < numFaces; f++)
        {
          for (unsigned int nt = 0; nt < polyhedronFaceTriangulations[f].size(); nt++)
            polyhedronFaceTriangulations[f][nt] = face2DsCCW[f][polyhedronFaceTriangulationsCCW[f][nt]];
        }

        result.Cell3DsFaces3DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces3DVertices[c],
                                                                                                            polyhedronFaceTriangulations);
        result.Cell3DsFaces2DTriangulations[c] = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(result.Cell3DsFaces2DVertices[c],
                                                                                                            polyhedronFaceTriangulations);
      }


      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdges3DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdges2DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdge3DTangents[c].resize(numFaces);
      result.Cell3DsFacesEdge2DTangents[c].resize(numFaces);

      for(unsigned int f = 0; f < numFaces; f++)
      {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d>& concaveCell2DTriangulationPoints = result.Cell3DsFaces2DTriangulations[c][f];
        const unsigned int& numConcaveCell2DTriangulation = concaveCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd concaveCell2DTriangulationAreas(numConcaveCell2DTriangulation);
        Eigen::MatrixXd concaveCell2DTriangulationCentroids(3, numConcaveCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConcaveCell2DTriangulation; cct++)
        {
          concaveCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(concaveCell2DTriangulationPoints[cct]);
          concaveCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(concaveCell2DTriangulationPoints[cct]);
        }

        result.Cell3DsFacesAreas[c][f] = concaveCell2DTriangulationAreas.sum();
        result.Cell3DsFaces2DCentroids[c][f] = geometryUtilities.PolygonCentroid(concaveCell2DTriangulationCentroids,
                                                                                 concaveCell2DTriangulationAreas,
                                                                                 result.Cell3DsFacesAreas[c][f]);
        result.Cell3DsFacesDiameters[c][f] = geometryUtilities.PolygonDiameter(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdgeLengths[c][f] = geometryUtilities.PolygonEdgeLengths(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdges3DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdges2DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge3DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);

        const double facesNormalOrientation = face2DCCWOrientation[f] ? +1.0 : -1.0;
        result.Cell3DsFacesEdge2DNormals[c][f] = facesNormalOrientation *
                                                 geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);

      }

      unsigned int cell3DTetrahedronsSize = 0;
      std::vector<std::vector<Eigen::MatrixXd>> convexCell3DTetrahedronsPoints(numConvexCell3Ds);
      Eigen::VectorXd convexCell3DsVolume(numConvexCell3Ds);
      Eigen::MatrixXd convexCell3DsCentroid(3, numConvexCell3Ds);
      std::vector<std::vector<Eigen::MatrixXd>> convexCell3DsFaces3DVertices(numConvexCell3Ds);
      std::vector<std::vector<std::vector<unsigned int>>> convexCell3DsFacesUnalignedVertices(numConvexCell3Ds);
      std::vector<std::vector<Eigen::Vector3d>> convexCell3DsNormal(numConvexCell3Ds);
      std::vector<std::vector<bool>> convexCell3DsNormalDirections(numConvexCell3Ds);

      for (unsigned int cc = 0; cc < numConvexCell3Ds; cc++)
      {
        const unsigned int convexCell3DIndex = convexCell3DIndices[cc];

        const GeometryUtilities::Polyhedron convexCell3DPolyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                                            convexCell3DIndex);

        convexCell3DsFaces3DVertices[cc] = geometryUtilities.PolyhedronFaceVertices(convexCell3DPolyhedron.Vertices,
                                                                                    convexCell3DPolyhedron.Faces);
        const std::vector<Eigen::Vector3d> convexCell3DFacesTranslation = geometryUtilities.PolyhedronFaceTranslations(convexCell3DsFaces3DVertices[cc]);
        convexCell3DsNormal[cc] = geometryUtilities.PolyhedronFaceNormals(convexCell3DsFaces3DVertices[cc]);
        const std::vector<Eigen::Matrix3d> convexCell3DFacesRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(convexCell3DsFaces3DVertices[cc],
                                                                                                                                convexCell3DsNormal[cc],
                                                                                                                                convexCell3DFacesTranslation);

        const std::vector<Eigen::MatrixXd> convexCell3DFaces2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(convexCell3DsFaces3DVertices[cc],
                                                                                                                         convexCell3DFacesTranslation,
                                                                                                                         convexCell3DFacesRotationMatrices);

        for (unsigned int ccf = 0; ccf < convexCell3DFaces2DVertices.size(); ccf++)
        {
          const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(convexCell3DFaces2DVertices[ccf],
                                                                                   false);

          if (geometryUtilities.PolygonOrientation(faceConvexHull) ==
              Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
          {
            const auto cell2DIndex = convexMesh.Cell3DFace(convexCell3DIndex,
                                                           ccf);
            std::cout<< "Cell3D "<< convexCell3DIndex<< " face "<< cell2DIndex<< std::endl;
            std::cout<< "Cell3D vertices "<< convexMesh.Cell3DVertices(convexCell3DIndex)<< std::endl;
            std::cout<< "Cell2D vertices "<< convexMesh.Cell2DVertices(cell2DIndex)<< std::endl;

            std::cout<< "Cell2D hull vertices";
            for (const auto faceConvexHullIndex : faceConvexHull)
              std::cout<< " "<< convexMesh.Cell2DVertex(cell2DIndex,
                                                        faceConvexHullIndex);
            std::cout<< std::endl;
          }

          Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                         Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
        }

        convexCell3DsFacesUnalignedVertices[cc] = geometryUtilities.PolyhedronFacesUnalignedVertices(convexCell3DFaces2DVertices);

        const std::vector<std::vector<unsigned int>> convexCell3DFacesTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(convexCell3DPolyhedron.Faces,
                                                                                                                                                   convexCell3DsFaces3DVertices[cc]);
        const std::vector<std::vector<Eigen::Matrix3d>> convexCell3DFaces2DTriangulations = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(convexCell3DFaces2DVertices,
                                                                                                                                                       convexCell3DFacesTriangulations);

        convexCell3DsNormalDirections[cc] = geometryUtilities.PolyhedronFaceNormalDirections(convexCell3DsFaces3DVertices[cc],
                                                                                             geometryUtilities.PolyhedronBarycenter(convexCell3DPolyhedron.Vertices),
                                                                                             convexCell3DsNormal[cc]);

        convexCell3DsVolume[cc] = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(convexCell3DFaces2DTriangulations,
                                                                                       convexCell3DsNormal[cc],
                                                                                       convexCell3DsNormalDirections[cc],
                                                                                       convexCell3DFacesTranslation,
                                                                                       convexCell3DFacesRotationMatrices);

        convexCell3DsCentroid.col(cc) = geometryUtilities.PolyhedronCentroid(convexCell3DFaces2DTriangulations,
                                                                             convexCell3DsNormal[cc],
                                                                             convexCell3DsNormalDirections[cc],
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

      const FindConcaveCell3DFacesConvexCell2DResult convexCell2Ds = FindConcaveCell3DFacesConvexCell2D(geometryUtilities,
                                                                                                        c,
                                                                                                        mesh,
                                                                                                        convexMesh,
                                                                                                        convexCell3DIndices,
                                                                                                        result.Cell3DsFaces3DVertices[c],
                                                                                                        face2DVerticesCCW,
                                                                                                        result.Cell3DsFacesTranslations[c],
                                                                                                        result.Cell3DsFacesRotationMatrices[c],
                                                                                                        result.Cell3DsFacesNormals[c],
                                                                                                        convexCell3DsFaces3DVertices,
                                                                                                        convexCell3DsFacesUnalignedVertices);

      result.Cell3DsFacesNormalDirections[c].resize(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
      {
        const Eigen::Vector3d& faceNormal = result.Cell3DsFacesNormals[c][f];

        const unsigned int cc = convexCell2Ds.ConcaveCell3DFacesConvexCell2D[f].ConvexCell3DIndex;
        const unsigned int ccf = convexCell2Ds.ConcaveCell3DFacesConvexCell2D[f].ConvexCell3DFaceIndex;

        const Eigen::Vector3d outgoingFaceNormal =
            convexCell3DsNormalDirections[cc][ccf] ?
              convexCell3DsNormal[cc][ccf] : -1.0 * convexCell3DsNormal[cc][ccf];

        result.Cell3DsFacesNormalDirections[c][f] = geometryUtilities.IsValuePositive(faceNormal.dot(outgoingFaceNormal),
                                                                                      geometryUtilities.Tolerance1DSquared());
      }

      result.Cell3DsFacesTangents[c] = geometryUtilities.PolyhedronFaceTangents(result.Cell3DsFaces3DVertices[c],
                                                                                result.Cell3DsFacesNormals[c],
                                                                                result.Cell3DsFacesNormalDirections[c]);
      result.Cell3DsFacesNormalGlobalDirection[c].resize(numFaces);
      result.Cell3DsFacesTangentsGlobalDirection[c].resize(numFaces);

      for (unsigned int f = 0; f < numFaces; ++f)
      {
        const unsigned int cell2D_index = mesh.Cell3DFace(c, f);

        Gedim::Output::Assert(mesh.Cell2DNumberNeighbourCell3D(cell2D_index) > 0);
        unsigned int first_cell3D_neigh_position = 0;
        for (unsigned int f_n = 0; f_n < mesh.Cell2DNumberNeighbourCell3D(cell2D_index); ++f_n)
        {
          if (!mesh.Cell2DHasNeighbourCell3D(cell2D_index, f_n))
            continue;

          first_cell3D_neigh_position = f_n;
          break;
        }

        result.Cell3DsFacesNormalGlobalDirection[c][f] =
            mesh.Cell2DNeighbourCell3D(cell2D_index,
                                             first_cell3D_neigh_position) == c;

        result.Cell3DsFacesTangentsGlobalDirection[c][f][0] =
            result.Cell3DsFacesEdgeDirections[c][f][0];
        result.Cell3DsFacesTangentsGlobalDirection[c][f][1] =
            result.Cell3DsFacesTangentsGlobalDirection[c][f][0] ==
            result.Cell3DsFacesNormalGlobalDirection[c][f];
      }
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData3D MeshUtilities::FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const std::vector<std::vector<Eigen::MatrixXd> >& cell3Ds_tetra_vertices,
                                                                            const std::vector<std::vector<Eigen::Matrix3d> >& cell2Ds_triangles_3D_vertices) const
  {
    MeshGeometricData3D result;

    result.Cell3DsVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdges.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsBoundingBox.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgesCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeLengths.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalGlobalDirection.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesTangentsGlobalDirection.resize(mesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces3DTriangulations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges3DCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdges2DCentroid.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge3DTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(mesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(mesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  mesh.Cell3DTotalNumber(); c++)
    {
      if (!mesh.Cell3DIsActive(c))
        continue;

      const GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh,
                                                                              c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsBoundingBox[c] = geometryUtilities.PointsBoundingBox(result.Cell3DsVertices[c]);

      result.Cell3DsEdgesCentroid[c] = geometryUtilities.PolyhedronEdgesCentroid(result.Cell3DsVertices[c],
                                                                                 result.Cell3DsEdges[c]);
      result.Cell3DsEdgeLengths[c] = geometryUtilities.PolyhedronEdgesLength(result.Cell3DsVertices[c],
                                                                             result.Cell3DsEdges[c]);

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
        const unsigned int cell2DIndex = mesh.Cell3DFace(c, f);

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
      std::vector<bool> face2DCCWOrientation(numFaces);
      std::vector<Eigen::MatrixXd> face2DVerticesCCW(numFaces);

      for (unsigned int f = 0; f < numFaces; f++)
      {
        const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(result.Cell3DsFaces2DVertices[c][f],
                                                                                 false);

        switch (geometryUtilities.PolygonOrientation(faceConvexHull))
        {
          case Gedim::GeometryUtilities::PolygonOrientations::Clockwise:
          {
            face2DCCWOrientation[f] = false;
            const auto new_orientation_points = geometryUtilities.ChangePolygonOrientation(result.Cell3DsFaces2DVertices[c][f].cols());
            face2DVerticesCCW[f] = geometryUtilities.ExtractPoints(result.Cell3DsFaces2DVertices[c][f],
                                                                   new_orientation_points);
          }
            break;
          case Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise:
          {
            face2DCCWOrientation[f] = true;
            face2DVerticesCCW[f] = result.Cell3DsFaces2DVertices[c][f];
          }
            break;
          default:
            throw std::runtime_error("Not managed face orientation case");
        }
      }

      for (unsigned int f = 0; f < numFaces; f++)
      {
        const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(face2DVerticesCCW[f],
                                                                                 false);

        Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                       Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
      }

      {
        result.Cell3DsFaces3DTriangulations[c].resize(numFaces);
        result.Cell3DsFaces2DTriangulations[c].resize(numFaces);
        for (unsigned int f = 0; f < numFaces; ++f)
        {
          const unsigned int c2D_index = mesh.Cell3DFace(c, f);

          const auto& face_rotation = result.Cell3DsFacesRotationMatrices[c][f];
          const auto& face_translation = result.Cell3DsFacesTranslations[c][f];


          result.Cell3DsFaces3DTriangulations[c][f] = cell2Ds_triangles_3D_vertices.at(c2D_index);
          const unsigned int num_triangles = result.Cell3DsFaces3DTriangulations[c][f].size();

          result.Cell3DsFaces2DTriangulations[c][f].resize(num_triangles);

          for (unsigned int t = 0; t < num_triangles; ++t)
          {
            result.Cell3DsFaces2DTriangulations[c][f][t] =
                geometryUtilities.RotatePointsFrom3DTo2D(result.Cell3DsFaces3DTriangulations[c][f][t],
                                                         face_rotation.transpose(),
                                                         face_translation);

            const vector<unsigned int> face_triangle_convex_hull = geometryUtilities.ConvexHull(result.Cell3DsFaces2DTriangulations[c][f][t],
                                                                                                false);

            switch (geometryUtilities.PolygonOrientation(face_triangle_convex_hull))
            {
              case Gedim::GeometryUtilities::PolygonOrientations::Clockwise:
              {
                const auto new_orientation = geometryUtilities.ChangePolygonOrientation(result.Cell3DsFaces2DTriangulations[c][f][t].cols());
                result.Cell3DsFaces2DTriangulations[c][f][t] = geometryUtilities.ExtractPoints(result.Cell3DsFaces2DTriangulations[c][f][t],
                                                                                               new_orientation);
              }
                break;
              case Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise:
                break;
              default:
                throw std::runtime_error("Not managed face orientation case");
            }
          }
        }
      }


      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdges3DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdges2DCentroid[c].resize(numFaces);
      result.Cell3DsFacesEdge3DTangents[c].resize(numFaces);
      result.Cell3DsFacesEdge2DTangents[c].resize(numFaces);

      for(unsigned int f = 0; f < numFaces; f++)
      {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d>& concaveCell2DTriangulationPoints = result.Cell3DsFaces2DTriangulations[c][f];
        const unsigned int& numConcaveCell2DTriangulation = concaveCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd concaveCell2DTriangulationAreas(numConcaveCell2DTriangulation);
        Eigen::MatrixXd concaveCell2DTriangulationCentroids(3, numConcaveCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConcaveCell2DTriangulation; cct++)
        {
          concaveCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(concaveCell2DTriangulationPoints[cct]);
          concaveCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(concaveCell2DTriangulationPoints[cct]);
        }

        result.Cell3DsFacesAreas[c][f] = concaveCell2DTriangulationAreas.sum();
        result.Cell3DsFaces2DCentroids[c][f] = geometryUtilities.PolygonCentroid(concaveCell2DTriangulationCentroids,
                                                                                 concaveCell2DTriangulationAreas,
                                                                                 result.Cell3DsFacesAreas[c][f]);
        result.Cell3DsFacesDiameters[c][f] = geometryUtilities.PolygonDiameter(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdgeLengths[c][f] = geometryUtilities.PolygonEdgeLengths(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdges3DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdges2DCentroid[c][f] = geometryUtilities.PolygonEdgesCentroid(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge3DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces3DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);

        const double facesNormalOrientation = face2DCCWOrientation[f] ? +1.0 : -1.0;
        result.Cell3DsFacesEdge2DNormals[c][f] = facesNormalOrientation *
                                                 geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);

      }

      result.Cell3DsTetrahedronPoints[c] = cell3Ds_tetra_vertices.at(c);
      const unsigned int num_cell3D_tetra = result.Cell3DsTetrahedronPoints.at(c).size();

      Eigen::VectorXd cell3D_tetra_volume(num_cell3D_tetra);
      Eigen::MatrixXd cell3D_tetra_centroid(3, num_cell3D_tetra);
      std::vector<std::vector<Eigen::MatrixXd>> cell3D_tetra_faces_3D_vertices(num_cell3D_tetra);
      std::vector<std::vector<std::vector<unsigned int>>> cell3D_tetra_faces_unaligned_vertices(num_cell3D_tetra);
      std::vector<std::vector<Eigen::Vector3d>> cell3D_tetra_normal(num_cell3D_tetra);
      std::vector<std::vector<bool>> cell3D_tetra_normal_directions(num_cell3D_tetra);

      for (unsigned int t = 0; t < num_cell3D_tetra; ++t)
      {
        const auto& tetra_vertices = result.Cell3DsTetrahedronPoints.at(c).at(t);

        const auto tetra = geometryUtilities.CreateTetrahedronWithVertices(tetra_vertices.col(0),
                                                                           tetra_vertices.col(1),
                                                                           tetra_vertices.col(2),
                                                                           tetra_vertices.col(3));

        cell3D_tetra_faces_3D_vertices[t] = geometryUtilities.PolyhedronFaceVertices(tetra.Vertices,
                                                                                     tetra.Faces);
        const std::vector<Eigen::Vector3d> cell3D_tetra_faces_translation = geometryUtilities.PolyhedronFaceTranslations(cell3D_tetra_faces_3D_vertices[t]);
        cell3D_tetra_normal[t] = geometryUtilities.PolyhedronFaceNormals(cell3D_tetra_faces_3D_vertices[t]);
        const std::vector<Eigen::Matrix3d> cell3D_tetra_faces_rotation = geometryUtilities.PolyhedronFaceRotationMatrices(cell3D_tetra_faces_3D_vertices[t],
                                                                                                                          cell3D_tetra_normal[t],
                                                                                                                          cell3D_tetra_faces_translation);

        const std::vector<Eigen::MatrixXd> cell3D_tetra_faces_2D_vertices = geometryUtilities.PolyhedronFaceRotatedVertices(cell3D_tetra_faces_3D_vertices[t],
                                                                                                                            cell3D_tetra_faces_translation,
                                                                                                                            cell3D_tetra_faces_rotation);

        for (unsigned int ccf = 0; ccf < cell3D_tetra_faces_2D_vertices.size(); ++ccf)
        {
          const vector<unsigned int> faceConvexHull = geometryUtilities.ConvexHull(cell3D_tetra_faces_2D_vertices[ccf],
                                                                                   false);

          if (geometryUtilities.PolygonOrientation(faceConvexHull) ==
              Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
          {
            std::cout<< "Cell3D "<< c<< " tetra "<< t<< " face "<< ccf<< " orientation wrong"<< std::endl;
          }

          Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                         Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
        }

        cell3D_tetra_faces_unaligned_vertices[t] = geometryUtilities.PolyhedronFacesUnalignedVertices(cell3D_tetra_faces_2D_vertices);

        const std::vector<std::vector<unsigned int>> cell3D_tetra_faces_triangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(tetra.Faces,
                                                                                                                                                     cell3D_tetra_faces_3D_vertices[t]);
        const std::vector<std::vector<Eigen::Matrix3d>> cell3D_tetra_faces_2D_triangulations = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(cell3D_tetra_faces_2D_vertices,
                                                                                                                                                          cell3D_tetra_faces_triangulations);

        cell3D_tetra_normal_directions[t] = geometryUtilities.PolyhedronFaceNormalDirections(cell3D_tetra_faces_3D_vertices[t],
                                                                                             geometryUtilities.PolyhedronBarycenter(tetra.Vertices),
                                                                                             cell3D_tetra_normal[t]);

        cell3D_tetra_volume[t] = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(cell3D_tetra_faces_2D_triangulations,
                                                                                      cell3D_tetra_normal[t],
                                                                                      cell3D_tetra_normal_directions[t],
                                                                                      cell3D_tetra_faces_translation,
                                                                                      cell3D_tetra_faces_rotation);

        cell3D_tetra_centroid.col(t) = geometryUtilities.PolyhedronCentroid(cell3D_tetra_faces_2D_triangulations,
                                                                            cell3D_tetra_normal[t],
                                                                            cell3D_tetra_normal_directions[t],
                                                                            cell3D_tetra_faces_translation,
                                                                            cell3D_tetra_faces_rotation,
                                                                            cell3D_tetra_volume[t]);
      }

      result.Cell3DsVolumes[c] = cell3D_tetra_volume.sum();
      result.Cell3DsCentroids[c] = cell3D_tetra_centroid * cell3D_tetra_volume / result.Cell3DsVolumes[c];

      const FindConcaveCell3DFacesConvexCell2DResult convexCell2Ds = FindConcaveCell3DFacesConvexCell2D(geometryUtilities,
                                                                                                        c,
                                                                                                        mesh,
                                                                                                        result.Cell3DsTetrahedronPoints.at(c),
                                                                                                        result.Cell3DsFaces2DTriangulations.at(c),
                                                                                                        result.Cell3DsFaces3DVertices.at(c),
                                                                                                        result.Cell3DsFacesTranslations.at(c),
                                                                                                        result.Cell3DsFacesRotationMatrices.at(c),
                                                                                                        result.Cell3DsFacesNormals.at(c),
                                                                                                        cell3D_tetra_faces_3D_vertices,
                                                                                                        cell3D_tetra_faces_unaligned_vertices);

      result.Cell3DsFacesNormalDirections[c].resize(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
      {
        const Eigen::Vector3d& faceNormal = result.Cell3DsFacesNormals[c][f];

        const unsigned int cc = convexCell2Ds.ConcaveCell3DFacesConvexCell2D[f].ConvexCell3DIndex;
        const unsigned int ccf = convexCell2Ds.ConcaveCell3DFacesConvexCell2D[f].ConvexCell3DFaceIndex;

        const Eigen::Vector3d outgoingFaceNormal =
            cell3D_tetra_normal_directions[cc][ccf] ?
              cell3D_tetra_normal[cc][ccf] : -1.0 * cell3D_tetra_normal[cc][ccf];

        result.Cell3DsFacesNormalDirections[c][f] = geometryUtilities.IsValuePositive(faceNormal.dot(outgoingFaceNormal),
                                                                                      geometryUtilities.Tolerance1DSquared());
      }

      result.Cell3DsFacesTangents[c] = geometryUtilities.PolyhedronFaceTangents(result.Cell3DsFaces3DVertices[c],
                                                                                result.Cell3DsFacesNormals[c],
                                                                                result.Cell3DsFacesNormalDirections[c]);

      result.Cell3DsFacesNormalGlobalDirection[c].resize(numFaces);
      result.Cell3DsFacesTangentsGlobalDirection[c].resize(numFaces);

      for (unsigned int f = 0; f < numFaces; ++f)
      {
        const unsigned int cell2D_index = mesh.Cell3DFace(c, f);

        Gedim::Output::Assert(mesh.Cell2DNumberNeighbourCell3D(cell2D_index) > 0);
        unsigned int first_cell3D_neigh_position = 0;
        for (unsigned int f_n = 0; f_n < mesh.Cell2DNumberNeighbourCell3D(cell2D_index); ++f_n)
        {
          if (!mesh.Cell2DHasNeighbourCell3D(cell2D_index, f_n))
            continue;

          first_cell3D_neigh_position = f_n;
          break;
        }

        result.Cell3DsFacesNormalGlobalDirection[c][f] =
            mesh.Cell2DNeighbourCell3D(cell2D_index,
                                             first_cell3D_neigh_position) == c;

        result.Cell3DsFacesTangentsGlobalDirection[c][f][0] =
            result.Cell3DsFacesEdgeDirections[c][f][0];
        result.Cell3DsFacesTangentsGlobalDirection[c][f][1] =
            result.Cell3DsFacesTangentsGlobalDirection[c][f][0] ==
            result.Cell3DsFacesNormalGlobalDirection[c][f];
      }
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::FindConcaveCell3DFacesConvexCell2DResult MeshUtilities::FindConcaveCell3DFacesConvexCell2D(const GeometryUtilities& geometryUtilities,
                                                                                                            const unsigned int& concaveCell3DIndex,
                                                                                                            const IMeshDAO& mesh,
                                                                                                            const IMeshDAO& convexMesh,
                                                                                                            const std::vector<unsigned int>& convexCell3DIndices,
                                                                                                            const std::vector<Eigen::MatrixXd>& concaveCell3DFaces3DVertices,
                                                                                                            const std::vector<Eigen::MatrixXd>& concaveCell3DFaces2DVertices,
                                                                                                            const std::vector<Eigen::Vector3d>& concaveCell3DFacesTranslation,
                                                                                                            const std::vector<Eigen::Matrix3d>& concaveCell3DFacesRotationMatrix,
                                                                                                            const std::vector<Eigen::Vector3d>& concaveCell3DFacesNormal,
                                                                                                            const std::vector<std::vector<Eigen::MatrixXd>>& convexCell3DsFaces3DVertices,
                                                                                                            const std::vector<std::vector<std::vector<unsigned int>>>& convexCell3DsFacesUnalignedVertices) const
  {
    FindConcaveCell3DFacesConvexCell2DResult result;

    const unsigned int& numConvexCell3Ds = convexCell3DIndices.size();
    const unsigned int numConcaveFaces = mesh.Cell3DNumberFaces(concaveCell3DIndex);

    result.ConcaveCell3DFacesConvexCell2D.resize(numConcaveFaces);

    for (unsigned int f = 0; f < numConcaveFaces; f++)
    {
      FindConcaveCell3DFacesConvexCell2DResult::ConvexCell2D& convexCell2D = result.ConcaveCell3DFacesConvexCell2D[f];

      const Eigen::MatrixXd& faceVertices = concaveCell3DFaces3DVertices[f];
      const Eigen::Vector3d& faceOrigin = faceVertices.col(0);
      const Eigen::MatrixXd& face2DVertices = concaveCell3DFaces2DVertices[f];
      const Eigen::Matrix3d& faceRotationMatrix = concaveCell3DFacesRotationMatrix[f];
      const Eigen::Vector3d& faceTranslation = concaveCell3DFacesTranslation[f];
      const Eigen::Vector3d& faceNormal = concaveCell3DFacesNormal[f];

      int convexCellFound = -1, convexFaceFound = -1;

      // find convex cell3D face parallel to concave face normal
      for (unsigned int cc = 0; cc < numConvexCell3Ds; cc++)
      {
        const unsigned int convexCell3DIndex = convexCell3DIndices[cc];

        const unsigned int numConvexFaces = convexMesh.Cell3DNumberFaces(convexCell3DIndex);

        for (unsigned int ccf = 0; ccf < numConvexFaces; ccf++)
        {
          // verify if the convex face is in the same plane of the concave face
          if (!geometryUtilities.IsPolygonCoplanar(faceNormal,
                                                   faceOrigin,
                                                   convexCell3DsFaces3DVertices[cc][ccf],
                                                   convexCell3DsFacesUnalignedVertices[cc][ccf]))
            continue;

          const Eigen::MatrixXd& convexFaceVertices3D = convexCell3DsFaces3DVertices[cc][ccf];
          const std::vector<unsigned int>& convexFaceUnalignedVertices = convexCell3DsFacesUnalignedVertices[cc][ccf];

          Eigen::Matrix3d convexFaceTriangle;
          convexFaceTriangle.col(0)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[0]);
          convexFaceTriangle.col(1)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[1]);
          convexFaceTriangle.col(2)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[2]);

          // rotate convex face to 2D with concave face matrices
          Eigen::MatrixXd convexFace2DTriangle = geometryUtilities.RotatePointsFrom3DTo2D(convexFaceTriangle,
                                                                                          faceRotationMatrix.transpose(),
                                                                                          faceTranslation);

          if (geometryUtilities.IsValueNegative(faceRotationMatrix.determinant(),
                                                geometryUtilities.Tolerance1D()))
            convexFace2DTriangle.block(0, 1, 3, convexFace2DTriangle.cols() - 1).rowwise().reverseInPlace();

          // check if convex face point is inside concave cell
          if (!geometryUtilities.IsPointInsidePolygon_RayCasting(convexFace2DTriangle.col(0),
                                                                 face2DVertices))
            continue;
          if (!geometryUtilities.IsPointInsidePolygon_RayCasting(convexFace2DTriangle.col(1),
                                                                 face2DVertices))
            continue;
          if (!geometryUtilities.IsPointInsidePolygon_RayCasting(convexFace2DTriangle.col(2),
                                                                 face2DVertices))
            continue;

          // check if convex cell centroid is inside concave cell
          if (!geometryUtilities.IsPointInsidePolygon_RayCasting(geometryUtilities.PolygonBarycenter(convexFace2DTriangle),
                                                                 face2DVertices))
            continue;

          convexFaceFound = ccf;

          break;
        }

        if (convexFaceFound >= 0)
        {
          convexCellFound = cc;
          break;
        }
      }

      if (convexCellFound < 0 || convexFaceFound < 0)
      {
        std::cerr<< "Concave Cell3D "<< concaveCell3DIndex<< " face "<< f<< " Cell2D "<< mesh.Cell3DFace(concaveCell3DIndex, f)<< " ";
        std::cerr<< "convex Cell3D "<< convexCellFound<< " face "<< convexFaceFound<< " ";
        std::cerr<< "convex Cell2D "<< ((convexCellFound >= 0 && convexFaceFound >= 0) ? convexMesh.Cell3DFace(convexCell3DIndices[convexCellFound],
                                                                                                               convexCellFound) : 0)<< std::endl;
        throw runtime_error("Convex Cell 3D face not found");
      }

      convexCell2D.ConvexCell3DIndex = convexCellFound;
      convexCell2D.ConvexCell3DFaceIndex = convexFaceFound;
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::FindConcaveCell3DFacesConvexCell2DResult MeshUtilities::FindConcaveCell3DFacesConvexCell2D(const GeometryUtilities& geometryUtilities,
                                                                                                            const unsigned int& concaveCell3DIndex,
                                                                                                            const IMeshDAO& mesh,
                                                                                                            const std::vector<Eigen::MatrixXd>& concaveCell3DTetra,
                                                                                                            const std::vector<std::vector<Eigen::Matrix3d>>& concaveCell3D_faces_2D_triangles,
                                                                                                            const std::vector<Eigen::MatrixXd>& concaveCell3DFaces3DVertices,
                                                                                                            const std::vector<Eigen::Vector3d>& concaveCell3DFacesTranslation,
                                                                                                            const std::vector<Eigen::Matrix3d>& concaveCell3DFacesRotationMatrix,
                                                                                                            const std::vector<Eigen::Vector3d>& concaveCell3DFacesNormal,
                                                                                                            const std::vector<std::vector<Eigen::MatrixXd>>& convexCell3DsFaces3DVertices,
                                                                                                            const std::vector<std::vector<std::vector<unsigned int>>>& convexCell3DsFacesUnalignedVertices) const
  {
    FindConcaveCell3DFacesConvexCell2DResult result;

    const unsigned int& numConvexCell3Ds = concaveCell3DTetra.size();
    const unsigned int numConcaveFaces = mesh.Cell3DNumberFaces(concaveCell3DIndex);

    result.ConcaveCell3DFacesConvexCell2D.resize(numConcaveFaces);

    auto det_2D = [](const Eigen::Vector3d& p1,
                  const Eigen::Vector3d& p2,
                  const Eigen::Vector3d& p3)
    {
      return
          p1.x() * (p2.y() - p3.y()) +
          p2.x() * (p3.y() - p1.y()) +
          p3.x() * (p1.y() - p2.y());
    };


    for (unsigned int f = 0; f < numConcaveFaces; f++)
    {
      FindConcaveCell3DFacesConvexCell2DResult::ConvexCell2D& convexCell2D = result.ConcaveCell3DFacesConvexCell2D[f];

      const Eigen::MatrixXd& faceVertices = concaveCell3DFaces3DVertices[f];
      const Eigen::Vector3d& faceOrigin = faceVertices.col(0);
      const Eigen::Matrix3d& faceRotationMatrix = concaveCell3DFacesRotationMatrix[f];
      const Eigen::Vector3d& faceTranslation = concaveCell3DFacesTranslation[f];
      const Eigen::Vector3d& faceNormal = concaveCell3DFacesNormal[f];

      Eigen::Matrix3d face_2D_triangle = concaveCell3D_faces_2D_triangles.at(f).at(0);

      int convexCellFound = -1, convexFaceFound = -1;

      if (geometryUtilities.IsValueNegative(det_2D(face_2D_triangle.col(0),
                                                   face_2D_triangle.col(1),
                                                   face_2D_triangle.col(2)),
                                            geometryUtilities.Tolerance1D()))
      {
        throw std::runtime_error("Face triangle not unclockwise " + std::to_string(concaveCell3DIndex) + " " + std::to_string(f));
        // std::cout<< "WARNING: concave face triangle not unclockwise "<< concaveCell3DIndex<< " "<< f<< std::endl;
        //face_2D_triangle.block(0, 1, 3, face_2D_triangle.cols() - 1).rowwise().reverseInPlace();
      }

      // find convex cell3D face parallel to concave face normal
      for (unsigned int cc = 0; cc < numConvexCell3Ds; cc++)
      {
        for (unsigned int ccf = 0; ccf < 4; ++ccf)
        {
          // verify if the convex face is in the same plane of the concave face
          if (!geometryUtilities.IsPolygonCoplanar(faceNormal,
                                                   faceOrigin,
                                                   convexCell3DsFaces3DVertices[cc][ccf],
                                                   convexCell3DsFacesUnalignedVertices[cc][ccf]))
            continue;

          const Eigen::MatrixXd& convexFaceVertices3D = convexCell3DsFaces3DVertices[cc][ccf];
          const std::vector<unsigned int>& convexFaceUnalignedVertices = convexCell3DsFacesUnalignedVertices[cc][ccf];

          Eigen::Matrix3d convexFaceTriangle;
          convexFaceTriangle.col(0)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[0]);
          convexFaceTriangle.col(1)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[1]);
          convexFaceTriangle.col(2)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[2]);

          // rotate convex face to 2D with concave face matrices
          Eigen::MatrixXd convexFace2DTriangle = geometryUtilities.RotatePointsFrom3DTo2D(convexFaceTriangle,
                                                                                          faceRotationMatrix.transpose(),
                                                                                          faceTranslation);


          if (geometryUtilities.IsValueNegative(det_2D(convexFace2DTriangle.col(0),
                                                       convexFace2DTriangle.col(1),
                                                       convexFace2DTriangle.col(2)),
                                                geometryUtilities.Tolerance1D()))
          {
            // std::cout<< "WARNING: convex face triangle not unclockwise "<< cc<< " "<< ccf<< std::endl;
            convexFace2DTriangle.block(0, 1, 3, convexFace2DTriangle.cols() - 1).rowwise().reverseInPlace();
          }

          if (!geometryUtilities.CheckTrianglesIntersection(face_2D_triangle,
                                                            convexFace2DTriangle,
                                                            false))
            continue;

          convexFaceFound = ccf;

          break;
        }

        if (convexFaceFound >= 0)
        {
          convexCellFound = cc;
          break;
        }
      }

      if (convexCellFound < 0 || convexFaceFound < 0)
      {
        std::cerr<< "Concave Cell3D "<< concaveCell3DIndex<< " face "<< f<< " Cell2D "<< mesh.Cell3DFace(concaveCell3DIndex, f)<< " ";
        std::cerr<< "convex Cell3D "<< convexCellFound<< " face "<< convexFaceFound<< " ";
        throw runtime_error("Convex Cell 3D face not found");
      }

      convexCell2D.ConvexCell3DIndex = convexCellFound;
      convexCell2D.ConvexCell3DFaceIndex = convexFaceFound;
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::FindPointMeshPositionResult MeshUtilities::FindPointMeshPosition(const MeshUtilities::FindPointCell3DResult& find_point_cell3D_result,
                                                                                  const IMeshDAO& mesh) const
  {
    if (find_point_cell3D_result.Cell3Ds_found.empty())
      return {};

    FindPointMeshPositionResult result;
    result.MeshPositions.resize(find_point_cell3D_result.Cell3Ds_found.size());

    for (unsigned int p = 0; p < find_point_cell3D_result.Cell3Ds_found.size(); p++)
    {
      const auto& cell3D_found = find_point_cell3D_result.Cell3Ds_found.at(p);

      switch (cell3D_found.Cell3D_Position.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
          result.MeshPositions[p] =
          {
            FindPointMeshPositionResult::PointMeshPosition::Types::Cell3D,
            cell3D_found.Cell3D_index
          };
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
          result.MeshPositions[p] =
          {
            FindPointMeshPositionResult::PointMeshPosition::Types::Cell2D,
            mesh.Cell3DFace(cell3D_found.Cell3D_index,
            cell3D_found.Cell3D_Position.BorderIndex)
          };
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
          result.MeshPositions[p] =
          {
            FindPointMeshPositionResult::PointMeshPosition::Types::Cell1D,
            mesh.Cell3DEdge(cell3D_found.Cell3D_index,
            cell3D_found.Cell3D_Position.BorderIndex)
          };
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
          result.MeshPositions[p] =
          {
            FindPointMeshPositionResult::PointMeshPosition::Types::Cell0D,
            mesh.Cell3DVertex(cell3D_found.Cell3D_index,
            cell3D_found.Cell3D_Position.BorderIndex)
          };
          break;
        default:
          throw std::runtime_error("Unknown PointPolyhedronPositionResult");
      }
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::FindPointCell3DResult MeshUtilities::FindPointCell3D(const GeometryUtilities& geometryUtilities,
                                                                      const Eigen::Vector3d& point,
                                                                      const IMeshDAO& mesh,
                                                                      const std::vector<std::vector<Eigen::MatrixXi> >& cell3DsFaces,
                                                                      const std::vector<std::vector<Eigen::MatrixXd> >& cell3DsFaceVertices,
                                                                      const std::vector<std::vector<Eigen::MatrixXd> >& cell3DsFaceRotatedVertices,
                                                                      const std::vector<std::vector<Eigen::Vector3d> >& cell3DsFaceNormals,
                                                                      const std::vector<std::vector<bool> >& cell3DsFaceNormalDirections,
                                                                      const std::vector<std::vector<Eigen::Vector3d> >& cell3DsFaceTranslations,
                                                                      const std::vector<std::vector<Eigen::Matrix3d> >& cell3DsFaceRotationMatrices,
                                                                      const std::vector<Eigen::MatrixXd>& cell3DsBoundingBox,
                                                                      const bool find_only_first_cell3D,
                                                                      const unsigned int starting_cell3D_index) const
  {
    std::list<FindPointCell3DResult::PointCell3DFound> cell3Ds_found;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      unsigned int c3D_index = (starting_cell3D_index + c) % mesh.Cell3DTotalNumber();

      if (!mesh.Cell3DIsActive(c3D_index))
        continue;

      if (!geometryUtilities.IsPointInBoundingBox(point,
                                                  cell3DsBoundingBox.at(c3D_index)))
        continue;

      const GeometryUtilities::PointPolyhedronPositionResult pointPosition = geometryUtilities.PointPolyhedronPosition(point,
                                                                                                                       cell3DsFaces.at(c3D_index),
                                                                                                                       cell3DsFaceVertices.at(c3D_index),
                                                                                                                       cell3DsFaceRotatedVertices.at(c3D_index),
                                                                                                                       cell3DsFaceNormals.at(c3D_index),
                                                                                                                       cell3DsFaceNormalDirections.at(c3D_index),
                                                                                                                       cell3DsFaceTranslations.at(c3D_index),
                                                                                                                       cell3DsFaceRotationMatrices.at(c3D_index));

      switch (pointPosition.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
          cell3Ds_found.push_back({ c3D_index, pointPosition });
          break;
        default:
          throw std::runtime_error("Unknown point polyhedron position");
      }

      if (find_only_first_cell3D &&
          cell3Ds_found.size() > 0)
        break;
    }

    return { std::vector<FindPointCell3DResult::PointCell3DFound>(cell3Ds_found.begin(),
                                                                  cell3Ds_found.end()) };
  }
  // ***************************************************************************
  MeshUtilities::FindPointCell3DResult MeshUtilities::FindPointCell3D(const GeometryUtilities& geometryUtilities,
                                                                      const Eigen::Vector3d& point,
                                                                      const IMeshDAO& mesh,
                                                                      const std::vector<std::vector<Eigen::MatrixXi> >& cell3DsFaces,
                                                                      const std::vector<std::vector<Eigen::MatrixXd> >& cell3DsFaceVertices,
                                                                      const std::vector<std::vector<Eigen::MatrixXd> >& cell3DsFaceRotatedVertices,
                                                                      const std::vector<std::vector<Eigen::Vector3d> >& cell3DsFaceNormals,
                                                                      const std::vector<std::vector<bool> >& cell3DsFaceNormalDirections,
                                                                      const std::vector<std::vector<Eigen::Vector3d> >& cell3DsFaceTranslations,
                                                                      const std::vector<std::vector<Eigen::Matrix3d> >& cell3DsFaceRotationMatrices,
                                                                      const std::vector<Eigen::MatrixXd>& cell3DsBoundingBox,
                                                                      const std::vector<std::vector<Eigen::MatrixXd>>& cell3DsTetrahedra,
                                                                      const bool find_only_first_cell3D,
                                                                      const unsigned int starting_cell3D_index) const
  {
    std::list<FindPointCell3DResult::PointCell3DFound> cell3Ds_found;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
      unsigned int c3D_index = (starting_cell3D_index + c) % mesh.Cell3DTotalNumber();

      if (!mesh.Cell3DIsActive(c3D_index))
        continue;

      if (!geometryUtilities.IsPointInBoundingBox(point,
                                                  cell3DsBoundingBox.at(c3D_index)))
        continue;

      const GeometryUtilities::PointPolyhedronPositionResult pointPosition = geometryUtilities.PointPolyhedronPosition(point,
                                                                                                                       cell3DsFaces.at(c3D_index),
                                                                                                                       cell3DsFaceVertices.at(c3D_index),
                                                                                                                       cell3DsFaceRotatedVertices.at(c3D_index),
                                                                                                                       cell3DsFaceNormals.at(c3D_index),
                                                                                                                       cell3DsFaceNormalDirections.at(c3D_index),
                                                                                                                       cell3DsFaceTranslations.at(c3D_index),
                                                                                                                       cell3DsFaceRotationMatrices.at(c3D_index),
                                                                                                                       cell3DsTetrahedra.at(c3D_index));

      switch (pointPosition.Type)
      {
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Outside:
          break;
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex:
        case GeometryUtilities::PointPolyhedronPositionResult::Types::Inside:
          cell3Ds_found.push_back({ c3D_index, pointPosition });
          break;
        default:
          throw std::runtime_error("Unknown point polyhedron position");
      }

      if (find_only_first_cell3D &&
          cell3Ds_found.size() > 0)
        break;
    }

    return { std::vector<FindPointCell3DResult::PointCell3DFound>(cell3Ds_found.begin(),
                                                                  cell3Ds_found.end()) };
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell1DCell3DNeighbours(IMeshDAO& mesh) const
  {
    // Compute Cell1D neighbours starting from cell3Ds
    std::vector<std::list<unsigned int>> cell1DsNeighbours(mesh.Cell1DTotalNumber());
    for (unsigned int c3D = 0; c3D < mesh.Cell3DTotalNumber(); c3D++)
    {
      if (!mesh.Cell3DIsActive(c3D))
        continue;

      const unsigned int numCell3DEdges = mesh.Cell3DNumberEdges(c3D);
      for (unsigned int e = 0; e < numCell3DEdges; e++)
      {
        const unsigned int cell1D = mesh.Cell3DEdge(c3D, e);
        cell1DsNeighbours[cell1D].push_back(c3D);
      }
    }

    std::vector<unsigned int> cell1DsNumNeighbours3D(mesh.Cell1DTotalNumber());
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
      cell1DsNumNeighbours3D[c1D] = cell1DsNeighbours[c1D].size();

    mesh.Cell1DsInitializeNeighbourCell3Ds(cell1DsNumNeighbours3D);
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
    {
      unsigned int n = 0;
      for (const auto& cell3DIndex : cell1DsNeighbours[c1D])
        mesh.Cell1DInsertNeighbourCell3D(c1D,
                                         n++,
                                         cell3DIndex);
    }
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell2DCell3DNeighbours(IMeshDAO& mesh) const
  {
    // Compute Cell2D neighbours starting from cell3Ds
    std::vector<std::list<unsigned int>> cell2DsNeighbours(mesh.Cell2DTotalNumber());
    for (unsigned int c3D = 0; c3D < mesh.Cell3DTotalNumber(); c3D++)
    {
      if (!mesh.Cell3DIsActive(c3D))
        continue;

      const unsigned int numCell3DFaces = mesh.Cell3DNumberFaces(c3D);
      for (unsigned int f = 0; f < numCell3DFaces; f++)
      {
        const unsigned int cell2D = mesh.Cell3DFace(c3D, f);
        cell2DsNeighbours[cell2D].push_back(c3D);
      }
    }

    std::vector<unsigned int> cell2DsNumNeighbours3D(mesh.Cell2DTotalNumber());
    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
      cell2DsNumNeighbours3D[c2D] = cell2DsNeighbours[c2D].size();

    mesh.Cell2DsInitializeNeighbourCell3Ds(cell2DsNumNeighbours3D);
    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
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
  std::vector<unsigned int> MeshUtilities::SplitCell3D(const unsigned int& cell3DIndex,
                                                       const std::vector<std::vector<unsigned int>>& subCell3DsVertices,
                                                       const std::vector<std::vector<unsigned int>>& subCell3DsEdges,
                                                       const std::vector<std::vector<unsigned int>>& subCell3DsFaces,
                                                       IMeshDAO& mesh) const
  {
    Gedim::Output::Assert(subCell3DsVertices.size() == subCell3DsEdges.size() &&
                          subCell3DsVertices.size() == subCell3DsFaces.size());

    const unsigned int numSubCells = subCell3DsVertices.size();
    unsigned int newCell3DsStartingIndex = mesh.Cell3DAppend(numSubCells);

    vector<unsigned int> newCell3DsIndex(numSubCells);

    mesh.Cell3DSetState(cell3DIndex, false);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell3DsIndex[c] = newCell3DsStartingIndex + c;

      const unsigned int& newCell3DIndex = newCell3DsIndex[c];
      mesh.Cell3DAddVertices(newCell3DIndex,
                             subCell3DsVertices[c]);
      mesh.Cell3DAddEdges(newCell3DIndex,
                          subCell3DsEdges[c]);
      mesh.Cell3DAddFaces(newCell3DIndex,
                          subCell3DsFaces[c]);

      mesh.Cell3DSetMarker(newCell3DIndex,
                           mesh.Cell3DMarker(cell3DIndex));
      mesh.Cell3DSetState(newCell3DIndex,
                          true);

      mesh.Cell3DInsertUpdatedCell3D(cell3DIndex,
                                     newCell3DIndex);

      for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(newCell3DIndex); e++)
      {
        const unsigned int cell1DIndex = mesh.Cell3DEdge(newCell3DIndex, e);

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell3D(cell1DIndex); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, n))
            continue;

          if (mesh.Cell1DNeighbourCell3D(cell1DIndex, n) == cell3DIndex)
            mesh.Cell1DInsertNeighbourCell3D(cell1DIndex,
                                             n,
                                             newCell3DIndex);
        }
      }

      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(newCell3DIndex); f++)
      {
        const unsigned int cell2DIndex = mesh.Cell3DFace(newCell3DIndex, f);

        for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(cell2DIndex); n++)
        {
          if (!mesh.Cell2DHasNeighbourCell3D(cell2DIndex, n))
            continue;

          if (mesh.Cell2DNeighbourCell3D(cell2DIndex, n) == cell3DIndex)
            mesh.Cell2DInsertNeighbourCell3D(cell2DIndex,
                                             n,
                                             newCell3DIndex);
        }
      }
    }

    return newCell3DsIndex;
  }
  // ***************************************************************************
  MeshUtilities::AgglomerateCell3DInformation MeshUtilities::AgglomerateCell3Ds(const std::unordered_set<unsigned int>& cell3DsIndex,
                                                                                const IMeshDAO& mesh) const
  {
    AgglomerateCell3DInformation result;

    if (cell3DsIndex.size() < 2)
      return result;

    std::unordered_set<unsigned int> agglomerated_cell2Ds;
    std::unordered_set<unsigned int> removed_cell2Ds;

    for (const unsigned int c3D_index : cell3DsIndex)
    {

      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c3D_index); f++)
      {
        const unsigned int c2D_index = mesh.Cell3DFace(c3D_index, f);

        if (removed_cell2Ds.find(c2D_index) != removed_cell2Ds.end())
          continue;

        bool to_remove = true;

        unsigned int num_neigh_cell3Ds = 0;
        for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(c2D_index); n++)
        {
          if (!mesh.Cell2DHasNeighbourCell3D(c2D_index, n))
          {
            to_remove = false;
            break;
          }

          const unsigned int c3D_neigh_index = mesh.Cell2DNeighbourCell3D(c2D_index,
                                                                          n);

          if (cell3DsIndex.find(c3D_neigh_index) == cell3DsIndex.end())
          {
            to_remove = false;
            break;
          }

          num_neigh_cell3Ds++;
        }

        if (num_neigh_cell3Ds < 2)
          to_remove = false;

        if (to_remove)
        {
          if (removed_cell2Ds.find(c2D_index) == removed_cell2Ds.end())
            removed_cell2Ds.insert(c2D_index);

          continue;
        }

        if (agglomerated_cell2Ds.find(c2D_index) != agglomerated_cell2Ds.end())
          continue;

        agglomerated_cell2Ds.insert(c2D_index);
      }
    }

    std::unordered_set<unsigned int> agglomerated_cell1Ds;
    std::unordered_set<unsigned int> removed_cell1Ds;

    for (const unsigned int c3D_index : cell3DsIndex)
    {

      for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(c3D_index); e++)
      {
        const unsigned int c1D_index = mesh.Cell3DEdge(c3D_index, e);

        if (removed_cell1Ds.find(c1D_index) != removed_cell1Ds.end())
          continue;

        bool to_remove = true;

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(c1D_index); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(c1D_index, n))
          {
            to_remove = false;
            break;
          }

          const unsigned int c2D_neigh_index = mesh.Cell1DNeighbourCell2D(c1D_index,
                                                                          n);

          if (removed_cell2Ds.find(c2D_neigh_index) == removed_cell2Ds.end())
          {
            to_remove = false;
            break;
          }
        }

        if (to_remove)
        {
          if (removed_cell1Ds.find(c1D_index) == removed_cell1Ds.end())
            removed_cell1Ds.insert(c1D_index);

          continue;
        }

        if (agglomerated_cell1Ds.find(c1D_index) != agglomerated_cell1Ds.end())
          continue;

        agglomerated_cell1Ds.insert(c1D_index);
      }
    }

    std::unordered_set<unsigned int> agglomerated_cell0Ds;
    std::unordered_set<unsigned int> removed_cell0Ds;

    for (const unsigned int c3D_index : cell3DsIndex)
    {

      for (unsigned int v = 0; v < mesh.Cell3DNumberVertices(c3D_index); v++)
      {
        const unsigned int c0D_index = mesh.Cell3DVertex(c3D_index, v);

        if (removed_cell0Ds.find(c0D_index) != removed_cell0Ds.end())
          continue;

        bool to_remove = true;

        for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell1D(c0D_index); n++)
        {
          if (!mesh.Cell0DHasNeighbourCell1D(c0D_index, n))
          {
            to_remove = false;
            break;
          }

          const unsigned int c1D_neigh_index = mesh.Cell0DNeighbourCell1D(c0D_index,
                                                                          n);

          if (removed_cell1Ds.find(c1D_neigh_index) == removed_cell1Ds.end())
          {
            to_remove = false;
            break;
          }
        }

        if (to_remove)
        {
          if (removed_cell0Ds.find(c0D_index) == removed_cell0Ds.end())
            removed_cell0Ds.insert(c0D_index);

          continue;
        }

        if (agglomerated_cell0Ds.find(c0D_index) != agglomerated_cell0Ds.end())
          continue;

        agglomerated_cell0Ds.insert(c0D_index);
      }
    }

    result.AgglomerateCell3DVertices = std::vector<unsigned int>(agglomerated_cell0Ds.begin(),
                                                                 agglomerated_cell0Ds.end());
    result.AgglomerateCell3DEdges = std::vector<unsigned int>(agglomerated_cell1Ds.begin(),
                                                              agglomerated_cell1Ds.end());
    result.AgglomerateCell3DFaces = std::vector<unsigned int>(agglomerated_cell2Ds.begin(),
                                                              agglomerated_cell2Ds.end());
    result.SubCell3DsRemovedVertices = std::vector<unsigned int>(removed_cell0Ds.begin(),
                                                                 removed_cell0Ds.end());
    result.SubCell3DsRemovedEdges = std::vector<unsigned int>(removed_cell1Ds.begin(),
                                                              removed_cell1Ds.end());
    result.SubCell3DsRemovedFaces = std::vector<unsigned int>(removed_cell2Ds.begin(),
                                                              removed_cell2Ds.end());


    return result;
  }
  // ***************************************************************************
  unsigned int MeshUtilities::AgglomerateCell3Ds(const std::unordered_set<unsigned int>& subCell3DsIndex,
                                                 const std::vector<unsigned int>& agglomerateCell3DVertices,
                                                 const std::vector<unsigned int>& agglomerateCell3DEdges,
                                                 const std::vector<unsigned int>& agglomerateCell3DFaces,
                                                 const std::vector<unsigned int>& subCell3DsRemovedCell0Ds,
                                                 const std::vector<unsigned int>& subCell3DsRemovedCell1Ds,
                                                 const std::vector<unsigned int>& subCell3DsRemovedCell2Ds,
                                                 IMeshDAO& mesh,
                                                 std::vector<std::vector<unsigned int> >& meshCell3DsOriginalCell3Ds) const
  {
    if (subCell3DsIndex.size() == 0)
      return mesh.Cell3DTotalNumber();

    if (subCell3DsIndex.size() == 1)
      return *subCell3DsIndex.begin();

    const unsigned int agglomeratedCell3DIndex = mesh.Cell3DAppend(1);

    unsigned int max_marker = 0;
    std::list<unsigned int> agglomeratedCell3DConvexCells;
    for (const auto subCell3DIndex : subCell3DsIndex)
    {
      mesh.Cell3DSetState(subCell3DIndex, false);

      if (mesh.Cell3DMarker(subCell3DIndex) > max_marker)
        max_marker = mesh.Cell3DMarker(subCell3DIndex);

      const auto& convexCells = meshCell3DsOriginalCell3Ds.at(subCell3DIndex);
      for (const auto& convexCell : convexCells)
        agglomeratedCell3DConvexCells.push_back(convexCell);
    }
    meshCell3DsOriginalCell3Ds.resize(meshCell3DsOriginalCell3Ds.size() + 1);
    meshCell3DsOriginalCell3Ds.at(agglomeratedCell3DIndex) =
        std::vector<unsigned int>(agglomeratedCell3DConvexCells.begin(),
                                  agglomeratedCell3DConvexCells.end());

    mesh.Cell3DAddVertices(agglomeratedCell3DIndex,
                           agglomerateCell3DVertices);
    mesh.Cell3DAddEdges(agglomeratedCell3DIndex,
                        agglomerateCell3DEdges);
    mesh.Cell3DAddFaces(agglomeratedCell3DIndex,
                        agglomerateCell3DFaces);

    mesh.Cell3DSetMarker(agglomeratedCell3DIndex,
                         max_marker);
    mesh.Cell3DSetState(agglomeratedCell3DIndex,
                        true);

    for (unsigned int v = 0; v < mesh.Cell3DNumberVertices(agglomeratedCell3DIndex); v++)
    {
      const unsigned int cell0DIndex = mesh.Cell3DVertex(agglomeratedCell3DIndex, v);

      std::list<unsigned int> neighCell3Ds;

      for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell3D(cell0DIndex); n++)
      {
        if (!mesh.Cell0DHasNeighbourCell3D(cell0DIndex, n))
          continue;

        const unsigned int neighCell3DIndex = mesh.Cell0DNeighbourCell3D(cell0DIndex, n);

        if (subCell3DsIndex.find(neighCell3DIndex) == subCell3DsIndex.end())
          neighCell3Ds.push_back(neighCell3DIndex);
      }

      if (neighCell3Ds.empty())
        continue;

      neighCell3Ds.push_back(agglomeratedCell3DIndex);
      mesh.Cell0DInitializeNeighbourCell3Ds(cell0DIndex,
                                            std::vector<unsigned int>(neighCell3Ds.begin(),
                                                                      neighCell3Ds.end()));
    }

    for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(agglomeratedCell3DIndex); e++)
    {
      const unsigned int cell1DIndex = mesh.Cell3DEdge(agglomeratedCell3DIndex, e);

      std::list<unsigned int> neighCell3Ds;

      for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell3D(cell1DIndex); n++)
      {
        if (!mesh.Cell1DHasNeighbourCell3D(cell1DIndex, n))
          continue;

        const unsigned int neighCell3DIndex = mesh.Cell1DNeighbourCell3D(cell1DIndex, n);

        if (subCell3DsIndex.find(neighCell3DIndex) == subCell3DsIndex.end())
          neighCell3Ds.push_back(neighCell3DIndex);
      }

      if (neighCell3Ds.empty())
        continue;

      neighCell3Ds.push_back(agglomeratedCell3DIndex);
      mesh.Cell1DInitializeNeighbourCell3Ds(cell1DIndex,
                                            std::vector<unsigned int>(neighCell3Ds.begin(),
                                                                      neighCell3Ds.end()));
    }

    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(agglomeratedCell3DIndex); f++)
    {
      const unsigned int cell2DIndex = mesh.Cell3DFace(agglomeratedCell3DIndex, f);

      for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(cell2DIndex); n++)
      {
        if (!mesh.Cell2DHasNeighbourCell3D(cell2DIndex, n))
          continue;

        const unsigned int neighCell3DIndex = mesh.Cell2DNeighbourCell3D(cell2DIndex, n);

        if (subCell3DsIndex.find(neighCell3DIndex) != subCell3DsIndex.end())
        {
          mesh.Cell2DInsertNeighbourCell3D(cell2DIndex,
                                           n,
                                           agglomeratedCell3DIndex);
        }
      }
    }

    for (const auto cell0DIndex : subCell3DsRemovedCell0Ds)
      mesh.Cell0DSetState(cell0DIndex, false);
    for (const auto cell1DIndex : subCell3DsRemovedCell1Ds)
      mesh.Cell1DSetState(cell1DIndex, false);
    for (const auto cell2DIndex : subCell3DsRemovedCell2Ds)
      mesh.Cell2DSetState(cell2DIndex, false);

    return agglomeratedCell3DIndex;
  }
  // ***************************************************************************
  MeshUtilities::FilterMeshData MeshUtilities::FilterMesh3D(const std::vector<unsigned int>& cell3DsFilter,
                                                            const IMeshDAO& mesh) const
  {
    list<unsigned int> cell3Ds;
    std::set<unsigned int> cell0Ds, cell1Ds, cell2Ds;

    for (const unsigned int cell3DIndex : cell3DsFilter)
    {
      if (!mesh.Cell3DIsActive(cell3DIndex))
        continue;

      cell3Ds.push_back(cell3DIndex);

      for (unsigned int v = 0; v < mesh.Cell3DNumberVertices(cell3DIndex); v++)
      {
        const unsigned int cell0DIndex = mesh.Cell3DVertex(cell3DIndex, v);

        if (cell0Ds.find(cell0DIndex) == cell0Ds.end())
          cell0Ds.insert(cell0DIndex);
      }

      for (unsigned int v = 0; v < mesh.Cell3DNumberEdges(cell3DIndex); v++)
      {
        const unsigned int cell1DIndex = mesh.Cell3DEdge(cell3DIndex, v);

        if (cell1Ds.find(cell1DIndex) == cell1Ds.end())
          cell1Ds.insert(cell1DIndex);
      }

      for (unsigned int v = 0; v < mesh.Cell3DNumberFaces(cell3DIndex); v++)
      {
        const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex, v);

        if (cell2Ds.find(cell2DIndex) == cell2Ds.end())
          cell2Ds.insert(cell2DIndex);
      }
    }

    MeshUtilities::FilterMeshData result;

    result.Cell0Ds = std::vector<unsigned int>(cell0Ds.begin(),
                                               cell0Ds.end());
    result.Cell1Ds = std::vector<unsigned int>(cell1Ds.begin(),
                                               cell1Ds.end());
    result.Cell2Ds = std::vector<unsigned int>(cell2Ds.begin(),
                                               cell2Ds.end());
    result.Cell3Ds = std::vector<unsigned int>(cell3Ds.begin(),
                                               cell3Ds.end());

    return result;
  }
  // ***************************************************************************
  MeshUtilities::ExtractMeshData MeshUtilities::ExtractMesh3D(const std::vector<unsigned int>& cell0DsFilter,
                                                              const std::vector<unsigned int>& cell1DsFilter,
                                                              const std::vector<unsigned int>& cell2DsFilter,
                                                              const std::vector<unsigned int>& cell3DsFilter,
                                                              const IMeshDAO& originalMesh,
                                                              IMeshDAO& mesh) const
  {
    ExtractMeshData result;
    result.NewCell0DToOldCell0D.resize(cell0DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.NewCell1DToOldCell1D.resize(cell1DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.NewCell2DToOldCell2D.resize(cell2DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.NewCell3DToOldCell3D.resize(cell3DsFilter.size(), std::numeric_limits<unsigned int>::max());
    result.OldCell0DToNewCell0D.resize(originalMesh.Cell0DTotalNumber(), std::numeric_limits<unsigned int>::max());
    result.OldCell1DToNewCell1D.resize(originalMesh.Cell1DTotalNumber(), std::numeric_limits<unsigned int>::max());
    result.OldCell2DToNewCell2D.resize(originalMesh.Cell2DTotalNumber(), std::numeric_limits<unsigned int>::max());
    result.OldCell3DToNewCell3D.resize(originalMesh.Cell3DTotalNumber(), std::numeric_limits<unsigned int>::max());

    Eigen::MatrixXd newCell0Ds(3, cell0DsFilter.size());
    for (unsigned int v = 0; v < cell0DsFilter.size(); v++)
    {
      const unsigned int oldCell0DIndex = cell0DsFilter[v];
      result.NewCell0DToOldCell0D[v] = oldCell0DIndex;
      result.OldCell0DToNewCell0D[oldCell0DIndex] = v;

      newCell0Ds.col(v) = originalMesh.Cell0DCoordinates(oldCell0DIndex);
    }

    Eigen::MatrixXi newCell1Ds(2, cell1DsFilter.size());
    for (unsigned int e = 0; e < cell1DsFilter.size(); e++)
    {
      const unsigned int oldCell1DIndex = cell1DsFilter[e];
      result.NewCell1DToOldCell1D[e] = oldCell1DIndex;
      result.OldCell1DToNewCell1D[oldCell1DIndex] = e;

      const Eigen::VectorXi cell1DExtremes = originalMesh.Cell1DExtremes(oldCell1DIndex);

      newCell1Ds(0, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[0]);
      newCell1Ds(1, e) = result.OldCell0DToNewCell0D.at(cell1DExtremes[1]);
    }

    std::vector<Eigen::MatrixXi> newCell2Ds(cell2DsFilter.size());
    for (unsigned int p = 0; p < cell2DsFilter.size(); p++)
    {
      const unsigned int oldCell2DIndex = cell2DsFilter[p];
      result.NewCell2DToOldCell2D[p] = oldCell2DIndex;
      result.OldCell2DToNewCell2D[oldCell2DIndex] = p;

      const std::vector<unsigned int> vertices = originalMesh.Cell2DVertices(oldCell2DIndex);
      const std::vector<unsigned int> edges = originalMesh.Cell2DEdges(oldCell2DIndex);

      Gedim::Output::Assert(vertices.size() == edges.size());

      newCell2Ds[p].resize(2, vertices.size());
      for (unsigned int v = 0; v < vertices.size(); v++)
      {
        newCell2Ds[p](0, v) = result.OldCell0DToNewCell0D.at(vertices[v]);
        newCell2Ds[p](1, v) = result.OldCell1DToNewCell1D.at(edges[v]);
      }
    }

    std::vector<Mesh3DPolyhedron> newCell3Ds(cell3DsFilter.size());
    for (unsigned int p = 0; p < cell3DsFilter.size(); p++)
    {
      const unsigned int oldCell3DIndex = cell3DsFilter[p];
      result.NewCell3DToOldCell3D[p] = oldCell3DIndex;
      result.OldCell3DToNewCell3D[oldCell3DIndex] = p;

      const std::vector<unsigned int> vertices = originalMesh.Cell3DVertices(oldCell3DIndex);
      const std::vector<unsigned int> edges = originalMesh.Cell3DEdges(oldCell3DIndex);
      const std::vector<unsigned int> faces = originalMesh.Cell3DFaces(oldCell3DIndex);

      Mesh3DPolyhedron& polyhedron = newCell3Ds.at(p);
      polyhedron.VerticesIndex.resize(vertices.size());
      polyhedron.EdgesIndex.resize(edges.size());
      polyhedron.FacesIndex.resize(faces.size());

      for (unsigned int v = 0; v < vertices.size(); v++)
        polyhedron.VerticesIndex[v] = result.OldCell0DToNewCell0D.at(vertices[v]);
      for (unsigned int v = 0; v < edges.size(); v++)
        polyhedron.EdgesIndex[v] = result.OldCell1DToNewCell1D.at(edges[v]);
      for (unsigned int v = 0; v < faces.size(); v++)
        polyhedron.FacesIndex[v] = result.OldCell2DToNewCell2D.at(faces[v]);
    }

    FillMesh3D(newCell0Ds,
               newCell1Ds,
               newCell2Ds,
               newCell3Ds,
               mesh);

    for (unsigned int v = 0; v < cell0DsFilter.size(); v++)
    {
      const unsigned int oldCell0DIndex = result.NewCell0DToOldCell0D[v];
      mesh.Cell0DSetMarker(v, originalMesh.Cell0DMarker(oldCell0DIndex));
      mesh.Cell0DSetState(v, originalMesh.Cell0DIsActive(oldCell0DIndex));
    }

    for (unsigned int e = 0; e < cell1DsFilter.size(); e++)
    {
      const unsigned int oldCell1DIndex = result.NewCell1DToOldCell1D[e];
      mesh.Cell1DSetMarker(e, originalMesh.Cell1DMarker(oldCell1DIndex));
      mesh.Cell1DSetState(e, originalMesh.Cell1DIsActive(oldCell1DIndex));
    }

    for (unsigned int f = 0; f < cell2DsFilter.size(); f++)
    {
      const unsigned int oldCell2DIndex = result.NewCell2DToOldCell2D[f];
      mesh.Cell2DSetMarker(f, originalMesh.Cell2DMarker(oldCell2DIndex));
      mesh.Cell2DSetState(f, originalMesh.Cell2DIsActive(oldCell2DIndex));
    }

    for (unsigned int p = 0; p < cell3DsFilter.size(); p++)
    {
      const unsigned int oldCell3DIndex = result.NewCell3DToOldCell3D[p];
      mesh.Cell3DSetMarker(p, originalMesh.Cell3DMarker(oldCell3DIndex));
      mesh.Cell3DSetState(p, originalMesh.Cell3DIsActive(oldCell3DIndex));
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::MakeMeshTriangularFaces(const std::vector<std::vector<unsigned int>>& faces_triangulation,
                                              IMeshDAO& mesh) const
  {
    const unsigned int num_faces = mesh.Cell2DTotalNumber();
    std::vector<std::vector<unsigned int>> cell2Ds_new_faces(num_faces);
    std::vector<std::vector<unsigned int>> cell2Ds_new_edges(num_faces);

    for (unsigned int f = 0; f < num_faces; f++)
    {
      if (!mesh.Cell2DIsActive(f))
        continue;

      const unsigned int num_face_vertices = mesh.Cell2DNumberVertices(f);

      if (num_face_vertices == 3)
      {
        cell2Ds_new_faces[f].push_back(f);
        continue;
      }

      const auto& face_cell0Ds_idex = mesh.Cell2DVertices(f);
      const auto& face_cell1Ds_idex = mesh.Cell2DEdges(f);
      const auto& face_triangles = faces_triangulation.at(f);
      Gedim::Output::Assert(3 * (num_face_vertices - 2) == face_triangles.size());

      const unsigned int numTriangles = face_triangles.size() / 3;

      // add new edges
      const unsigned int new_cell1D_starting = mesh.Cell1DAppend(numTriangles - 1);
      cell2Ds_new_edges[f].resize(numTriangles - 1);

      for (unsigned int t = 0; t < numTriangles - 1; t++)
      {
        const unsigned int new_cell1D_index = new_cell1D_starting + t;

        cell2Ds_new_edges[f][t] = new_cell1D_index;
        mesh.Cell1DSetState(new_cell1D_index, true);
        mesh.Cell1DSetMarker(new_cell1D_index,
                             mesh.Cell2DMarker(f));
        mesh.Cell1DInsertExtremes(new_cell1D_index,
                                  face_cell0Ds_idex.at(face_triangles.at(3 * t + 1)),
                                  face_cell0Ds_idex.at(face_triangles.at(3 * t + 2)));
      }

      // add new faces
      const unsigned int new_cell2D_starting = mesh.Cell2DAppend(numTriangles);
      cell2Ds_new_faces[f].resize(numTriangles);

      for (unsigned int t = 0; t < numTriangles; t++)
      {
        const unsigned int new_cell2D_index = new_cell2D_starting + t;

        cell2Ds_new_faces[f][t] = new_cell2D_index;
        mesh.Cell2DSetMarker(new_cell2D_index,
                             mesh.Cell2DMarker(f));
        mesh.Cell2DSetState(new_cell2D_index, true);



        Eigen::MatrixXi vertices_edges(2, 3);
        vertices_edges.row(0)<<
                                face_cell0Ds_idex.at(face_triangles.at(3 * t)),
            face_cell0Ds_idex.at(face_triangles.at(3 * t + 1)),
            face_cell0Ds_idex.at(face_triangles.at(3 * t + 2));

        const unsigned int new_cell2D_first_edge =
            (t == 0) ? face_cell1Ds_idex.at(0) :
                       new_cell1D_starting + (t - 1);
        const unsigned int new_cell2D_third_edge =
            (t == numTriangles - 1) ? face_cell1Ds_idex.at(num_face_vertices - 1) :
                                      new_cell1D_starting + t;
        vertices_edges.row(1)<< new_cell2D_first_edge,
            face_cell1Ds_idex.at(t + 1),
            new_cell2D_third_edge;

        mesh.Cell2DAddVerticesAndEdges(new_cell2D_index,
                                       vertices_edges);
        mesh.Cell2DInsertUpdatedCell2D(f, new_cell2D_index);
      }

      mesh.Cell2DSetState(f, false);
    }

    const unsigned int num_cell3Ds = mesh.Cell3DTotalNumber();
    for (unsigned int c = 0; c < num_cell3Ds; c++)
    {
      if (!mesh.Cell3DIsActive(c))
        continue;

      const auto& vertices = mesh.Cell3DVertices(c);
      const auto& edges = mesh.Cell3DEdges(c);
      const auto& faces = mesh.Cell3DFaces(c);

      std::list<unsigned int> cell3D_new_edges;
      std::list<unsigned int> cell3D_new_faces;
      for (const auto& cell2D_index : faces)
      {
        for (const auto& cell1D_index : cell2Ds_new_edges.at(cell2D_index))
          cell3D_new_edges.push_back(cell1D_index);

        for (const auto& new_cell2D_index : cell2Ds_new_faces.at(cell2D_index))
          cell3D_new_faces.push_back(new_cell2D_index);
      }

      if (cell3D_new_faces.size() == faces.size())
        continue;

      const unsigned int new_cell3D_index = mesh.Cell3DAppend(1);
      mesh.Cell3DSetMarker(new_cell3D_index,
                           mesh.Cell3DMarker(c));
      mesh.Cell3DSetState(new_cell3D_index,
                          true);

      mesh.Cell3DInitializeVertices(new_cell3D_index,
                                    vertices.size());
      mesh.Cell3DInitializeEdges(new_cell3D_index,
                                 edges.size() + cell3D_new_edges.size());
      mesh.Cell3DInitializeFaces(new_cell3D_index,
                                 cell3D_new_faces.size());

      for (unsigned int v = 0; v < vertices.size(); v++)
        mesh.Cell3DInsertVertex(new_cell3D_index,
                                v,
                                vertices.at(v));

      for (unsigned int e = 0; e < edges.size(); e++)
        mesh.Cell3DInsertEdge(new_cell3D_index,
                              e,
                              edges.at(e));

      unsigned int e = 0;
      for (const auto& cell1D_index : cell3D_new_edges)
        mesh.Cell3DInsertEdge(new_cell3D_index,
                              edges.size() + e++,
                              cell1D_index);

      unsigned int f = 0;
      for (const unsigned int& cell2D_index : cell3D_new_faces)
        mesh.Cell3DInsertFace(new_cell3D_index,
                              f++,
                              cell2D_index);

      mesh.Cell3DSetState(c, false);
    }
  }
  // ***************************************************************************
  void MeshUtilities::ChangePolyhedronMeshMarkers(const GeometryUtilities& geometryUtilities,
                                                  const Eigen::MatrixXd& polyhedron_vertices,
                                                  const Eigen::MatrixXi& polyhedron_edges,
                                                  const std::vector<Eigen::MatrixXi>& polyhedron_faces,
                                                  const Eigen::MatrixXd& polyhedron_edges_tangent,
                                                  const Eigen::VectorXd& polyhedron_edges_length,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_normal,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_vertices,
                                                  const std::vector<Eigen::MatrixXd>& polyhedron_faces_vertices_2D,
                                                  const std::vector<Eigen::Vector3d>& polyhedron_faces_translation,
                                                  const std::vector<Eigen::Matrix3d>& polyhedron_faces_rotation_matrix,
                                                  const std::vector<unsigned int>& polyhedron_vertices_marker,
                                                  const std::vector<unsigned int>& polyhedron_edges_marker,
                                                  const std::vector<unsigned int>& polyhedron_faces_marker,
                                                  const std::vector<Eigen::Vector3d>& cell1Ds_centroid,
                                                  const std::vector<Eigen::Vector3d>& cell2Ds_centroid,
                                                  IMeshDAO& mesh) const
  {
    const unsigned int num_polyhedron_edges = polyhedron_edges.cols();
    const unsigned int num_polyhedron_faces = polyhedron_faces.size();

    for (unsigned int f = 0; f < num_polyhedron_faces; ++f)
    {
      const auto& face_vertices_2d = polyhedron_faces_vertices_2D.at(f);
      const auto& face_plane_normal = polyhedron_faces_normal.at(f);
      const Eigen::Vector3d face_plane_origin = polyhedron_faces_vertices.at(f).col(0);
      const auto& face_translation = polyhedron_faces_translation.at(f);
      const auto& face_rotation_matrix = polyhedron_faces_rotation_matrix.at(f);
      const unsigned int face_marker = polyhedron_faces_marker.at(f);

      SetMeshMarkersOnPolygon(geometryUtilities,
                              face_plane_normal,
                              face_plane_origin,
                              face_vertices_2d,
                              face_translation,
                              face_rotation_matrix,
                              cell1Ds_centroid,
                              cell2Ds_centroid,
                              face_marker,
                              mesh);
    }

    for (unsigned int e = 0; e < num_polyhedron_edges; ++e)
    {
      const auto edge_line_tangent = polyhedron_edges_tangent.col(e);
      const auto edge_length = polyhedron_edges_length[e];
      const Eigen::VectorXd edge_line_origin = polyhedron_vertices.col(polyhedron_edges(0, e));
      const unsigned int edge_marker = polyhedron_edges_marker.at(e);

      SetMeshMarkersOnSegment(geometryUtilities,
                              edge_line_origin,
                              edge_line_tangent,
                              edge_length * edge_length,
                              edge_marker,
                              mesh);
    }


    for (unsigned int c_0D = 0; c_0D < mesh.Cell0DTotalNumber(); ++c_0D)
    {
      const Eigen::Vector3d cell0D_coordinate = mesh.Cell0DCoordinates(c_0D);

      const auto vertex_found = geometryUtilities.FindPointInPoints(polyhedron_vertices,
                                                                    cell0D_coordinate);

      for (const auto& vertex_index : vertex_found)
      {
        mesh.Cell0DSetMarker(c_0D,
                             polyhedron_vertices_marker.at(vertex_index));
      }
    }
  }
  // ***************************************************************************
}
