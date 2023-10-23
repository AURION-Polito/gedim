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
    result.Cell3DsFaces3DTriangulations.resize(convexMesh.Cell3DTotalNumber());
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
    result.Cell3DsFaces3DTriangulations.resize(mesh.Cell3DTotalNumber());
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

        if (geometryUtilities.PolygonOrientation(faceConvexHull) ==
            Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
        {
          face2DCCWOrientation[f] = false;
          face2DsCCW[f] = geometryUtilities.ChangePolygonOrientation(result.Cell3DsFaces2DVertices[c][f].cols());
          face2DVerticesCCW[f] = geometryUtilities.ExtractPoints(result.Cell3DsFaces2DVertices[c][f],
                                                                 face2DsCCW[f]);
        }
        else
        {
          face2DCCWOrientation[f] = true;
          face2DVerticesCCW[f] = result.Cell3DsFaces2DVertices[c][f];

          face2DsCCW[f].resize(result.Cell3DsFaces2DVertices[c][f].cols());
          std::iota(face2DsCCW[f].begin(), face2DsCCW[f].end(), 0);
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

          Output::Assert(geometryUtilities.PolygonOrientation(faceConvexHull) ==
                         Gedim::GeometryUtilities::PolygonOrientations::CounterClockwise);
        }

        convexCell3DsFacesUnalignedVertices[cc].resize(convexCell3DFaces2DVertices.size());
        for (unsigned int ccf = 0; ccf < convexCell3DsFacesUnalignedVertices[cc].size(); ccf++)
          convexCell3DsFacesUnalignedVertices[cc][ccf] = geometryUtilities.UnalignedPoints(convexCell3DFaces2DVertices[ccf]);

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
          const Eigen::MatrixXd& convexFaceVertices3D = convexCell3DsFaces3DVertices[cc][ccf];
          const std::vector<unsigned int>& convexFaceUnalignedVertices = convexCell3DsFacesUnalignedVertices[cc][ccf];

          Eigen::Matrix3d convexFaceTriangle;
          convexFaceTriangle.col(0)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[0]);
          convexFaceTriangle.col(1)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[1]);
          convexFaceTriangle.col(2)<< convexFaceVertices3D.col(convexFaceUnalignedVertices[2]);

          // verify if the convex face is in the same plane of the concave face
          if (!geometryUtilities.IsPointOnPlane(convexFaceTriangle.col(0),
                                                faceNormal,
                                                faceOrigin))
            continue;

          if (!geometryUtilities.IsPointOnPlane(convexFaceTriangle.col(1),
                                                faceNormal,
                                                faceOrigin))
            continue;

          if (!geometryUtilities.IsPointOnPlane(convexFaceTriangle.col(2),
                                                faceNormal,
                                                faceOrigin))
            continue;

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
}
