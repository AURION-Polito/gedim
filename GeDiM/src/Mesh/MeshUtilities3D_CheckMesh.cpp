// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "MeshUtilities.hpp"

#include "IOStream.hpp"
#include <cassert>
#include <fstream>
#include <numeric>

namespace Gedim
{
// ***************************************************************************
void MeshUtilities::CheckMesh3D(const CheckMesh3DConfiguration &configuration,
                                const GeometryUtilities &geometryUtilities,
                                const Gedim::IMeshDAO &mesh) const
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

                if (geometryUtilities.PointsAreCoincident(mesh.Cell0DCoordinates(p1), mesh.Cell0DCoordinates(p2)))
                {
                    throw std::runtime_error("Cell0D " + std::to_string(p1) + " and " + std::to_string(p2) + " are coincident");
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

            if (mesh.Cell1DByExtremes(mesh.Cell1DOrigin(e1), mesh.Cell1DEnd(e1)) != e1)
            {
                throw std::runtime_error("Cell1D " + std::to_string(e1) + " has wrong extremes");
            }

            if (mesh.Cell1DByExtremes(mesh.Cell1DEnd(e1), mesh.Cell1DOrigin(e1)) != mesh.Cell1DTotalNumber())
            {
                throw std::runtime_error("Cell1D " + std::to_string(e1) + " has wrong extremes");
            }

            for (unsigned int e2 = e1 + 1; e2 < mesh.Cell1DTotalNumber(); e2++)
            {
                if (!mesh.Cell1DIsActive(e2))
                    continue;

                if ((mesh.Cell1DOrigin(e1) == mesh.Cell1DOrigin(e2) && mesh.Cell1DEnd(e1) == mesh.Cell1DEnd(e2)))
                {
                    throw std::runtime_error("Cell1D " + std::to_string(e1) + " and " + std::to_string(e2) + " are duplicated");
                }
                if ((mesh.Cell1DEnd(e1) == mesh.Cell1DOrigin(e2) && mesh.Cell1DOrigin(e1) == mesh.Cell1DEnd(e2)))
                {
                    throw std::runtime_error("Cell1D " + std::to_string(e1) + " and " + std::to_string(e2) + " are duplicated");
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
                    geometryUtilities.SegmentLength(mesh.Cell1DOriginCoordinates(e), mesh.Cell1DEndCoordinates(e)),
                    geometryUtilities.Tolerance1D()))
            {
                throw std::runtime_error("Cell1D " + std::to_string(e) + " is zero");
            }
        }
    }

    if (configuration.Cell2D_CheckEdges)
    {
        for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
        {
            if (!mesh.Cell2DIsActive(p))
                continue;

            const unsigned int cell2DNumVertices = mesh.Cell2DNumberVertices(p);
            const unsigned int cell2DNumEdges = mesh.Cell2DNumberEdges(p);

            for (unsigned int v = 0; v < cell2DNumEdges; v++)
            {
                const unsigned int eO = mesh.Cell2DVertex(p, v);
                const unsigned int eE = mesh.Cell2DVertex(p, (v + 1) % cell2DNumEdges);

                const unsigned int edgeFromVerticesOE = mesh.Cell2DFindEdgeByExtremes(p, eO, eE);
                const unsigned int edgeFromVerticesEO = mesh.Cell2DFindEdgeByExtremes(p, eE, eO);

                if (!((edgeFromVerticesOE < cell2DNumEdges && edgeFromVerticesOE == v) ||
                      (edgeFromVerticesEO < cell2DNumEdges && edgeFromVerticesEO == v)))
                {
                    throw std::runtime_error("Cell2D " + std::to_string(p) + " has wrong edge " + std::to_string(v));
                }

                if (mesh.Cell2DFindVertex(p, eO) == cell2DNumVertices)
                    throw std::runtime_error("Cell2D " + std::to_string(p) + " vertex " + std::to_string(eO) + " not found");

                if (mesh.Cell2DFindVertex(p, eE) == cell2DNumVertices)
                    throw std::runtime_error("Cell2D " + std::to_string(p) + " vertex " + std::to_string(eE) + " not found");
            }
        }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
        for (unsigned int p1 = 0; p1 < mesh.Cell2DTotalNumber(); p1++)
        {
            if (!mesh.Cell2DIsActive(p1))
                continue;

            std::vector<unsigned int> cell2D1Vertices = mesh.Cell2DVertices(p1);
            sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
            std::vector<unsigned int> cell2D1Edges = mesh.Cell2DEdges(p1);
            sort(cell2D1Edges.begin(), cell2D1Edges.end());

            for (unsigned int p2 = p1 + 1; p2 < mesh.Cell2DTotalNumber(); p2++)
            {
                if (!mesh.Cell2DIsActive(p2))
                    continue;

                std::vector<unsigned int> cell2D2Vertices = mesh.Cell2DVertices(p2);
                sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
                std::vector<unsigned int> cell2D2Edges = mesh.Cell2DEdges(p2);
                sort(cell2D2Edges.begin(), cell2D2Edges.end());

                if (!(cell2D1Vertices.size() != cell2D2Vertices.size() ||
                      !equal(cell2D1Vertices.begin(), cell2D1Vertices.end(), cell2D2Vertices.begin())))
                {
                    throw std::runtime_error("Cell2D " + std::to_string(p1) + " and " + std::to_string(p2) + " are duplicated");
                }

                if (!(cell2D1Edges.size() != cell2D2Edges.size() ||
                      !equal(cell2D1Edges.begin(), cell2D1Edges.end(), cell2D2Edges.begin())))
                {
                    throw std::runtime_error("Cell2D " + std::to_string(p1) + " and " + std::to_string(p2) + " are duplicated");
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
            const Eigen::Matrix3d cell2DRotationMatrix =
                geometryUtilities.PolygonRotationMatrix(cell2DVertices3D, cell2DNormal, cell2DTranslation);
            const Eigen::MatrixXd cell2DVertices2D =
                geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D, cell2DRotationMatrix.transpose(), cell2DTranslation);
            const std::vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices2D);
            const Eigen::MatrixXd convexCell2DUnalignedVertices =
                geometryUtilities.ExtractPoints(cell2DVertices2D, convexCell2DUnalignedVerticesFilter);
            const std::vector<unsigned int> convexHull = geometryUtilities.ConvexHull(convexCell2DUnalignedVertices);
            const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(convexCell2DUnalignedVertices, convexHull);

            if (!geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices, convexHullVertices))
            {
                throw std::runtime_error("Cell2D " + std::to_string(p) + " is not convex ");
            }

            if (geometryUtilities.PolygonOrientation(convexHull) == Gedim::GeometryUtilities::PolygonOrientations::Clockwise)
            {
                throw std::runtime_error("Cell2D " + std::to_string(p) + " is not CounterClockwise");
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
            const Eigen::Matrix3d cell2DRotationMatrix =
                geometryUtilities.PolygonRotationMatrix(cell2DVertices3D, cell2DNormal, cell2DTranslation);
            const Eigen::MatrixXd cell2DVertices2D =
                geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D, cell2DRotationMatrix.transpose(), cell2DTranslation);

            if (geometryUtilities.IsValueZero(geometryUtilities.PolygonArea(cell2DVertices2D), geometryUtilities.Tolerance2D()))
            {
                throw std::runtime_error("Cell2D " + std::to_string(p) + " is zero");
            }
        }
    }

    if (configuration.Cell3D_CheckDuplications)
    {
        for (unsigned int p1 = 0; p1 < mesh.Cell3DTotalNumber(); p1++)
        {
            if (!mesh.Cell3DIsActive(p1))
                continue;

            std::vector<unsigned int> cell3D1Vertices = mesh.Cell3DVertices(p1);
            sort(cell3D1Vertices.begin(), cell3D1Vertices.end());
            std::vector<unsigned int> cell3D1Edges = mesh.Cell3DEdges(p1);
            sort(cell3D1Edges.begin(), cell3D1Edges.end());
            std::vector<unsigned int> cell3D1Faces = mesh.Cell3DFaces(p1);
            sort(cell3D1Faces.begin(), cell3D1Faces.end());

            for (unsigned int p2 = p1 + 1; p2 < mesh.Cell3DTotalNumber(); p2++)
            {
                if (!mesh.Cell3DIsActive(p2))
                    continue;

                std::vector<unsigned int> cell3D2Vertices = mesh.Cell3DVertices(p2);
                sort(cell3D2Vertices.begin(), cell3D2Vertices.end());
                std::vector<unsigned int> cell3D2Edges = mesh.Cell3DEdges(p2);
                sort(cell3D2Edges.begin(), cell3D2Edges.end());
                std::vector<unsigned int> cell3D2Faces = mesh.Cell3DFaces(p2);
                sort(cell3D2Faces.begin(), cell3D2Faces.end());

                if (!(cell3D1Vertices.size() != cell3D2Vertices.size() ||
                      !equal(cell3D1Vertices.begin(), cell3D1Vertices.end(), cell3D2Vertices.begin())))
                {
                    throw std::runtime_error("Cell3D " + std::to_string(p1) + " and " + std::to_string(p2) + " are duplicated");
                }

                if (!(cell3D1Edges.size() != cell3D2Edges.size() ||
                      !equal(cell3D1Edges.begin(), cell3D1Edges.end(), cell3D2Edges.begin())))
                {
                    throw std::runtime_error("Cell3D " + std::to_string(p1) + " and " + std::to_string(p2) + " are duplicated");
                }

                if (!(cell3D1Faces.size() != cell3D2Faces.size() ||
                      !equal(cell3D1Faces.begin(), cell3D1Faces.end(), cell3D2Faces.begin())))
                {
                    throw std::runtime_error("Cell3D " + std::to_string(p1) + " and " + std::to_string(p2) + " are duplicated");
                }
            }
        }
    }

    if (configuration.Cell3D_CheckEdgesAreActive)
    {
        for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
        {
            if (!mesh.Cell3DIsActive(p))
                continue;

            const unsigned int cell3DNumEdges = mesh.Cell3DNumberEdges(p);
            for (unsigned int e = 0; e < cell3DNumEdges; e++)
            {
                if (!mesh.Cell1DIsActive(mesh.Cell3DEdge(p, e)))
                {
                    throw std::runtime_error("Cell3D " + std::to_string(p) + " has edge " + std::to_string(e) + " inactive");
                }
            }
        }
    }

    if (configuration.Cell3D_CheckEdges)
    {
        for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
        {
            if (!mesh.Cell3DIsActive(p))
                continue;

            const unsigned int c_num_edges = mesh.Cell3DNumberEdges(p);
            for (unsigned int e = 0; e < c_num_edges; ++e)
            {
                const unsigned int c_1D = mesh.Cell3DEdge(p, e);

                const unsigned int c_1D_origin = mesh.Cell1DOrigin(c_1D);
                if (mesh.Cell3DFindVertex(p, c_1D_origin) == c_num_edges)
                    throw std::runtime_error("Cell3D " + std::to_string(p) + " vertex " + std::to_string(c_1D_origin) + " not found");

                const unsigned int c_1D_end = mesh.Cell1DEnd(c_1D);
                if (mesh.Cell3DFindVertex(p, c_1D_end) == c_num_edges)
                    throw std::runtime_error("Cell3D " + std::to_string(p) + " vertex " + std::to_string(c_1D_end) + " not found");
            }

            std::vector<unsigned int> cell3DEdges = mesh.Cell3DEdges(p);
            sort(cell3DEdges.begin(), cell3DEdges.end());

            const unsigned int cell3DNumFaces = mesh.Cell3DNumberFaces(p);
            std::set<unsigned int> cell3DFaceEdges;
            for (unsigned int f = 0; f < cell3DNumFaces; f++)
            {
                const unsigned int faceIndex = mesh.Cell3DFace(p, f);
                for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(faceIndex); e++)
                    cell3DFaceEdges.insert(mesh.Cell2DEdge(faceIndex, e));
            }

            std::vector<unsigned int> commonEdges;
            std::set_intersection(cell3DEdges.begin(),
                                  cell3DEdges.end(),
                                  cell3DFaceEdges.begin(),
                                  cell3DFaceEdges.end(),
                                  std::back_inserter(commonEdges));

            if (commonEdges.size() != cell3DEdges.size() || commonEdges.size() != cell3DFaceEdges.size())
            {
                throw std::runtime_error("Cell3D " + std::to_string(p) + " has wrong edges");
            }
        }
    }

    if (configuration.Cell3D_CheckConvexity)
    {
        for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
        {
            if (!mesh.Cell3DIsActive(p))
                continue;

            GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh, p);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const std::vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const std::vector<Eigen::Vector3d> polyhedronFaceBarycenters =
                geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

            const std::vector<Eigen::Vector3d> polyhedronFaceNormals =
                geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const std::vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const std::vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const std::vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const std::vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            const GeometryUtilities::PointPolyhedronPositionResult polyhedronBarycenterPosition =
                geometryUtilities.PointPolyhedronPosition(polyhedronBarycenter,
                                                          polyhedron.Faces,
                                                          polyhedronFace3DVertices,
                                                          polyhedronFace2DVertices,
                                                          polyhedronFaceNormals,
                                                          polyhedronFaceNormalDirections,
                                                          polyhedronFaceTranslations,
                                                          polyhedronFaceRotationMatrices);

            Output::Assert(polyhedronBarycenterPosition.Type == GeometryUtilities::PointPolyhedronPositionResult::Types::Inside);

            if (!geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
                                                      polyhedronFaceBarycenters,
                                                      polyhedronFaceNormals,
                                                      polyhedronFaceNormalDirections,
                                                      polyhedronBarycenter))
            {
                throw std::runtime_error("Cell3D " + std::to_string(p) + " is not convex");
            }
        }
    }

    if (configuration.Cell3D_CheckMeasure)
    {
        for (unsigned int p = 0; p < mesh.Cell3DTotalNumber(); p++)
        {
            if (!mesh.Cell3DIsActive(p))
                continue;

            GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(mesh, p);

            const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
            const std::vector<Eigen::MatrixXd> polyhedronFace3DVertices =
                geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
            const std::vector<std::vector<unsigned int>> polyhedronFaceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, polyhedronFace3DVertices);

            const std::vector<Eigen::Vector3d> polyhedronFaceNormals =
                geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
            const std::vector<bool> polyhedronFaceNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices, polyhedronBarycenter, polyhedronFaceNormals);
            const std::vector<Eigen::Vector3d> polyhedronFaceTranslations =
                geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
            const std::vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices, polyhedronFaceNormals, polyhedronFaceTranslations);

            const std::vector<Eigen::MatrixXd> polyhedronFace2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices, polyhedronFaceTranslations, polyhedronFaceRotationMatrices);

            const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(polyhedronFace2DVertices, polyhedronFaceTriangulations);

            if (geometryUtilities.IsValueZero(geometryUtilities.PolyhedronVolumeByBoundaryIntegral(polyhedronFace2DTriangulationPoints,
                                                                                                   polyhedronFaceNormals,
                                                                                                   polyhedronFaceNormalDirections,
                                                                                                   polyhedronFaceTranslations,
                                                                                                   polyhedronFaceRotationMatrices),
                                              geometryUtilities.Tolerance3D()))
            {
                throw std::runtime_error("Cell3D " + std::to_string(p) + " is zero");
            }
        }
    }

    if (configuration.Cell3D_CheckFaces)
    {
        for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
        {
            if (!mesh.Cell3DIsActive(c))
                continue;

            for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c); ++f)
            {
                const unsigned int c2D_index = mesh.Cell3DFace(c, f);

                if (!mesh.Cell2DIsActive(c2D_index))
                    throw std::runtime_error("Cell3D " + std::to_string(c) + " face " + std::to_string(c2D_index) + " is not active");

                const unsigned int c_n_vertices = mesh.Cell3DNumberVertices(c);
                for (unsigned int f_v = 0; f_v < mesh.Cell2DNumberVertices(c2D_index); ++f_v)
                {
                    const unsigned int c0D_index = mesh.Cell2DVertex(c2D_index, f_v);

                    if (mesh.Cell3DFindVertex(c, c0D_index) == c_n_vertices)
                        throw std::runtime_error("Cell3D " + std::to_string(c) + " vertex " + std::to_string(c0D_index) + " not found");
                }

                const unsigned int c_n_edges = mesh.Cell3DNumberEdges(c);
                for (unsigned int f_e = 0; f_e < mesh.Cell2DNumberEdges(c2D_index); ++f_e)
                {
                    const unsigned int c1D_index = mesh.Cell2DEdge(c2D_index, f_e);

                    if (mesh.Cell3DFindEdge(c, c1D_index) == c_n_edges)
                        throw std::runtime_error("Cell3D " + std::to_string(c) + " edge " + std::to_string(c1D_index) + " not found");
                }
            }
        }
    }
}
// ***************************************************************************
} // namespace Gedim
