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

#ifndef __test_mesh_bulk_face_H
#define __test_mesh_bulk_face_H

#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "SphereMeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numeric>

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestBulkFaceMesh, TestCreateCircle)
{

    std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateCircle";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    {
        const Eigen::Vector3d center = Eigen::Vector3d::Zero();
        const double diameter = 1.0;
        const unsigned int num_points = 4;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

        vector<unsigned int> vertexMarkers(vertices.cols());
        vector<unsigned int> edgeMarkers(vertices.cols());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), vertices.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        mesh_utilities.Mesh2DFromPolygon(vertices, vertexMarkers, edgeMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolygon(vertices);
            vtkUtilities.Export(exportFolder + "/unit_circle.vtu", Gedim::VTKUtilities::Ascii);
        }
    }

    {
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        center << 1.5, 2.4, 0.0;
        const double diameter = 2.0;
        const unsigned int num_points = 30;
        Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

        vector<unsigned int> vertexMarkers(vertices.cols());
        vector<unsigned int> edgeMarkers(vertices.cols());

        std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
        std::iota(edgeMarkers.begin(), edgeMarkers.end(), vertices.cols() + 1);

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        mesh_utilities.Mesh2DFromPolygon(vertices, vertexMarkers, edgeMarkers, mesh);

        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

        // Export to VTK
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolygon(vertices);
            vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
        }
    }
}

TEST(TestBulkFaceMesh, TestCreateTriangularMeshCircle)
{

    std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateTriangularMeshCircle";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const double diameter = 1.0;
    const unsigned int num_points = 100;
    Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

#if ENABLE_TRIANGLE == 0
    throw std::runtime_error("Triangle library not active");
#endif
    const double max_relative_area = 0.001;
    const double max_cell_area = 2.0 * std::numbers::pi * max_relative_area;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.CreateTriangularMesh(vertices, max_cell_area, mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "CircleDomain");
}

void CreateProjectedMesh(const Gedim::GeometryUtilities &geometryUtilities,
                         const Eigen::MatrixXd &boundary_vertices,
                         const Gedim::IMeshDAO &original_mesh,
                         Gedim::IMeshDAO &mesh)
{
    mesh.InitializeDimension(2);
    const unsigned int marker = 1;

    std::map<unsigned int, unsigned int> old_to_new_vertices;
    std::map<unsigned int, unsigned int> old_to_new_edges;

    std::vector<Gedim::GeometryUtilities::PointPolygonPositionResult> result(original_mesh.Cell0DTotalNumber());
    for (unsigned int c = 0; c < original_mesh.Cell0DTotalNumber(); c++)
    {
        result[c] = geometryUtilities.PointPolygonPosition(original_mesh.Cell0DCoordinates(c), boundary_vertices);

        switch (result[c].Type)
        {
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside: {
            const unsigned int new_id = mesh.Cell0DAppend(1);
            mesh.Cell0DInsertCoordinates(new_id, original_mesh.Cell0DCoordinates(c));
            mesh.Cell0DSetMarker(new_id, 0);
            mesh.Cell0DSetState(new_id, true);

            old_to_new_vertices.insert({c, new_id});
        }
        break;
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge:
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex: {
            const unsigned int new_id = mesh.Cell0DAppend(1);
            mesh.Cell0DInsertCoordinates(new_id, original_mesh.Cell0DCoordinates(c));
            mesh.Cell0DSetMarker(new_id, marker);
            mesh.Cell0DSetState(new_id, true);

            old_to_new_vertices.insert({c, new_id});
        }
        break;
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside:
            break;
        case Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Unknown:
        default:
            throw std::runtime_error("not valid position");
        }
    }

    const unsigned int num_boundary_vertices = boundary_vertices.cols();
    for (unsigned int c = 0; c < original_mesh.Cell1DTotalNumber(); c++)
    {
        const unsigned int id_origin = original_mesh.Cell1DOrigin(c);
        const unsigned int id_end = original_mesh.Cell1DEnd(c);

        if (old_to_new_vertices.find(id_origin) != old_to_new_vertices.end() &&
            old_to_new_vertices.find(id_end) != old_to_new_vertices.end())
            continue;

        if (result[id_origin].Type == Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside &&
            result[id_end].Type == Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside)
        {
            const Eigen::Vector3d origin = original_mesh.Cell0DCoordinates(id_origin);
            const Eigen::Vector3d end = original_mesh.Cell0DCoordinates(id_end);

            bool found_edge = false;
            for (unsigned int vd = 0; vd < num_boundary_vertices; vd++)
            {
                Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                    geometryUtilities.IntersectionSegmentSegment(origin,
                                                                 end,
                                                                 boundary_vertices.col(vd),
                                                                 boundary_vertices.col((vd + 1) % num_boundary_vertices));

                switch (result.IntersectionSegmentsType)
                {
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection: {
                    const unsigned int new_end = mesh.Cell0DAppend(1);
                    Eigen::Vector3d new_end_coordinates =
                        origin + result.FirstSegmentIntersections.at(0).CurvilinearCoordinate * (end - origin);
                    mesh.Cell0DInsertCoordinates(new_end, new_end_coordinates);
                    mesh.Cell0DSetMarker(new_end, marker);
                    mesh.Cell0DSetState(new_end, true);

                    old_to_new_vertices.insert({id_end, new_end});

                    found_edge = true;
                }
                break;
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection:
                    break;
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections:
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::Unknown:
                default:
                    throw std::runtime_error("not valid case");
                }

                if (found_edge)
                    break;
            }
        }
        else if (result[id_origin].Type == Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside &&
                 result[id_end].Type == Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside)
        {
            const Eigen::Vector3d origin = original_mesh.Cell0DCoordinates(id_origin);
            const Eigen::Vector3d end = original_mesh.Cell0DCoordinates(id_end);

            bool found_edge = false;
            for (unsigned int vd = 0; vd < num_boundary_vertices; vd++)
            {
                Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                    geometryUtilities.IntersectionSegmentSegment(origin,
                                                                 end,
                                                                 boundary_vertices.col(vd),
                                                                 boundary_vertices.col((vd + 1) % num_boundary_vertices));

                switch (result.IntersectionSegmentsType)
                {
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection: {
                    const unsigned int new_origin = mesh.Cell0DAppend(1);
                    Eigen::Vector3d new_origin_coordinates =
                        origin + result.FirstSegmentIntersections.at(0).CurvilinearCoordinate * (end - origin);
                    mesh.Cell0DInsertCoordinates(new_origin, new_origin_coordinates);
                    mesh.Cell0DSetMarker(new_origin, marker);
                    mesh.Cell0DSetState(new_origin, true);

                    old_to_new_vertices.insert({id_origin, new_origin});

                    found_edge = true;
                }
                break;
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection:
                    break;
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections:
                case Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::Unknown:
                default:
                    throw std::runtime_error("not valid case");
                }

                if (found_edge)
                    break;
            }
        }
        else
        {
            // non lo inserisco
        }
    }

    for (unsigned int c = 0; c < original_mesh.Cell1DTotalNumber(); c++)
    {
        const unsigned int id_origin = original_mesh.Cell1DOrigin(c);
        const unsigned int id_end = original_mesh.Cell1DEnd(c);

        const bool id_origin_exist = !(old_to_new_vertices.find(id_origin) == old_to_new_vertices.end());
        const bool id_end_exist = !(old_to_new_vertices.find(id_end) == old_to_new_vertices.end());

        if (id_origin_exist && id_end_exist)
        {
            const unsigned int new_id = mesh.Cell1DAppend(1);
            const unsigned int new_origin = old_to_new_vertices[id_origin];
            const unsigned int new_end = old_to_new_vertices[id_end];
            mesh.Cell1DInsertExtremes(new_id, new_origin, new_end);
            mesh.Cell1DSetMarker(new_id, std::max(mesh.Cell0DMarker(new_end), mesh.Cell0DMarker(new_origin)));
            mesh.Cell1DSetState(new_id, true);

            old_to_new_edges.insert({c, new_id});
        }
    }

    for (unsigned int c = 0; c < original_mesh.Cell2DTotalNumber(); c++)
    {
        unsigned int count = 0;
        const unsigned int num_vertices = original_mesh.Cell2DNumberVertices(c);
        for (unsigned int v = 0; v < num_vertices; v++)
        {
            const unsigned int id_vertex = original_mesh.Cell2DVertex(c, v);
            if (old_to_new_vertices.find(id_vertex) != old_to_new_vertices.end())
                count++;
        }

        if (count == num_vertices)
        {
            const unsigned int new_id = mesh.Cell2DAppend(1);
            mesh.Cell2DInitializeVertices(new_id, num_vertices);
            mesh.Cell2DInitializeEdges(new_id, num_vertices);

            for (unsigned int v = 0; v < num_vertices; v++)
            {
                mesh.Cell2DInsertVertex(new_id, v, old_to_new_vertices[original_mesh.Cell2DVertex(c, v)]);
                mesh.Cell2DInsertEdge(new_id, v, old_to_new_edges[original_mesh.Cell2DEdge(c, v)]);
            }
            mesh.Cell2DSetMarker(new_id, 0);
            mesh.Cell2DSetState(new_id, true);
        }
        else if (count >= 3)
        {
            const unsigned int new_id = mesh.Cell2DAppend(1);
            mesh.Cell2DInitializeVertices(new_id, count);
            mesh.Cell2DInitializeEdges(new_id, count);

            unsigned int id = 0;
            for (unsigned int v = 0; v < num_vertices; v++)
            {
                const unsigned int id_vertex = original_mesh.Cell2DVertex(c, v);

                if (old_to_new_vertices.find(id_vertex) == old_to_new_vertices.end())
                    continue;

                unsigned int id_next = original_mesh.Cell2DVertex(c, (v + 1) % num_vertices);
                if (old_to_new_vertices.find(id_next) != old_to_new_vertices.end())
                {
                    mesh.Cell2DInsertVertex(new_id, id, old_to_new_vertices[id_vertex]);
                    mesh.Cell2DInsertEdge(new_id, id, old_to_new_edges[original_mesh.Cell2DEdge(c, v)]);
                    id++;
                }
                else
                {
                    do
                    {
                        v++;
                        id_next = original_mesh.Cell2DVertex(c, (v + 1) % num_vertices);
                    } while (old_to_new_vertices.find(id_next) == old_to_new_vertices.end());

                    const unsigned int new_id_edge = mesh.Cell1DAppend(1);
                    mesh.Cell1DInsertExtremes(new_id_edge, old_to_new_vertices[id_vertex], old_to_new_vertices[id_next]);
                    mesh.Cell1DSetMarker(new_id_edge, 1);
                    mesh.Cell1DSetState(new_id_edge, true);

                    mesh.Cell2DInsertVertex(new_id, id, old_to_new_vertices[id_vertex]);
                    mesh.Cell2DInsertEdge(new_id, id, new_id_edge);
                    id++;
                }
            }
        }
    }
}

TEST(TestBulkFaceMesh, TestCreateProjectedCartesianMeshCircle)
{

    std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateProjectedCartesianMeshCircle";
    Gedim::Output::CreateFolder(exportFolder);

    const Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::MeshUtilities mesh_utilities;
    const Gedim::SphereMeshUtilities sphere_mesh_utilities(geometryUtilities, mesh_utilities);

    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const double diameter = 1.0;
    const unsigned int num_points = 100;
    Eigen::MatrixXd vertices = sphere_mesh_utilities.circle(center, diameter, num_points);

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    const Eigen::MatrixXd bounding_box = geometryUtilities.PointsBoundingBox(vertices);
    const Eigen::Vector3d first_vertex = (Eigen::Vector3d() << bounding_box(0, 0), bounding_box(1, 0), 0.0).finished();
    const Eigen::Vector3d second_vertex = (Eigen::Vector3d() << bounding_box(0, 1), bounding_box(1, 0), 0.0).finished();
    const Eigen::Vector3d fourth_vertex = (Eigen::Vector3d() << bounding_box(0, 0), bounding_box(1, 1), 0.0).finished();

    const vector<double> baseMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(8, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(8, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(first_vertex,
                                       second_vertex - first_vertex,
                                       fourth_vertex - first_vertex,
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       original_mesh);

    mesh_utilities.ExportMeshToVTU(original_mesh, exportFolder, "original_mesh");

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    CreateProjectedMesh(geometryUtilities, vertices, original_mesh, mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    mesh_utilities.CheckMesh2D(config, geometryUtilities, mesh);
}

void CreateIntersectionMesh(const Gedim::GeometryUtilities &geometry_utilities,
                            const Eigen::MatrixXd &interface_vertices,
                            Gedim::IMeshDAO &mesh)
{
    const unsigned int num_domains_1D = interface_vertices.cols();

    struct cell2Ds_intersected
    {
        std::vector<unsigned int> internal_points;
        std::vector<unsigned int> internal_edges;
    };

    struct cell2Ds_updated
    {
        std::vector<unsigned int> vertices;
        std::vector<unsigned int> edges;
    };

    struct cell1Ds_intersected
    {
        std::vector<unsigned int> internal_points;
        std::vector<unsigned int> internal_edges;
    };

    std::map<unsigned int, cell2Ds_updated> cell2Ds_update;
    std::map<unsigned int, cell2Ds_intersected> cell2Ds_intersect;
    std::map<unsigned int, cell1Ds_intersected> cell1Ds_intersect;

    for (unsigned int i = 0; i < num_domains_1D; ++i)
    {
        const unsigned int next = (i + 1) % num_domains_1D;
        const auto domain_1D_length = geometry_utilities.SegmentLength(interface_vertices.col(i), interface_vertices.col(next));
        const auto domain_1D_tangent = geometry_utilities.SegmentTangent(interface_vertices.col(0), interface_vertices.col(next));
        const auto domain_1D_barycenter =
            geometry_utilities.SegmentBarycenter(interface_vertices.col(0), interface_vertices.col(next));

        Gedim::IntersectorMesh2DSegment intersector(mesh, geometry_utilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh intersection_mesh;
        intersector.CreateIntersectionMesh(interface_vertices.col(i),
                                           interface_vertices.col(next),
                                           domain_1D_tangent,
                                           domain_1D_barycenter,
                                           domain_1D_length,
                                           intersection_mesh);

        std::map<double, unsigned int> id_points;
        for (const auto &point : intersection_mesh.Points)
        {
            if (point.second.Cell0DIds.size() == 0)
            {
                const unsigned int id = mesh.Cell0DAppend(1);
                Eigen::Vector3d coordinates;
                for (unsigned int d = 0; d < 3; d++)
                    coordinates(d) = interface_vertices(d, i) + point.first * domain_1D_tangent[d];
                mesh.Cell0DInsertCoordinates(id, coordinates);
                mesh.Cell0DSetMarker(id, 0);
                mesh.Cell0DSetState(id, true);
                id_points.insert({point.first, id});

                if (point.second.Cell2DIds.size() == 1)
                {
                    const unsigned int cell2D = *point.second.Cell2DIds.begin();
                    if (cell2Ds_intersect.find(cell2D) == cell2Ds_intersect.end())
                    {
                        cell2Ds_intersected obj;
                        obj.internal_points = {id};
                        cell2Ds_intersect.insert({cell2D, obj});
                    }
                    else
                        cell2Ds_intersect[cell2D].internal_points.push_back(id);
                }

                if (point.second.Cell1DIds.size() == 1)
                {
                    const unsigned int cell1D = *point.second.Cell1DIds.begin();
                    if (cell1Ds_intersect.find(cell1D) == cell1Ds_intersect.end())
                    {
                        cell1Ds_intersected obj;
                        obj.internal_points = {id};
                        cell1Ds_intersect.insert({cell1D, obj});
                    }
                    else
                        cell1Ds_intersect[cell1D].internal_points.push_back(id);
                }
            }
            else if (point.second.Cell0DIds.size() == 1)
            {
                const unsigned int id = *point.second.Cell0DIds.begin();
                id_points.insert({point.first, id});
            }
            else
                throw std::runtime_error("not valid original mesh");
        }

        for (const auto &segment : intersection_mesh.Segments)
        {
            if (segment.Cell2DIds.size() == 1)
            {
                const unsigned int id = mesh.Cell1DAppend(1);
                mesh.Cell1DInsertExtremes(id, id_points[segment.Points[0]], id_points[segment.Points[1]]);
                mesh.Cell1DSetMarker(id, 0);
                mesh.Cell1DSetState(id, true);

                const unsigned int cell2D = *segment.Cell2DIds.begin();
                if (cell2Ds_intersect.find(cell2D) == cell2Ds_intersect.end())
                {
                    cell2Ds_intersected obj;
                    obj.internal_edges = {id};
                    cell2Ds_intersect.insert({cell2D, obj});
                }
                else
                    cell2Ds_intersect[cell2D].internal_edges.push_back(id);
            }
            else if (segment.Cell1DIds.size() == 1)
            {
                const unsigned int cell1D = *segment.Cell1DIds.begin();

                if ((mesh.Cell1DOrigin(cell1D) == id_points[segment.Points[0]] &&
                     mesh.Cell1DEnd(cell1D) == id_points[segment.Points[1]]) ||
                    (mesh.Cell1DOrigin(cell1D) == id_points[segment.Points[1]] &&
                     mesh.Cell1DEnd(cell1D) == id_points[segment.Points[0]]))
                    continue;

                const unsigned int id = mesh.Cell1DAppend(1);
                mesh.Cell1DInsertExtremes(id, id_points[segment.Points[0]], id_points[segment.Points[1]]);
                mesh.Cell1DSetMarker(id, mesh.Cell1DMarker(cell1D));
                mesh.Cell1DSetState(id, true);

                if (cell1Ds_intersect.find(cell1D) == cell1Ds_intersect.end())
                {
                    cell1Ds_intersected obj;
                    obj.internal_edges = {id};
                    cell1Ds_intersect.insert({cell1D, obj});
                }
                else
                    cell1Ds_intersect[cell1D].internal_edges.push_back(id);
            }
            else
                throw std::runtime_error("not valid original mesh");
        }
    }

    Gedim::MeshUtilities mesh_utilities;
    for (const auto &cell1D : cell1Ds_intersect)
    {
        if (cell1D.second.internal_edges.size() < 2)
            throw std::runtime_error("not valid intersection");

        const unsigned int origin = mesh.Cell1DOrigin(cell1D.first);
        unsigned int first_internal_edge = std::numeric_limits<unsigned int>::max();
        const unsigned int num_internal_edge = cell1D.second.internal_edges.size();
        for (unsigned int e = 0; e < num_internal_edge; e++)
        {
            const unsigned int ee = cell1D.second.internal_edges[e];
            if (mesh.Cell1DOrigin(ee) == origin || mesh.Cell1DEnd(ee) == origin)
            {
                first_internal_edge = e;
                break;
            }
        }

        if (first_internal_edge == std::numeric_limits<unsigned int>::max())
            throw std::runtime_error("not valid intersection");

        const auto neigh_cell2D = mesh.Cell1DNeighbourCell2Ds(cell1D.first);

        for (unsigned int n = 0; n < neigh_cell2D.size(); n++)
        {
            if (neigh_cell2D[n] == std::numeric_limits<unsigned int>::max())
                continue;

            std::vector<unsigned int> old_vertices;
            std::vector<unsigned int> old_edges;

            if (cell2Ds_update.find(neigh_cell2D[n]) == cell2Ds_update.end())
            {
                old_vertices = mesh.Cell2DVertices(neigh_cell2D[n]);
                old_edges = mesh.Cell2DEdges(neigh_cell2D[n]);
            }
            else
            {
                old_vertices = cell2Ds_update[neigh_cell2D[n]].vertices;
                old_edges = cell2Ds_update[neigh_cell2D[n]].edges;
            }

            const unsigned int num_vertices = old_vertices.size();
            unsigned int id_local = std::numeric_limits<unsigned int>::max();
            for (unsigned int e = 0; e < num_vertices; e++)
            {
                if (old_edges[e] == cell1D.first)
                {
                    id_local = e;
                    break;
                }
            }

            if (id_local == std::numeric_limits<unsigned int>::max())
                throw std::runtime_error("not valid mesh");

            std::vector<unsigned int> new_vertices;
            std::vector<unsigned int> new_edges;
            for (unsigned int e = 0; e < id_local; e++)
            {
                new_edges.push_back(old_edges[e]);
                new_vertices.push_back(old_vertices[e]);
            }

            const unsigned int last_vertex = old_vertices[id_local];

            if (last_vertex == origin)
            {
                unsigned int local_origin = last_vertex;
                unsigned int local_id = first_internal_edge;
                for (unsigned int e = 0; e < num_internal_edge; e++)
                {
                    new_edges.push_back(cell1D.second.internal_edges[local_id]);

                    if (mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]) == local_origin)
                    {
                        new_vertices.push_back(mesh.Cell1DEnd(cell1D.second.internal_edges[local_id]));
                        local_origin = mesh.Cell1DEnd(cell1D.second.internal_edges[local_id]);
                    }
                    else
                    {
                        new_vertices.push_back(mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]));
                        local_origin = mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]);
                    }

                    local_id = (local_id + 1) % num_internal_edge;
                }
            }
            else if (last_vertex == mesh.Cell1DEnd(cell1D.first))
            {
                unsigned int local_origin = last_vertex;
                unsigned int local_id = (first_internal_edge == 0) ? num_internal_edge - 1 : first_internal_edge - 1;
                for (unsigned int e = 0; e < num_internal_edge; e++)
                {
                    new_edges.push_back(cell1D.second.internal_edges[local_id]);

                    if (mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]) == local_origin)
                    {
                        new_vertices.push_back(mesh.Cell1DEnd(cell1D.second.internal_edges[local_id]));
                        local_origin = mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]);
                    }
                    else
                    {
                        new_vertices.push_back(mesh.Cell1DOrigin(cell1D.second.internal_edges[local_id]));
                        local_origin = mesh.Cell1DEnd(cell1D.second.internal_edges[local_id]);
                    }

                    local_id = (local_id == 0) ? num_internal_edge - 1 : local_id - 1;
                }
            }
            else
                throw std::runtime_error("not valid intersection");

            for (unsigned int e = id_local + 1; e < num_vertices; e++)
            {
                new_edges.push_back(old_edges[e]);
                new_vertices.push_back(old_vertices[e]);
            }

            if (cell2Ds_update.find(neigh_cell2D[n]) == cell2Ds_update.end())
            {
                cell2Ds_updated obj;
                obj.edges = new_edges;
                obj.vertices = new_vertices;
                cell2Ds_update.insert({neigh_cell2D[n], obj});
            }
            else
            {
                cell2Ds_update[neigh_cell2D[n]].edges = new_edges;
                cell2Ds_update[neigh_cell2D[n]].vertices = new_vertices;
            }
        }

        mesh.Cell1DSetState(cell1D.first, false);
    }

    mesh_utilities.ComputeCell0DCell1DNeighbours(mesh);
    for (const auto &cell2D : cell2Ds_intersect)
    {
        std::vector<unsigned int> old_vertices;
        std::vector<unsigned int> old_edges;

        if (cell2Ds_update.find(cell2D.first) == cell2Ds_update.end())
        {
            old_vertices = mesh.Cell2DVertices(cell2D.first);
            old_edges = mesh.Cell2DEdges(cell2D.first);
        }
        else
        {
            old_vertices = cell2Ds_update[cell2D.first].vertices;
            old_edges = cell2Ds_update[cell2D.first].edges;
        }

        std::vector<unsigned int> vertices_1;
        std::vector<unsigned int> edges_1;

        std::vector<unsigned int> vertices_2;
        std::vector<unsigned int> edges_2;

        const unsigned int num_vertices = old_vertices.size();
        unsigned int intersect_point_1 = std::numeric_limits<unsigned int>::max();
        unsigned int edge_internal_1 = std::numeric_limits<unsigned int>::max();
        bool found = false;
        for (unsigned int v = 0; v < num_vertices; v++)
        {
            unsigned int id_vert = old_vertices.at(v);
            std::vector<unsigned int> neigh_edges = mesh.Cell0DNeighbourCell1Ds(id_vert);

            for (unsigned int n = 0; n < neigh_edges.size(); n++)
            {
                const auto it = find(cell2D.second.internal_edges.begin(), cell2D.second.internal_edges.end(), neigh_edges[n]);
                if (it != cell2D.second.internal_edges.end())
                {
                    found = true;
                    edge_internal_1 = it - cell2D.second.internal_edges.begin();
                    break;
                }
            }

            if (found)
            {
                intersect_point_1 = v;
                break;
            }
        }

        if (!found || intersect_point_1 == std::numeric_limits<unsigned int>::max())
            throw std::runtime_error("not valid intersections");

        unsigned int intersect_point_2 = std::numeric_limits<unsigned int>::max();
        unsigned int edge_internal_2 = std::numeric_limits<unsigned int>::max();
        found = false;
        for (unsigned int v = intersect_point_1 + 1; v < num_vertices; v++)
        {
            unsigned int id_vert = old_vertices.at(v);
            std::vector<unsigned int> neigh_edges = mesh.Cell0DNeighbourCell1Ds(id_vert);

            for (unsigned int n = 0; n < neigh_edges.size(); n++)
            {
                const auto it = find(cell2D.second.internal_edges.begin(), cell2D.second.internal_edges.end(), neigh_edges[n]);
                if (it != cell2D.second.internal_edges.end())
                {
                    found = true;
                    edge_internal_2 = it - cell2D.second.internal_edges.begin();
                    break;
                }
            }

            if (found)
            {
                intersect_point_2 = v;
                break;
            }
        }

        if (!found || intersect_point_2 == std::numeric_limits<unsigned int>::max())
            throw std::runtime_error("not valid intersections");

        unsigned int v = intersect_point_1;
        while (v != intersect_point_2)
        {
            const unsigned int id_vert = old_vertices.at(v);
            vertices_1.push_back(id_vert);

            const unsigned int id_edge = old_edges.at(v);
            edges_1.push_back(id_edge);

            v = (v + 1) % num_vertices;
        }
        unsigned int id_vert = old_vertices.at(v);
        vertices_1.push_back(id_vert);

        unsigned int origin = id_vert;
        const unsigned int num_internal_edges = cell2D.second.internal_edges.size();
        unsigned int e = edge_internal_2;
        for (unsigned int c = 0; c < num_internal_edges; c++)
        {
            const unsigned int id_edge = cell2D.second.internal_edges[e];
            edges_1.push_back(id_edge);

            unsigned int end = 0;
            if (mesh.Cell1DOrigin(id_edge) == origin)
                end = mesh.Cell1DEnd(id_edge);
            else
                end = mesh.Cell1DOrigin(id_edge);

            if (end == vertices_1[0] && c != num_internal_edges - 1)
                throw std::runtime_error("not valid");

            if (c < num_internal_edges - 1)
                vertices_1.push_back(end);

            origin = end;

            e = (e + 1) % num_internal_edges;
        }

        while (v != intersect_point_1)
        {
            const unsigned int id_vert = old_vertices.at(v);
            vertices_2.push_back(id_vert);

            const unsigned int id_edge = old_edges.at(v);
            edges_2.push_back(id_edge);

            v = (v + 1) % num_vertices;
        }

        id_vert = old_vertices.at(v);
        vertices_2.push_back(id_vert);

        origin = id_vert;
        e = edge_internal_1;
        for (unsigned int c = 0; c < num_internal_edges; c++)
        {
            const unsigned int id_edge = cell2D.second.internal_edges[e];
            edges_2.push_back(id_edge);

            unsigned int end = 0;
            if (mesh.Cell1DOrigin(id_edge) == origin)
                end = mesh.Cell1DEnd(id_edge);
            else
                end = mesh.Cell1DOrigin(id_edge);

            if (end == vertices_2[0] && c != num_internal_edges - 1)
                throw std::runtime_error("not valid");

            if (c < num_internal_edges - 1)
                vertices_2.push_back(end);

            origin = end;

            e = (e + 1) % num_internal_edges;
        }

        mesh.Cell2DSetState(cell2D.first, false);

        unsigned int id = mesh.Cell2DAppend(2);

        mesh.Cell2DInitializeEdges(id, edges_1.size());
        mesh.Cell2DInitializeVertices(id, vertices_1.size());
        mesh.Cell2DSetMarker(id, mesh.Cell2DMarker(cell2D.first));
        mesh.Cell2DSetState(id, true);
        mesh.Cell2DInsertVertices(id, vertices_1);
        mesh.Cell2DInsertEdges(id, edges_1);

        id++;
        mesh.Cell2DInitializeEdges(id, edges_2.size());
        mesh.Cell2DInitializeVertices(id, vertices_2.size());
        mesh.Cell2DSetMarker(id, mesh.Cell2DMarker(cell2D.first));
        mesh.Cell2DSetState(id, true);
        mesh.Cell2DInsertVertices(id, vertices_2);
        mesh.Cell2DInsertEdges(id, edges_2);
    }

    for (const auto &cell2D : cell2Ds_update)
    {
        if (cell2Ds_intersect.find(cell2D.first) != cell2Ds_intersect.end())
            continue;

        unsigned int id = mesh.Cell2DAppend(1);

        mesh.Cell2DInitializeEdges(id, cell2D.second.edges.size());
        mesh.Cell2DInitializeVertices(id, cell2D.second.vertices.size());
        mesh.Cell2DSetMarker(id, mesh.Cell2DMarker(cell2D.first));
        mesh.Cell2DSetState(id, true);
        mesh.Cell2DInsertVertices(id, cell2D.second.vertices);
        mesh.Cell2DInsertEdges(id, cell2D.second.edges);

        mesh.Cell2DSetState(cell2D.first, false);
    }
}

// TEST(TestBulkFaceMesh, TestCreateIntersectionMesh)
// {
//     Gedim::GeometryUtilitiesConfig geometry_utilities_config;
//     geometry_utilities_config.Tolerance1D = 1e-8;
//     Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

//     const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh";
//     Gedim::Output::CreateFolder(exportFolder);

//     const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
//     const double domain_2D_area = 1.0;

//     {
//         std::vector<double> domain_2D_id(1, 1.0);
//         Gedim::VTKUtilities exporter;
//         exporter.AddPolygon(domain_2D,
//                             {{"domain_id",
//                               Gedim::VTPProperty::Formats::Cells,
//                               static_cast<unsigned int>(domain_2D_id.size()),
//                               domain_2D_id.data()}});
//         exporter.Export(exportFolder + "/domain_2D.vtu");
//     }

//     const unsigned int num_domains_1D = 5;
//     Eigen::MatrixXd vertices(3, num_domains_1D);
//     vertices.col(0) << 1.0 / 2.0, 1.0 / 4.0, 0.0;
//     vertices.col(1) << 5.0 / 8.0, 3.0 / 8.0, 0.0;
//     vertices.col(2) << 3.0 / 4.0, 1.0 / 2.0, 0.0;
//     vertices.col(3) << 1.0 / 2.0, 3.0 / 4.0, 0.0;
//     vertices.col(4) << 1.0 / 4.0, 1.0 / 2.0, 0.0;

//     {
//         Gedim::VTKUtilities vtkUtilities;
//         vtkUtilities.AddPolygon(vertices);
//         vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
//     }

//     Gedim::MeshUtilities mesh_utilities;
//     // working meshes
//     // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
//     const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

//     Gedim::MeshMatrices domain_2D_mesh_data;
//     Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

//     const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0,
//     true); const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1,
//     0.0, 1.0, true);

//     Gedim::MeshMatrices original_mesh_data;
//     Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
//     mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
//                                        domain_2D.col(1) - domain_2D.col(0),
//                                        domain_2D.col(3) - domain_2D.col(0),
//                                        baseMeshCurvilinearCoordinates,
//                                        heightMeshCurvilinearCoordinates,
//                                        domain_2D_mesh);

//     // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
//     mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

//     {
//         mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
//     }

//     CreateIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

//     Gedim::MeshMatrices mesh_data;
//     Gedim::MeshMatricesDAO mesh(mesh_data);

//     const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
//     mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh,
//     mesh);

//     mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

//     {
//         mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
//     }

//     Gedim::MeshUtilities::CheckMesh2DConfiguration config;
//     mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
// }

// TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_1)
// {
//     Gedim::GeometryUtilitiesConfig geometry_utilities_config;
//     geometry_utilities_config.Tolerance1D = 1e-8;
//     Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

//     const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_1";
//     Gedim::Output::CreateFolder(exportFolder);

//     const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
//     const double domain_2D_area = 1.0;

//     {
//         std::vector<double> domain_2D_id(1, 1.0);
//         Gedim::VTKUtilities exporter;
//         exporter.AddPolygon(domain_2D,
//                             {{"domain_id",
//                               Gedim::VTPProperty::Formats::Cells,
//                               static_cast<unsigned int>(domain_2D_id.size()),
//                               domain_2D_id.data()}});
//         exporter.Export(exportFolder + "/domain_2D.vtu");
//     }

//     const unsigned int num_domains_1D = 5;
//     Eigen::MatrixXd vertices(3, num_domains_1D);
//     vertices.col(0) << 5.0 / 8.0, 3.0 / 8.0, 0.0;
//     vertices.col(1) << 3.0 / 4.0, 1.0 / 2.0, 0.0;
//     vertices.col(2) << 1.0 / 2.0, 3.0 / 4.0, 0.0;
//     vertices.col(3) << 1.0 / 4.0, 1.0 / 2.0, 0.0;
//     vertices.col(4) << 1.0 / 2.0, 1.0 / 4.0, 0.0;

//     {
//         Gedim::VTKUtilities vtkUtilities;
//         vtkUtilities.AddPolygon(vertices);
//         vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
//     }

//     Gedim::MeshUtilities mesh_utilities;
//     // working meshes
//     // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
//     const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

//     Gedim::MeshMatrices domain_2D_mesh_data;
//     Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

//     const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0,
//     true); const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1,
//     0.0, 1.0, true);

//     Gedim::MeshMatrices original_mesh_data;
//     Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
//     mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
//                                        domain_2D.col(1) - domain_2D.col(0),
//                                        domain_2D.col(3) - domain_2D.col(0),
//                                        baseMeshCurvilinearCoordinates,
//                                        heightMeshCurvilinearCoordinates,
//                                        domain_2D_mesh);

//     // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
//     mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

//     {
//         mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
//     }

//     CreateIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

//     Gedim::MeshMatrices mesh_data;
//     Gedim::MeshMatricesDAO mesh(mesh_data);

//     const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
//     mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh,
//     mesh);

//     mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

//     {
//         mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
//     }

//     Gedim::MeshUtilities::CheckMesh2DConfiguration config;
//     mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
// }

TEST(TestBulkFaceMesh, TestCreateIntersectionMesh_2)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "./Export/TestBulkFaceMesh/TestCreateIntersectionMesh_2";
    Gedim::Output::CreateFolder(exportFolder);

    const auto domain_2D = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(), 1.0);
    const double domain_2D_area = 1.0;

    {
        std::vector<double> domain_2D_id(1, 1.0);
        Gedim::VTKUtilities exporter;
        exporter.AddPolygon(domain_2D,
                            {{"domain_id",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(domain_2D_id.size()),
                              domain_2D_id.data()}});
        exporter.Export(exportFolder + "/domain_2D.vtu");
    }

    const unsigned int num_domains_1D = 5;
    Eigen::MatrixXd vertices(3, num_domains_1D);
    vertices.col(0) << 1.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(1) << 3.0 / 8.0, 1.0 / 4.0, 0.0;
    vertices.col(2) << 3.0 / 4.0, 1.0 / 4.0, 0.0;
    vertices.col(3) << 3.0 / 4.0, 3.0 / 4.0, 0.0;
    vertices.col(4) << 1.0 / 4.0, 3.0 / 4.0, 0.0;

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(vertices);
        vtkUtilities.Export(exportFolder + "/circle.vtu", Gedim::VTKUtilities::Ascii);
    }

    Gedim::MeshUtilities mesh_utilities;
    // working meshes
    // 0.05, 0.01, 0.005, 0.001 (a = 0.25 OK, a = 0.1 NO)
    const double domain_2D_max_area = 0.005 * domain_2D_area; // 0.005

    Gedim::MeshMatrices domain_2D_mesh_data;
    Gedim::MeshMatricesDAO domain_2D_mesh(domain_2D_mesh_data);

    const vector<double> baseMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);
    const vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4 + 1, 0.0, 1.0, true);

    Gedim::MeshMatrices original_mesh_data;
    Gedim::MeshMatricesDAO original_mesh(original_mesh_data);
    mesh_utilities.CreateRectangleMesh(domain_2D.col(0),
                                       domain_2D.col(1) - domain_2D.col(0),
                                       domain_2D.col(3) - domain_2D.col(0),
                                       baseMeshCurvilinearCoordinates,
                                       heightMeshCurvilinearCoordinates,
                                       domain_2D_mesh);

    // mesh_utilities.CreateTriangularMesh(domain_2D, domain_2D_max_area, domain_2D_mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(domain_2D_mesh);

    {
        mesh_utilities.ExportMeshToVTU(domain_2D_mesh, exportFolder, "original_mesh");
    }

    CreateIntersectionMesh(geometry_utilities, vertices, domain_2D_mesh);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    const auto filter_data = mesh_utilities.FilterActiveMesh(domain_2D_mesh);
    mesh_utilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, domain_2D_mesh, mesh);

    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    {
        mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "mesh");
    }

    // Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    // mesh_utilities.CheckMesh2D(config, geometry_utilities, mesh);
}

} // namespace GedimUnitTesting
#endif
