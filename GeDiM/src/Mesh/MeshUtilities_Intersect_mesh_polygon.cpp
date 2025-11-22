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

#include "IntersectorMesh2DSegment.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
void Gedim::MeshUtilities::CreatePolygonIntersectionMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                                         const Eigen::MatrixXd &interface_vertices,
                                                         Gedim::IMeshDAO &mesh) const
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
        const auto domain_1D_tangent = geometry_utilities.SegmentTangent(interface_vertices.col(i), interface_vertices.col(next));
        const auto domain_1D_barycenter =
            geometry_utilities.SegmentBarycenter(interface_vertices.col(i), interface_vertices.col(next));

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

    for (auto &cell1D : cell1Ds_intersect)
    {
        if (!(cell1D.second.internal_edges.size() == 0 && cell1D.second.internal_points.size() == 1))
            continue;

        unsigned int id = mesh.Cell1DAppend(2);

        mesh.Cell1DInsertExtremes(id, mesh.Cell1DOrigin(cell1D.first), cell1D.second.internal_points[0]);
        mesh.Cell1DSetMarker(id, mesh.Cell1DMarker(cell1D.first));
        mesh.Cell1DSetState(id, true);
        cell1D.second.internal_edges.push_back(id);

        id++;
        mesh.Cell1DInsertExtremes(id, cell1D.second.internal_points[0], mesh.Cell1DEnd(cell1D.first));
        mesh.Cell1DSetMarker(id, mesh.Cell1DMarker(cell1D.first));
        mesh.Cell1DSetState(id, true);
        cell1D.second.internal_edges.push_back(id);
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
            new_vertices.push_back(last_vertex);

            if (last_vertex == origin)
            {
                unsigned int local_origin = last_vertex;
                unsigned int local_id = first_internal_edge;
                for (unsigned int e = 0; e < num_internal_edge; e++)
                {
                    new_edges.push_back(cell1D.second.internal_edges[local_id]);

                    if (e < num_internal_edge - 1)
                    {
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

                    if (e < num_internal_edge - 1)
                    {
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
        bool direction = edge_internal_2 == (edge_internal_1 + 1) % num_internal_edges;
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

            if (direction)
                e = (e + 1) % num_internal_edges;
            else
                e = (e == 0) ? num_internal_edges - 1 : (e - 1);
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

            if (!direction)
                e = (e + 1) % num_internal_edges;
            else
                e = (e == 0) ? num_internal_edges - 1 : (e - 1);
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
// ***************************************************************************
} // namespace Gedim
