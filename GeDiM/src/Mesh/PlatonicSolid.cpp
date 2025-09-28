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

#include "PlatonicSolid.hpp"

#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include <numeric>

namespace Gedim
{
/// \brief MeshUtilities
/// \copyright See top level LICENSE file for details.
///
/// https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html

Gedim::GeometryUtilities::Polyhedron PlatonicSolid::dual_polyhedron(const Gedim::GeometryUtilities::Polyhedron &polyhedron) const
{

    std::vector<unsigned int> vertexMarkers(polyhedron.Vertices.cols());
    std::vector<unsigned int> edgeMarkers(polyhedron.Edges.cols());
    std::vector<unsigned int> faceMarkers(polyhedron.Faces.size());

    std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
    std::iota(edgeMarkers.begin(), edgeMarkers.end(), polyhedron.Vertices.cols() + 1);
    std::iota(faceMarkers.begin(), faceMarkers.end(), polyhedron.Vertices.cols() + polyhedron.Edges.cols() + 1);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    meshUtilities.Mesh3DFromPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);

    meshUtilities.ComputeCell0DCell2DNeighbours(mesh);
    meshUtilities.ComputeCell1DCell2DNeighbours(mesh);

    // the new dual mesh
    Gedim::GeometryUtilities::Polyhedron dual;
    dual.Vertices.resize(3, mesh.Cell2DTotalNumber());

    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
        dual.Vertices.col(f) = geometryUtilities.SimplexBarycenter(mesh.Cell2DVerticesCoordinates(f));

    std::vector<Eigen::VectorXi> faces(mesh.Cell0DTotalNumber());

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
        std::vector<unsigned int> neighbourhood = mesh.Cell0DNeighbourCell2Ds(v);
        Eigen::VectorXi face = -1 * Eigen::VectorXi::Ones(neighbourhood.size());

        unsigned int p = 0;
        unsigned int neighbour = neighbourhood[neighbourhood.size() - 1];
        neighbourhood.pop_back();
        while (!neighbourhood.empty())
        {
            face(p++) = neighbour;
            std::vector<unsigned int> edges = mesh.Cell2DEdges(neighbour);
            for (unsigned int e = 0; e < edges.size(); e++)
            {
                const std::vector<unsigned int> neighbourhood_edge = mesh.Cell1DNeighbourCell2Ds(edges[e]);
                std::vector<unsigned int> tmp;
                std::set_intersection(neighbourhood_edge.begin(),
                                      neighbourhood_edge.end(),
                                      neighbourhood.begin(),
                                      neighbourhood.end(),
                                      back_inserter(tmp));

                assert(tmp.size() <= 1);

                if (tmp.size() == 1)
                {
                    neighbour = tmp[0];
                    auto it = find(neighbourhood.begin(), neighbourhood.end(), neighbour);
                    if (it != neighbourhood.end())
                        neighbourhood.erase(it);
                    else
                        throw std::runtime_error("Not valid configuration in dual polyhedron");
                    break;
                }
            }

            if (neighbourhood.empty())
                face(p++) = neighbour;
        }

        faces[v] = face;
    }

    Gedim::MeshUtilities::ComputeMesh2DCell1DsResult result = meshUtilities.ComputeMesh2DCell1Ds(dual.Vertices, faces);

    dual.Edges = result.Cell1Ds;
    dual.Faces = result.Cell2Ds;

    return dual;
}

void PlatonicSolid::first_class_geodesic_polyhedron(const GeometryUtilities::Polyhedron &starting_polyhedron,
                                                    const unsigned int &frequency,
                                                    Gedim::MeshMatricesDAO &filter_mesh) const
{

    Output::Assert(frequency >= 1);

    std::vector<unsigned int> vertexMarkers(starting_polyhedron.Vertices.cols());
    std::vector<unsigned int> edgeMarkers(starting_polyhedron.Edges.cols());
    std::vector<unsigned int> faceMarkers(starting_polyhedron.Faces.size());

    std::iota(vertexMarkers.begin(), vertexMarkers.end(), 1);
    std::iota(edgeMarkers.begin(), edgeMarkers.end(), starting_polyhedron.Vertices.cols() + 1);
    std::iota(faceMarkers.begin(), faceMarkers.end(), starting_polyhedron.Vertices.cols() + starting_polyhedron.Edges.cols() + 1);

    if (frequency == 1)
    {
        meshUtilities.Mesh3DFromPolyhedron(starting_polyhedron.Vertices,
                                           starting_polyhedron.Edges,
                                           starting_polyhedron.Faces,
                                           vertexMarkers,
                                           edgeMarkers,
                                           faceMarkers,
                                           filter_mesh);
        return;
    }

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    meshUtilities.Mesh3DFromPolyhedron(starting_polyhedron.Vertices,
                                       starting_polyhedron.Edges,
                                       starting_polyhedron.Faces,
                                       vertexMarkers,
                                       edgeMarkers,
                                       faceMarkers,
                                       mesh);

    meshUtilities.ComputeCell1DCell2DNeighbours(mesh);

    const Eigen::VectorXd abscissa = Eigen::VectorXd::LinSpaced(frequency + 1, 0.0, 1.0);
    const unsigned int numCell1Ds = mesh.Cell1DTotalNumber();
    const unsigned int numCell2Ds = mesh.Cell2DTotalNumber();
    unsigned int id_vertices =
        mesh.Cell0DAppend((frequency - 1) * numCell1Ds + (0.5 * (frequency - 1) * (frequency - 2)) * numCell2Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
        const Eigen::Vector3d origin = mesh.Cell1DOriginCoordinates(e);
        const Eigen::Vector3d end = mesh.Cell1DEndCoordinates(e);
        const Eigen::Vector3d tangent = end - origin;

        unsigned int origin_id = mesh.Cell1DOrigin(e);
        Eigen::MatrixXi new_edges(2, frequency);
        for (unsigned p = 1; p < frequency; p++)
        {
            mesh.Cell0DInsertCoordinates(id_vertices, origin + abscissa(p) * tangent);
            mesh.Cell0DSetMarker(id_vertices, mesh.Cell1DMarker(e));
            mesh.Cell0DSetState(id_vertices, true);

            new_edges.col(p - 1) << origin_id, id_vertices;
            origin_id = id_vertices;

            id_vertices++;
        }

        new_edges.col(frequency - 1) << id_vertices - 1, mesh.Cell1DEnd(e);

        meshUtilities.SplitCell1D(e, new_edges, mesh);
    }

    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
        assert(mesh.Cell2DNumberVertices(f) == 3);

        const unsigned int second_vertex = mesh.Cell2DVertex(f, 1);

        const unsigned int first_edge = mesh.Cell2DEdge(f, 0);

        std::list<unsigned int> first_tmp;
        mesh.Cell1DUpdatedCell1Ds(first_edge, first_tmp);

        std::vector<unsigned int> first_edge_subcells(first_tmp.size());
        if (mesh.Cell1DOrigin(first_edge) == second_vertex)
        {
            unsigned int c = first_tmp.size() - 1;
            for (unsigned int cell : first_tmp)
                first_edge_subcells[c--] = cell;
        }
        else
        {

            unsigned int c = 0;
            for (unsigned int cell : first_tmp)
                first_edge_subcells[c++] = cell;
        }

        const unsigned int second_edge = mesh.Cell2DEdge(f, 1);

        std::list<unsigned int> second_tmp;
        mesh.Cell1DUpdatedCell1Ds(second_edge, second_tmp);

        std::vector<unsigned int> second_edge_subcells(second_tmp.size());
        if (mesh.Cell1DOrigin(second_edge) == second_vertex)
        {
            unsigned int c = second_tmp.size() - 1;
            for (unsigned int cell : second_tmp)
                second_edge_subcells[c--] = cell;
        }
        else
        {

            unsigned int c = 0;
            for (unsigned int cell : second_tmp)
                second_edge_subcells[c++] = cell;
        }

        std::vector<Eigen::MatrixXi> subCells;

        // first level
        unsigned int v2_first;
        unsigned int v2_second;
        unsigned int id_edges = mesh.Cell1DAppend(1.5 * frequency * (frequency + 1) - 2);
        {
            v2_first = mesh.Cell1DOrigin(first_edge_subcells[0]) == second_vertex
                           ? mesh.Cell1DEnd(first_edge_subcells[0])
                           : mesh.Cell1DOrigin(first_edge_subcells[0]);
            v2_second = mesh.Cell1DOrigin(second_edge_subcells[0]) == second_vertex
                            ? mesh.Cell1DEnd(second_edge_subcells[0])
                            : mesh.Cell1DOrigin(second_edge_subcells[0]);

            mesh.Cell1DInsertExtremes(id_edges, v2_first, v2_second);
            mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
            mesh.Cell1DSetState(id_edges, true);
            id_edges++;

            Eigen::MatrixXi face(2, 3);
            face << v2_first, second_vertex, v2_second, first_edge_subcells[0], second_edge_subcells[0], id_edges - 1;

            subCells.push_back(face);
        }

        unsigned int id_previous_horizontal_edge;
        for (unsigned int level = 1; level < frequency - 1; level++)
        {
            const unsigned int v1_first = v2_first;
            const unsigned int v1_second = v2_second;

            v2_first = mesh.Cell1DOrigin(first_edge_subcells[level]) == v1_first
                           ? mesh.Cell1DEnd(first_edge_subcells[level])
                           : mesh.Cell1DOrigin(first_edge_subcells[level]);
            v2_second = mesh.Cell1DOrigin(second_edge_subcells[level]) == v1_second
                            ? mesh.Cell1DEnd(second_edge_subcells[level])
                            : mesh.Cell1DOrigin(second_edge_subcells[level]);

            const Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(level + 2, 0.0, 1.0);

            const Eigen::Vector3d origin = mesh.Cell0DCoordinates(v2_first);
            const Eigen::Vector3d end = mesh.Cell0DCoordinates(v2_second);
            const Eigen::Vector3d tangent = end - origin;

            unsigned int origin_id = v2_first;
            Eigen::MatrixXi new_edges(2, frequency);
            for (unsigned p = 1; p < level + 1; p++)
            {
                mesh.Cell0DInsertCoordinates(id_vertices, origin + s(p) * tangent);
                mesh.Cell0DSetMarker(id_vertices, mesh.Cell2DMarker(f));
                mesh.Cell0DSetState(id_vertices, true);

                mesh.Cell1DInsertExtremes(id_edges, origin_id, id_vertices);
                mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                mesh.Cell1DSetState(id_edges, true);
                id_edges++;

                if (p == 1)
                {
                    mesh.Cell1DInsertExtremes(id_edges, id_vertices, v1_first);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    id_previous_horizontal_edge = id_edges - (3 * (level - 1) + 1) - 2;

                    Eigen::MatrixXi face(2, 3);
                    face << v2_first, v1_first, id_vertices, first_edge_subcells[level], id_edges - 1, id_edges - 2;
                    subCells.push_back(face);
                }
                else
                {
                    unsigned int jj = id_vertices - (level - 1) - 1;
                    mesh.Cell1DInsertExtremes(id_edges, id_vertices, jj);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face(2, 3);
                    face << jj, id_vertices, id_vertices - 1, id_edges - 1, id_edges - 2, id_edges - 3;
                    subCells.push_back(face);
                }

                if (p == level)
                {
                    mesh.Cell1DInsertExtremes(id_edges, id_vertices, v2_second);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    mesh.Cell1DInsertExtremes(id_edges, id_vertices, v1_second);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face1(2, 3);
                    face1 << id_vertices, v1_second, v2_second, id_edges - 1, second_edge_subcells[level], id_edges - 2;
                    subCells.push_back(face1);

                    unsigned int id_horizontal_edge;
                    if (level == 1)
                        id_horizontal_edge = id_edges - 5;
                    else
                        id_horizontal_edge = id_edges - (3 * level + 1) - 2;
                    Eigen::MatrixXi face2(2, 3);
                    face2 << v1_second, id_vertices, mesh.Cell1DOrigin(id_horizontal_edge), id_edges - 1, id_edges - 3,
                        id_horizontal_edge;
                    subCells.push_back(face2);
                }
                else
                {
                    unsigned int jj = id_vertices - (level - 1);
                    mesh.Cell1DInsertExtremes(id_edges, id_vertices, jj);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face(2, 3);
                    face << jj, id_vertices, mesh.Cell1DOrigin(id_previous_horizontal_edge + 3 * (p - 1)), id_edges - 1,
                        id_edges - 2, id_previous_horizontal_edge + 3 * (p - 1);
                    subCells.push_back(face);
                }

                origin_id = id_vertices;
                id_vertices++;
            }
        }

        // last level
        {

            const unsigned int level = frequency - 1;
            const unsigned int v1_first = v2_first;
            const unsigned int v1_second = v2_second;

            v2_first = mesh.Cell1DOrigin(first_edge_subcells[level]) == v1_first
                           ? mesh.Cell1DEnd(first_edge_subcells[level])
                           : mesh.Cell1DOrigin(first_edge_subcells[level]);
            v2_second = mesh.Cell1DOrigin(second_edge_subcells[level]) == v1_second
                            ? mesh.Cell1DEnd(second_edge_subcells[level])
                            : mesh.Cell1DOrigin(second_edge_subcells[level]);

            const unsigned int third_edge = mesh.Cell2DEdge(f, 2);

            std::list<unsigned int> third_tmp;
            mesh.Cell1DUpdatedCell1Ds(third_edge, third_tmp);

            std::vector<unsigned int> third_edge_subcells(third_tmp.size());
            if (mesh.Cell1DOrigin(third_edge) == v2_first)
            {
                unsigned int c = third_tmp.size() - 1;
                for (unsigned int cell : third_tmp)
                    third_edge_subcells[c--] = cell;
            }
            else
            {
                unsigned int c = 0;
                for (unsigned int cell : third_tmp)
                    third_edge_subcells[c++] = cell;
            }

            unsigned int origin_id = v2_first;
            for (unsigned p = 1; p < level + 1; p++)
            {

                const unsigned int last_id_vertices = mesh.Cell1DOrigin(third_edge_subcells[p - 1]) == origin_id
                                                          ? mesh.Cell1DEnd(third_edge_subcells[p - 1])
                                                          : mesh.Cell1DOrigin(third_edge_subcells[p - 1]);

                if (p == 1)
                {
                    mesh.Cell1DInsertExtremes(id_edges, last_id_vertices, v1_first);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    id_previous_horizontal_edge = id_edges - (3 * (level - 1) + 1) - 1;

                    Eigen::MatrixXi face(2, 3);
                    face << v2_first, v1_first, last_id_vertices, first_edge_subcells[level], id_edges - 1,
                        third_edge_subcells[p - 1];
                    subCells.push_back(face);
                }
                else
                {
                    unsigned int jj = id_vertices - (level - 1) + (p - 2);
                    mesh.Cell1DInsertExtremes(id_edges, last_id_vertices, jj);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face(2, 3);
                    face << jj, last_id_vertices, origin_id, id_edges - 1, third_edge_subcells[p - 1], id_edges - 2;
                    subCells.push_back(face);
                }

                if (p == level)
                {

                    mesh.Cell1DInsertExtremes(id_edges, last_id_vertices, v1_second);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face1(2, 3);
                    face1 << last_id_vertices, v1_second, v2_second, id_edges - 1, second_edge_subcells[level],
                        third_edge_subcells[p];
                    subCells.push_back(face1);

                    unsigned int id_horizontal_edge;
                    if (frequency == 2)
                        id_horizontal_edge = id_edges - 3;
                    else
                        id_horizontal_edge = id_edges - (2 * level) - 2;

                    Eigen::MatrixXi face2(2, 3);
                    face2 << v1_second, last_id_vertices, mesh.Cell1DOrigin(id_horizontal_edge), id_edges - 1,
                        id_edges - 2, id_horizontal_edge;
                    subCells.push_back(face2);
                }
                else
                {
                    unsigned int jj = id_vertices - (level - 1) + (p - 1);
                    mesh.Cell1DInsertExtremes(id_edges, last_id_vertices, jj);
                    mesh.Cell1DSetMarker(id_edges, mesh.Cell2DMarker(f));
                    mesh.Cell1DSetState(id_edges, true);
                    id_edges++;

                    Eigen::MatrixXi face(2, 3);
                    face << jj, last_id_vertices, mesh.Cell1DOrigin(id_previous_horizontal_edge + 3 * (p - 1)),
                        id_edges - 1, id_edges - 2, id_previous_horizontal_edge + 3 * (p - 1);
                    subCells.push_back(face);
                }

                origin_id = last_id_vertices;
            }
        }

        meshUtilities.SplitCell2D(f, subCells, mesh);
    }

    const auto filter_active_mesh = meshUtilities.FilterActiveMesh(mesh);
    const auto extract_mesh_data =
        meshUtilities.ExtractMesh2D(filter_active_mesh.Cell0Ds, filter_active_mesh.Cell1Ds, filter_active_mesh.Cell2Ds, mesh, filter_mesh);

    filter_mesh.InitializeDimension(3);
}

void PlatonicSolid::second_class_geodesic_polyhedron(const GeometryUtilities::Polyhedron &starting_polyhedron,
                                                     const unsigned int &frequency,
                                                     Gedim::MeshMatricesDAO &mesh) const
{

    Output::Assert(frequency >= 1);

    Gedim::MeshMatrices mesh_first_class_data;
    Gedim::MeshMatricesDAO mesh_first_class(mesh_first_class_data);
    first_class_geodesic_polyhedron(starting_polyhedron, frequency, mesh_first_class);

    meshUtilities.ComputeCell1DCell2DNeighbours(mesh_first_class);

    const unsigned int num_faces = mesh_first_class.Cell2DTotalNumber();
    const unsigned int num_edges = mesh_first_class.Cell1DTotalNumber();
    const unsigned int num_vertices = mesh_first_class.Cell0DTotalNumber();
    const unsigned int lower_bound = starting_polyhedron.Vertices.cols() + 1;
    const unsigned int upper_bound = starting_polyhedron.Vertices.cols() + starting_polyhedron.Edges.cols();

    mesh.InitializeDimension(3);
    mesh.Cell0DAppend(mesh_first_class.Cell0DTotalNumber());
    mesh.Cell0DsInsertCoordinates(mesh_first_class.Cell0DsCoordinates());
    for (unsigned int v = 0; v < num_vertices; v++)
    {
        mesh.Cell0DSetState(v, true);
        mesh.Cell0DSetMarker(v, mesh_first_class.Cell0DMarker(v));
    }

    unsigned int id_vertices = mesh.Cell0DAppend(starting_polyhedron.Edges.cols() * frequency + num_faces);
    unsigned int id_edges = mesh.Cell1DAppend(starting_polyhedron.Edges.cols() * 2 * frequency +
                                              starting_polyhedron.Faces.size() * 1.5 * (3 * frequency * frequency + frequency));
    unsigned int id_faces = mesh.Cell2DAppend(starting_polyhedron.Faces.size() * 3.0 * (frequency * frequency + frequency));

    const unsigned int offset_centroid = id_vertices;
    for (unsigned int f = 0; f < num_faces; f++)
    {
        mesh.Cell0DInsertCoordinates(id_vertices,
                                     geometryUtilities.PolygonBarycenter(mesh_first_class.Cell2DVerticesCoordinates(f)));
        mesh.Cell0DSetMarker(id_vertices, mesh_first_class.Cell2DMarker(f));
        mesh.Cell0DSetState(id_vertices, true);
        id_vertices++;

        for (unsigned int i = 0; i < 3; i++)
        {
            mesh.Cell1DInsertExtremes(id_edges, id_vertices - 1, mesh_first_class.Cell2DVertex(f, i));
            mesh.Cell1DSetMarker(id_edges, mesh_first_class.Cell2DMarker(f));
            mesh.Cell1DSetState(id_edges, true);
            id_edges++;
        }
    }

    for (unsigned int e = 0; e < num_edges; e++)
    {
        const unsigned int marker = mesh_first_class.Cell1DMarker(e);
        if (marker >= lower_bound && marker <= upper_bound)
        {
            const Eigen::Vector3d origin = mesh_first_class.Cell1DOriginCoordinates(e);
            const Eigen::Vector3d tangent = mesh_first_class.Cell1DEndCoordinates(e) - origin;

            mesh.Cell0DInsertCoordinates(id_vertices, origin + 0.5 * tangent);
            mesh.Cell0DSetMarker(id_vertices, mesh_first_class.Cell1DMarker(e));
            mesh.Cell0DSetState(id_vertices, true);
            id_vertices++;

            mesh.Cell1DInsertExtremes(id_edges, mesh_first_class.Cell1DOrigin(e), id_vertices - 1);
            mesh.Cell1DSetMarker(id_edges, mesh_first_class.Cell1DMarker(e));
            mesh.Cell1DSetState(id_edges, true);
            id_edges++;

            mesh.Cell1DInsertExtremes(id_edges, id_vertices - 1, mesh_first_class.Cell1DEnd(e));
            mesh.Cell1DSetMarker(id_edges, mesh_first_class.Cell1DMarker(e));
            mesh.Cell1DSetState(id_edges, true);
            id_edges++;

            const auto neighbourhood = mesh_first_class.Cell1DNeighbourCell2Ds(e);

            if (neighbourhood.size() != 2)
                throw std::runtime_error("not valid edge");

            for (unsigned int n = 0; n < 2; n++)
            {
                mesh.Cell1DInsertExtremes(id_edges, id_vertices - 1, offset_centroid + neighbourhood[n]);
                mesh.Cell1DSetMarker(id_edges, mesh_first_class.Cell1DMarker(e));
                mesh.Cell1DSetState(id_edges, true);
                id_edges++;

                const unsigned int id_local_edge_origin =
                    mesh_first_class.Cell2DFindVertex(neighbourhood[n], mesh_first_class.Cell1DOrigin(e));
                const unsigned int id_local_edge_end =
                    mesh_first_class.Cell2DFindVertex(neighbourhood[n], mesh_first_class.Cell1DEnd(e));

                std::vector<unsigned int> face_vertices_1 = {mesh_first_class.Cell1DOrigin(e),
                                                             id_vertices - 1,
                                                             offset_centroid + neighbourhood[n]};
                std::vector<unsigned int> face_edges_1 = {id_edges - 3 - n, id_edges - 1, 3 * neighbourhood[n] + id_local_edge_origin};
                mesh.Cell2DInitializeVertices(id_faces, 3);
                mesh.Cell2DInitializeEdges(id_faces, 3);
                mesh.Cell2DInsertVertices(id_faces, face_vertices_1);
                mesh.Cell2DInsertEdges(id_faces, face_edges_1);
                mesh.Cell2DSetState(id_faces, true);
                mesh.Cell2DSetMarker(id_faces, mesh_first_class.Cell2DMarker(neighbourhood[n]));
                id_faces++;

                std::vector<unsigned int> face_vertices_2 = {mesh_first_class.Cell1DEnd(e),
                                                             id_vertices - 1,
                                                             offset_centroid + neighbourhood[n]};
                std::vector<unsigned int> face_edges_2 = {id_edges - 2 - n, id_edges - 1, 3 * neighbourhood[n] + id_local_edge_end};
                mesh.Cell2DInitializeVertices(id_faces, 3);
                mesh.Cell2DInitializeEdges(id_faces, 3);
                mesh.Cell2DInsertVertices(id_faces, face_vertices_2);
                mesh.Cell2DInsertEdges(id_faces, face_edges_2);
                mesh.Cell2DSetState(id_faces, true);
                mesh.Cell2DSetMarker(id_faces, mesh_first_class.Cell2DMarker(neighbourhood[n]));
                id_faces++;
            }
        }
        else
        {
            const auto neighbourhood = mesh_first_class.Cell1DNeighbourCell2Ds(e);

            if (neighbourhood.size() != 2)
                throw std::runtime_error("not valid edge");

            mesh.Cell1DInsertExtremes(id_edges, offset_centroid + neighbourhood[0], offset_centroid + neighbourhood[1]);
            mesh.Cell1DSetMarker(id_edges, mesh_first_class.Cell2DMarker(neighbourhood[0]));
            mesh.Cell1DSetState(id_edges, true);
            id_edges++;

            const unsigned int id_local_origin_n0 =
                mesh_first_class.Cell2DFindVertex(neighbourhood[0], mesh_first_class.Cell1DOrigin(e));
            const unsigned int id_local_origin_n1 =
                mesh_first_class.Cell2DFindVertex(neighbourhood[1], mesh_first_class.Cell1DOrigin(e));

            const unsigned int id_local_end_n0 =
                mesh_first_class.Cell2DFindVertex(neighbourhood[0], mesh_first_class.Cell1DEnd(e));
            const unsigned int id_local_end_n1 =
                mesh_first_class.Cell2DFindVertex(neighbourhood[1], mesh_first_class.Cell1DEnd(e));

            std::vector<unsigned int> face_vertices_1 = {offset_centroid + neighbourhood[0],
                                                         offset_centroid + neighbourhood[1],
                                                         mesh_first_class.Cell1DOrigin(e)};
            std::vector<unsigned int> face_edges_1 = {id_edges - 1,
                                                      3 * neighbourhood[1] + id_local_origin_n1,
                                                      3 * neighbourhood[0] + id_local_origin_n0};
            mesh.Cell2DInitializeVertices(id_faces, 3);
            mesh.Cell2DInitializeEdges(id_faces, 3);
            mesh.Cell2DInsertVertices(id_faces, face_vertices_1);
            mesh.Cell2DInsertEdges(id_faces, face_edges_1);
            mesh.Cell2DSetState(id_faces, true);
            mesh.Cell2DSetMarker(id_faces, mesh_first_class.Cell2DMarker(neighbourhood[0]));
            id_faces++;

            std::vector<unsigned int> face_vertices_2 = {offset_centroid + neighbourhood[0],
                                                         offset_centroid + neighbourhood[1],
                                                         mesh_first_class.Cell1DEnd(e)};
            std::vector<unsigned int> face_edges_2 = {id_edges - 1,
                                                      3 * neighbourhood[1] + id_local_end_n1,
                                                      3 * neighbourhood[0] + id_local_end_n0};
            mesh.Cell2DInitializeVertices(id_faces, 3);
            mesh.Cell2DInitializeEdges(id_faces, 3);
            mesh.Cell2DInsertVertices(id_faces, face_vertices_2);
            mesh.Cell2DInsertEdges(id_faces, face_edges_2);
            mesh.Cell2DSetState(id_faces, true);
            mesh.Cell2DSetMarker(id_faces, mesh_first_class.Cell2DMarker(neighbourhood[0]));
            id_faces++;
        }
    }
}

} // namespace Gedim
