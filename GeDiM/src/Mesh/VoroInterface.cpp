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

#include "VoroInterface.hpp"

#include "Eigen/Eigen"
#include "Gedim_Macro.hpp"
#include "MeshUtilities.hpp"
#include <numeric>
#include <unistd.h>

#include "VTKUtilities.hpp"

using namespace std;

namespace Gedim
{
// ************************************************************************* //
VoroInterface::VoroInterface(const Gedim::GeometryUtilities &geometryUtilities) : geometryUtilities(geometryUtilities)
{
}
// ************************************************************************* //
#if ENABLE_VORO == 0
void VoroInterface::GenerateVoronoiTassellations3D(const Eigen::MatrixXd &polyhedronVertices,
                                                   const Eigen::MatrixXi &polyhedronEdges,
                                                   const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                                   const unsigned int &numPoints,
                                                   const unsigned int &numIterations,
                                                   Gedim::IMeshDAO &mesh,
                                                   const unsigned int random_seed)
{
    Gedim::Utilities::Unused(polyhedronVertices);
    Gedim::Utilities::Unused(polyhedronEdges);
    Gedim::Utilities::Unused(polyhedronFaces);
    Gedim::Utilities::Unused(numPoints);
    Gedim::Utilities::Unused(numIterations);
    Gedim::Utilities::Unused(mesh);
    Gedim::Utilities::Unused(random_seed);
    throw std::runtime_error("Not active module VORO");
}

void VoroInterface::GenerateVoronoiTassellations2D(const Eigen::MatrixXd &polygonVertices,
                                                   const unsigned int &numPoints,
                                                   const unsigned int &numIterations,
                                                   Gedim::IMeshDAO &mesh,
                                                   const unsigned int random_seed)
{
    Gedim::Utilities::Unused(polygonVertices);
    Gedim::Utilities::Unused(numPoints);
    Gedim::Utilities::Unused(numIterations);
    Gedim::Utilities::Unused(mesh);
    Gedim::Utilities::Unused(random_seed);
    throw std::runtime_error("Not active module VORO");
}
#else
// ************************************************************************* //
unsigned int VoroInterface::InsertNewPoint(const Cell0D &cell0D, std::map<unsigned int, Cell0D> &cell0Ds)
{
    const unsigned int numCell0Ds = cell0Ds.size();
    if (numCell0Ds == 0)
    {
        cell0Ds.insert({0, cell0D});
        return 0;
    }
    else
    {
        for (const auto &other_cell0D : cell0Ds)
        {
            if (geometryUtilities.PointsAreCoincident(other_cell0D.second.coordinates, cell0D.coordinates))
                return other_cell0D.first;
        }

        cell0Ds.insert({numCell0Ds, cell0D});

        return numCell0Ds;
    }
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations2D(const Eigen::MatrixXd &polygonVertices,
                                                   const unsigned int &numPoints,
                                                   const unsigned int &numIterations,
                                                   Gedim::IMeshDAO &mesh,
                                                   const unsigned int random_seed)
{
    Eigen::MatrixXd VoronoiPoints = GenerateRandomPoints(polygonVertices, numPoints, random_seed);

    VoroInterface::GenerateVoronoiTassellations2D(polygonVertices, numIterations, VoronoiPoints, mesh);
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations2D(const Eigen::MatrixXd &domain_vertices,
                                                   const unsigned int &num_iterations,
                                                   Eigen::MatrixXd &VoronoiPoints,
                                                   Gedim::IMeshDAO &mesh)
{

    const unsigned int num_points = VoronoiPoints.cols();

    if (num_points == 0)
    {
        std::vector<unsigned int> vertex_markers(domain_vertices.cols());
        std::iota(vertex_markers.begin(), vertex_markers.end(), 1);
        std::vector<unsigned int> edge_markers(domain_vertices.cols());
        std::iota(edge_markers.begin(), edge_markers.end(), domain_vertices.cols() + 1);

        Gedim::MeshUtilities mesh_utilities;
        mesh_utilities.Mesh2DFromPolygon(domain_vertices, vertex_markers, edge_markers, mesh);

        VoronoiPoints = geometryUtilities.PolygonBarycenter(domain_vertices);

        return;
    }

    const unsigned int numVerticesDomain = domain_vertices.cols();
    const unsigned int numEdgesDomain = numVerticesDomain;

    vector<Eigen::Vector3d> edgeDomainTangent(numEdgesDomain);
    vector<double> edgeDomainTangentSquaredNorm(numEdgesDomain);
    vector<Eigen::Vector3d> edgeDomainLineOrigin(numEdgesDomain);
    for (unsigned int e = 0; e < numEdgesDomain; e++)
    {
        edgeDomainLineOrigin[e] = domain_vertices.col(e);
        const Eigen::Vector3d end = domain_vertices.col((e + 1) % numEdgesDomain);

        Eigen::Vector3d tangent = end - edgeDomainLineOrigin[e];

        edgeDomainTangent[e] = tangent;
        edgeDomainTangentSquaredNorm[e] = tangent.squaredNorm();
    }

    //    vector<unsigned int> map_marker(numFacesDomain, 0);

    const double x_min = domain_vertices.row(0).minCoeff(), x_max = domain_vertices.row(0).maxCoeff();
    const double y_min = domain_vertices.row(1).minCoeff(), y_max = domain_vertices.row(1).maxCoeff();
    const double cvol = (x_max - x_min) * (y_max - y_min);

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = n_x, n_z = n_x;

    // Loop on Voronoi cells
    for (unsigned int it = 0; it < num_iterations; it++)
    {
        // Initialize Voronoi mesh
        voro::container con(x_min, x_max, y_min, y_max, -0.5, 0.5, n_x, n_y, n_z, false, false, false, 8);

        for (unsigned int i = 0; i < VoronoiPoints.cols(); i++)
            con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));

        unsigned int countPoints = 0;

        voro::c_loop_all vl(con);
        voro::voronoicell_neighbor c;
        if (vl.start())
        {
            do
            {
                if (con.compute_cell(c, vl))
                {
                    const unsigned int ijk = vl.ijk;
                    const unsigned int q = vl.q;
                    const double *pp = con.p[ijk] + con.ps * q;

                    // Cell barycenter
                    const double cx = pp[0];
                    const double cy = pp[1];
                    const double cz = pp[2];

                    // Get cell generator
                    std::vector<double> vertices;
                    c.vertices(cx, cy, cz, vertices);

                    // Get cell face cardinality
                    std::vector<int> faces_num_vertices;
                    c.face_orders(faces_num_vertices);

                    // Get face-vertices connectivity
                    std::vector<int> faces_vertices;
                    c.face_vertices(faces_vertices);

                    // Get cell neighbors IDs
                    std::vector<int> faces_neighs;
                    c.neighbors(faces_neighs);

                    const unsigned int num_cell_vertices = vertices.size() / 3;
                    std::vector<bool> consider_vertices(num_cell_vertices, true);
                    Eigen::MatrixXd cell_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_cell_vertices);
                    for (unsigned int v = 0; v < num_cell_vertices; v++)
                    {
                        cell_vertices_coordinates.col(v) << vertices[3 * v], vertices[3 * v + 1], 0.0;

                        if (!geometryUtilities.IsValueZero(vertices[3 * v + 2] + 0.5, geometryUtilities.Tolerance1D()))
                            consider_vertices[v] = false;
                    }

                    // cell faces
                    const unsigned int num_faces = faces_num_vertices.size();

                    unsigned int count = 0;
                    for (unsigned int i = 0; i < num_faces; i++)
                    {
                        const unsigned num_face_vertices = faces_num_vertices[i];
                        Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, num_face_vertices);

                        count++;
                        unsigned int num_non_considered_vertices = 0;
                        for (unsigned int j = 0; j < num_face_vertices; j++)
                        {
                            if (!consider_vertices[faces_vertices[count]])
                                num_non_considered_vertices++;
                            else
                                faceVertCoordinates.col(j) = cell_vertices_coordinates.col(faces_vertices[count]);

                            count++;
                        }

                        if (num_non_considered_vertices != 0)
                            continue;

                        if (faces_neighs[i] >= 0) // if it is a boundary face
                            throw runtime_error("wrong face");

                        VoronoiPoints.col(countPoints++) = geometryUtilities.PolygonBarycenter(faceVertCoordinates);
                    }
                }
            } while (vl.inc());
        }
    }

    double vvol = 0.0;

    // Get cell number
    std::map<unsigned int, Cell0D> cell0Ds;
    std::map<unsigned int, Cell2D> cell2Ds;
    // key: origin-end, value: id
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> origin_end_cell1Ds;
    std::map<unsigned int, Cell1D> cell1Ds;

    // Initialize Voronoi mesh
    voro::container con(x_min, x_max, y_min, y_max, -0.5, 0.5, n_x, n_y, n_z, false, false, false, 8);

    for (unsigned int i = 0; i < num_points; i++)
        con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do
        {
            if (con.compute_cell(c, vl))
            {
                const unsigned int ijk = vl.ijk;
                const unsigned int q = vl.q;
                const double *pp = con.p[ijk] + con.ps * q;

                // Cell barycenter
                const double cx = pp[0];
                const double cy = pp[1];
                const double cz = pp[2];

                const unsigned int cell2DIndex = con.id[ijk][q];

                // Get cell generator
                std::vector<double> vertices;
                c.vertices(cx, cy, cz, vertices);

                // Get cell face cardinality
                std::vector<int> faces_num_vertices;
                c.face_orders(faces_num_vertices);

                // Get face-vertices connectivity
                std::vector<int> faces_vertices;
                c.face_vertices(faces_vertices);

                // Get cell neighbors IDs
                std::vector<int> faces_neighs;
                c.neighbors(faces_neighs);

                // cell vertices
                const unsigned int num_cell_vertices = vertices.size() / 3;
                std::vector<bool> consider_vertices(num_cell_vertices, true);
                Eigen::MatrixXd cell_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_cell_vertices);
                std::vector<unsigned int> cell_vertices(num_cell_vertices, std::numeric_limits<unsigned int>::max());
                for (unsigned int v = 0; v < num_cell_vertices; v++)
                {
                    cell_vertices_coordinates.col(v) << vertices[3 * v], vertices[3 * v + 1], 0.0;

                    if (!geometryUtilities.IsValueZero(vertices[3 * v + 2] + 0.5, geometryUtilities.Tolerance1D()))
                        consider_vertices[v] = false;
                    else
                    {
                        Cell0D cell0D(cell_vertices_coordinates.col(v));
                        cell_vertices[v] = InsertNewPoint(cell0D, cell0Ds);
                    }
                }

                // cell faces
                const unsigned int num_faces = faces_num_vertices.size();

                unsigned int count = 0;
                bool found_face = false;

                Cell2D cell2D;
                cell2D.marker = 0;
#ifdef TEST_VORO
                unsigned int num_border_edge = 0;
#endif
                for (unsigned int f = 0; f < num_faces; f++)
                {
                    const unsigned num_face_vertices = faces_num_vertices[f];
                    std::vector<unsigned int> face_vertices(num_face_vertices);
                    Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, num_face_vertices);

                    count++;
                    unsigned int num_non_considered_vertices = 0;
                    for (unsigned int j = 0; j < num_face_vertices; j++)
                    {
                        if (!consider_vertices[faces_vertices[count]])
                            num_non_considered_vertices++;
                        else
                        {
                            face_vertices[j] = cell_vertices[faces_vertices[count]];
                            faceVertCoordinates.col(j) = cell_vertices_coordinates.col(faces_vertices[count]);
                        }
                        count++;
                    }

                    if (num_non_considered_vertices > 0)
                    {
                        if (faces_neighs[f] >= 0)
                            cell2D.neighbors_of_related_3D_cells.push_back(faces_neighs[f]);

                        continue;
                    }

                    if (faces_neighs[f] >= 0) // if it is a boundary face
                        throw runtime_error("wrong face");

                    cell2D.vertices = face_vertices;

                    cell2D.edges.resize(num_face_vertices);
                    for (unsigned int v = 0; v < num_face_vertices; v++)
                    {
                        const unsigned int id_min = min(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);
                        const unsigned int id_max = max(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);

                        if (id_min == id_max)
                            throw std::runtime_error("a degenerated edge is created");

                        const auto it = origin_end_cell1Ds.find({id_min, id_max});
                        if (it != origin_end_cell1Ds.end())
                        {
                            const unsigned int edgeId = it->second;
                            cell2D.edges[v] = edgeId;

#ifdef TEST_VORO
                            if (cell1Ds[edgeId].marker != 0)
                                num_border_edge++;
#endif
                        }
                        else
                        {
                            unsigned int marker = 0;

                            for (unsigned int e = 0; e < numEdgesDomain; e++)
                            {

                                if (geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                    edgeDomainLineOrigin[e],
                                                                    edgeDomainTangent[e],
                                                                    edgeDomainTangentSquaredNorm[e]) &&
                                    geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v + 1) % num_face_vertices),
                                                                    edgeDomainLineOrigin[e],
                                                                    edgeDomainTangent[e],
                                                                    edgeDomainTangentSquaredNorm[e]))
                                {
                                    marker = numVerticesDomain + e + 1;
                                    break;
                                }
                            }

                            const unsigned int cell1DIndex = cell1Ds.size();
                            Cell1D cell1D;
                            cell1D.origin = id_min;
                            cell1D.end = id_max;
                            cell1D.marker = marker;
                            origin_end_cell1Ds.insert({{id_min, id_max}, cell1DIndex});
                            cell1Ds.insert({cell1DIndex, cell1D});
                            cell2D.edges[v] = cell1DIndex;
#ifdef TEST_VORO
                            if (marker != 0)
                                num_border_edge++;
#endif
                        }
                    }

                    found_face = true;
                }

                if (!found_face)
                    throw std::runtime_error("face - cell not uniquely identified");
#ifdef TEST_VORO
                if (cell2D.neighbors_of_related_3D_cells.size() != cell2D.edges.size() - num_border_edge)
                    throw std::runtime_error("not considered case");
#endif

                cell2Ds.insert({cell2DIndex, cell2D});
                vvol += c.volume();
            }
        } while (vl.inc());
    }

    // Compute Cell1D neighbours starting from cell2Ds
    for (const auto &cell2D : cell2Ds)
    {
        const unsigned int c2D = cell2D.first;

        const unsigned int numCell2DEdges = cell2D.second.edges.size();
        for (unsigned int e = 0; e < numCell2DEdges; e++)
        {
            const unsigned int cell1DIndex = cell2D.second.edges[e];
            cell1Ds[cell1DIndex].neighbors_2D.push_back(c2D);
        }

        const unsigned int numCell2DVertices = cell2D.second.vertices.size();
        for (unsigned int v = 0; v < numCell2DVertices; v++)
        {
            const unsigned int cell0DIndex = cell2D.second.vertices[v];
            cell0Ds.at(cell0DIndex).neighbors_2D.push_back(c2D);
        }
    }

    for (const auto &cell1D : cell1Ds)
    {
        cell0Ds.at(cell1D.second.origin).neighbors_1D.push_back(cell1D.first);
        cell0Ds.at(cell1D.second.end).neighbors_1D.push_back(cell1D.first);
    }

    struct CollapsingEdge
    {
        unsigned int id = std::numeric_limits<unsigned int>::max();
        double dist = std::numeric_limits<double>::max();
        unsigned int local_edge_id;
        unsigned int cell_2D;
        unsigned int corresponding_origin;
        unsigned int corresponding_end;
    };

    unsigned int numberOfPointsMesh = cell0Ds.size();
    unsigned int numberOfEdgesMesh = cell1Ds.size();
    unsigned int numberOfCellsMesh = cell2Ds.size();

    for (auto &cell1D : cell1Ds)
    {
        const unsigned int num_neighbors_2D = cell1D.second.neighbors_2D.size();

        Gedim::Output::Assert(num_neighbors_2D > 0);

        if (!cell1D.second.active)
            continue;

        if (cell1D.second.marker != 0 && num_neighbors_2D == 1)
            continue;

        if (cell1D.second.marker == 0 && num_neighbors_2D == 2)
            continue;

        if (cell1D.second.marker == 0 && num_neighbors_2D == 1)
        {
            const unsigned int cell2DIndex = cell1D.second.neighbors_2D[0];
            const Cell2D &cell2D = cell2Ds.at(cell2DIndex);
            const unsigned int num_cell2D_neighbors_3D = cell2D.neighbors_of_related_3D_cells.size();

            CollapsingEdge collapsing_edge;
            for (unsigned int n_2D = 0; n_2D < num_cell2D_neighbors_3D; n_2D++)
            {
                // Check which no-border edge of each neigh has only one neighbor -> these are the candidate equal edges
                const Cell2D &neighCell2D = cell2Ds.at(cell2D.neighbors_of_related_3D_cells[n_2D]);

                for (unsigned int e = 0; e < neighCell2D.edges.size(); e++)
                {
                    if (cell1Ds.at(neighCell2D.edges[e]).marker != 0)
                        continue;

                    if (cell1Ds.at(neighCell2D.edges[e]).neighbors_2D.size() == 2)
                        continue;

                    const Cell1D &candidateCell1D = cell1Ds.at(neighCell2D.edges[e]);

                    const double d_origin_origin =
                        (cell0Ds.at(cell1D.second.origin).coordinates - cell0Ds.at(candidateCell1D.origin).coordinates).squaredNorm();
                    const double d_origin_end =
                        (cell0Ds.at(cell1D.second.origin).coordinates - cell0Ds.at(candidateCell1D.end).coordinates).squaredNorm();

                    double dist = 0.0;
                    unsigned int corresponding_origin = std::numeric_limits<unsigned int>::max();
                    unsigned int corresponding_end = std::numeric_limits<unsigned int>::max();
                    if (d_origin_origin > d_origin_end)
                    {
                        dist =
                            d_origin_end +
                            (cell0Ds.at(cell1D.second.end).coordinates - cell0Ds.at(candidateCell1D.origin).coordinates).squaredNorm();
                        corresponding_origin = candidateCell1D.end;
                        corresponding_end = candidateCell1D.origin;
                    }
                    else
                    {
                        dist = d_origin_origin +
                               (cell0Ds.at(cell1D.second.end).coordinates - cell0Ds.at(candidateCell1D.end).coordinates).squaredNorm();
                        corresponding_origin = candidateCell1D.origin;
                        corresponding_end = candidateCell1D.end;
                    }

                    if (dist < collapsing_edge.dist)
                    {
                        collapsing_edge.id = neighCell2D.edges[e];
                        collapsing_edge.dist = dist;
                        collapsing_edge.local_edge_id = e;
                        collapsing_edge.cell_2D = cell2D.neighbors_of_related_3D_cells[n_2D];
                        collapsing_edge.corresponding_origin = corresponding_origin;
                        collapsing_edge.corresponding_end = corresponding_end;
                    }
                }
            }

            if (collapsing_edge.id == std::numeric_limits<unsigned int>::max())
                throw std::runtime_error("not considered case");

            // Turn off the edge
            cell1Ds.at(collapsing_edge.id).active = false;
            numberOfEdgesMesh--;
            cell1D.second.neighbors_2D.push_back(collapsing_edge.cell_2D);

            // Propagate information to cell2D
            cell2Ds.at(collapsing_edge.cell_2D).edges.at(collapsing_edge.local_edge_id) = cell1D.first;

            if (collapsing_edge.corresponding_origin != cell1D.second.origin)
            {
                Cell0D &correspondingCell0D = cell0Ds.at(collapsing_edge.corresponding_origin);
                for (unsigned int c2D = 0; c2D < correspondingCell0D.neighbors_2D.size(); c2D++)
                {
                    for (unsigned int v = 0; v < cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices.size(); v++)
                    {
                        if (cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices[v] == collapsing_edge.corresponding_origin)
                            cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices[v] = cell1D.second.origin;
                    }
                }

                for (unsigned int c1D = 0; c1D < correspondingCell0D.neighbors_1D.size(); c1D++)
                {

                    if (cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).origin == collapsing_edge.corresponding_origin)
                        cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).origin = cell1D.second.origin;
                    else if (cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).end == collapsing_edge.corresponding_origin)
                        cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).end = cell1D.second.origin;
                    else
                        throw std::runtime_error("not valid case");
                }

                numberOfPointsMesh--;
                correspondingCell0D.active = false;
            }

            if (collapsing_edge.corresponding_end != cell1D.second.end)
            {
                Cell0D &correspondingCell0D = cell0Ds.at(collapsing_edge.corresponding_end);
                for (unsigned int c2D = 0; c2D < correspondingCell0D.neighbors_2D.size(); c2D++)
                {
                    for (unsigned int v = 0; v < cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices.size(); v++)
                    {
                        if (cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices[v] == collapsing_edge.corresponding_end)
                            cell2Ds.at(correspondingCell0D.neighbors_2D[c2D]).vertices[v] = cell1D.second.end;
                    }
                }

                for (unsigned int c1D = 0; c1D < correspondingCell0D.neighbors_1D.size(); c1D++)
                {

                    if (cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).origin == collapsing_edge.corresponding_end)
                        cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).origin = cell1D.second.end;
                    else if (cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).end == collapsing_edge.corresponding_end)
                        cell1Ds.at(correspondingCell0D.neighbors_1D[c1D]).end = cell1D.second.end;
                    else
                        throw std::runtime_error("not valid case");
                }

                numberOfPointsMesh--;
                correspondingCell0D.active = false;
            }

            continue;
        }

        throw std::runtime_error("not considered case");
    }

    /// <li> Set Cell0Ds
    mesh.InitializeDimension(2);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    mesh.Cell2DsInitialize(numberOfCellsMesh);
    mesh.Cell3DsInitialize(0);

    unsigned int cell0DIndex = 0;
    std::vector<unsigned int> map_cell0D_id(cell0Ds.size(), std::numeric_limits<unsigned int>::max());
    for (const auto &cell0D : cell0Ds)
    {
        if (!cell0D.second.active)
            continue;

        map_cell0D_id[cell0D.first] = cell0DIndex;

        mesh.Cell0DSetState(cell0DIndex, true);
        mesh.Cell0DInsertCoordinates(cell0DIndex, cell0D.second.coordinates);

        const bool is_xmin =
            (geometryUtilities.IsValueZero(cell0D.second.coordinates(0) - x_min, geometryUtilities.Tolerance1D())); // l == 0
        const bool is_xmax =
            (geometryUtilities.IsValueZero(cell0D.second.coordinates(0) - x_max, geometryUtilities.Tolerance1D())); // l == (numLengthPoints - 1)
        const bool is_ymin =
            (geometryUtilities.IsValueZero(cell0D.second.coordinates(1) - y_min, geometryUtilities.Tolerance1D())); // h == 0
        const bool is_ymax =
            (geometryUtilities.IsValueZero(cell0D.second.coordinates(1) - y_max, geometryUtilities.Tolerance1D())); // h == (numHeightPoints - 1)

        const unsigned int marker = 1 * (is_xmin && is_ymin) + 2 * (is_xmax && is_ymin) + 4 * (is_xmin && is_ymax) +
                                    3 * (is_xmax && is_ymax) + 5 * (is_ymin && !is_xmin && !is_xmax) +
                                    7 * (is_ymax && !is_xmin && !is_xmax) + 8 * (is_xmin && !is_ymin && !is_ymax) +
                                    6 * (is_xmax && !is_ymin && !is_ymax);

        mesh.Cell0DSetMarker(cell0DIndex, marker);
        cell0DIndex++;
    }

    /// <li> Set Edges
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    unsigned int cell1DIndex = 0;
    std::vector<unsigned int> map_cell1D_id(cell1Ds.size(), std::numeric_limits<unsigned int>::max());
    for (const auto &cell1D : cell1Ds)
    {
        if (!cell1D.second.active)
            continue;

        map_cell1D_id[cell1D.first] = cell1DIndex;
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, cell1D.second.marker);
        mesh.Cell1DInsertExtremes(cell1DIndex, map_cell0D_id[cell1D.second.origin], map_cell0D_id[cell1D.second.end]);

        cell1DIndex++;
    }

    /// <li> Set cells
    double gvol = 0.0;
    unsigned int cell2DIndex = 0;
    for (const auto &cell2D : cell2Ds)
    {
        const unsigned int f = cell2DIndex;

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, cell2D.second.marker);

        const unsigned int numCellVertices = cell2D.second.vertices.size();
        const unsigned int numCellEdges = cell2D.second.edges.size();

        Gedim::Output::Assert(numCellVertices == numCellEdges);

        bool oriented = false;
        for (unsigned int v = 0; v < numCellVertices; v++)
        {
            Eigen::MatrixXd vert = Eigen::MatrixXd::Zero(3, 3);
            vert.col(0) = mesh.Cell0DCoordinates(map_cell0D_id[cell2D.second.vertices[v]]);
            vert.col(1) = mesh.Cell0DCoordinates(map_cell0D_id[cell2D.second.vertices[(v + 1) % numCellVertices]]);
            vert.col(2) = mesh.Cell0DCoordinates(map_cell0D_id[cell2D.second.vertices[(v + 2) % numCellVertices]]);

            oriented = geometryUtilities.IsValuePositive(geometryUtilities.PolygonArea(vert), geometryUtilities.Tolerance2D());
            if (oriented)
                break;
        }

        mesh.Cell2DInitializeVertices(f, numCellVertices);
        mesh.Cell2DInitializeEdges(f, numCellEdges);

        if (oriented)
        {
            for (unsigned int v = 0; v < numCellVertices; v++)
                mesh.Cell2DInsertVertex(f, v, map_cell0D_id[cell2D.second.vertices[v]]);

            for (unsigned int e = 0; e < numCellEdges; e++)
                mesh.Cell2DInsertEdge(f, e, map_cell1D_id[cell2D.second.edges[e]]);
        }
        else
        {
            for (unsigned int v = 0; v < numCellVertices; v++)
                mesh.Cell2DInsertVertex(f, v, map_cell0D_id[cell2D.second.vertices[(numCellVertices - 1) - v]]);

            for (unsigned int e = 0; e < numCellEdges - 1; e++)
                mesh.Cell2DInsertEdge(f, e, map_cell1D_id[cell2D.second.edges[(numCellEdges - 2) - e]]);

            mesh.Cell2DInsertEdge(f, numCellEdges - 1, map_cell1D_id[cell2D.second.edges[numCellEdges - 1]]);
        }

        gvol += geometryUtilities.PolygonArea(mesh.Cell2DVerticesCoordinates(f));
        cell2DIndex++;
    }

    // if (!geometryUtilities.AreValuesEqual(cvol, vvol, geometryUtilities.Tolerance2D()))
    // {
    //     std::cout.precision(16);
    //     std::cout << std::scientific << cvol << " vs " << vvol << " " << gvol << std::endl;
    //     // throw runtime_error("Error generating Voronoi cells: areas do not mathc each other");
    // }

    if (!geometryUtilities.AreValuesEqual(cvol, gvol, geometryUtilities.Tolerance2D()))
        throw runtime_error("Error generating Voronoi cells: areas do not match each other");
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations3D(const Eigen::MatrixXd &domain_vertices,
                                                   const Eigen::MatrixXi &domain_edges,
                                                   const std::vector<Eigen::MatrixXi> &domain_faces,
                                                   const unsigned int &num_points,
                                                   const unsigned int &num_iterations,
                                                   Gedim::IMeshDAO &mesh,
                                                   const unsigned int random_seed)
{
    Eigen::MatrixXd VoronoiPoints = GenerateRandomPoints(domain_vertices, num_points, random_seed);
    VoroInterface::GenerateVoronoiTassellations3D(domain_vertices, domain_edges, domain_faces, num_iterations, VoronoiPoints, mesh);
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations3D(const Eigen::MatrixXd &domain_vertices,
                                                   const Eigen::MatrixXi &domain_edges,
                                                   const std::vector<Eigen::MatrixXi> &domain_faces,
                                                   const unsigned int &num_iterations,
                                                   Eigen::MatrixXd &VoronoiPoints,
                                                   Gedim::IMeshDAO &mesh)
{
    const unsigned int num_points = VoronoiPoints.cols();
    if (num_points == 0)
    {
        std::vector<unsigned int> vertex_markers(domain_vertices.cols());
        std::iota(vertex_markers.begin(), vertex_markers.end(), 1);
        std::vector<unsigned int> edge_markers(domain_edges.cols());
        std::iota(edge_markers.begin(), edge_markers.end(), domain_vertices.cols() + 1);
        std::vector<unsigned int> face_markers(domain_faces.size());
        std::iota(face_markers.begin(), face_markers.end(), domain_vertices.cols() + domain_edges.cols() + 1);

        Gedim::MeshUtilities mesh_utilities;
        mesh_utilities.Mesh3DFromPolyhedron(domain_vertices, domain_edges, domain_faces, vertex_markers, edge_markers, face_markers, mesh);

        VoronoiPoints.resize(3, 1);
        VoronoiPoints.col(0) = geometryUtilities.PolyhedronBarycenter(domain_vertices);

        return;
    }

    const unsigned int num_vertices_domain = domain_vertices.cols();
    const unsigned int num_edges_domain = domain_edges.cols();
    const unsigned int num_faces_domain = domain_faces.size();

    vector<Eigen::Vector3d> edgeDomainTangent(num_edges_domain);
    vector<double> edgeDomainTangentSquaredNorm(num_edges_domain);
    vector<Eigen::Vector3d> edgeDomainLineOrigin(num_edges_domain);
    for (unsigned int e = 0; e < num_edges_domain; e++)
    {
        edgeDomainLineOrigin[e] = domain_vertices.col(domain_edges(0, e));
        const Eigen::Vector3d end = domain_vertices.col(domain_edges(1, e));

        Eigen::Vector3d tangent = end - edgeDomainLineOrigin[e];

        edgeDomainTangent[e] = tangent;
        edgeDomainTangentSquaredNorm[e] = tangent.squaredNorm();
    }

    vector<Eigen::Vector3d> faceDomainNormal(num_faces_domain);
    for (unsigned int f = 0; f < num_faces_domain; f++)
    {
        const unsigned int numVerticesFaceDomain = domain_faces[f].cols();
        Eigen::MatrixXd vert = Eigen::MatrixXd::Zero(3, numVerticesFaceDomain);

        for (unsigned int v = 0; v < numVerticesFaceDomain; v++)
            vert.col(v) = domain_vertices.col(domain_faces[f](0, v));

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(vert);

        faceDomainNormal[f] = normal;
    }

    std::vector<unsigned int> map_marker(num_faces_domain, 0);

    const double x_min = domain_vertices.row(0).minCoeff(), x_max = domain_vertices.row(0).maxCoeff();
    const double y_min = domain_vertices.row(1).minCoeff(), y_max = domain_vertices.row(1).maxCoeff();
    const double z_min = domain_vertices.row(2).minCoeff(), z_max = domain_vertices.row(2).maxCoeff();
    const double cvol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = n_x, n_z = n_x;

    // Lyod-iterations
    for (unsigned int it = 0; it < num_iterations; it++)
    {
        // Initialize Voronoi mesh
        voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);

        // Set Voronoi generators
        for (unsigned int i = 0; i < VoronoiPoints.cols(); i++)
            con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));

        unsigned int countPoints = 0;

        // Create an interator to visit all the cells
        voro::c_loop_all vl(con);
        voro::voronoicell_neighbor c;
        if (vl.start())
        {
            do
            {
                if (con.compute_cell(c, vl))
                {
                    const unsigned int ijk = vl.ijk; // id of the block -> 0,..., n_x, n_y, n_z
                    const unsigned int q = vl.q;     // id of the generator inside the block

                    // con.ps denotes how much memory the generator use
                    const double *pp = con.p[ijk] + con.ps * q;

                    // Cell barycenter
                    const double cx = pp[0];
                    const double cy = pp[1];
                    const double cz = pp[2];

                    // Get cell generator
                    std::vector<double> vertices;
                    c.vertices(cx, cy, cz, vertices);

                    // Get cell face cardinality
                    std::vector<int> faces_num_vertices;
                    c.face_orders(faces_num_vertices);

                    // Get face-vertices connectivity
                    std::vector<int> faces_vertices;
                    c.face_vertices(faces_vertices);

                    // cell vertices
                    const unsigned int num_cell_vertices = vertices.size() / 3;
                    Eigen::MatrixXd cell_vertices = Eigen::MatrixXd::Zero(3, num_cell_vertices);
                    for (unsigned int v = 0; v < num_cell_vertices; v++)
                        cell_vertices.col(v) << vertices[3 * v], vertices[3 * v + 1], vertices[3 * v + 2];

                    // cell faces
                    const unsigned int num_cell_faces = faces_num_vertices.size();
                    vector<Eigen::MatrixXi> Cell3DsFaces(num_cell_faces);

                    unsigned int count = 0;
                    for (unsigned int i = 0; i < num_cell_faces; i++)
                    {
                        const unsigned num_face_vertices = faces_num_vertices[i];
                        Cell3DsFaces[i] = Eigen::MatrixXi::Zero(2, num_face_vertices);

                        count++;
                        for (unsigned int j = 0; j < num_face_vertices; j++)
                        {
                            Cell3DsFaces[i](0, j) = faces_vertices[count];
                            count++;
                        }
                    }

                    VoronoiPoints.col(countPoints++) = geometryUtilities.PolyhedronBarycenter(cell_vertices);
                }
            } while (vl.inc());
        }
    }

    double vvol = 0.0;

    // Initialize Voronoi mesh
    voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8);

    for (unsigned int i = 0; i < VoronoiPoints.cols(); i++)
        con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));

    // Get cell number
    std::map<unsigned int, Cell0D> cell0Ds;
    std::map<unsigned int, Cell2D> cell2Ds;
    std::map<unsigned int, Cell3D> cell3Ds;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> origin_end_cell1Ds;
    std::map<unsigned int, Cell1D> cell1Ds;

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do
        {
            if (con.compute_cell(c, vl))
            {
                const unsigned int ijk = vl.ijk; // id of the block -> 0,..., n_x, n_y, n_z
                const unsigned int q = vl.q;     // id of the generator inside the block

                // con.ps denotes how much memory the generator use
                const double *pp = con.p[ijk] + con.ps * q;

                // Cell barycenter
                const double cx = pp[0];
                const double cy = pp[1];
                const double cz = pp[2];

                // Get cell3D ID
                const unsigned int cell3DIndex = con.id[ijk][q];

                // Get cell vertices
                std::vector<double> vertices;
                c.vertices(cx, cy, cz, vertices);

                // Get cell face cardinality
                std::vector<int> faces_num_vertices;
                c.face_orders(faces_num_vertices);

                // Get face-vertices connectivity
                std::vector<int> faces_vertices;
                c.face_vertices(faces_vertices);

                // Get cell neighbors IDs
                std::vector<int> faces_neighs;
                c.neighbors(faces_neighs);

                Cell3D cell3D;
                cell3D.marker = 0;

                // // Cell 3D vertices
                const unsigned int num_cell_vertices = vertices.size() / 3;
                Eigen::MatrixXd cell_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_cell_vertices);
                for (unsigned int v = 0; v < num_cell_vertices; v++)
                    cell_vertices_coordinates.col(v) =
                        (Eigen::Vector3d() << vertices[3 * v], vertices[3 * v + 1], vertices[3 * v + 2]).finished();

                const unsigned int num_cell_faces = faces_num_vertices.size();
                unsigned int count_v = 0;
                cell3D.vertices.resize(num_cell_vertices, std::numeric_limits<unsigned int>::max());
                cell3D.faces.resize(num_cell_faces);
                cell3D.neighbors.resize(num_cell_faces);
                for (unsigned int f = 0; f < num_cell_faces; f++)
                {
                    const unsigned num_face_vertices = faces_num_vertices[f];
                    Eigen::MatrixXd face_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_face_vertices);
                    std::vector<unsigned int> local_id_vertices(num_face_vertices);
                    cell3D.neighbors[f] = faces_neighs[f];

                    count_v++;
                    for (unsigned int j = 0; j < num_face_vertices; j++)
                    {
                        face_vertices_coordinates.col(j) = cell_vertices_coordinates.col(faces_vertices[count_v]);
                        local_id_vertices[j] = faces_vertices[count_v];
                        count_v++;
                    }

                    if (faces_neighs[f] >= 0)
                    {
                        const auto it = cell3Ds.find(static_cast<unsigned int>(faces_neighs[f]));
                        if (it != cell3Ds.end())
                        {
                            const unsigned int neighbors_index = it->first;
                            const Cell3D neighbor = cell3Ds[neighbors_index];

                            bool found_face = false;
                            for (unsigned int n_f = 0; n_f < neighbor.faces.size(); n_f++)
                            {
                                if (neighbor.neighbors[n_f] == static_cast<int>(cell3DIndex))
                                {
                                    // Insert face
                                    cell3D.faces[f] = neighbor.faces[n_f];

                                    // Insert edges
                                    const Cell2D face = cell2Ds[cell3D.faces[f]];

                                    if (face.vertices.size() != num_face_vertices)
                                    {

                                        // {
                                        //     std::string dbg_export_folder = "./TEST_DBG";
                                        //     Gedim::Output::CreateFolder(dbg_export_folder);

                                        //     {
                                        //         Gedim::VTKUtilities exporter;
                                        //         exporter.AddPoints(cell_vertices_coordinates);
                                        //         exporter.Export(dbg_export_folder + "/cell3D_points.vtu");
                                        //     }

                                        //     {
                                        //         Eigen::MatrixXd neigh_cell_vertices_coordinates =
                                        //         Eigen::MatrixXd::Zero(3, neighbor.vertices.size()); for (unsigned int
                                        //         j = 0; j < neighbor.vertices.size(); j++)
                                        //         {
                                        //             const auto c0D = cell0Ds.at(neighbor.vertices.at(j));
                                        //             neigh_cell_vertices_coordinates.col(j) = c0D.coordinates;
                                        //         }

                                        //         Gedim::VTKUtilities exporter;
                                        //         exporter.AddPoints(neigh_cell_vertices_coordinates);
                                        //         exporter.Export(dbg_export_folder + "/neigh_cell3D_points.vtu");
                                        //     }

                                        //     {
                                        //         Gedim::VTKUtilities exporter;
                                        //         exporter.AddPoints(face_vertices_coordinates);
                                        //         exporter.Export(dbg_export_folder +
                                        //                         "/cell3D_face_" +
                                        //                         std::to_string(f) +
                                        //                         "_points.vtu");
                                        //     }

                                        //     {
                                        //         Eigen::MatrixXd neigh_face_vertices_coordinates =
                                        //         Eigen::MatrixXd::Zero(3, face.vertices.size()); for (unsigned int j =
                                        //         0; j < face.vertices.size(); j++)
                                        //         {
                                        //             const auto c0D = cell0Ds.at(face.vertices.at(j));
                                        //             neigh_face_vertices_coordinates.col(j) = c0D.coordinates;
                                        //         }

                                        //         {
                                        //             Gedim::VTKUtilities exporter;
                                        //             exporter.AddPoints(neigh_face_vertices_coordinates);
                                        //             exporter.Export(dbg_export_folder +
                                        //                             "/cell3D_neigh_face_" +
                                        //                             std::to_string(n_f) +
                                        //                             "_points.vtu");
                                        //         }
                                        //     }
                                        // }

                                        throw std::runtime_error("not valid neigh face");
                                    }

                                    for (unsigned int e = 0; e < face.edges.size(); e++)
                                        cell3D.edges.insert(face.edges[e]);

                                    // Find vertices
                                    for (unsigned int v_f = 0; v_f < num_face_vertices; v_f++)
                                    {
                                        bool found_vertex = false;
                                        std::pair<unsigned int, double> candidate_point = {
                                            std::numeric_limits<unsigned int>::max(),
                                            std::numeric_limits<double>::max()};
                                        for (unsigned int v_c = 0; v_c < num_face_vertices; v_c++)
                                        {
                                            if (cell3D.vertices[local_id_vertices[v_c]] == std::numeric_limits<unsigned int>::max())
                                            {
                                                if (geometryUtilities.PointsAreCoincident(cell0Ds.at(face.vertices[v_f]).coordinates,
                                                                                          face_vertices_coordinates.col(v_c)))
                                                {
                                                    cell3D.vertices[local_id_vertices[v_c]] = face.vertices[v_f];
                                                    found_vertex = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    const double dist = (cell0Ds.at(face.vertices[v_f]).coordinates -
                                                                         face_vertices_coordinates.col(v_c))
                                                                            .norm();
                                                    if (candidate_point.second > dist)
                                                    {
                                                        candidate_point.second = dist;
                                                        candidate_point.first = local_id_vertices[v_c];
                                                    }
                                                }
                                            }
                                            else if (cell3D.vertices[local_id_vertices[v_c]] == face.vertices[v_f])
                                            {
                                                found_vertex = true;
                                                break;
                                            }
                                        }

                                        if (!found_vertex)
                                        {
                                            if (candidate_point.first != std::numeric_limits<unsigned int>::max())
                                            {
                                                cell3D.vertices[candidate_point.first] = face.vertices[v_f];
                                            }
                                            else
                                                throw runtime_error("ID vertex not found");
                                        }
                                    }

                                    found_face = true;
                                    break;
                                }
                            }

                            if (!found_face)
                                throw runtime_error("ID face not found");
                        }
                    }
                }

                // Cell 3D vertices
                for (unsigned int v = 0; v < num_cell_vertices; v++)
                {
                    if (cell3D.vertices[v] == std::numeric_limits<unsigned int>::max())
                    {
                        const unsigned int numCell0Ds = cell0Ds.size();
                        Cell0D cell0D(cell_vertices_coordinates.col(v));
                        cell0Ds.insert({numCell0Ds, cell0D});

                        cell3D.vertices[v] = numCell0Ds;
                    }
                }

                // Cell 3D faces
                count_v = 0;
                for (unsigned int f = 0; f < num_cell_faces; f++)
                {
                    const unsigned num_face_vertices = faces_num_vertices[f];
                    std::vector<unsigned int> face_vertices(num_face_vertices);
                    Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, num_face_vertices);
                    cell3D.neighbors[f] = faces_neighs[f];

                    count_v++;
                    for (unsigned int j = 0; j < num_face_vertices; j++)
                    {
                        face_vertices[j] = cell3D.vertices[faces_vertices[count_v]];
                        const auto it = cell0Ds.find(face_vertices[j]);

                        if (it == cell0Ds.end())
                            throw std::runtime_error("not valid face");

                        faceVertCoordinates.col(j) = it->second.coordinates;
                        count_v++;
                    }

                    if (faces_neighs[f] < 0) // if it is a boundary face
                    {
                        // prima volta che incontro la faccia  -  nuovi lati possibili - faccia di bordo

                        const unsigned int cell2DIndex = cell2Ds.size();

                        Cell2D cell2D;
                        cell2D.vertices = face_vertices;
                        cell2D.marker = 1;

                        if (map_marker[std::abs(faces_neighs[f]) - 1] == 0)
                        {
                            const Eigen::Vector3d barycenter = faceVertCoordinates.rowwise().mean();

                            for (unsigned int df = 0; df < num_faces_domain; df++)
                            {

                                if (geometryUtilities.IsPointOnPlane(barycenter,
                                                                     faceDomainNormal[df],
                                                                     domain_vertices.col(domain_faces[df](0, 0))))
                                {
                                    map_marker[std::abs(faces_neighs[f]) - 1] = num_vertices_domain + num_edges_domain + df + 1;
                                    cell2D.marker = map_marker[std::abs(faces_neighs[f]) - 1];
                                    break;
                                }
                            }

                            if (map_marker[std::abs(faces_neighs[f]) - 1] == 0)
                                throw runtime_error("Marker face wrong");
                        }
                        else
                            cell2D.marker = map_marker[std::abs(faces_neighs[f]) - 1];

                        cell2D.edges.resize(num_face_vertices);
                        for (unsigned int v = 0; v < num_face_vertices; v++)
                        {
                            const unsigned int id_min = min(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);
                            const unsigned int id_max = max(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);

                            const auto it = origin_end_cell1Ds.find({id_min, id_max});

                            if (it != origin_end_cell1Ds.end())
                                cell2D.edges[v] = it->second;
                            else
                            {
                                const unsigned int cell1DIndex = cell1Ds.size();
                                origin_end_cell1Ds.insert({{id_min, id_max}, cell1DIndex});

                                Cell1D cell1D;
                                cell1D.origin = id_min;
                                cell1D.end = id_max;
                                cell1D.marker = 1;
                                cell1Ds.insert({cell1DIndex, cell1D});

                                cell2D.edges[v] = cell1DIndex;
                            }
                        }

                        cell2Ds.insert({cell2DIndex, cell2D});

                        cell3D.faces[f] = cell2DIndex;
                        for (unsigned int e = 0; e < cell2D.edges.size(); e++)
                            cell3D.edges.insert(cell2D.edges[e]);
                    }
                    else
                    {
                        const auto it = cell3Ds.find(static_cast<unsigned int>(faces_neighs[f]));
                        if (it == cell3Ds.end())
                        {
                            // prima volta che incontro la faccia -  nuovi lati possibili
                            const unsigned int cell2DIndex = cell2Ds.size();

                            Cell2D cell2D;
                            cell2D.marker = 0;
                            cell2D.vertices = face_vertices;

                            cell2D.edges.resize(num_face_vertices);
                            for (unsigned int v = 0; v < num_face_vertices; v++)
                            {
                                const unsigned int id_min = min(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);
                                const unsigned int id_max = max(face_vertices[v], face_vertices[(v + 1) % num_face_vertices]);

                                const auto it = origin_end_cell1Ds.find({id_min, id_max});

                                if (it != origin_end_cell1Ds.end())
                                    cell2D.edges[v] = it->second;
                                else
                                {
                                    const unsigned int cell1DIndex = cell1Ds.size();
                                    origin_end_cell1Ds.insert({{id_min, id_max}, cell1DIndex});

                                    Cell1D cell1D;
                                    cell1D.origin = id_min;
                                    cell1D.end = id_max;
                                    cell1D.marker = 0;
                                    cell1Ds.insert({cell1DIndex, cell1D});

                                    cell2D.edges[v] = cell1DIndex;
                                }
                            }

                            cell2Ds.insert({cell2DIndex, cell2D});

                            cell3D.faces[f] = cell2DIndex;
                            for (unsigned int e = 0; e < cell2D.edges.size(); e++)
                                cell3D.edges.insert(cell2D.edges[e]);
                        }
                    }
                }

                cell3Ds.insert({cell3DIndex, cell3D});

                vvol += c.volume();
            }
        } while (vl.inc());
    }

    /// <li> Set Cell0Ds
    const unsigned int numberOfPointsMesh = cell0Ds.size();
    const unsigned int numberOfEdgesMesh = cell1Ds.size();
    const unsigned int numberOfFacesMesh = cell2Ds.size();
    const unsigned int numberOfCellsMesh = cell3Ds.size();

    mesh.InitializeDimension(3);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    mesh.Cell2DsInitialize(numberOfFacesMesh);
    mesh.Cell3DsInitialize(numberOfCellsMesh);

    for (const auto &cell0D : cell0Ds)
    {
        const unsigned int id = cell0D.first;

        mesh.Cell0DSetState(id, true);
        mesh.Cell0DInsertCoordinates(id, cell0D.second.coordinates);

        const Eigen::Vector3d coordinates = cell0D.second.coordinates;

        const bool is_xmin = (geometryUtilities.IsValueZero(coordinates(0) - x_min, geometryUtilities.Tolerance1D())); // l ==
        // 0
        const bool is_xmax = (geometryUtilities.IsValueZero(coordinates(0) - x_max, geometryUtilities.Tolerance1D())); // l ==
        // (numLengthPoints
        // - 1)
        const bool is_ymin = (geometryUtilities.IsValueZero(coordinates(1) - y_min, geometryUtilities.Tolerance1D())); // h ==
        // 0
        const bool is_ymax = (geometryUtilities.IsValueZero(coordinates(1) - y_max, geometryUtilities.Tolerance1D())); // h ==
        // (numHeightPoints
        // - 1)
        const bool is_zmin = (geometryUtilities.IsValueZero(coordinates(2) - z_min, geometryUtilities.Tolerance1D())); // w ==
        // 0
        const bool is_zmax = (geometryUtilities.IsValueZero(coordinates(2) - z_max, geometryUtilities.Tolerance1D())); // w ==
        // (numWidthPoints
        // - 1)

        const unsigned int marker =
            1 * (is_xmin && is_ymin && is_zmin) + 2 * (is_xmax && is_ymin && is_zmin) +
            3 * (is_xmax && is_ymax && is_zmin) + 4 * (is_xmin && is_ymax && is_zmin) + 5 * (is_xmin && is_ymin && is_zmax) +
            6 * (is_xmax && is_ymin && is_zmax) + 7 * (is_xmax && is_ymax && is_zmax) + 8 * (is_xmin && is_ymax && is_zmax) +
            9 * (!is_xmin && !is_xmax && is_ymin && is_zmin) + 10 * (is_xmax && !is_ymin && !is_ymax && is_zmin) +
            11 * (!is_xmin && !is_xmax && is_ymax && is_zmin) + 12 * (is_xmin && !is_ymin && !is_ymax && is_zmin) +
            13 * (!is_xmin && !is_xmax && is_ymin && is_zmax) + 14 * (is_xmax && !is_ymin && !is_ymax && is_zmax) +
            15 * (!is_xmin && !is_xmax && is_ymax && is_zmax) + 16 * (is_xmin && !is_ymin && !is_ymax && is_zmax) +
            17 * (is_xmin && is_ymin && !is_zmin && !is_zmax) + 18 * (is_xmax && is_ymin && !is_zmin && !is_zmax) +
            19 * (is_xmax && is_ymax && !is_zmin && !is_zmax) + 20 * (is_xmin && is_ymax && !is_zmin && !is_zmax) +
            21 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmin) +
            22 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmax) +
            23 * (is_xmin && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
            24 * (is_xmax && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
            25 * (!is_xmin && !is_xmax && is_ymin && !is_zmin && !is_zmax) +
            26 * (!is_xmin && !is_xmax && is_ymax && !is_zmin && !is_zmax);

        mesh.Cell0DSetMarker(id, marker);
    }

    /// <li> Set Edges
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    for (const auto &cell1D : cell1Ds)
    {
        const unsigned int e = cell1D.first;

        mesh.Cell1DSetState(e, true);
        mesh.Cell1DInsertExtremes(e, cell1D.second.origin, cell1D.second.end);

        const unsigned int marker_origin = mesh.Cell0DMarker(cell1D.second.origin);
        const unsigned int marker_end = mesh.Cell0DMarker(cell1D.second.end);
        if (marker_origin == 0 || marker_end == 0)
            mesh.Cell1DSetMarker(e, 0);
        else if (marker_origin > 8 && marker_end > 8)
            mesh.Cell1DSetMarker(e, std::max(marker_origin, marker_end));
        else
        {
            unsigned int marker = 0;
            for (unsigned int e = 0; e < num_edges_domain; e++)
            {

                if (geometryUtilities.IsPointOnLine(mesh.Cell0DCoordinates(cell1D.second.origin),
                                                    edgeDomainLineOrigin[e],
                                                    edgeDomainTangent[e],
                                                    edgeDomainTangentSquaredNorm[e]) &&
                    geometryUtilities.IsPointOnLine(mesh.Cell0DCoordinates(cell1D.second.end),
                                                    edgeDomainLineOrigin[e],
                                                    edgeDomainTangent[e],
                                                    edgeDomainTangentSquaredNorm[e]))
                {
                    marker = num_vertices_domain + e + 1;
                    break;
                }
            }

            if (marker == 0)
                throw std::runtime_error("not considered case");

            mesh.Cell1DSetMarker(e, marker);
        }
    }

    /// <li> Set Faces
    for (const auto &cell2D : cell2Ds)
    {
        const unsigned int f = cell2D.first;

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, cell2D.second.marker);

        const unsigned int numCellVertices = cell2D.second.vertices.size();
        const unsigned int numCellEdges = cell2D.second.edges.size();

        Gedim::Output::Assert(numCellVertices == numCellEdges);

        mesh.Cell2DInitializeVertices(f, numCellVertices);
        mesh.Cell2DInitializeEdges(f, numCellEdges);

        for (unsigned int v = 0; v < numCellVertices; v++)
            mesh.Cell2DInsertVertex(f, v, cell2D.second.vertices[v]);

        for (unsigned int e = 0; e < numCellEdges; e++)
            mesh.Cell2DInsertEdge(f, e, cell2D.second.edges[e]);
    }

    /// <li> Set Cells
    unsigned int cell3DIndex = 0;
    double gvol = 0.0;
    for (const auto &cell3D : cell3Ds)
    {
        const unsigned int numVertices = cell3D.second.vertices.size();
        const unsigned int numEdges = cell3D.second.edges.size();
        const unsigned int numFaces = cell3D.second.faces.size();

        if (numVertices == 0)
            continue;

        mesh.Cell3DSetState(cell3DIndex, true);
        mesh.Cell3DSetMarker(cell3DIndex, 0);

        mesh.Cell3DInitializeVertices(cell3DIndex, numVertices);
        mesh.Cell3DInitializeEdges(cell3DIndex, numEdges);
        mesh.Cell3DInitializeFaces(cell3DIndex, numFaces);

        for (unsigned int v = 0; v < numVertices; v++)
            mesh.Cell3DInsertVertex(cell3DIndex, v, cell3D.second.vertices[v]);

        unsigned int e = 0;
        for (const auto &edge_id : cell3D.second.edges)
        {
            mesh.Cell3DInsertEdge(cell3DIndex, e, edge_id);
            e++;
        }

        for (unsigned int f = 0; f < numFaces; f++)
            mesh.Cell3DInsertFace(cell3DIndex, f, cell3D.second.faces[f]);

        {
            const Eigen::MatrixXd Cell3DsVertices = mesh.Cell3DVerticesCoordinates(cell3DIndex);
            std::vector<Eigen::MatrixXi> Cell3DsFaces(numFaces);
            for (unsigned int f = 0; f < numFaces; f++)
            {
                const unsigned int faceId = mesh.Cell3DFace(cell3DIndex, f);
                Cell3DsFaces[f] = Eigen::MatrixXi::Zero(2, mesh.Cell2DNumberEdges(faceId));
                for (unsigned int v = 0; v < mesh.Cell2DNumberEdges(faceId); v++)
                {
                    Cell3DsFaces[f](0, v) = mesh.Cell3DFindVertex(cell3DIndex, mesh.Cell2DVertex(faceId, v));
                    Cell3DsFaces[f](1, v) = mesh.Cell3DFindEdge(cell3DIndex, mesh.Cell2DEdge(faceId, v));
                }
            }

            const std::vector<Eigen::MatrixXd> Cell3DsFaces3DVertices =
                geometryUtilities.PolyhedronFaceVertices(Cell3DsVertices, Cell3DsFaces);
            const std::vector<Eigen::Vector3d> Cell3DsFacesTranslations =
                geometryUtilities.PolyhedronFaceTranslations(Cell3DsFaces3DVertices);
            const std::vector<Eigen::Vector3d> Cell3DsFacesNormals = geometryUtilities.PolyhedronFaceNormals(Cell3DsFaces3DVertices);
            const std::vector<Eigen::Matrix3d> Cell3DsFacesRotationMatrices =
                geometryUtilities.PolyhedronFaceRotationMatrices(Cell3DsFaces3DVertices, Cell3DsFacesNormals, Cell3DsFacesTranslations);

            const vector<vector<unsigned int>> polyhedronFaceTriangulations =
                geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(Cell3DsFaces, Cell3DsFaces3DVertices);

            const std::vector<Eigen::MatrixXd> Cell3DsFaces2DVertices =
                geometryUtilities.PolyhedronFaceRotatedVertices(Cell3DsFaces3DVertices, Cell3DsFacesTranslations, Cell3DsFacesRotationMatrices);

            const std::vector<std::vector<Eigen::Matrix3d>> Cell3DsFaces2DTriangulations =
                geometryUtilities.PolyhedronFaceExtractTriangulationPoints(Cell3DsFaces2DVertices, polyhedronFaceTriangulations);

            const std::vector<bool> Cell3DsFacesNormalDirections =
                geometryUtilities.PolyhedronFaceNormalDirections(Cell3DsFaces3DVertices,
                                                                 geometryUtilities.PolyhedronBarycenter(Cell3DsVertices),
                                                                 Cell3DsFacesNormals);

            gvol += geometryUtilities.PolyhedronVolumeByBoundaryIntegral(Cell3DsFaces2DTriangulations,
                                                                         Cell3DsFacesNormals,
                                                                         Cell3DsFacesNormalDirections,
                                                                         Cell3DsFacesTranslations,
                                                                         Cell3DsFacesRotationMatrices);
        }

        cell3DIndex++;
    }

    // if (!geometryUtilities.AreValuesEqual(cvol, gvol, geometryUtilities.Tolerance3D()) ||
    //     !geometryUtilities.AreValuesEqual(cvol, vvol, geometryUtilities.Tolerance3D()))
    // {
    //     std::cout.precision(16);
    //     std::cout << std::scientific << cvol << " " << vvol << " " << gvol << std::endl;

    //     // throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");
    // }

    if (!geometryUtilities.AreValuesEqual(cvol, gvol, geometryUtilities.Tolerance3D()))
        throw runtime_error("Error generating Voronoi cells: volumes do not match each other");
}
// ************************************************************************* //
Eigen::MatrixXd VoroInterface::GenerateRandomPoints(const Eigen::MatrixXd &domainVertices,
                                                    const unsigned int &numPoints,
                                                    const unsigned int random_seed)
{
    const double x_min = domainVertices.row(0).minCoeff(), x_max = domainVertices.row(0).maxCoeff();
    const double y_min = domainVertices.row(1).minCoeff(), y_max = domainVertices.row(1).maxCoeff();
    const double z_min = domainVertices.row(2).minCoeff(), z_max = domainVertices.row(2).maxCoeff();

    srand(random_seed);

    Eigen::MatrixXd VoronoiPoints = Eigen::MatrixXd::Zero(3, numPoints);
    for (unsigned int i = 0; i < numPoints; i++)
    {
        VoronoiPoints(0, i) = x_min + rnd() * (x_max - x_min);
        VoronoiPoints(1, i) = y_min + rnd() * (y_max - y_min);
        VoronoiPoints(2, i) = z_min + rnd() * (z_max - z_min);
    }

    return VoronoiPoints;
}
// ************************************************************************* //
void VoroInterface::GenerateCartesianPoints3D(const Eigen::MatrixXd &polyhedronVertices, const unsigned int &numPoints, voro::container &con)
{
    Eigen::Vector3d parallelepipedLengthTangent = polyhedronVertices.col(1) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedHeightTangent = polyhedronVertices.col(3) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedWidthTangent = polyhedronVertices.col(4) - polyhedronVertices.col(0);

    vector<double> lengthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1, 0.0, 1.0, 1);
    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1, 0.0, 1.0, 1);
    vector<double> widthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1, 0.0, 1.0, 1);

    const unsigned int numLengthPoints = lengthMeshCurvilinearCoordinates.size();
    const unsigned int numHeightPoints = heightMeshCurvilinearCoordinates.size();
    const unsigned int numWidthPoints = widthMeshCurvilinearCoordinates.size();

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
        for (unsigned int h = 0; h < numHeightPoints; h++)
        {
            for (unsigned int l = 0; l < numLengthPoints; l++)
            {
                const Eigen::Vector3d coordinate = polyhedronVertices.col(0) +
                                                   lengthMeshCurvilinearCoordinates[l] * parallelepipedLengthTangent +
                                                   heightMeshCurvilinearCoordinates[h] * parallelepipedHeightTangent +
                                                   widthMeshCurvilinearCoordinates[w] * parallelepipedWidthTangent;
                con.put(cell0DIndex, coordinate(0), coordinate(1), coordinate(2));
                cell0DIndex++;
            }
        }
    }
}
#endif
// ************************************************************************* //
} // namespace Gedim
