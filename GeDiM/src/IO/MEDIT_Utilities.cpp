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

#include "MEDIT_Utilities.hpp"
#include <fstream>

using namespace std;

namespace Gedim
{
// ***************************************************************************
void MEDIT_Utilities::ExportPoints(const std::string &filePath,
                                   const Eigen::MatrixXd &points,
                                   const std::vector<unsigned int> &points_reference_ids) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportMEDIT_Ascii(points, points_reference_ids, {}, {}, filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void MEDIT_Utilities::ExportSegments(const std::string &filePath,
                                     const Eigen::MatrixXd &points,
                                     const Eigen::MatrixXi &cells,
                                     const std::vector<unsigned int> &points_reference_ids,
                                     const std::vector<unsigned int> &cells_reference_ids) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportMEDIT_Ascii(points, points_reference_ids, {}, CreateEdgeCells(cells, cells_reference_ids), filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void MEDIT_Utilities::ExportPolygons(const std::string &filePath,
                                     const Eigen::MatrixXd &points,
                                     const std::vector<std::vector<unsigned int>> &cells_vertices,
                                     const std::vector<unsigned int> &points_reference_ids,
                                     const std::vector<unsigned int> &cells_reference_ids,
                                     const unsigned int &order) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportMEDIT_Ascii(points, points_reference_ids, {}, CreatePolygonCells(cells_vertices, cells_reference_ids, order), filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void MEDIT_Utilities::ExportPolyhedra(const std::string &filePath,
                                      const Eigen::MatrixXd &points,
                                      const std::vector<std::vector<unsigned int>> &faces_vertices,
                                      const std::vector<std::vector<unsigned int>> &cells_vertices,
                                      const std::vector<unsigned int> &points_reference_ids,
                                      const std::vector<unsigned int> &faces_reference_ids,
                                      const std::vector<unsigned int> &cells_reference_ids,
                                      const unsigned int &order) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportMEDIT_Ascii(points,
                          points_reference_ids,
                          CreatePolygonCells(faces_vertices, faces_reference_ids, order),
                          CreatePolyhedraCells(cells_vertices, cells_reference_ids, order),
                          filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
std::vector<MEDIT_Cell> MEDIT_Utilities::CreateVertexCells(const Eigen::MatrixXd &points,
                                                           const std::vector<unsigned int> &reference_ids) const
{
    std::vector<MEDIT_Cell> cells;
    cells.reserve(points.cols());

    for (unsigned int p = 0; p < points.cols(); p++)
    {
        cells.push_back(MEDIT_Cell(MEDIT_Cell::Types::Vertex,
                                   {p},
                                   reference_ids.size() == static_cast<unsigned int>(points.cols()) ? reference_ids[p] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<MEDIT_Cell> MEDIT_Utilities::CreateEdgeCells(const Eigen::MatrixXi &edges, const std::vector<unsigned int> &reference_ids) const
{
    std::vector<MEDIT_Cell> cells;
    cells.reserve(edges.cols());

    for (unsigned int l = 0; l < edges.cols(); l++)
    {
        cells.push_back(MEDIT_Cell(MEDIT_Cell::Types::Edge,
                                   {static_cast<unsigned int>(edges(0, l)), static_cast<unsigned int>(edges(1, l))},
                                   reference_ids.size() == static_cast<unsigned int>(edges.cols()) ? reference_ids[l] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<MEDIT_Cell> MEDIT_Utilities::CreatePolygonCells(const std::vector<std::vector<unsigned int>> &polygons_vertices,
                                                            const std::vector<unsigned int> &reference_ids,
                                                            const unsigned int &order) const
{
    std::vector<MEDIT_Cell> cells;
    cells.reserve(polygons_vertices.size());

    for (unsigned int l = 0; l < polygons_vertices.size(); l++)
    {
        MEDIT_Cell::Types polygon_type = MEDIT_Cell::Types::Unknown;

        if (polygons_vertices[l].size() == 3 * order)
            polygon_type = MEDIT_Cell::Types::Triangle;
        else if (polygons_vertices[l].size() == 4 * order)
            polygon_type = MEDIT_Cell::Types::Quadrilateral;
        else
            throw std::runtime_error("Polygon type not supported");

        cells.push_back(
            MEDIT_Cell(polygon_type,
                       polygons_vertices[l],
                       reference_ids.size() == static_cast<unsigned int>(polygons_vertices.size()) ? reference_ids[l] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<MEDIT_Cell> MEDIT_Utilities::CreatePolyhedraCells(const std::vector<std::vector<unsigned int>> &polyhedra_vertices,
                                                              const std::vector<unsigned int> &reference_ids,
                                                              const unsigned int &order) const
{
    std::vector<MEDIT_Cell> cells;
    cells.reserve(polyhedra_vertices.size());

    for (unsigned int l = 0; l < polyhedra_vertices.size(); l++)
    {
        MEDIT_Cell::Types polyhedra_type = MEDIT_Cell::Types::Unknown;

        if ((polyhedra_vertices[l].size() == 4 && order == 1) || (polyhedra_vertices[l].size() == 10 && order == 2))
            polyhedra_type = MEDIT_Cell::Types::Tetrahedron;
        else if ((polyhedra_vertices[l].size() == 8 && order == 1) || (polyhedra_vertices[l].size() == 20 && order == 2))
            polyhedra_type = MEDIT_Cell::Types::Hexahedron;
        else
            throw std::runtime_error("Polyhedron type not supported");

        cells.push_back(
            MEDIT_Cell(polyhedra_type,
                       polyhedra_vertices.at(l),
                       reference_ids.size() == static_cast<unsigned int>(polyhedra_vertices.size()) ? reference_ids[l] : 0));
    }

    return cells;
}
// ***************************************************************************
void MEDIT_Utilities::ExportMEDIT_Ascii(const Eigen::MatrixXd &points,
                                        const std::vector<unsigned int> &points_reference_ids,
                                        const std::vector<MEDIT_Cell> &faces,
                                        const std::vector<MEDIT_Cell> &cells,
                                        const std::string &filePath) const
{
    std::ofstream file;

    file.open(filePath.c_str());

    if (file.fail())
        throw runtime_error("File '" + filePath + "' cannot be opened");

    const char sep = ' ';
    file.precision(16);

    // header
    file << "MeshVersionFormatted 2" << std::endl; // file writing with double precision
    file << "Dimension" << std::endl;
    file << "3" << std::endl;

    // export vertices
    if (points.cols() > 0)
    {
        file << MEDIT_Cell::CellLabel(MEDIT_Cell::Types::Vertex) << std::endl;
        file << points.cols() << std::endl;

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            file << std::scientific << points(0, p) << sep;
            file << std::scientific << points(1, p) << sep;
            file << std::scientific << points(2, p) << sep;
            file << std::scientific << points_reference_ids.at(p) << std::endl;
        }
    }

    // export faces
    if (faces.size() > 0)
    {
        file << MEDIT_Cell::CellLabel(faces.at(0).Type) << std::endl;
        file << faces.size() << std::endl;
        for (unsigned int f = 0; f < faces.size(); f++)
        {
            const auto &face = faces.at(f);
            for (unsigned int p = 0; p < face.Points_id.size(); p++)
                file << (face.Points_id[p] + 1) << sep;
            file << face.Reference_id << std::endl;
        }
    }

    // export cells
    if (cells.size() > 0)
    {
        file << MEDIT_Cell::CellLabel(cells.at(0).Type) << std::endl;
        file << cells.size() << std::endl;
        for (unsigned int c = 0; c < cells.size(); c++)
        {
            const auto &cell = cells.at(c);
            for (unsigned int p = 0; p < cell.Points_id.size(); p++)
                file << (cell.Points_id[p] + 1) << sep;
            file << cell.Reference_id << std::endl;
        }
    }

    file << "End" << std::endl;
    file.close();
}
// ***************************************************************************
std::string MEDIT_Cell::CellLabel(const MEDIT_Cell::Types type)
{
    switch (type)
    {
    case MEDIT_Cell::Types::Vertex:
        return "Vertices";
    case MEDIT_Cell::Types::Edge:
        return "Edges";
    case MEDIT_Cell::Types::Triangle:
        return "Triangles";
    case MEDIT_Cell::Types::Quadrilateral:
        return "Quadrilaterals";
    case MEDIT_Cell::Types::Tetrahedron:
        return "Tetrahedra";
    case MEDIT_Cell::Types::Hexahedron:
        return "Hexahedra";
    default:
        throw std::runtime_error("Type not supported");
    }
}

// ***************************************************************************
} // namespace Gedim
