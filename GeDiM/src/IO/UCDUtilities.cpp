#include "UCDUtilities.hpp"
#include <fstream>

using namespace std;

namespace Gedim
{
// ***************************************************************************
void UCDUtilities::ExportPoints(const std::string &filePath,
                                const Eigen::MatrixXd &points,
                                const std::vector<UCDProperty<double>> &points_properties,
                                const Eigen::VectorXi &materials) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportUCDAscii(points, {}, CreatePointCells(points, materials), points_properties, filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void UCDUtilities::ExportSegments(const std::string &filePath,
                                  const Eigen::MatrixXd &points,
                                  const Eigen::MatrixXi &segments,
                                  const std::vector<UCDProperty<double>> &points_properties,
                                  const std::vector<UCDProperty<double>> &segmnents_properties,
                                  const Eigen::VectorXi &materials) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportUCDAscii(points, points_properties, CreateLineCells(segments, materials), segmnents_properties, filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void UCDUtilities::ExportPolygons(const std::string &filePath,
                                  const Eigen::MatrixXd &points,
                                  const std::vector<std::vector<unsigned int>> &polygons_vertices,
                                  const std::vector<UCDProperty<double>> &points_properties,
                                  const std::vector<UCDProperty<double>> &polygons_properties,
                                  const Eigen::VectorXi &materials) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportUCDAscii(points, points_properties, CreatePolygonCells(polygons_vertices, materials), polygons_properties, filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
void UCDUtilities::ExportPolyhedra(const std::string &filePath,
                                   const Eigen::MatrixXd &points,
                                   const std::vector<std::vector<unsigned int>> &polyhedra_vertices,
                                   const std::vector<UCDProperty<double>> &points_properties,
                                   const std::vector<UCDProperty<double>> &polyhedra_properties,
                                   const Eigen::VectorXi &materials) const
{
    ExportFormats format = ExportFormats::Ascii;

    switch (format)
    {
    case ExportFormats::Ascii:
        ExportUCDAscii(points, points_properties, CreatePolyhedraCells(polyhedra_vertices, materials), polyhedra_properties, filePath);
        break;
    default:
        throw std::runtime_error("Unknown format");
    }
}
// ***************************************************************************
std::vector<UCDCell> UCDUtilities::CreatePointCells(const Eigen::MatrixXd &points, const Eigen::VectorXi &materials) const
{
    std::vector<UCDCell> cells;
    cells.reserve(points.cols());

    for (unsigned int p = 0; p < points.cols(); p++)
    {
        cells.push_back(UCDCell(UCDCell::Types::Point, {p}, materials.size() == points.cols() ? materials[p] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<UCDCell> UCDUtilities::CreateLineCells(const Eigen::MatrixXi &lines, const Eigen::VectorXi &materials) const
{
    std::vector<UCDCell> cells;
    cells.reserve(lines.cols());

    for (unsigned int l = 0; l < lines.cols(); l++)
    {
        cells.push_back(UCDCell(UCDCell::Types::Line,
                                {static_cast<unsigned int>(lines(0, l)), static_cast<unsigned int>(lines(1, l))},
                                materials.size() == lines.cols() ? materials[l] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<UCDCell> UCDUtilities::CreatePolygonCells(const std::vector<std::vector<unsigned int>> &polygons_vertices,
                                                      const Eigen::VectorXi &materials) const
{
    std::vector<UCDCell> cells;
    cells.reserve(polygons_vertices.size());

    for (unsigned int l = 0; l < polygons_vertices.size(); l++)
    {
        UCDCell::Types polygon_type = UCDCell::Types::Unknown;

        if (polygons_vertices[l].size() == 3)
            polygon_type = UCDCell::Types::Triangle;
        else if (polygons_vertices[l].size() == 4)
            polygon_type = UCDCell::Types::Quadrilateral;
        else
            throw std::runtime_error("Polygon type not supported");

        cells.push_back(UCDCell(polygon_type,
                                polygons_vertices[l],
                                materials.size() == static_cast<unsigned int>(polygons_vertices.size()) ? materials[l] : 0));
    }

    return cells;
}
// ***************************************************************************
std::vector<UCDCell> UCDUtilities::CreatePolyhedraCells(const std::vector<std::vector<unsigned int>> &polyhedra_vertices,
                                                        const Eigen::VectorXi &materials) const
{
    std::vector<UCDCell> cells;
    cells.reserve(polyhedra_vertices.size());

    for (unsigned int l = 0; l < polyhedra_vertices.size(); l++)
    {
        UCDCell::Types polyhedra_type = UCDCell::Types::Unknown;

        if (polyhedra_vertices[l].size() == 4)
            polyhedra_type = UCDCell::Types::Tetrahedron;
        else
            throw std::runtime_error("Polygon type not supported");

        cells.push_back(UCDCell(polyhedra_type,
                                polyhedra_vertices[l],
                                materials.size() == static_cast<unsigned int>(polyhedra_vertices.size()) ? materials[l] : 0));
    }

    return cells;
}
// ***************************************************************************
void UCDUtilities::ExportUCDAscii(const Eigen::MatrixXd &points,
                                  const std::vector<UCDProperty<double>> &point_properties,
                                  const std::vector<UCDCell> &cells,
                                  const std::vector<UCDProperty<double>> &cell_properties,
                                  const std::string &filePath) const
{
    std::ofstream file;

    file.open(filePath.c_str());

    if (file.fail())
        throw runtime_error("File '" + filePath + "' cannot be opened");

    const char sep = ' ';
    file.precision(16);

    // export full info
    file << points.cols() << sep;
    file << cells.size() << sep;
    file << point_properties.size() << sep;
    file << cell_properties.size() << sep;
    file << 0 << std::endl; // model not supported

    // export points
    for (unsigned int p = 0; p < points.cols(); p++)
    {
        file << (p + 1) << sep;
        file << std::scientific << points(0, p) << sep;
        file << std::scientific << points(1, p) << sep;
        file << std::scientific << points(2, p) << std::endl;
    }

    // export cells
    for (unsigned int c = 0; c < cells.size(); c++)
    {
        file << (c + 1) << sep;
        file << std::scientific << cells[c].MaterialId << sep;
        file << std::scientific << cells[c].CellLabel(cells[c].Type);
        for (unsigned int p = 0; p < cells[c].PointIds.size(); p++)
        {
            file << sep << (cells[c].PointIds[p] + 1);
        }
        file << std::endl;
    }

    // export points properties
    if (point_properties.size() > 0)
    {
        file << point_properties.size();
        for (unsigned int pr = 0; pr < point_properties.size(); pr++)
        {
            file << sep << point_properties[pr].NumComponents;
        }
        file << std::endl;

        for (unsigned int pr = 0; pr < point_properties.size(); pr++)
        {
            file << point_properties[pr].Label << ',' << sep;
            file << point_properties[pr].UnitLabel << std::endl;
        }

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            file << (p + 1);
            for (unsigned int pr = 0; pr < point_properties.size(); pr++)
            {
                for (unsigned int cp = 0; cp < point_properties[pr].NumComponents; cp++)
                {
                    file << sep << point_properties[pr].Data[point_properties[pr].NumComponents * p + cp];
                }
            }
            file << std::endl;
        }
    }

    // export cells properties
    if (cell_properties.size() > 0)
    {
        file << cell_properties.size();
        for (unsigned int pr = 0; pr < cell_properties.size(); pr++)
        {
            file << sep << cell_properties[pr].NumComponents;
        }
        file << std::endl;

        for (unsigned int pr = 0; pr < cell_properties.size(); pr++)
        {
            file << cell_properties[pr].Label << ',' << sep;
            file << cell_properties[pr].UnitLabel << std::endl;
        }

        for (unsigned int c = 0; c < cells.size(); c++)
        {
            file << (c + 1);
            for (unsigned int pr = 0; pr < cell_properties.size(); pr++)
            {
                for (unsigned int cp = 0; cp < cell_properties[pr].NumComponents; cp++)
                {
                    file << sep << cell_properties[pr].Data[cell_properties[pr].NumComponents * c + cp];
                }
            }
            file << std::endl;
        }
    }

    file.close();
}
// ***************************************************************************
const string UCDCell::CellLabel(const UCDCell::Types type) const
{
    switch (type)
    {
    case UCDCell::Types::Line:
        return "line";
    case UCDCell::Types::Triangle:
        return "tri";
    case UCDCell::Types::Quadrilateral:
        return "quad";
    case UCDCell::Types::Hexahedron:
        return "hex";
    case UCDCell::Types::Prism:
        return "prism";
    case UCDCell::Types::Tetrahedron:
        return "tet";
    case UCDCell::Types::Pyramid:
        return "pyr";
    case UCDCell::Types::Point:
        return "pt";
    default:
        throw std::runtime_error("Type not supported");
    }
}

// ***************************************************************************
} // namespace Gedim
