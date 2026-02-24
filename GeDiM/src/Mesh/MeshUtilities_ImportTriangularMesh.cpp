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

#include "FileTextReader.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
// ***************************************************************************
void Gedim::MeshUtilities::ImportTriangularMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                                const std::string &cell0Ds_file_path,
                                                const std::string &cell2Ds_file_path,
                                                const std::string &marker_file_path,
                                                const char separator,
                                                IMeshDAO &mesh) const
{
    Eigen::MatrixXd cell0Ds;
    std::vector<Eigen::VectorXi> originalCell2Ds;
    std::vector<Eigen::VectorXi> cell2Ds;

    {
        std::vector<std::string> cell0DsLines;
        Gedim::FileReader csvFileReader(cell0Ds_file_path);

        if (!csvFileReader.Open())
            throw std::runtime_error("File cell0Ds not found");

        csvFileReader.GetAllLines(cell0DsLines);
        csvFileReader.Close();

        unsigned int numCell0Ds = cell0DsLines.size() - 1;
        if (numCell0Ds == 0)
            throw std::runtime_error("File cell0Ds empty");

        cell0Ds.setZero(3, numCell0Ds);

        char temp;
        unsigned int id;
        for (unsigned int v = 0; v < numCell0Ds; v++)
        {
            std::istringstream converter(cell0DsLines[v + 1]);

            converter >> id;
            if (separator != ' ')
                converter >> temp;
            converter >> cell0Ds(0, v);
            if (separator != ' ')
                converter >> temp;
            converter >> cell0Ds(1, v);
        }
    }

    {
        std::vector<std::string> cell2DsLines;
        Gedim::FileReader csvFileReader(cell2Ds_file_path);

        if (!csvFileReader.Open())
            throw std::runtime_error("File cell2Ds not found");

        csvFileReader.GetAllLines(cell2DsLines);
        csvFileReader.Close();

        unsigned int numCell2Ds = cell2DsLines.size() - 1;
        if (numCell2Ds == 0)
            throw std::runtime_error("File cell2Ds empty");

        cell2Ds.resize(numCell2Ds, Eigen::VectorXi::Zero(3));
        originalCell2Ds.resize(numCell2Ds, Eigen::VectorXi::Zero(3));

        char temp;
        unsigned int id;
        for (unsigned int t = 0; t < numCell2Ds; t++)
        {
            std::istringstream converter(cell2DsLines[t + 1]);

            converter >> id;
            if (separator != ' ')
                converter >> temp;
            converter >> originalCell2Ds[t](0);
            if (separator != ' ')
                converter >> temp;
            converter >> originalCell2Ds[t](1);
            if (separator != ' ')
                converter >> temp;
            converter >> originalCell2Ds[t](2);

            Eigen::Matrix3d points;
            points.col(0) << cell0Ds.col(originalCell2Ds.at(t)(0));
            points.col(1) << cell0Ds.col(originalCell2Ds.at(t)(1));
            points.col(2) << cell0Ds.col(originalCell2Ds.at(t)(2));

            std::vector<unsigned int> convexPoints = geometry_utilities.ConvexHull(points);
            cell2Ds.at(t)(0) = originalCell2Ds.at(t)(convexPoints.at(0));
            cell2Ds.at(t)(1) = originalCell2Ds.at(t)(convexPoints.at(1));
            cell2Ds.at(t)(2) = originalCell2Ds.at(t)(convexPoints.at(2));
        }
    }

    const auto cell1Ds = ComputeMesh2DCell1Ds(cell0Ds, cell2Ds);

    FillMesh2D(cell0Ds, cell1Ds.Cell1Ds, cell1Ds.Cell2Ds, mesh);

    {
        std::vector<std::string> cell2DsLines;
        Gedim::FileReader csvFileReader(marker_file_path);

        if (!csvFileReader.Open())
            throw std::runtime_error("File cell2Ds_marker not found");

        csvFileReader.GetAllLines(cell2DsLines);
        csvFileReader.Close();

        unsigned int numCell2Ds = cell2DsLines.size() - 1;
        if (numCell2Ds == 0)
            throw std::runtime_error("File cell2DsMarker empty");

        char temp;
        unsigned int cell2DIndex, vertexIndex, marker;
        for (unsigned int t = 0; t < numCell2Ds; t++)
        {
            std::istringstream converter(cell2DsLines[t + 1]);

            converter >> cell2DIndex;
            if (separator != ' ')
                converter >> temp;
            converter >> vertexIndex;
            if (separator != ' ')
                converter >> temp;
            converter >> marker;

            if (marker != 0)
            {
                const unsigned int correctedOriginIndex = originalCell2Ds.at(cell2DIndex)[(vertexIndex + 1) % 3];
                const unsigned int correctedEndIndex = originalCell2Ds.at(cell2DIndex)[(vertexIndex + 2) % 3];
                unsigned int cell1DIndex = mesh.Cell1DByExtremes(correctedOriginIndex, correctedEndIndex);
                if (cell1DIndex == mesh.Cell1DTotalNumber())
                {
                    cell1DIndex = mesh.Cell1DByExtremes(correctedEndIndex, correctedOriginIndex);
                }

                mesh.Cell0DSetMarker(correctedOriginIndex, marker);
                mesh.Cell0DSetMarker(correctedEndIndex, marker);
                mesh.Cell1DSetMarker(cell1DIndex, marker);
            }
        }
    }
}
// ***************************************************************************
} // namespace Gedim
