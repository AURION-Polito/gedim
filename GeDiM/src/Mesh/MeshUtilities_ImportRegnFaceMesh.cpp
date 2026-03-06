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
void Gedim::MeshUtilities::ImportRegnFaceMesh(const Gedim::GeometryUtilities &geometry_utilities,
                                                const std::string &node_file_path,
                                                const std::string &ele_file_path,
                                                const char separator,
                                                IMeshDAO &mesh) const
{
    Eigen::MatrixXd cell0Ds;
    std::vector<Eigen::VectorXi> originalCell2Ds;
    std::vector<Eigen::VectorXi> cell2Ds;

    {
        std::vector<std::string> cell0DsLines;
        Gedim::FileReader csvFileReader(node_file_path);

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
}
// ***************************************************************************
} // namespace Gedim
