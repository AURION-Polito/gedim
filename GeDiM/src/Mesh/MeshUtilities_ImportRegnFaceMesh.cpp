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
                                                IMeshDAO &mesh) const
  {
    Eigen::MatrixXd cell0Ds;
    std::vector<std::vector<unsigned int>> cell2Ds_vertices;
    std::vector<std::vector<unsigned int>> cell3Ds_faces;

    {
      std::vector<std::string> cell0DsLines;
      Gedim::FileReader csvFileReader(node_file_path);

      if (!csvFileReader.Open())
        throw std::runtime_error("File cell0Ds not found");

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0Ds = 0;
      {
        std::istringstream converter(cell0DsLines[1]);
        converter >> numCell0Ds;
      }

      if (numCell0Ds == 0)
        throw std::runtime_error("File cell0Ds empty");

      cell0Ds.setZero(3, numCell0Ds);

      unsigned int id;
      for (unsigned int v = 0; v < numCell0Ds; v++)
      {
        std::istringstream converter(cell0DsLines[v + 2]);

        converter >> id;
        converter >> cell0Ds(0, v);
        converter >> cell0Ds(1, v);
        converter >> cell0Ds(2, v);
      }
    }

    FillMesh3D(cell0Ds,
               cell2Ds_vertices,
               cell3Ds_faces,
               mesh);
  }
  // ***************************************************************************
  void Gedim::MeshUtilities::FillMesh3D(const Eigen::MatrixXd &cell0Ds,
                                        const std::vector<std::vector<unsigned int>> &cell2Ds_vertices,
                                        const std::vector<std::vector<unsigned int>> &cell3Ds_faces,
                                        Gedim::IMeshDAO &mesh) const
  {

  }
  // ***************************************************************************
} // namespace Gedim
