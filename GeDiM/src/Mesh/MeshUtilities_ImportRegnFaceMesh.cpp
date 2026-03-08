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
    std::vector<std::vector<std::vector<unsigned int>>> cell3Ds_faces_vertices;

    {
      std::vector<std::string> cell0DsLines;
      Gedim::FileReader csvFileReader(node_file_path);

      if (!csvFileReader.Open())
        throw std::runtime_error("File node not found");

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0Ds = 0;
      {
        std::istringstream converter(cell0DsLines[1]);
        converter >> numCell0Ds;
      }

      if (numCell0Ds == 0)
        throw std::runtime_error("File node empty");

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

    {
      std::vector<std::string> cell3DsLines;
      Gedim::FileReader csvFileReader(ele_file_path);

      if (!csvFileReader.Open())
        throw std::runtime_error("File ele not found");

      csvFileReader.GetAllLines(cell3DsLines);
      csvFileReader.Close();

      unsigned int numCell3Ds = 0;
      {
        std::istringstream converter(cell3DsLines[1]);
        converter >> numCell3Ds;
      }

      if (numCell3Ds == 0)
        throw std::runtime_error("File ele empty");

      cell3Ds_faces_vertices.resize(numCell3Ds);

      unsigned int num_line = 2;
      unsigned int id, num_faces;
      for (unsigned int c = 0; c < numCell3Ds; c++)
      {
        std::istringstream converter(cell3DsLines.at(num_line++));

        converter >> id;
        converter >> num_faces;

        cell3Ds_faces_vertices.at(c).resize(num_faces);

        unsigned int num_face_vertices;
        for (unsigned int f = 0; f < num_faces; ++f)
        {
          std::istringstream converter_face(cell3DsLines.at(num_line++));
          converter_face >> id;
          converter_face >> num_face_vertices;

          cell3Ds_faces_vertices.at(c).at(f).resize(num_face_vertices);


          unsigned int ver_id;
          for (unsigned int f_v = 0; f_v < num_face_vertices; ++f_v)
          {
            converter_face >> ver_id;
            cell3Ds_faces_vertices.at(c).at(f).at(f_v) = ver_id;
          }
        }
      }
    }

    FillMesh3D(cell0Ds,
               cell3Ds_faces_vertices,
               mesh);
  }
  // ***************************************************************************
  void Gedim::MeshUtilities::FillMesh3D(const Eigen::MatrixXd &cell0Ds,
                                        const std::vector<std::vector<std::vector<unsigned int>>> &cell3Ds_faces_vertices,
                                        Gedim::IMeshDAO &mesh) const
  {    
    FillMesh3D(cell0Ds,
               {},
               {},
               {},
               {},
               {},
               {},
               mesh);
  }
  // ***************************************************************************
} // namespace Gedim
