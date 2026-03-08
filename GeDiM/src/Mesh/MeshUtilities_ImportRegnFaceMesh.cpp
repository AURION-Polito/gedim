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
    const unsigned int num_cell3Ds = cell3Ds_faces_vertices.size();

    std::vector<std::pair<unsigned int, unsigned int>> cell1Ds_vertices;
    std::vector<std::vector<unsigned int>> cell2Ds_vertices;
    std::vector<std::vector<unsigned int>> cell2Ds_edges;
    std::vector<std::vector<unsigned int>> cell3Ds_vertices(num_cell3Ds);
    std::vector<std::vector<unsigned int>> cell3Ds_edges(num_cell3Ds);
    std::vector<std::vector<unsigned int>> cell3Ds_faces(num_cell3Ds);

    {
      std::map<std::array<unsigned int, 3>, unsigned int> cell2Ds_id;
      std::list<std::vector<unsigned int>> list_cell2Ds_vertices;

      auto make_cell2D_id = [](const unsigned int v_0, const unsigned int v_1, const unsigned int v_2)
      {
        std::array<unsigned int, 3> cell2D_id = { v_0, v_1, v_2};
        std::sort(cell2D_id.begin(),
                  cell2D_id.end());

        return cell2D_id;
      };

      for (unsigned int c = 0; c < num_cell3Ds; ++c)
      {
        const auto& cell3D_faces_vertices = cell3Ds_faces_vertices.at(c);
        unsigned int cell3D_num_face = cell3D_faces_vertices.size();

        cell3Ds_faces.at(c).resize(cell3D_num_face);

        for (unsigned int c_f = 0; c_f < cell3D_num_face; ++c_f)
        {
          const auto& cell3D_face = cell3D_faces_vertices.at(c_f);

          int face_id = -1;
          const unsigned int n_vertices = cell3D_face.size();
          for (unsigned int v = 0; v < n_vertices; ++v)
          {
            const auto cell2d_id = make_cell2D_id(cell3D_face.at(v),
                                                  cell3D_face.at((v + 1) % n_vertices),
                                                  cell3D_face.at((v + 2) % n_vertices));

            if (cell2Ds_id.contains(cell2d_id))
            {
              face_id = cell2Ds_id.at(cell2d_id);
              break;
            }
          }

          if (face_id == -1)
          {
            const auto cell2d_id = make_cell2D_id(cell3D_face.at(0),
                                                  cell3D_face.at(1),
                                                  cell3D_face.at(2));

            face_id = cell2Ds_id.size();
            cell2Ds_id.insert(std::make_pair(cell2d_id, face_id));

            list_cell2Ds_vertices.push_back(cell3D_face);
          }

          cell3Ds_faces.at(c).at(c_f) = face_id;
        }
      }

      cell2Ds_vertices = std::vector<std::vector<unsigned int>>(list_cell2Ds_vertices.begin(),
                                                                list_cell2Ds_vertices.end());
    }

    {
      std::map<std::pair<unsigned int, unsigned int>, unsigned int> cell1Ds_id;
      std::list<std::pair<unsigned int, unsigned int>> list_cell1Ds_vertices;

      cell2Ds_edges.resize(cell2Ds_vertices.size());

      auto make_cell1D_id = [](const unsigned int v_0, const unsigned int v_1)
      {
        return v_0 > v_1 ?
              std::make_pair(v_1, v_0) :
              std::make_pair(v_0, v_1);
      };

      for (unsigned int f = 0; f < cell2Ds_vertices.size(); ++f)
      {
        const auto cell2d_vertices = cell2Ds_vertices.at(f);
        const unsigned int n_vertices = cell2d_vertices.size();

        cell2Ds_edges.at(f).resize(n_vertices);

        for (unsigned int f_v = 0; f_v < n_vertices; ++f_v)
        {
          int edge_id = -1;

          const auto cell1d_id = make_cell1D_id(cell2d_vertices.at(f_v),
                                                cell2d_vertices.at((f_v + 1) % n_vertices));

          if (!cell1Ds_id.contains(cell1d_id))
          {
            edge_id = cell1Ds_id.size();
            cell1Ds_id.insert(std::make_pair(cell1d_id, edge_id));

            list_cell1Ds_vertices.push_back(cell1d_id);
          }
          else
            edge_id = cell1Ds_id.at(cell1d_id);

          cell2Ds_edges.at(f).at(f_v) = edge_id;
        }
      }

      for (unsigned int c = 0; c < cell3Ds_faces.size(); ++c)
      {
        std::set<unsigned int> cell3D_vertices;
        std::set<unsigned int> cell3D_edges;

        for (unsigned int c_f = 0; c_f < cell3Ds_faces.at(c).size(); ++c_f)
        {
          const auto& face_vertices = cell2Ds_vertices.at(c_f);

          for (unsigned int f_v = 0; f_v < face_vertices.size(); ++f_v)
            cell3D_vertices.insert(face_vertices.at(f_v));

          const auto& face_edges = cell2Ds_edges.at(c_f);

          for (unsigned int f_e = 0; f_e < face_edges.size(); ++f_e)
            cell3D_edges.insert(face_edges.at(f_e));
        }

        cell3Ds_vertices.at(c) = std::vector<unsigned int>(cell3D_vertices.begin(),
                                                           cell3D_vertices.end());

        cell3Ds_edges.at(c) = std::vector<unsigned int>(cell3D_edges.begin(),
                                                        cell3D_edges.end());
      }

      cell1Ds_vertices = std::vector<std::pair<unsigned int, unsigned int>>(list_cell1Ds_vertices.begin(),
                                                                            list_cell1Ds_vertices.end());
    }

    FillMesh3D(cell0Ds,
               cell1Ds_vertices,
               cell2Ds_vertices,
               cell2Ds_edges,
               cell3Ds_vertices,
               cell3Ds_edges,
               cell3Ds_faces,
               mesh);
  }
  // ***************************************************************************
} // namespace Gedim
