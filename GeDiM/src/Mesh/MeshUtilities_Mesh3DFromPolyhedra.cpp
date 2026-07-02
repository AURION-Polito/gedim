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

#include "MeshUtilities.hpp"

namespace Gedim
{
  // ***************************************************************************
  void Gedim::MeshUtilities::Mesh3DFromPolyhedra(const Gedim::GeometryUtilities &geometry_utilities,
                                                 const std::vector<Gedim::GeometryUtilities::Polyhedron>& polyhedra,
                                                 Gedim::IMeshDAO &mesh) const
  {
    if (polyhedra.empty())
      return;

    struct cell0D_filter
    {
        Eigen::Vector3d coordinate;
        unsigned int idx;
    };

    struct cell2D_filter
    {
        std::vector<unsigned int> vertices;
        std::vector<unsigned int> edges;
        unsigned int idx;
    };

    std::vector<std::vector<unsigned int>> polyhedra_cell0Ds(polyhedra.size());
    std::vector<std::vector<unsigned int>> polyhedra_cell1Ds(polyhedra.size());
    std::vector<std::vector<unsigned int>> polyhedra_cell2Ds(polyhedra.size());
    std::list<cell0D_filter> cell0Ds_filter;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> cell1Ds_filter;
    std::map<std::array<unsigned int, 3>, cell2D_filter> cell2Ds_filter;

    {
      unsigned int p = 0;
      for (const auto& polyhedron : polyhedra)
      {
        auto& vertices = polyhedra_cell0Ds.at(p);
        vertices.resize(polyhedron.Vertices.cols());

        for (unsigned int p_v = 0; p_v < polyhedron.Vertices.cols(); ++p_v)
        {
          const Eigen::Vector3d point = polyhedron.Vertices.col(p_v);

          auto found = cell0Ds_filter.end();
          for (auto cell0D_it = cell0Ds_filter.begin();
               cell0D_it != cell0Ds_filter.end();
               cell0D_it++)
          {
            if (geometry_utilities.PointsAreCoincident(point,
                                                       cell0D_it->coordinate))
            {
              found = cell0D_it;
              break;
            }
          }

          if (found == cell0Ds_filter.end())
          {
            vertices.at(p_v) = cell0Ds_filter.size();

            cell0Ds_filter.push_back({
                                       point,
                                       static_cast<unsigned int>(cell0Ds_filter.size())
                                     });
          }
          else
          {
            vertices.at(p_v) = found->idx;
          }
        }

        auto& edges = polyhedra_cell1Ds.at(p);
        edges.resize(polyhedron.Edges.cols());

        for (unsigned int p_e = 0; p_e < polyhedron.Edges.cols(); ++p_e)
        {
          const unsigned int edge_org = vertices.at(polyhedron.Edges(0, p_e));
          const unsigned int edge_end = vertices.at(polyhedron.Edges(1, p_e));

          std::pair<unsigned int, unsigned int> edge_extr =
          { std::min(edge_org, edge_end), std::max(edge_org, edge_end) };

          auto found = cell1Ds_filter.find(edge_extr);

          if (found == cell1Ds_filter.end())
          {
            edges.at(p_e) = cell1Ds_filter.size();

            cell1Ds_filter.insert(std::make_pair(edge_extr,
                                                 static_cast<unsigned int>(cell1Ds_filter.size())));
          }
          else
          {
            edges.at(p_e) = found->second;
          }
        }

        auto& faces = polyhedra_cell2Ds.at(p);
        faces.resize(polyhedron.Faces.size());

        for (unsigned int p_f = 0; p_f < polyhedron.Faces.size(); ++p_f)
        {
          std::vector<unsigned int> face_vertices(polyhedron.Faces.at(p_f).cols());
          std::vector<unsigned int> face_edges(polyhedron.Faces.at(p_f).cols());

          unsigned int prev_p_f_v_max = 0;
          unsigned int p_f_v_max = 0;
          unsigned int next_p_f_v_max = 0;
          unsigned int vertex_max = 0;

          for (unsigned int p_f_v = 0; p_f_v < face_vertices.size(); ++p_f_v)
          {
            face_vertices.at(p_f_v) = vertices.at(polyhedron.Faces.at(p_f)(0, p_f_v));
            face_edges.at(p_f_v) = edges.at(polyhedron.Faces.at(p_f)(1, p_f_v));

            if (vertex_max < face_vertices.at(p_f_v))
            {
              vertex_max = face_vertices.at(p_f_v);
              p_f_v_max = p_f_v;
              prev_p_f_v_max = (p_f_v == 0) ? face_vertices.size() - 1 :
                                              p_f_v - 1;
              next_p_f_v_max = (p_f_v + 1) % face_vertices.size();
            }
          }

          std::array<unsigned int, 3> face_key = {
            face_vertices.at(prev_p_f_v_max),
            face_vertices.at(p_f_v_max),
            face_vertices.at(next_p_f_v_max)
          };
          std::sort(face_key.begin(),
                    face_key.end());

          auto found = cell2Ds_filter.find(face_key);

          if (found == cell2Ds_filter.end())
          {
            faces.at(p_f) = cell2Ds_filter.size();

            cell2Ds_filter.insert(std::make_pair(face_key,
                                                 cell2D_filter({
                                                                 face_vertices,
                                                                 face_edges,
                                                                 static_cast<unsigned int>(cell2Ds_filter.size())
                                                               })));
          }
          else
          {
            faces.at(p_f) = found->second.idx;
          }
        }

        p++;
      }
    }

    Eigen::MatrixXd cell0Ds(3, cell0Ds_filter.size());
    for (const auto& cell0D_filter : cell0Ds_filter)
      cell0Ds.col(cell0D_filter.idx) = cell0D_filter.coordinate;

    std::vector<std::pair<unsigned int, unsigned int>> cell1Ds(cell1Ds_filter.size());
    for (const auto& cell1D_filter : cell1Ds_filter)
      cell1Ds.at(cell1D_filter.second) = cell1D_filter.first;

    std::vector<std::vector<unsigned int>> cell2Ds_vertices(cell2Ds_filter.size());
    std::vector<std::vector<unsigned int>> cell2Ds_edges(cell2Ds_filter.size());
    for (const auto& cell2D_filter : cell2Ds_filter)
    {
      cell2Ds_vertices.at(cell2D_filter.second.idx) = cell2D_filter.second.vertices;
      cell2Ds_edges.at(cell2D_filter.second.idx) = cell2D_filter.second.edges;
    }

    FillMesh3D(cell0Ds,
               cell1Ds,
               cell2Ds_vertices,
               cell2Ds_edges,
               polyhedra_cell0Ds,
               polyhedra_cell1Ds,
               polyhedra_cell2Ds,
               mesh);

    std::cout<< "Num vertices "<< cell0Ds_filter.size()<< std::endl;
  }
  // ***************************************************************************
} // namespace Gedim
