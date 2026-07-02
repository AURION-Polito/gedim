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

    std::vector<std::vector<unsigned int>> polyhedra_cell0Ds(polyhedra.size());
    std::list<cell0D_filter> cell0Ds_filter;

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

        p++;
      }
    }

    std::cout<< "Num vertices "<< cell0Ds_filter.size()<< std::endl;
  }
  // ***************************************************************************
} // namespace Gedim
