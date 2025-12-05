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

#ifndef __MEDIT_Utilities_H
#define __MEDIT_Utilities_H

#include "Eigen/Eigen"
#include <list>
#include <string>

namespace Gedim
{
  class MEDIT_Cell final
  {
    public:
      enum struct Types
      {
        Unknown = -1,
        Vertex = 0,
        Edge = 1,
        Triangle = 2,
        Quadrilateral = 3,
        Tetrahedron = 4,
        Hexahedron = 5
      };

      const Types Type;
      const std::vector<unsigned int> Points_id;
      const unsigned int Reference_id;

      MEDIT_Cell(const Types type,
                 const std::vector<unsigned int>& points_id,
                 const unsigned int reference_id)
        : Type(type), Points_id(points_id), Reference_id(reference_id)
      {
      }

      static std::string CellLabel(const MEDIT_Cell::Types type);
  };

  ///
  /// \brief The MEDIT_Utilities class for interface with MEDIT
  /// @see https://github.com/ISCDtoolbox/Medit
  /// @see https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf
  ///
  class MEDIT_Utilities final
  {
    public:
      enum ExportFormats
      {
        Ascii = 0
      };

    private:
      std::vector<MEDIT_Cell> CreateVertexCells(const Eigen::MatrixXd &points,
                                                const std::vector<unsigned int> &reference_ids) const;
      std::vector<MEDIT_Cell> CreateEdgeCells(const Eigen::MatrixXi &edges,
                                              const std::vector<unsigned int> &reference_ids) const;
      std::vector<MEDIT_Cell> CreatePolygonCells(const std::vector<std::vector<unsigned int>> &polygons_vertices,
                                                 const std::vector<unsigned int> &reference_ids) const;
      std::vector<MEDIT_Cell> CreatePolyhedraCells(const std::vector<std::vector<unsigned int>> &polyhedra_vertices,
                                                   const std::vector<unsigned int> &reference_ids) const;

      void ExportMEDIT_Ascii(const Eigen::MatrixXd &points,
                             const std::vector<unsigned int>& points_reference_ids,
                             const std::vector<MEDIT_Cell> &faces,
                             const std::vector<MEDIT_Cell> &cells,
                             const std::string &filePath) const;

    public:
      MEDIT_Utilities()
      {
      }
      virtual ~MEDIT_Utilities()
      {
      }

      void ExportPoints(const std::string &filePath,
                        const Eigen::MatrixXd &points,
                        const std::vector<unsigned int> &points_reference_ids) const;

      void ExportSegments(const std::string &filePath,
                          const Eigen::MatrixXd &points,
                          const Eigen::MatrixXi &cells,
                          const std::vector<unsigned int> &points_reference_ids,
                          const std::vector<unsigned int> &cells_reference_ids) const;

      void ExportPolygons(const std::string &filePath,
                          const Eigen::MatrixXd &points,
                          const std::vector<std::vector<unsigned int>> &cells_vertices,
                          const std::vector<unsigned int> &points_reference_ids,
                          const std::vector<unsigned int> &cells_reference_ids) const;

      void ExportPolyhedra(const std::string &filePath,
                           const Eigen::MatrixXd &points,
                           const std::vector<std::vector<unsigned int>> &faces_vertices,
                           const std::vector<std::vector<unsigned int>> &cells_vertices,
                           const std::vector<unsigned int> &points_reference_ids,
                           const std::vector<unsigned int> &faces_reference_ids,
                           const std::vector<unsigned int> &cells_reference_ids) const;
  };
} // namespace Gedim

#endif // __MEDIT_Utilities_H
