#ifndef __UCDUtilities_H
#define __UCDUtilities_H

#include "Eigen/Eigen"
#include <list>
#include <string>

namespace Gedim
{
  template <typename T>
  struct UCDProperty final
  {
      enum Formats
      {
        Points = 0,
        Cells = 1
      };

      std::string Label;
      std::string UnitLabel;
      Formats Format;

      unsigned int Size;
      unsigned int NumComponents;
      const T* Data;
  };


  class UCDCell final
  {
    public:
      enum struct Types
      {
        Point = 0,
        Line = 1,
        Triangle = 2,
        Quadrilateral = 3,
        Hexahedron = 4,
        Prism = 5,
        Tetrahedron = 6,
        Pyramid = 7
      };

      const Types Type;
      const std::vector<unsigned int> PointIds;
      const unsigned int MaterialId;

      UCDCell(const Types type,
              const std::vector<unsigned int> pointIds,
              const unsigned int materialId) :
        Type(type),
        PointIds(pointIds),
        MaterialId(materialId)
      { }

      const std::string CellLabel(const UCDCell::Types type) const;
  };

  class UCDUtilities final
  {
    public:
      enum ExportFormats
      {
        Ascii = 0
      };

    private:
      std::vector<UCDCell> CreatePointCells(const Eigen::MatrixXd& points,
                                            const Eigen::VectorXi& points_material) const;

      void ExportUCDAscii(const Eigen::MatrixXd& points,
                          const std::vector<UCDProperty<double>>& point_properties,
                          const std::vector<UCDCell>& cells,
                          const std::vector<UCDProperty<double>>& cell_properties,
                          const std::string& filePath) const;

    public:
      UCDUtilities() { }
      virtual ~UCDUtilities() { }

      void ExportPoints(const std::string& filePath,
                        const Eigen::MatrixXd& points,
                        const std::vector<UCDProperty<double>>& points_properties = {},
                        const Eigen::VectorXi& points_material = {}) const;
  };
}

#endif // __UCDUtilities_H
