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

  template<typename pointCollection>
  struct UCDPoints final
  {
      const pointCollection& Points;

      UCDPoints(const pointCollection& points) :
        Points(points)
      { }
  };

  class UCDUtilities final
  {
    public:
      enum ExportFormats
      {
        Ascii = 0
      };

    private:
      void ExportUCDAscii(const Eigen::MatrixXd& points,
                          const std::vector<UCDProperty<double>>& point_properties,
                          const std::string& filePath) const;

    public:
      UCDUtilities() { }
      virtual ~UCDUtilities() { }

      void ExportPoints(const Eigen::MatrixXd& points,
                        const std::vector<UCDProperty<double>>& points_properties,
                        const std::string& filePath,
                        const ExportFormats& format = ExportFormats::Ascii) const;
  };
}

#endif // __UCDUtilities_H
