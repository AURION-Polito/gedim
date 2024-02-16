#include "UCDUtilities.hpp"
#include <fstream>

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  void UCDUtilities::ExportPoints(const Eigen::MatrixXd& points,
                                  const std::vector<UCDProperty<double>>& points_properties,
                                  const std::string& filePath,
                                  const ExportFormats& format) const
  {
    switch (format)
    {
      case ExportFormats::Ascii:
        ExportUCDAscii(points,
                       points_properties,
                       filePath);
        break;
      default:
        throw std::runtime_error("Unknown format");
    }
  }
  // ***************************************************************************
  void UCDUtilities::ExportUCDAscii(const Eigen::MatrixXd& points,
                                    const std::vector<UCDProperty<double> >& point_properties,
                                    const std::string& filePath) const
  {
    std::ofstream file;

    file.open(filePath.c_str());

    if (file.fail())
      throw runtime_error("File '" + filePath + "' cannot be opened");

    const char sep = ' ';
    file.precision(16);

    // export full info
    file<< points.cols()<< sep;
    file<< 0<< sep;
    file<< point_properties.size()<< sep;
    file<< 0<< std::endl;

    // export points
    for (unsigned int p = 0; p < points.cols(); p++)
    {
      file<< p<< sep;
      file<< std::scientific<< points(0, p)<< sep;
      file<< std::scientific<< points(0, p)<< sep;
      file<< std::scientific<< points(1, p)<< std::endl;
    }

    // export cells

    // export points properties
    file<< point_properties.size();
    for (unsigned int pr = 0; pr < point_properties.size(); pr++)
    {
      file<< sep<< point_properties[pr].NumComponents;
    }
    file<< std::endl;

    for (unsigned int pr = 0; pr < point_properties.size(); pr++)
    {
      file<< point_properties[pr].Label<< sep;
      file<< point_properties[pr].UnitLabel<< std::endl;
    }

    for (unsigned int p = 0; p < points.cols(); p++)
    {
      file<< p;
      for (unsigned int pr = 0; pr < point_properties.size(); pr++)
      {
        for (unsigned int cp = 0; cp < point_properties[pr].NumComponents; cp++)
        {
          file<< sep<< point_properties[pr].Data[point_properties[pr].NumComponents * p + cp];
        }
      }
    }

    // export cells properties

    file.close();
  }
  // ***************************************************************************
}
