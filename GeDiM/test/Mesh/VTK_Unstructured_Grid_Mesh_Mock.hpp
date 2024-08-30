#ifndef __VTK_Unstructured_Grid_Mesh_Mock_H
#define __VTK_Unstructured_Grid_Mesh_Mock_H

#include <fstream>
#include <string>
#include <vector>

namespace GedimUnitTesting
{
  class VTK_Unstructured_Grid_Mesh_Mock final
  {
    public:
      static void ExportFile(const std::string& file_path)
      {
        const std::vector<std::string> file_lines {
          "# vtk DataFile Version 3.1",
          "MCVE VTK file",
          "ASCII",
          "DATASET UNSTRUCTURED_GRID",
          "POINTS      16 double",
          " 0.   0.   0.",
          " 0.   0.   3.",
          " 0.   2.   0.",
          " 0.   2.   3.",
          " 4.   0.   0.",
          " 4.   0.   3.",
          " 4.   2.   0.",
          " 4.   2.   3.",
          " 5.   0.   0.",
          " 5.   0.   3.",
          " 5.   2.   0.",
          " 5.   2.   3.",
          "13.   0.   0.",
          "13.   0.   3.",
          "13.   2.   0.",
          "13.   2.   3.",
          "",
          "CELLS        3     27",
          " 8    0   1   3   2   4   5   7   6",
          " 8    4   5   7   6   8   9  11  10",
          " 8    8   9  11  10  12  13  15  14",
          "",
          "CELL_TYPES        3",
          "          12          12          12",
          "",
          "CELL_DATA        3",
          "SCALARS elem_val float",
          "LOOKUP_TABLE default",
          " 1",
          " 2",
          " 3"
        };

        std::ofstream exportFile(file_path);
        if (!exportFile.is_open())
          throw std::runtime_error("Unable to export file " + file_path);

        for (const std::string& line : file_lines)
          exportFile << line<< std::endl;
        exportFile.close();
      }
  };
}

#endif // __VTK_Unstructured_Grid_Mesh_Mock_H
