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
    static void ExportFile(const std::string &file_path)
    {
        const std::vector<std::string> file_lines{"# vtk DataFile Version 3.0",
                                                  "VTK file with 2 tetrahedrons",
                                                  "ASCII",
                                                  "DATASET UNSTRUCTURED_GRID",
                                                  "",
                                                  "POINTS 5 float",
                                                  "0 0 0",
                                                  "1 0 0",
                                                  "0 1 0",
                                                  "0 0 1",
                                                  "1 1 1",
                                                  "",
                                                  "CELLS 2 10",
                                                  "4 0 1 2 3",
                                                  "4 4 1 2 3",
                                                  "",
                                                  "CELL_TYPES 2",
                                                  "10",
                                                  "10",
                                                  "",
                                                  "POINT_DATA 5",
                                                  "SCALARS PointScalars float 1",
                                                  "LOOKUP_TABLE default",
                                                  "0.0",
                                                  "1.0",
                                                  "2.0",
                                                  "3.0",
                                                  "4.0"};

        std::ofstream exportFile(file_path);
        if (!exportFile.is_open())
            throw std::runtime_error("Unable to export file " + file_path);

        for (const std::string &line : file_lines)
            exportFile << line << std::endl;
        exportFile.close();
    }
};
} // namespace GedimUnitTesting

#endif // __VTK_Unstructured_Grid_Mesh_Mock_H
