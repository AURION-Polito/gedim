#include "VtkMeshInterface.hpp"

#include <vtkUnstructuredGrid.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  VtkMeshInterface::VtkMeshInterface()
  {
  }
  VtkMeshInterface::~VtkMeshInterface()
  {
  }
  // ***************************************************************************
  void VtkMeshInterface::VtkMeshToMeshDAO(const VtkMesh& originalMesh,
                                          IMeshDAO& convertedMesh) const
  {
    convertedMesh.InitializeDimension(3);

    // Create Cell0Ds
    convertedMesh.Cell0DsInitialize(originalMesh.NumCell0Ds);
    convertedMesh.Cell0DsInsertCoordinates(originalMesh.Cell0Ds);
    for (unsigned int v = 0; v < originalMesh.NumCell0Ds; v++)
      convertedMesh.Cell0DSetState(v, true);
  }
  // ***************************************************************************
  void VtkMeshInterface::ImportMeshFromFile(const std::string& vtkFilePath,
                                            IMeshDAO& mesh) const
  {
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
        vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(vtkFilePath.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if (reader->IsFileUnstructuredGrid())
    {
      VtkMesh vtk_mesh;
      const vtkSmartPointer<vtkUnstructuredGrid> output = reader->GetUnstructuredGridOutput();
      const vtkIdType num_points = output->GetNumberOfPoints();

      if (num_points == 0)
        return;

      vtk_mesh.NumCell0Ds = num_points;
      vtk_mesh.Cell0Ds.resize(3, num_points);

      const vtkSmartPointer<vtkPoints> points_data = output->GetPoints();
      for (vtkIdType p = 0; p < num_points; p++)
      {
        std::array<double, 3> point;
        points_data->GetPoint(p, point.data());

        vtk_mesh.Cell0Ds.col(p)<< point.at(0), point.at(1), point.at(2);
      }

      const unsigned int num_cells = output->GetNumberOfCells();

      vtk_mesh.NumCell3Ds = num_cells;
      vtk_mesh.Cell3Ds.resize(num_cells);

      for (vtkIdType c = 0; c < num_cells; c++)
      {
        const vtkSmartPointer<vtkCell> cell = output->GetCell(c);
        const vtkIdType cell_type = cell->GetCellType();
        const vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();

        auto& cell3D = vtk_mesh.Cell3Ds.at(c);
        switch (cell_type)
        {
          case 10:
            cell3D.Type = VtkMesh::Cell3D::Types::Tetrahedron;
            break;
          case 11:
            cell3D.Type = VtkMesh::Cell3D::Types::Hexahedron;
            break;
          default:
            throw std::runtime_error("VTK mesh Cell type not supported");
        }

        const vtkIdType num_cell_points_id = pointIds->GetNumberOfIds();
        cell3D.Cell0D_Id.resize(num_cell_points_id);
        for (vtkIdType p = 0; p < num_cell_points_id; ++p)
        {
          const vtkIdType point_id = pointIds->GetId(p);
          cell3D.Cell0D_Id[p] = point_id;
        }
      }

      VtkMeshToMeshDAO(vtk_mesh,
                       mesh);
    }
  }
  // ***************************************************************************
}
