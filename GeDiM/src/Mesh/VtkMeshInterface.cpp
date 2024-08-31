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
  VtkMeshInterface::VtkMesh3D VtkMeshInterface::ComputeVtkMesh3D(VtkMesh& originalMesh) const
  {
    VtkMesh3D result;

    result.Cell0Ds = originalMesh.Cell0Ds;

    // Create Cell1Ds
    Eigen::SparseMatrix<unsigned int> edges;
    edges.resize(originalMesh.NumCell0Ds,
                 originalMesh.NumCell0Ds);

    std::list<Eigen::Triplet<unsigned int>> triplets;
    for (unsigned int c = 0; c < originalMesh.NumCell3Ds; c++)
    {
      const auto& cell3D = originalMesh.Cell3Ds.at(c);

      switch (cell3D.Type)
      {
        case VtkMesh::Cell3D::Types::Tetrahedron:
          break;
        default:
          throw std::runtime_error("VTK mesh Cell type not supported");
      }

      const auto& cell3D_vertices = cell3D.Cell0D_Id;
      Gedim::Output::Assert(cell3D_vertices.size() == 4);

      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(0),
                                                          cell3D_vertices.at(1)),
                                                      MAX(cell3D_vertices.at(0),
                                                          cell3D_vertices.at(1)),
                                                      1));
      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(1),
                                                          cell3D_vertices.at(2)),
                                                      MAX(cell3D_vertices.at(1),
                                                          cell3D_vertices.at(2)),
                                                      1));
      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(2),
                                                          cell3D_vertices.at(0)),
                                                      MAX(cell3D_vertices.at(2),
                                                          cell3D_vertices.at(0)),
                                                      1));
      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(0),
                                                          cell3D_vertices.at(3)),
                                                      MAX(cell3D_vertices.at(0),
                                                          cell3D_vertices.at(3)),
                                                      1));
      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(1),
                                                          cell3D_vertices.at(3)),
                                                      MAX(cell3D_vertices.at(1),
                                                          cell3D_vertices.at(3)),
                                                      1));
      triplets.push_back(Eigen::Triplet<unsigned int>(MIN(cell3D_vertices.at(2),
                                                          cell3D_vertices.at(3)),
                                                      MAX(cell3D_vertices.at(2),
                                                          cell3D_vertices.at(3)),
                                                      1));
    }

    edges.setFromTriplets(triplets.begin(), triplets.end());
    edges.makeCompressed();

    unsigned int num_cell1Ds = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
        it.valueRef() = 1 + num_cell1Ds++;
    }

    result.Cell1Ds.resize(2, num_cell1Ds);

    num_cell1Ds = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
        result.Cell1Ds.col(num_cell1Ds++)<< it.row(), it.col();
    }

    return result;
  }
  // ***************************************************************************
  VtkMeshInterface::VtkMesh3D VtkMeshInterface::ImportMesh3DFromFile(const std::string& vtkFilePath) const
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
        return {};

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
          case 12:
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

      return ComputeVtkMesh3D(vtk_mesh);
    }

    return {};
  }
  // ***************************************************************************
}
