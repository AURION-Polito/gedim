#include "VtkMeshInterface.hpp"

#include "Macro.hpp"

#include "CommonUtilities.hpp"

#if ENABLE_VTK == 1
#include <vtkUnstructuredGrid.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#endif

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
  unsigned int VtkMeshInterface::InsertNewEdge(const unsigned int origin,
                                               const unsigned int end,
                                               std::map<std::pair<unsigned int, unsigned int>, unsigned int>& edges) const
  {
    const auto edge = std::make_pair(MIN(origin,
                                         end),
                                     MAX(origin,
                                         end));

    const auto edge_id = edges.find(edge);
    if (edge_id != edges.end())
      return edge_id->second;

    const unsigned int new_edge_id = edges.size();
    edges.insert(std::make_pair(edge,
                                new_edge_id));

    return new_edge_id;
  }
  // ***************************************************************************
  unsigned int VtkMeshInterface::InsertNewFace(const std::vector<unsigned int>& face_vertices,
                                               const std::vector<unsigned int>& face_edges,
                                               std::map<std::pair<unsigned int, unsigned int>, Mesh_Face>& faces) const
  {
    for (unsigned int v = 0; v < 3; v++)
    {
      const auto face = std::make_pair(face_edges.at(v),
                                       face_vertices.at((v + 2) % 3));

      const auto face_id = faces.find(face);
      if (face_id != faces.end())
      {
        face_id->second.Marker = 0;
        return face_id->second.Index;
      }
    }

    Mesh_Face new_face =
    {
      Eigen::MatrixXi(2, 3),
      static_cast<unsigned int>(faces.size()),
      1
    };
    new_face.Extremes.row(0)<< face_vertices.at(0), face_vertices.at(1), face_vertices.at(2);
    new_face.Extremes.row(1)<< face_edges.at(0), face_edges.at(1), face_edges.at(2);

    faces.insert(std::make_pair(std::make_pair(face_edges.at(0),
                                               face_vertices.at(2)),
                                new_face));

    return new_face.Index;
  }
  // ***************************************************************************
  VtkMeshInterface::VtkMesh3D VtkMeshInterface::ComputeVtkMesh3D(VtkMesh& originalMesh) const
  {
    VtkMesh3D result;

    result.Cell0Ds = originalMesh.Cell0Ds;
    result.Markers[0].resize(originalMesh.NumCell0Ds, 0);

    result.Cell3Ds.resize(originalMesh.NumCell3Ds);
    result.Markers[3].resize(originalMesh.NumCell3Ds, 0);

    std::map<std::pair<unsigned int, unsigned int>, unsigned int> cell1Ds;
    std::map<std::pair<unsigned int, unsigned int>, Mesh_Face> cell2Ds;

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

      Gedim::Output::Assert(cell3D.Cell0D_Id.size() == 4);

      auto& new_cell3D = result.Cell3Ds.at(c);

      new_cell3D.VerticesIndex = cell3D.Cell0D_Id;
      new_cell3D.EdgesIndex.resize(6);
      new_cell3D.FacesIndex.resize(4);

      const auto& cell3D_vertices = new_cell3D.VerticesIndex;
      auto& cell3D_edges = new_cell3D.EdgesIndex;
      auto& cell3D_faces = new_cell3D.FacesIndex;

      cell3D_edges[0] = InsertNewEdge(cell3D_vertices.at(0),
                                      cell3D_vertices.at(1),
                                      cell1Ds);
      cell3D_edges[1] = InsertNewEdge(cell3D_vertices.at(1),
                                      cell3D_vertices.at(2),
                                      cell1Ds);
      cell3D_edges[2] = InsertNewEdge(cell3D_vertices.at(2),
                                      cell3D_vertices.at(0),
                                      cell1Ds);
      cell3D_edges[3] = InsertNewEdge(cell3D_vertices.at(0),
                                      cell3D_vertices.at(3),
                                      cell1Ds);
      cell3D_edges[4] = InsertNewEdge(cell3D_vertices.at(1),
                                      cell3D_vertices.at(3),
                                      cell1Ds);
      cell3D_edges[5] = InsertNewEdge(cell3D_vertices.at(2),
                                      cell3D_vertices.at(3),
                                      cell1Ds);

      cell3D_faces[0] = InsertNewFace({
                                        cell3D_vertices.at(0),
                                        cell3D_vertices.at(1),
                                        cell3D_vertices.at(2)
                                      },
                                      {
                                        cell3D_edges.at(0),
                                        cell3D_edges.at(1),
                                        cell3D_edges.at(2),
                                      },
                                      cell2Ds);

      cell3D_faces[1] = InsertNewFace({
                                        cell3D_vertices.at(0),
                                        cell3D_vertices.at(1),
                                        cell3D_vertices.at(3)
                                      },
                                      {
                                        cell3D_edges.at(0),
                                        cell3D_edges.at(4),
                                        cell3D_edges.at(3),
                                      },
                                      cell2Ds);

      cell3D_faces[2] = InsertNewFace({
                                        cell3D_vertices.at(1),
                                        cell3D_vertices.at(2),
                                        cell3D_vertices.at(3)
                                      },
                                      {
                                        cell3D_edges.at(1),
                                        cell3D_edges.at(5),
                                        cell3D_edges.at(4),
                                      },
                                      cell2Ds);

      cell3D_faces[3] = InsertNewFace({
                                        cell3D_vertices.at(2),
                                        cell3D_vertices.at(3),
                                        cell3D_vertices.at(0)
                                      },
                                      {
                                        cell3D_edges.at(5),
                                        cell3D_edges.at(3),
                                        cell3D_edges.at(2),
                                      },
                                      cell2Ds);
    }

    result.Cell1Ds.resize(2, cell1Ds.size());
    result.Markers[1].resize(cell1Ds.size(), 0);
    for (const auto& cell1D : cell1Ds)
    {
      result.Cell1Ds.col(cell1D.second)<< cell1D.first.first, cell1D.first.second;
    }

    result.Cell2Ds.resize(cell2Ds.size());
    result.Markers[2].resize(cell2Ds.size(), 0);
    for (const auto& cell2D : cell2Ds)
    {
      const unsigned int cell2D_index = cell2D.second.Index;
      auto& cell2D_data = result.Cell2Ds.at(cell2D_index);
      cell2D_data = cell2D.second.Extremes;
      result.Markers[2][cell2D_index] = cell2D.second.Marker;

      if (cell2D.second.Marker > 0)
      {
        for (unsigned int v = 0; v < cell2D.second.Extremes.cols(); v++)
        {
          result.Markers[0][cell2D.second.Extremes(0, v)] = cell2D.second.Marker;
          result.Markers[1][cell2D.second.Extremes(1, v)] = cell2D.second.Marker;
        }
      }
    }

    return result;
  }
  // ***************************************************************************
  VtkMeshInterface::VtkMesh3D VtkMeshInterface::ImportMesh3DFromFile(const std::string& vtkFilePath) const
  {
#if ENABLE_VTK == 1
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
#else
    Gedim::Utilities::Unused(vtkFilePath);
#endif

    return {};
  }
  // ***************************************************************************
}
