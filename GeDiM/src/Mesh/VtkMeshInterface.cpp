#include "FileTextReader.hpp"
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
  VtkMeshInterface::VtkMesh VtkMeshInterface::StringsToVtkMesh(const std::vector<std::string>& fileLines) const
  {
    VtkMesh mesh;

    unsigned int lineCounter = 0;

    // convert vertices
    lineCounter += 2;
    Gedim::Output::Assert(fileLines.size() >= lineCounter);
    mesh.NumCell0Ds = std::atoi(fileLines[lineCounter++].c_str());

    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell0Ds);

    mesh.Cell0Ds.resize(3, mesh.NumCell0Ds);
    for (unsigned int v = 0; v < mesh.NumCell0Ds; v++)
    {
      istringstream converter(fileLines[lineCounter++]);

      converter >> mesh.Cell0Ds(0, v)>> mesh.Cell0Ds(1, v)>> mesh.Cell0Ds(2, v);
    }

    // convert edges
    lineCounter += 1;
    Gedim::Output::Assert(fileLines.size() >= lineCounter);
    mesh.NumCell1Ds = std::atoi(fileLines[lineCounter++].c_str());

    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell1Ds);

    mesh.Cell1Ds.resize(2, 2 * mesh.NumCell1Ds);
    for (unsigned int e = 0; e < mesh.NumCell1Ds; e++)
    {
      istringstream converter(fileLines[lineCounter++]);

      unsigned int edgeOrigin, edgeEnd;
      converter >> edgeOrigin>> edgeEnd;

      mesh.Cell1Ds.col(2 * e)<< edgeOrigin, edgeEnd;
      mesh.Cell1Ds.col(2 * e + 1)<< edgeEnd, edgeOrigin;
    }

    // convert faces
    lineCounter += 1;
    Gedim::Output::Assert(fileLines.size() >= lineCounter);
    mesh.NumCell2Ds = std::atoi(fileLines[lineCounter++].c_str());

    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell2Ds);

    mesh.Cell2Ds.resize(mesh.NumCell2Ds);
    for (unsigned int f = 0; f < mesh.NumCell2Ds; f++)
    {
      istringstream converter(fileLines[lineCounter++]);

      unsigned int numFaceEdges = 0;
      converter >> numFaceEdges;
      mesh.Cell2Ds[f].resize(2, numFaceEdges);
      for (unsigned int fv = 0; fv < numFaceEdges; fv++)
      {
        unsigned int edgeFileIndex = 0;
        converter >> edgeFileIndex;
        const unsigned int edgeOriginIndex = mesh.Cell1Ds(0, edgeFileIndex);

        mesh.Cell2Ds[f].col(fv)<< edgeOriginIndex, edgeFileIndex;
      }
    }

    // convert polyhedra
    lineCounter += 1;
    Gedim::Output::Assert(fileLines.size() >= lineCounter);
    mesh.NumCell3Ds = std::atoi(fileLines[lineCounter++].c_str());

    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell3Ds);

    mesh.Cell3Ds.resize(mesh.NumCell3Ds);
    for (unsigned int p = 0; p < mesh.NumCell3Ds; p++)
    {
      istringstream converter(fileLines[lineCounter++]);

      unsigned int numPolyhedronVertices = 0;
      converter >> numPolyhedronVertices;
      mesh.Cell3Ds[p].FacesIndex.resize(numPolyhedronVertices);
      for (unsigned int pv = 0; pv < numPolyhedronVertices; pv++)
      {
        converter >> mesh.Cell3Ds[p].FacesIndex[pv];
      }
    }

    return mesh;
  }
  // ***************************************************************************
  std::vector<string> VtkMeshInterface::VtkMeshToStrings(const VtkMesh& mesh) const
  {
    list<string> lines;

    lines.push_back("Vtk ASCII");

    lines.push_back("Vertices");
    lines.push_back(to_string(mesh.NumCell0Ds));
    for (unsigned int v = 0; v < mesh.NumCell0Ds; v++)
    {
      ostringstream stream;
      stream.precision(16);
      stream<< scientific<< mesh.Cell0Ds(0, v)<< " ";
      stream<< scientific<< mesh.Cell0Ds(1, v)<< " ";
      stream<< scientific<< mesh.Cell0Ds(2, v);

      lines.push_back(stream.str());
    }

    lines.push_back("Edges");
    lines.push_back(to_string(mesh.NumCell1Ds));
    for (unsigned int e = 0; e < mesh.NumCell1Ds; e++)
    {
      ostringstream stream;
      stream<< mesh.Cell1Ds(0, 2 * e)<< " ";
      stream<< mesh.Cell1Ds(1, 2 * e);

      lines.push_back(stream.str());
    }

    lines.push_back("Faces");
    lines.push_back(to_string(mesh.NumCell2Ds));
    for (unsigned int f = 0; f < mesh.NumCell2Ds; f++)
    {
      ostringstream stream;
      stream<< mesh.Cell2Ds[f].cols();
      for (unsigned int fv = 0; fv < mesh.Cell2Ds[f].cols(); fv++)
        stream<< " "<< mesh.Cell2Ds[f](1, fv);

      lines.push_back(stream.str());
    }

    lines.push_back("Polyhedra");
    lines.push_back(to_string(mesh.NumCell3Ds));
    for (unsigned int p = 0; p < mesh.NumCell3Ds; p++)
    {
      ostringstream stream;
      stream<< mesh.Cell3Ds[p].FacesIndex.size();
      for (unsigned int pf = 0; pf < mesh.Cell3Ds[p].FacesIndex.size(); pf++)
        stream<< " "<< mesh.Cell3Ds[p].FacesIndex[pf];

      lines.push_back(stream.str());
    }

    return std::vector<string>(lines.begin(),
                               lines.end());
  }
  // ***************************************************************************
  VtkMeshInterface::VtkMesh VtkMeshInterface::MeshDAOToVtkMesh(const IMeshDAO& originalMesh) const
  {
    VtkMeshInterface::VtkMesh convertedMesh;

    convertedMesh.NumCell0Ds = originalMesh.Cell0DTotalNumber();
    convertedMesh.NumCell1Ds = originalMesh.Cell1DTotalNumber();
    convertedMesh.NumCell2Ds = originalMesh.Cell2DTotalNumber();
    convertedMesh.NumCell3Ds = originalMesh.Cell3DTotalNumber();

    convertedMesh.Cell0Ds = originalMesh.Cell0DsCoordinates();

    convertedMesh.Cell1Ds.resize(2, 2 * convertedMesh.NumCell1Ds);
    for (unsigned int e = 0; e < convertedMesh.NumCell1Ds; e++)
    {
      const unsigned int edgeOrigin = originalMesh.Cell1DOrigin(e);
      const unsigned int edgeEnd = originalMesh.Cell1DEnd(e);

      convertedMesh.Cell1Ds.col(2 * e)<< edgeOrigin, edgeEnd;
      convertedMesh.Cell1Ds.col(2 * e + 1)<< edgeEnd, edgeOrigin;
    }

    convertedMesh.Cell2Ds.resize(convertedMesh.NumCell2Ds);
    for (unsigned int f = 0; f < convertedMesh.NumCell2Ds; f++)
    {
      const std::vector<unsigned int> vertices = originalMesh.Cell2DVertices(f);
      const std::vector<unsigned int> edges = originalMesh.Cell2DEdges(f);
      convertedMesh.Cell2Ds[f].resize(2, vertices.size());
      for (unsigned int v = 0; v < vertices.size(); v++)
      {
        convertedMesh.Cell2Ds[f](0, v) = vertices[v];
        const bool edgeDirection = (originalMesh.Cell1DOrigin(edges[v]) == vertices[v]);
        convertedMesh.Cell2Ds[f](1, v) = edgeDirection ? 2 * edges[v] :
                                                         2 * edges[v] + 1;
      }
    }

    convertedMesh.Cell3Ds.resize(convertedMesh.NumCell3Ds);
    for (unsigned int p = 0; p < convertedMesh.NumCell3Ds; p++)
    {
      const std::vector<unsigned int> faces = originalMesh.Cell3DFaces(p);

      convertedMesh.Cell3Ds[p].FacesIndex.resize(faces.size());
      for (unsigned int fp = 0; fp < faces.size(); fp++)
      {
        const bool faceOrientation = true;
        convertedMesh.Cell3Ds[p].FacesIndex[fp] = faceOrientation ? 2 * faces[fp] :
                                                                    2 * faces[fp] + 1;
      }
    }

    return convertedMesh;
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

    // Create Cell1Ds
    convertedMesh.Cell1DsInitialize(originalMesh.NumCell1Ds);
    for (unsigned int e = 0; e < originalMesh.NumCell1Ds; e++)
    {
      const unsigned int Vtk_edgeIndex = 2 * e;
      convertedMesh.Cell1DInsertExtremes(e,
                                         originalMesh.Cell1Ds(0, Vtk_edgeIndex),
                                         originalMesh.Cell1Ds(1, Vtk_edgeIndex));
      convertedMesh.Cell1DSetState(e, true);
    }

    // Create Cell2Ds
    {
      std::vector<unsigned int> numberVertices(originalMesh.NumCell2Ds), numberEdges(originalMesh.NumCell2Ds);
      for (unsigned int f = 0; f < originalMesh.NumCell2Ds; f++)
      {
        const Eigen::MatrixXi& Vtk_face = originalMesh.Cell2Ds[f];

        numberVertices[f] = Vtk_face.cols();
        numberEdges[f] = Vtk_face.cols();
      }

      convertedMesh.Cell2DsInitialize(originalMesh.NumCell2Ds);
      convertedMesh.Cell2DsInitializeVertices(numberVertices);
      convertedMesh.Cell2DsInitializeEdges(numberEdges);
      for (unsigned int f = 0; f < originalMesh.NumCell2Ds; f++)
      {
        const Eigen::MatrixXi& Vtk_face = originalMesh.Cell2Ds[f];
        const unsigned int numFaceVertices = Vtk_face.cols();

        for (unsigned int fv = 0; fv < numFaceVertices; fv++)
        {
          const unsigned int& Vtk_edgeIndex = Vtk_face(1, fv);
          const unsigned int edgeIndex = Vtk_edgeIndex % 2 == 0 ? Vtk_edgeIndex / 2 :
                                                                  (Vtk_edgeIndex - 1) / 2;

          convertedMesh.Cell2DInsertVertex(f, fv, Vtk_face(0, fv));
          convertedMesh.Cell2DInsertEdge(f, fv, edgeIndex);
        }

        convertedMesh.Cell2DSetState(f, true);
      }
    }

    // Create Cell3Ds
    {
      std::vector<std::unordered_set<unsigned int>> cell3DsVertices(originalMesh.NumCell3Ds);
      std::vector<std::unordered_set<unsigned int>> cell3DsEdges(originalMesh.NumCell3Ds);
      std::vector<unsigned int> numberVertices(originalMesh.NumCell3Ds);
      std::vector<unsigned int> numberEdges(originalMesh.NumCell3Ds);
      std::vector<unsigned int> numberFaces(originalMesh.NumCell3Ds);

      for (unsigned int c = 0; c < originalMesh.NumCell3Ds; c++)
      {
        const VtkMesh::Cell3D& Vtk_faces = originalMesh.Cell3Ds[c];
        const unsigned int& numPolyhedronFaces = Vtk_faces.FacesIndex.size();

        std::unordered_set<unsigned int>& cell3DVertices = cell3DsVertices[c];
        std::unordered_set<unsigned int>& cell3DEdges = cell3DsEdges[c];

        for (unsigned pf = 0; pf < numPolyhedronFaces; pf++)
        {
          const unsigned int& Vtk_faceIndex = Vtk_faces.FacesIndex[pf];
          const unsigned int faceIndex = Vtk_faceIndex / 2;
          for (unsigned int fv = 0; fv < convertedMesh.Cell2DNumberVertices(faceIndex); fv++)
          {
            const unsigned int faceVertex = convertedMesh.Cell2DVertex(faceIndex, fv);
            const unsigned int faceEdge = convertedMesh.Cell2DEdge(faceIndex, fv);

            if (cell3DVertices.find(faceVertex) == cell3DVertices.end())
              cell3DVertices.insert(faceVertex);

            if (cell3DEdges.find(faceEdge) == cell3DEdges.end())
              cell3DEdges.insert(faceEdge);
          }
        }

        numberVertices[c] = cell3DVertices.size();
        numberEdges[c] = cell3DEdges.size();
        numberFaces[c] = Vtk_faces.FacesIndex.size();
      }

      convertedMesh.Cell3DsInitialize(originalMesh.NumCell3Ds);
      convertedMesh.Cell3DsInitializeVertices(numberVertices);
      convertedMesh.Cell3DsInitializeEdges(numberEdges);
      convertedMesh.Cell3DsInitializeFaces(numberFaces);

      for (unsigned int p = 0; p < originalMesh.NumCell3Ds; p++)
      {
        const VtkMesh::Cell3D& Vtk_faces = originalMesh.Cell3Ds[p];
        const unsigned int& numPolyhedronFaces = Vtk_faces.FacesIndex.size();

        const std::unordered_set<unsigned int>& cell3DVertices = cell3DsVertices[p];
        const std::unordered_set<unsigned int>& cell3DEdges = cell3DsEdges[p];

        for (unsigned pf = 0; pf < numPolyhedronFaces; pf++)
        {
          const unsigned int& Vtk_faceIndex = Vtk_faces.FacesIndex[pf];
          const unsigned int faceIndex = Vtk_faceIndex / 2;
          convertedMesh.Cell3DInsertFace(p, pf, faceIndex);
        }

        unsigned int v = 0;
        for (const unsigned int& vertexIndex : cell3DVertices)
          convertedMesh.Cell3DInsertVertex(p, v++, vertexIndex);

        unsigned int e = 0;
        for (const unsigned int& edgeIndex : cell3DEdges)
          convertedMesh.Cell3DInsertEdge(p, e++, edgeIndex);

        convertedMesh.Cell3DSetState(p, true);
      }
    }
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
      std::cout << "output is a ugrid" << std::endl;
      vtkUnstructuredGrid* output = reader->GetUnstructuredGridOutput();
      std::cout << "output has " << output->GetNumberOfPoints() << " points." << std::endl;
      std::cout << "output has " << output->GetNumberOfCells() << " cells." << std::endl;
      std::cout << "first cell type has " << output->GetCellType(0) << "." << std::endl;

      const unsigned int num_points = output->GetNumberOfPoints();
      const auto point_data = output->GetPoints();
      for (unsigned int p = 0; p < num_points; p++)
      {
        const auto& point = point_data[p];
        std::cout<< "Point p "<< p<< std::endl;
      }

        const unsigned int num_cells = output->GetNumberOfCells();
      const auto cell_data = output->GetCells();

    }

    //    vector<string> fileLines;

    //    {
    //      FileReader fileReader(VtkFilePath);

    //      if (!fileReader.Open())
    //        throw runtime_error("File " + VtkFilePath + " not found");

    //      fileReader.GetAllLines(fileLines);
    //      fileReader.Close();
    //    }

    //    const VtkMesh meshImported = StringsToVtkMesh(fileLines);
    //    VtkMeshToMeshDAO(meshImported,
    //                     mesh);
  }
  // ***************************************************************************
  void VtkMeshInterface::ExportMeshToFile(const IMeshDAO& mesh,
                                          const std::string& VtkFilePath) const
  {
    const VtkMesh meshToExport = MeshDAOToVtkMesh(mesh);
    const vector<string> fileLines = VtkMeshToStrings(meshToExport);

    ofstream exportFile(VtkFilePath);
    if (!exportFile.is_open())
      throw runtime_error("Unable to export file " + VtkFilePath);

    for (const string& line : fileLines)
      exportFile << line<< endl;
    exportFile.close();
  }
  // ***************************************************************************
}
