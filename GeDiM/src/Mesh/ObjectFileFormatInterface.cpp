#include "FileTextReader.hpp"
#include "ObjectFileFormatInterface.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  ObjectFileFormatInterface::ObjectFileFormatInterface()
  {
  }
  ObjectFileFormatInterface::~ObjectFileFormatInterface()
  {
  }
  // ***************************************************************************
  ObjectFileFormatInterface::OFFMesh ObjectFileFormatInterface::StringsToOFFMesh(const std::vector<std::string>& fileLines) const
  {
    OFFMesh mesh;

    unsigned int lineCounter = 0;

    // convert vertices
    lineCounter += 1;
    Gedim::Output::Assert(fileLines.size() >= lineCounter);

    {
      istringstream converter(fileLines[lineCounter++]);

      converter >> mesh.NumCell0Ds>> mesh.NumCell2Ds>> mesh.NumCell1Ds;
    }

    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell0Ds);

    mesh.Cell0Ds.resize(3, mesh.NumCell0Ds);
    for (unsigned int v = 0; v < mesh.NumCell0Ds; v++)
    {
      istringstream converter(fileLines[lineCounter++]);

      converter >> mesh.Cell0Ds(0, v)>> mesh.Cell0Ds(1, v)>> mesh.Cell0Ds(2, v);
    }

    // convert faces
    Gedim::Output::Assert(fileLines.size() >= lineCounter + mesh.NumCell2Ds);

    mesh.Cell2Ds.resize(mesh.NumCell2Ds);
    for (unsigned int f = 0; f < mesh.NumCell2Ds; f++)
    {
      istringstream converter(fileLines[lineCounter++]);

      unsigned int numFaceVertices = 0;
      converter >> numFaceVertices;
      mesh.Cell2Ds[f].resize(numFaceVertices);
      for (unsigned int fv = 0; fv < numFaceVertices; fv++)
        converter >> mesh.Cell2Ds[f][fv];
    }

    return mesh;
  }
  // ***************************************************************************
  std::vector<string> ObjectFileFormatInterface::OFFMeshToStrings(const OFFMesh& mesh) const
  {
    list<string> lines;

    lines.push_back("OFF");

    lines.push_back(std::to_string(mesh.NumCell0Ds) + " " +
                    std::to_string(mesh.NumCell2Ds) + " " +
                    std::to_string(mesh.NumCell1Ds));

    for (unsigned int v = 0; v < mesh.NumCell0Ds; v++)
    {
      ostringstream stream;
      stream.precision(16);
      stream<< scientific<< mesh.Cell0Ds(0, v)<< " ";
      stream<< scientific<< mesh.Cell0Ds(1, v)<< " ";
      stream<< scientific<< mesh.Cell0Ds(2, v);

      lines.push_back(stream.str());
    }

    for (unsigned int f = 0; f < mesh.NumCell2Ds; f++)
    {
      ostringstream stream;
      stream<< mesh.Cell2Ds[f].size();
      for (unsigned int fv = 0; fv < mesh.Cell2Ds[f].size(); fv++)
        stream<< " "<< mesh.Cell2Ds[f][fv];

      lines.push_back(stream.str());
    }

    return std::vector<string>(lines.begin(),
                               lines.end());
  }
  // ***************************************************************************
  ObjectFileFormatInterface::OFFMesh ObjectFileFormatInterface::MeshDAOToOFFMesh(const IMeshDAO& originalMesh) const
  {
    ObjectFileFormatInterface::OFFMesh convertedMesh;

    convertedMesh.NumCell0Ds = originalMesh.Cell0DTotalNumber();
    convertedMesh.NumCell1Ds = originalMesh.Cell1DTotalNumber();
    convertedMesh.NumCell2Ds = originalMesh.Cell2DTotalNumber();

    convertedMesh.Cell0Ds = originalMesh.Cell0DsCoordinates();

    convertedMesh.Cell2Ds.resize(convertedMesh.NumCell2Ds);
    for (unsigned int f = 0; f < convertedMesh.NumCell2Ds; f++)
    {
      const std::vector<unsigned int> vertices = originalMesh.Cell2DVertices(f);
      convertedMesh.Cell2Ds[f].resize(vertices.size());
      for (unsigned int v = 0; v < vertices.size(); v++)
        convertedMesh.Cell2Ds[f][v] = vertices[v];
    }

    return convertedMesh;
  }
  // ***************************************************************************
  void ObjectFileFormatInterface::OFFMeshToMeshDAO(const OFFMesh& originalMesh,
                                                   const MeshUtilities& meshUtilities,
                                                   IMeshDAO& convertedMesh) const
  {
    convertedMesh.InitializeDimension(2);

    // Create Cell0Ds
    convertedMesh.Cell0DsInitialize(originalMesh.NumCell0Ds);
    convertedMesh.Cell0DsInsertCoordinates(originalMesh.Cell0Ds);
    for (unsigned int v = 0; v < originalMesh.NumCell0Ds; v++)
      convertedMesh.Cell0DSetState(v, true);

    // Compute Cell1Ds
    const MeshUtilities::ComputeMesh2DCell1DsResult cell1Ds = meshUtilities.ComputeMesh2DCell1Ds(originalMesh.Cell0Ds,
                                                                                                 originalMesh.Cell2Ds);

    // Create Cell1Ds
    convertedMesh.Cell1DsInitialize(cell1Ds.Cell1Ds.cols());
    convertedMesh.Cell1DsInsertExtremes(cell1Ds.Cell1Ds);
    for (unsigned int e = 0; e < cell1Ds.Cell1Ds.cols(); e++)
      convertedMesh.Cell1DSetState(e, true);

    // Create Cell2Ds
    {
      std::vector<unsigned int> facesNumberVertices(originalMesh.NumCell2Ds);
      for (unsigned int f = 0; f < originalMesh.NumCell2Ds; f++)
      {
        const Eigen::MatrixXi& off_face = cell1Ds.Cell2Ds[f];
        facesNumberVertices[f] = off_face.cols();
      }

      convertedMesh.Cell2DsInitialize(originalMesh.NumCell2Ds);
      convertedMesh.Cell2DsInitializeVertices(facesNumberVertices);
      convertedMesh.Cell2DsInitializeEdges(facesNumberVertices);
      for (unsigned int f = 0; f < originalMesh.NumCell2Ds; f++)
      {
        const Eigen::MatrixXi& off_face = cell1Ds.Cell2Ds[f];
        const unsigned int numFaceVertices = off_face.cols();

        for (unsigned int fv = 0; fv < numFaceVertices; fv++)
        {
          convertedMesh.Cell2DInsertVertex(f, fv, off_face(0, fv));
          convertedMesh.Cell2DInsertEdge(f, fv, off_face(1, fv));
        }

        convertedMesh.Cell2DSetState(f, true);
      }
    }
  }
  // ***************************************************************************
  void ObjectFileFormatInterface::ImportMeshFromFile(const std::string& offFilePath,
                                                     const MeshUtilities& meshUtilities,
                                                     IMeshDAO& mesh) const
  {
    vector<string> fileLines;

    {
      FileReader fileReader(offFilePath);

      if (!fileReader.Open())
        throw runtime_error("File " + offFilePath + " not found");

      fileReader.GetAllLines(fileLines);
      fileReader.Close();
    }

    const OFFMesh meshImported = StringsToOFFMesh(fileLines);
    OFFMeshToMeshDAO(meshImported,
                     meshUtilities,
                     mesh);
  }
  // ***************************************************************************
  void ObjectFileFormatInterface::ExportMeshToFile(const IMeshDAO& mesh,
                                                   const std::string& offFilePath) const
  {
    const OFFMesh meshToExport = MeshDAOToOFFMesh(mesh);
    const vector<string> fileLines = OFFMeshToStrings(meshToExport);

    ofstream exportFile(offFilePath);
    if (!exportFile.is_open())
      throw runtime_error("Unable to export file " + offFilePath);

    for (const string& line : fileLines)
      exportFile << line<< endl;
    exportFile.close();
  }
  // ***************************************************************************
}
