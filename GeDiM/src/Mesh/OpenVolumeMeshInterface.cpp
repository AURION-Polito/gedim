#include "FileTextReader.hpp"
#include "OpenVolumeMeshInterface.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  OpenVolumeMeshInterface::OpenVolumeMeshInterface()
  {
  }
  OpenVolumeMeshInterface::~OpenVolumeMeshInterface()
  {
  }
  // ***************************************************************************
  OpenVolumeMeshInterface::OVMMesh OpenVolumeMeshInterface::ConvertOVMMesh(const std::vector<std::string>& fileLines) const
  {
    OVMMesh mesh;

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
      mesh.Cell3Ds[p].resize(numPolyhedronVertices);
      for (unsigned int pv = 0; pv < numPolyhedronVertices; pv++)
        converter >> mesh.Cell3Ds[p][pv];
    }

    return mesh;
  }
  // ***************************************************************************
  void OpenVolumeMeshInterface::ConvertGedimMesh(const OVMMesh& meshImported,
                                                 IMeshDAO& convertedMesh) const
  {
    convertedMesh.InitializeDimension(3);

    // Create Cell0Ds
    convertedMesh.Cell0DsInitialize(meshImported.NumCell0Ds);
    convertedMesh.Cell0DsInsertCoordinates(meshImported.Cell0Ds);
    for (unsigned int v = 0; v < meshImported.NumCell0Ds; v++)
      convertedMesh.Cell0DSetState(v, true);

    // Create Cell1Ds
    convertedMesh.Cell1DsInitialize(meshImported.NumCell1Ds);
    for (unsigned int e = 0; e < meshImported.NumCell1Ds; e++)
    {
      const unsigned int ovm_edgeIndex = 2 * e;
      convertedMesh.Cell1DInsertExtremes(e,
                                         meshImported.Cell1Ds(0, ovm_edgeIndex),
                                         meshImported.Cell1Ds(1, ovm_edgeIndex));
      convertedMesh.Cell1DSetState(e, true);
    }

    // Create Cell2Ds
    convertedMesh.Cell2DsInitialize(meshImported.NumCell2Ds);
    for (unsigned int f = 0; f < meshImported.NumCell2Ds; f++)
    {
      const Eigen::MatrixXi& ovm_face = meshImported.Cell2Ds[f];
      const unsigned int& numFaceVertices = ovm_face.cols();
      convertedMesh.Cell2DInitializeVertices(f, numFaceVertices);
      convertedMesh.Cell2DInitializeEdges(f, numFaceVertices);

      for (unsigned int fv = 0; fv < numFaceVertices; fv++)
      {
        const unsigned int& ovm_edgeIndex = ovm_face(1, fv);
        const unsigned int edgeIndex = ovm_edgeIndex % 2 == 0 ? ovm_edgeIndex / 2 :
                                                                (ovm_edgeIndex - 1) / 2;

        convertedMesh.Cell2DInsertVertex(f, fv, ovm_face(0, fv));
        convertedMesh.Cell2DInsertEdge(f, fv, edgeIndex);
      }

      convertedMesh.Cell2DSetState(f, true);
    }

    // Create Cell3Ds
    convertedMesh.Cell3DsInitialize(meshImported.NumCell3Ds);
    for (unsigned int p = 0; p < meshImported.NumCell3Ds; p++)
    {
      const vector<unsigned int>& ovm_faces = meshImported.Cell3Ds[p];
      const unsigned int& numPolyhedronFaces = ovm_faces.size();

      std::unordered_set<unsigned int> cell3DVertices, cell3DEdges;

      convertedMesh.Cell3DInitializeFaces(p, numPolyhedronFaces);
      for (unsigned pf = 0; pf < numPolyhedronFaces; pf++)
      {
        const unsigned int& ovm_faceIndex = ovm_faces[pf];
        const unsigned int faceIndex = ovm_faceIndex % 2 == 0 ? ovm_faceIndex / 2 :
                                                                (ovm_faceIndex - 1) / 2;

        for (unsigned int fv = 0; fv < convertedMesh.Cell2DNumberVertices(faceIndex); fv++)
        {
          const unsigned int faceVertex = convertedMesh.Cell2DVertex(faceIndex, fv);
          const unsigned int faceEdge = convertedMesh.Cell2DEdge(faceIndex, fv);

          if (cell3DVertices.find(faceVertex) == cell3DVertices.end())
            cell3DVertices.insert(faceVertex);

          if (cell3DEdges.find(faceEdge) == cell3DEdges.end())
            cell3DEdges.insert(faceEdge);
        }

        convertedMesh.Cell3DInsertFace(p, pf, faceIndex);
      }

      convertedMesh.Cell3DInitializeVertices(p, cell3DVertices.size());
      unsigned int v = 0;
      for (const unsigned int& vertexIndex : cell3DVertices)
        convertedMesh.Cell3DInsertVertex(p, v++, vertexIndex);

      convertedMesh.Cell3DInitializeEdges(p, cell3DEdges.size());
      unsigned int e = 0;
      for (const unsigned int& edgeIndex : cell3DEdges)
        convertedMesh.Cell3DInsertEdge(p, e++, edgeIndex);

      convertedMesh.Cell3DSetState(p, true);
    }
  }
  // ***************************************************************************
  void OpenVolumeMeshInterface::ImportMesh(IMeshDAO& mesh,
                                           const std::string& ovmFilePath) const
  {
    vector<string> fileLines;

    {
      FileReader fileReader(ovmFilePath);

      if (!fileReader.Open())
        throw runtime_error("File " + ovmFilePath + " not found");

      fileReader.GetAllLines(fileLines);
      fileReader.Close();
    }

    const OVMMesh meshImported = ConvertOVMMesh(fileLines);
    ConvertGedimMesh(meshImported,
                     mesh);
  }
  // ***************************************************************************
}
