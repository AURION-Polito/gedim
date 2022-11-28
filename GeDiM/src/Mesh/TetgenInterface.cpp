#include "TetgenInterface.hpp"
#include <unordered_set>

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  TetgenInterface::TetgenInterface()
  {
  }
  TetgenInterface::~TetgenInterface()
  {
  }
  // ***************************************************************************
  void TetgenInterface::DeleteTetgenStructure(tetgenio& tetgenInput,
                                              tetgenio& ) const
  {
    delete[] tetgenInput.pointlist; tetgenInput.pointlist = NULL;
    delete[] tetgenInput.pointmarkerlist; tetgenInput.pointmarkerlist = NULL;

    delete[] tetgenInput.edgelist; tetgenInput.edgelist = NULL;
    delete[] tetgenInput.edgemarkerlist; tetgenInput.edgemarkerlist = NULL;

    for (int f = 0; f < tetgenInput.numberoffacets; f++)
    {
      tetgenio::facet* tetgenFace = &tetgenInput.facetlist[f];

      tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];
      delete[] tetgenPolygon->vertexlist; tetgenPolygon->vertexlist = NULL;

      delete[] tetgenFace->polygonlist; tetgenFace->polygonlist = NULL;
    }

    delete[] tetgenInput.facetlist; tetgenInput.facetlist = NULL;
  }
  // ***************************************************************************
  void TetgenInterface::CreateMesh(const Eigen::MatrixXd& polyhedronVertices,
                                   const Eigen::MatrixXi& polyhedronEdges,
                                   const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                   const double& maxTetrahedronVolume,
                                   IMeshDAO& mesh) const
  {
    tetgenio* tetgenInput = new tetgenio();
    tetgenio* tetgenOutput = new tetgenio();

    CreateTetgenInput(polyhedronVertices,
                      polyhedronEdges,
                      polyhedronFaces,
                      *tetgenInput);
    CreateTetgenOutput(maxTetrahedronVolume,
                       *tetgenInput,
                       *tetgenOutput);

    ConvertTetgenOutputToMeshDAO(*tetgenOutput,
                                 mesh);

    DeleteTetgenStructure(*tetgenInput,
                          *tetgenOutput);
    delete tetgenInput;
    delete tetgenOutput;
  }
  // ***************************************************************************
  void TetgenInterface::CreateTetgenInput(const Eigen::MatrixXd& polyhedronVertices,
                                          const Eigen::MatrixXi& polyhedronEdges,
                                          const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                          tetgenio& tetgenInput,
                                          const MatrixXd& constrainedPoints,
                                          const std::vector<Eigen::VectorXi>& constrainedFaces) const
  {
    const unsigned int& numberOfVertices = polyhedronVertices.cols();
    const unsigned int numberOfConstrainedPoints = constrainedPoints.cols();

    const unsigned int& numberOfEdges = polyhedronEdges.cols();
    const unsigned int& numberOfFaces = polyhedronFaces.size();
    const unsigned int numberOfConstrainedFaces = constrainedFaces.size();

    Output::Assert(numberOfVertices > 0 && numberOfEdges > 0 && numberOfFaces > 0);

    tetgenInput.firstnumber = 0;
    tetgenInput.numberofpoints = numberOfVertices + numberOfConstrainedPoints;
    tetgenInput.pointlist = new REAL[(numberOfVertices + numberOfConstrainedPoints) * 3];
    tetgenInput.pointmarkerlist = new int[numberOfVertices + numberOfConstrainedPoints];

    tetgenInput.numberofedges = numberOfEdges;
    tetgenInput.edgelist = new int[(numberOfEdges)* 2];
    tetgenInput.edgemarkerlist = new int[(numberOfEdges )];

    tetgenInput.numberoffacets = numberOfFaces + numberOfConstrainedFaces;
    tetgenInput.facetlist = new tetgenio::facet[numberOfFaces + numberOfConstrainedFaces];
    tetgenInput.facetmarkerlist = new int[numberOfFaces + numberOfConstrainedFaces];

    double* point_list = tetgenInput.pointlist;
    int* point_markerlist = tetgenInput.pointmarkerlist;

    int* edge_list = tetgenInput.edgelist;
    int* edge_markerlist = tetgenInput.edgemarkerlist;

    tetgenio::facet* face_list = tetgenInput.facetlist;
    int* face_markerlist = tetgenInput.facetmarkerlist;

    for (unsigned int v = 0; v < numberOfVertices; v++)
    {
      const Eigen::Vector3d& point = polyhedronVertices.col(v);
      point_list[3 * v] = point(0);
      point_list[3 * v + 1] = point(1);
      point_list[3 * v + 2] = point(2);

      point_markerlist[v] = v + 1;
    }

    if(numberOfConstrainedPoints > 0)
    {
      unsigned int localOffset = numberOfVertices;
      for (unsigned int j = 0; j < constrainedPoints.cols(); j++)
      {
        point_list[3 * (localOffset + j)] = constrainedPoints(0, j);
        point_list[3 * (localOffset + j) + 1] = constrainedPoints(1, j);
        point_list[3 * (localOffset + j) + 2] = constrainedPoints(2, j);

        point_markerlist[(localOffset + j)] = numberOfVertices +
                                              j + 1;
      }
    }

    for (unsigned int e = 0; e < numberOfEdges; e++)
    {
      const unsigned int& vId1 = polyhedronEdges(0, e);
      const unsigned int& vId2 = polyhedronEdges(1, e);

      Output::Assert(vId1 < numberOfVertices && vId2 < numberOfVertices);

      edge_list[2 * e] = vId1;
      edge_list[2 * e + 1] = vId2;

      edge_markerlist[e]= numberOfVertices +
                          numberOfConstrainedPoints +
                          e + 1;
    }

    for (unsigned int f = 0; f < numberOfFaces; f++)
    {
      tetgenio::facet* tetgenFace = &face_list[f];

      tetgenFace->numberofpolygons = 1;
      tetgenFace->polygonlist = new tetgenio::polygon[1];
      tetgenFace->numberofholes = 0;
      tetgenFace->holelist = NULL;

      tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];

      const MatrixXi& face = polyhedronFaces[f];
      const size_t numberFacePoints = face.cols();
      tetgenPolygon->numberofvertices = numberFacePoints;
      tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

      for (unsigned int v = 0; v < numberFacePoints; v++)
      {
        const unsigned int& vId =  face(0, v);
        Output::Assert(vId < numberOfVertices);

        tetgenPolygon->vertexlist[v] = vId;
      }

      face_markerlist[f] = numberOfVertices +
                           numberOfConstrainedPoints +
                           numberOfEdges +
                           f + 1;
    }

    if(constrainedFaces.size() > 0)
    {
      unsigned int localOffset = numberOfFaces;
      unsigned int localOffsetPoints = numberOfVertices;
      for(unsigned int numFac = 0; numFac < constrainedFaces.size(); numFac++)
      {
        const VectorXi& localIds = constrainedFaces[numFac];
        tetgenio::facet* tetgenFace = &face_list[localOffset + numFac];

        tetgenFace->numberofpolygons = 1;
        tetgenFace->polygonlist = new tetgenio::polygon[1];
        tetgenFace->numberofholes = 0;
        tetgenFace->holelist = NULL;

        tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];

        const size_t numberFacePoints = localIds.size();
        tetgenPolygon->numberofvertices = numberFacePoints;
        tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

        for (unsigned int v = 0; v < numberFacePoints; v++)
        {
          unsigned int vId = localIds(v) + localOffsetPoints;
          tetgenPolygon->vertexlist[v] = vId;
        }
        face_markerlist[localOffset + numFac] = numberOfVertices +
                                                numberOfConstrainedPoints +
                                                numberOfEdges +
                                                numberOfFaces +
                                                numFac + 1;
      }
    }
  }
  // ***************************************************************************
  void TetgenInterface::CreateTetgenOutput(const double& maxTetrahedronArea,
                                           tetgenio& tetgenInput,
                                           tetgenio& tetgenOutput,
                                           const std::string& tetgenOptions) const
  {
    Output::Assert(maxTetrahedronArea > 0.0);

    tetgenbehavior b;

    ostringstream options;
    options.precision(16);
    options<< tetgenOptions;
    options<< maxTetrahedronArea;
    size_t sizeOptions = options.str().size();
    char* optionPointer = new char[sizeOptions + 1];
    options.str().copy(optionPointer, sizeOptions);
    optionPointer[sizeOptions] = '\0';

    b.parse_commandline(optionPointer);

    tetrahedralize(&b, &tetgenInput, &tetgenOutput);

    delete[] optionPointer;
  }
  // ***************************************************************************
  void TetgenInterface::ConvertTetgenOutputToMeshDAO(const tetgenio& tetgenOutput,
                                                     IMeshDAO& mesh) const
  {
    /// <li>	Fill mesh structures
    unsigned int numberOfCellsMesh = tetgenOutput.numberoftetrahedra;
    unsigned int numberOfFacesMesh = tetgenOutput.numberoftrifaces;
    unsigned int numberOfEgdesMesh = tetgenOutput.numberofedges;
    unsigned int numberOfPointsMesh = tetgenOutput.numberofpoints;

    mesh.InitializeDimension(3);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEgdesMesh);
    mesh.Cell2DsInitialize(numberOfFacesMesh);
    mesh.Cell3DsInitialize(numberOfCellsMesh);

    /// <li> Set Cell0Ds
    for (unsigned int p = 0; p < numberOfPointsMesh; p++)
    {
      mesh.Cell0DSetState(p, true);
      mesh.Cell0DInsertCoordinates(p,
                                   Vector3d(tetgenOutput.pointlist[3 * p],
                                   tetgenOutput.pointlist[3 * p + 1],
          tetgenOutput.pointlist[3 * p + 2]));
      mesh.Cell0DSetMarker(p, tetgenOutput.pointmarkerlist[p]);
    }

    /// <li> Set Cell1Ds
    for(unsigned int e = 0; e < numberOfEgdesMesh; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                tetgenOutput.edgelist[2 * e + 0],
          tetgenOutput.edgelist[2 * e + 1]);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, tetgenOutput.edgemarkerlist[e]);
    }

    /// <li> Set Faces
    vector<unsigned int> vertices(3);
    vector<unsigned int> edgeEndPoints(2);

    unsigned long estimatedValue = 2* numberOfPointsMesh * numberOfPointsMesh + 2*numberOfPointsMesh + 1;
    SparseMatrix<int> connectivityPointsFaces(numberOfPointsMesh,estimatedValue);
    vector< Triplet<int, unsigned long> > tripletsFaces;
    tripletsFaces.reserve(numberOfFacesMesh);
    for(unsigned int f = 0; f < numberOfFacesMesh; f++)
    {
      const unsigned int numCell2DVertices = 3;
      mesh.Cell2DInitializeVertices(f, numCell2DVertices);
      mesh.Cell2DInitializeEdges(f, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int vertexId = tetgenOutput.trifacelist[3 * f + v];
        const unsigned int nextVertexId = tetgenOutput.trifacelist[3 * f + (v + 1) % 3];
        const unsigned int edgeId = mesh.Cell1DExists(vertexId,
                                                      nextVertexId) ?
                                      mesh.Cell1DByExtremes(vertexId,
                                                            nextVertexId) :
                                      mesh.Cell1DByExtremes(nextVertexId,
                                                            vertexId);
        vertices[v] = vertexId;

        mesh.Cell2DInsertVertex(f, v, vertexId);
        mesh.Cell2DInsertEdge(f, v, edgeId);
        if (mesh.Cell1DMarker(edgeId) == 0)
          mesh.Cell1DSetMarker(edgeId, tetgenOutput.trifacemarkerlist[f]);
      }

      mesh.Cell2DSetState(f, true);
      mesh.Cell2DSetMarker(f,
                           tetgenOutput.trifacemarkerlist[f]);

      sort(vertices.begin(), vertices.end());
      unsigned long indexI = vertices[0];
      unsigned long indexJK = (vertices[1] + vertices[2]) * (vertices[1] + vertices[2] + 1) * 0.5 + vertices[2] + 1 ;
      tripletsFaces.push_back(Triplet<int, unsigned long>(indexI,
                                                          indexJK,
                                                          f + 1));

    }

    connectivityPointsFaces.setFromTriplets(tripletsFaces.begin(), tripletsFaces.end());
    connectivityPointsFaces.makeCompressed();

    /// <li> Set Cells
    vector<unsigned int> faceVertices(3);
    for (unsigned int c = 0; c < numberOfCellsMesh; c++)
    {
      const unsigned int numVertices = 4;
      const unsigned int numEdges = 6;
      const unsigned int numFaces = 4;

      mesh.Cell3DInitializeVertices(c, numVertices);
      mesh.Cell3DInitializeEdges(c, numEdges);
      mesh.Cell3DInitializeFaces(c, numFaces);

      for (unsigned int v = 0; v < numVertices; v++)
        mesh.Cell3DInsertVertex(c, v, tetgenOutput.tetrahedronlist[tetgenOutput.numberofcorners * c + v]);

      unordered_set<unsigned int> cell3DEdges;      // find cell edges and faces
      for (unsigned int j = 0; j < numFaces; j++)
      {
        for (unsigned int k = 0; k < 3; k++)
        {
          const unsigned int faceVertexId =(mesh.Cell3DVertex(c, (j + k) % 4));
          faceVertices[k] = faceVertexId;
        }

        sort(faceVertices.begin(),faceVertices.end());
        const unsigned int indexI = faceVertices[0];
        const unsigned int indexJK = (faceVertices[1] + faceVertices[2]) * (faceVertices[1] + faceVertices[2] + 1) * 0.5 + faceVertices[2] + 1 ;

        const int faceId = connectivityPointsFaces.coeff(indexI, indexJK) - 1;
        mesh.Cell3DInsertFace(c, j, faceId);

        if (j < 3)
        {
          for (unsigned int e = 0; e < 3; e++)
          {
            const unsigned int cell1D = mesh.Cell2DEdge(faceId, e);
            if (cell3DEdges.find(cell1D) != cell3DEdges.end())
              continue;

            mesh.Cell3DInsertEdge(c, cell3DEdges.size(), cell1D);
            cell3DEdges.insert(cell1D);
          }
        }
      }

      mesh.Cell3DSetState(c, true);
      mesh.Cell3DSetMarker(c, 0);
    }
  }
  // ***************************************************************************
  void TetgenInterface::ExportTetgenOutput(const string& nameFolder,
                                           const string& nameFile,
                                           tetgenio& tetgenOutput) const
  {

    ostringstream nameFolderStream, nameFileStream;

    nameFolderStream<< nameFolder<< "/";
    nameFolderStream<< "Tetgen/";

    Output::CreateFolder(nameFolderStream.str());

    nameFileStream<< nameFolderStream.str()<< nameFile;

    Output::CreateFolder(nameFolderStream.str());

    tetgenOutput.firstnumber = 0;
    tetgenOutput.save_nodes((char*)nameFileStream.str().c_str());
    tetgenOutput.save_elements((char*)nameFileStream.str().c_str());
    tetgenOutput.save_faces((char*)nameFileStream.str().c_str());
    tetgenOutput.save_edges((char*)nameFileStream.str().c_str());
  }
  // ***************************************************************************
}
