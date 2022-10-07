#include "MeshCreator3DTetgen.hpp"
#include "Output.hpp"

using namespace MainApplication;

namespace GeDiM
{
  // ***************************************************************************
  MeshCreator3DTetgen::MeshCreator3DTetgen() : MeshCreator()
  {
    inputMeshPointer = NULL;
    outputMeshPointer = NULL;
    tetgenOptions = "Qpqfezna";
  }
  MeshCreator3DTetgen::~MeshCreator3DTetgen()
  {
    if (inputMeshPointer != NULL)
    {
      delete[] inputMeshPointer->pointlist; inputMeshPointer->pointlist = NULL;
      delete[] inputMeshPointer->pointmarkerlist; inputMeshPointer->pointmarkerlist = NULL;

      delete[] inputMeshPointer->edgelist; inputMeshPointer->edgelist = NULL;
      delete[] inputMeshPointer->edgemarkerlist; inputMeshPointer->edgemarkerlist = NULL;

      for (int f = 0; f < inputMeshPointer->numberoffacets; f++)
      {
        tetgenio::facet* tetgenFace = &inputMeshPointer->facetlist[f];

        tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];
        delete[] tetgenPolygon->vertexlist; tetgenPolygon->vertexlist = NULL;

        delete[] tetgenFace->polygonlist; tetgenFace->polygonlist = NULL;
      }

      delete[] inputMeshPointer->facetlist; inputMeshPointer->facetlist = NULL;
    }

    delete inputMeshPointer; inputMeshPointer = NULL;
    delete outputMeshPointer; outputMeshPointer = NULL;
  }
  // ***************************************************************************
  Output::ExitCodes MeshCreator3DTetgen::CreateTetgenInput(const Polyhedron& domain)
  {
    delete inputMeshPointer; inputMeshPointer = NULL;

    const unsigned int& numberOfVertices = domain.NumberOfVertices();
    unsigned int numberOfConstrainedPoints = 0;
    for(unsigned int numConst = 0; numConst < constrainedPoints.size(); numConst++)
      numberOfConstrainedPoints += constrainedPoints[numConst].rows();

    const unsigned int& numberOfEdges = domain.NumberOfEdges();
    const unsigned int& numberOfFaces = domain.NumberOfFaces();
    unsigned int numberOfConstrainedFaces = 0;
    for(unsigned int numConst = 0; numConst < constrainedFaces.size(); numConst++)
      numberOfConstrainedFaces += constrainedFaces[numConst].size();

    if (numberOfVertices == 0 || numberOfEdges == 0 || numberOfFaces == 0)
    {
      Output::PrintErrorMessage("Wrong initialization of the domain %d, no vertices or faces", false, domain.GlobalId());
      return Output::GenericError;
    }

    inputMeshPointer = new tetgenio();

    inputMeshPointer->firstnumber = 0;
    inputMeshPointer->numberofpoints = numberOfVertices + numberOfConstrainedPoints;
    inputMeshPointer->pointlist = new REAL[(numberOfVertices + numberOfConstrainedPoints) * 3];
    inputMeshPointer->pointmarkerlist = new int[numberOfVertices + numberOfConstrainedPoints];

    inputMeshPointer->numberofedges = numberOfEdges;
    inputMeshPointer->edgelist = new int[(numberOfEdges)* 2];
    inputMeshPointer->edgemarkerlist = new int[(numberOfEdges )];

    inputMeshPointer->numberoffacets = numberOfFaces + numberOfConstrainedFaces;
    inputMeshPointer->facetlist = new tetgenio::facet[numberOfFaces + numberOfConstrainedFaces];
    inputMeshPointer->facetmarkerlist = new int[numberOfFaces + numberOfConstrainedFaces];

    double* point_list = inputMeshPointer->pointlist;
    int* point_markerlist = inputMeshPointer->pointmarkerlist;

    int* edge_list = inputMeshPointer->edgelist;
    int* edge_markerlist = inputMeshPointer->edgemarkerlist;

    tetgenio::facet* face_list = inputMeshPointer->facetlist;
    int* face_markerlist = inputMeshPointer->facetmarkerlist;

    for (unsigned int v = 0; v < numberOfVertices; v++)
    {
      const Point& point = (*domain.Vertex(v));
      point_list[3 * v] = point(0);
      point_list[3 * v + 1] = point(1);
      point_list[3 * v + 2] = point(2);

      point_markerlist[v] = 0;
    }

    for (unsigned int e = 0; e < numberOfEdges; e++)
    {
      const unsigned int& vId1 = domain.PositionPoint(*domain.Edge(e)->Vertex(0));
      const unsigned int& vId2 = domain.PositionPoint(*domain.Edge(e)->Vertex(1));

      if (vId1 >= numberOfVertices)
      {
        Output::PrintErrorMessage("%s: Vertex id  %dnot correct", false, __func__, vId1);
        return Output::GenericError;
      }

      if (vId2 >= numberOfVertices)
      {
        Output::PrintErrorMessage("%s: Vertex id  %dnot correct", false, __func__, vId2);
        return Output::GenericError;
      }

      edge_list[2 * e] = vId1;
      edge_list[2 * e + 1] = vId2;

      edge_markerlist[e]= 0;
    }

    for (unsigned int f = 0; f < numberOfFaces; f++)
    {
      tetgenio::facet* tetgenFace = &face_list[f];

      tetgenFace->numberofpolygons = 1;
      tetgenFace->polygonlist = new tetgenio::polygon[1];
      tetgenFace->numberofholes = 0;
      tetgenFace->holelist = NULL;

      tetgenio::polygon* tetgenPolygon =  &tetgenFace->polygonlist[0];

      const Polygon& face = *domain.Face(f);
      const size_t numberFacePoints = face.NumberOfVertices();
      tetgenPolygon->numberofvertices = numberFacePoints;
      tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

      for (unsigned int v = 0; v < numberFacePoints; v++)
      {
        const unsigned int& vId =  domain.PositionPoint(*face.Vertex(v));

        if (vId >= numberOfVertices)
        {
          Output::PrintErrorMessage("%s: Vertex id  %dnot correct", false, __func__, vId);
          return Output::GenericError;
        }

        tetgenPolygon->vertexlist[v] = vId;
      }

      face_markerlist[f] = 0;
    }

    if(numberOfConstrainedPoints > 0)
    {
      unsigned int localOffset = numberOfVertices;
      for(unsigned int numConst = 0; numConst < constrainedPoints.size(); numConst++)
      {
        MatrixXd& localConstrainedPoints = constrainedPoints[numConst];

        for (unsigned int j = 0; j < localConstrainedPoints.rows(); j++)
        {
          point_list[3 * (localOffset + j)] = localConstrainedPoints(j, 0);
          point_list[3 * (localOffset + j) + 1] = localConstrainedPoints(j, 1);
          point_list[3 * (localOffset + j) + 2] = localConstrainedPoints(j, 2);

          point_markerlist[(localOffset + j)] = 0;
        }

        localOffset += localConstrainedPoints.rows();
      }
    }

    if(constrainedFaces.size() > 0)
    {
      unsigned int localOffset = numberOfFaces;
      unsigned int localOffsetPoints = numberOfVertices;
      for(unsigned int numConst = 0; numConst < constrainedFaces.size(); numConst++)
      {
        vector<VectorXi>& localConstrainedFaces = constrainedFaces[numConst];

        for(unsigned int numFac = 0; numFac < localConstrainedFaces.size(); numFac++)
        {
          VectorXi& localIds = localConstrainedFaces[numFac];
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
          face_markerlist[localOffset + numFac] = 0;
        }

        localOffset += localConstrainedFaces.size();
        localOffsetPoints += constrainedPoints[numConst].rows();
      }
    }

    if (!unidimensionalVertexMarkers.empty())
      memcpy(point_markerlist, unidimensionalVertexMarkers.data(), (numberOfVertices) * sizeof(int));
    if (!unidimensionalEdgeMarkers.empty())
      memcpy(edge_markerlist, unidimensionalEdgeMarkers.data(), (numberOfEdges) * sizeof(int));
    if (!unidimensionalFaceMarkers.empty())
      memcpy(face_markerlist, unidimensionalFaceMarkers.data(), (numberOfFaces) * sizeof(int));

    return Output::Success;
  }
  // ***************************************************************************
  Output::ExitCodes MeshCreator3DTetgen::CreateTetgenOutput(const Polyhedron& domain)
  {
    if (minimumNumberOfCells == 0 && maximumCellSize <= 0)
    {
      Output::PrintErrorMessage("Wrong initialization of the minimumNumberOfCells or minimumCellSize", false);
      return Output::GenericError;
    }

    if (inputMeshPointer == NULL)
    {
      Output::PrintErrorMessage("No Tetgen input in domain %d", false, domain.GlobalId());
      return Output::GenericError;
    }

    if (minimumNumberOfCells > 0 && domain.Measure() <= 0)
    {
      Output::PrintErrorMessage("%s: Wrong initialization of the domain %d, no measure computed", false, __func__, domain.GlobalId());
      return Output::GenericError;
    }

    const double& domainVolume = domain.Measure();
    double cellVolume = minimumNumberOfCells == 0 ? maximumCellSize : domainVolume / (double)minimumNumberOfCells;

    tetgenbehavior b;

    ostringstream options;
    options.precision(16);
    options<< tetgenOptions;
    options<< cellVolume;
    size_t sizeOptions = options.str().size();
    char* optionPointer = new char[sizeOptions + 1];
    options.str().copy(optionPointer, sizeOptions);
    optionPointer[sizeOptions] = '\0';

    b.parse_commandline(optionPointer);

    delete outputMeshPointer; outputMeshPointer = NULL;
    outputMeshPointer = new tetgenio();

    tetrahedralize(&b, inputMeshPointer, outputMeshPointer);

    delete[] optionPointer;

    return Output::Success;
  }
  // ***************************************************************************
  Output::ExitCodes MeshCreator3DTetgen::CreateMesh(const IGeometricObject& domain,
                                                    IMesh& mesh)
  {
    /// <ul>

    const Polyhedron& polyhedron = static_cast<const Polyhedron&>(domain);
    CreateUniDimMarkers();
    CreateMappedMarkers();
    CreateTetgenInput(polyhedron);
    CreateTetgenOutput(polyhedron);

    if (outputMeshPointer == NULL)
    {
      Output::PrintErrorMessage("No Tetgen ouput in domain %d", false, polyhedron.GlobalId());
      return Output::GenericError;
    }

    const tetgenio& tetgenMesh = *outputMeshPointer;

    /// <li>	Fill mesh structures
    unsigned int numberOfCellsMesh = tetgenMesh.numberoftetrahedra;
    unsigned int numberOfFacesMesh = tetgenMesh.numberoftrifaces;
    unsigned int numberOfEgdesMesh = tetgenMesh.numberofedges;
    unsigned int numberOfPointsMesh = tetgenMesh.numberofpoints;

    mesh.SetDimension(3);
    mesh.InitializeCells3D(numberOfCellsMesh);
    mesh.InitializeCells2D(numberOfFacesMesh);
    mesh.InitializeCells1D(numberOfEgdesMesh);
    mesh.InitializeCells0D(numberOfPointsMesh);
    vector<Cell3D*> cells(numberOfCellsMesh);
    vector<Cell2D*> faces(numberOfFacesMesh);
    vector<Cell1D*> edges(numberOfEgdesMesh);
    vector<Cell0D*> points(numberOfPointsMesh);

    /// <li> Set Points
    for (unsigned int p = 0; p < numberOfPointsMesh; p++)
    {
      points[p] = mesh.CreateCell0D();

      Cell0D* point = points[p];

      Vector3d point3d(tetgenMesh.pointlist[3 * p], tetgenMesh.pointlist[3 * p + 1], tetgenMesh.pointlist[3 * p + 2]);

      point->SetCoordinates(point3d);
      if(!unidimensionalVertexMarkers.empty())
      {
        unsigned int position = tetgenMesh.pointmarkerlist[p];

        if(position == 0)
          point->SetMarkers(vector<unsigned int>(markerDimension, 0));
        else if (position <= mappedMarkers.size())
          point->SetMarkers(mappedMarkers[position - 1]);
        else
          point->SetMarkers(vector<unsigned int>(markerDimension, position));
      }
      else
      {
        point->SetMarkers(vector<unsigned int>(markerDimension, 0));
      }
      point->InitializeNeighbourhood3D(4);
      point->InitializeNeighbourhood2D(4);
      point->InitializeNeighbourhood1D(4);


      mesh.AddCell0D(point);
    }

    /// <li> Set Edges
    SparseMatrix<int> connectivityPointsEdges(numberOfPointsMesh,numberOfPointsMesh);
    vector< Triplet<int> > triplets;
    triplets.reserve(numberOfEgdesMesh);

    for(unsigned int ed = 0; ed < numberOfEgdesMesh; ed++)
    {
      edges[ed] = mesh.CreateCell1D();

      Cell1D* edge = edges[ed];
      edge->InitializeVertices(2);
      if(!unidimensionalEdgeMarkers.empty())
      {
        unsigned int position = tetgenMesh.edgemarkerlist[ed];

        if(position == 0)
          edge->SetMarkers(vector<unsigned int>(markerDimension, 0));
        else if (position <= mappedMarkers.size())
          edge->SetMarkers(mappedMarkers[position - 1]);
        else
          edge->SetMarkers(vector<unsigned int>(markerDimension, position));
      }
      else
      {
        edge->SetMarkers(vector<unsigned int>(markerDimension, 0));
      }
      edge->InitializeNeighbourhood3D(10);
      edge->InitializeNeighbourhood2D(10);

      for(int i = 0; i < 2; i++)
      {
        Cell0D* point = points[tetgenMesh.edgelist[2 * ed + i]];
        edge->AddVertex(*point);
        point->AddNeighCell1D(*edge);
      }

      triplets.push_back(Triplet<int>(tetgenMesh.edgelist[2 * ed], tetgenMesh.edgelist[2 * ed + 1], ed + 1));
      mesh.AddCell1D(edge);
    }
    connectivityPointsEdges.setFromTriplets(triplets.begin(), triplets.end());

    /// <li> Set Faces
    vector<unsigned int> vertices(3);
    vector<unsigned int> edgeEndPoints(2);

    unsigned long estimatedValue = 2* numberOfPointsMesh * numberOfPointsMesh + 2*numberOfPointsMesh + 1;
    SparseMatrix<int> connectivityPointsFaces(numberOfPointsMesh,estimatedValue);
    vector< Triplet<int, unsigned long> > tripletsFaces;
    tripletsFaces.reserve(numberOfFacesMesh);
    for(unsigned int f = 0; f < numberOfFacesMesh; f++)
    {
      faces[f] = mesh.CreateCell2D();

      Cell2D* face = faces[f];

      face->AllocateNeighCell3D(2);
      face->InitializeEdges(3);
      face->AllocateVertices(3);
      face->SetType(Polygon::Triangle);

      if(!unidimensionalFaceMarkers.empty())
      {
        unsigned int position = tetgenMesh.trifacemarkerlist[f];

        if(position == 0)
          face->SetMarkers(vector<unsigned int>(markerDimension, 0));
        else if (position <= mappedMarkers.size())
          face->SetMarkers(mappedMarkers[position - 1]);
        else
          face->SetMarkers(vector<unsigned int>(markerDimension, position));
      }
      else
      {
        face->SetMarkers(vector<unsigned int>(markerDimension, 0));
      }

      for (unsigned int j = 0; j < 3; j++)
      {
        vertices[j] = tetgenMesh.trifacelist[3 * f + j];
        Cell0D* point = points[vertices[j]];
        face->InsertVertex(*point, j);
        point->AddNeighCell2D(*face);
      }

      face->ComputePlane();

      sort(vertices.begin(),vertices.end());
      unsigned long indexI = vertices[0];
      unsigned long indexJK = (vertices[1] + vertices[2]) * (vertices[1] + vertices[2] + 1) * 0.5 + vertices[2] + 1 ;
      tripletsFaces.push_back(Triplet<int, unsigned long>(indexI, indexJK, f + 1));

      for (unsigned int j = 0; j < 3; j++)
      {
        const Point* firstPoint = face->Vertex(j);
        const Point* secondPoint = face->Vertex((j + 1) % 3);
        edgeEndPoints[0] = firstPoint->Id();
        edgeEndPoints[1] = secondPoint->Id();

        int value = 0;
        value = connectivityPointsEdges.coeff(edgeEndPoints[0], edgeEndPoints[1]) - 1;
        if(value < 0)
          value = connectivityPointsEdges.coeff(edgeEndPoints[1], edgeEndPoints[0]) - 1;

        face->AddEdge(*edges[value]);

        if(tetgenMesh.trifacemarkerlist[f] != 0)
        {
          if(tetgenMesh.edgemarkerlist[value] == 0)
          {
            Cell1D* edge = edges[value];
            unsigned int position = tetgenMesh.trifacemarkerlist[f];

            if(position == 0)
              edge->SetMarkers(vector<unsigned int>(markerDimension, 0));
            else if (position <= mappedMarkers.size())
              edge->SetMarkers(mappedMarkers[position - 1]);
            else
              edge->SetMarkers(vector<unsigned int>(markerDimension, position));
          }
        }
        Cell1D& segment = *edges[value];
        segment.AddNeighCell2D(*face);
      }
      mesh.AddCell2D(face);
    }

    connectivityPointsFaces.setFromTriplets(tripletsFaces.begin(), tripletsFaces.end());
    connectivityPointsFaces.makeCompressed();

    /// <li> Set Cells
    vector<unsigned int> faceVertices(3);
    for (unsigned int c = 0; c < numberOfCellsMesh; c++)
    {
      cells[c] = mesh.CreateCell3D();

      Cell3D* cell = cells[c];

      cell->SetMarkers(vector<unsigned int>(markerDimension, 0));
      cell->AllocateVertices(4);
      cell->AllocateNeighCell3D(4);
      cell->AllocateFaces(4);
      cell->InitializeEdges(6);
      cell->SetType(Polyhedron::Tetrahedron);
      cell->AllocateNormalSign();

      for (int i = 0; i < 4; i++)
      {
        Cell0D* point = points[tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners * c + i]];
        cell->InsertVertex(*point, i);
        point->AddNeighCell3D(*cell);
      }
      Point barycenter;
      cell->ComputeBarycenter(barycenter);

      for (unsigned int j = 0; j < 4; j++)
      {
        for(unsigned int k = 0; k < 3; k++)
        {
          const Point* point =(cell->Vertex((j + k) % 4));
          faceVertices[k] = point->Id();
        }

        sort(faceVertices.begin(),faceVertices.end());
        unsigned int indexI = faceVertices[0];
        unsigned int indexJK = (faceVertices[1] + faceVertices[2]) * (faceVertices[1] + faceVertices[2] + 1) * 0.5 + faceVertices[2] + 1 ;

        int faceId = connectivityPointsFaces.coeff(indexI, indexJK) - 1;
        Cell2D& face = *faces[faceId];
        cell->InsertFace(*faces[faceId], j);

        if(face.PointInPlane(barycenter) == Polygon::Negative)
        {
          face.InsertNeighCell3D(*cell, 1);
          cell->SetNormalSign(true, j);
        }
        else
        {
          face.InsertNeighCell3D(*cell, 0);
          cell->SetNormalSign(false, j);
        }

        if (j < 3)
        {
          for(unsigned int k = 0; k < 3; k++)
          {
            bool flag = true;
            unsigned int pos = 0;

            while (pos < cell->NumberOfEdges() && flag)
            {
              const Polygon* poly = cell->Face(j);
              flag = poly->Edge(k) != cell->Edge(pos++);
            }

            if(flag)
            {
              cell->AddEdge(*(cell->Face(j)->Edge(k)));
              const Cell1D* edge = static_cast<const Cell1D*>(cell->Face(j)->Edge(k));
              edges[edge->Id()]->AddNeighCell3D(*cell);
            }
          }
        }
      }

      mesh.AddCell3D(cell);
    }
    if ( outputMeshPointer->neighborlist != NULL)
    {
      for (unsigned int c = 0; c < numberOfCellsMesh; c++)
      {
        Cell3D* cell = cells[c];

        for (int i = 0; i < 4; i++)
        {
          if (outputMeshPointer->neighborlist[tetgenMesh.numberofcorners * c + i] > -1)
          {
            Cell3D* cellNeigh = cells[tetgenMesh.neighborlist[4 * c + i]];
            cell->InsertNeighCell3D(*cellNeigh, (i+1)%4);
          }
        }
      }
    }
    return Output::Success;

  }
  // ***************************************************************************
  Output::ExitCodes MeshCreator3DTetgen::ExportTetgenMesh(const string& nameFolder, const string& nameFile) const
  {
    if (outputMeshPointer == NULL)
      return Output::GenericError;

    ostringstream nameFolderStream, nameFileStream;

    nameFolderStream<< nameFolder<< "/";
    nameFolderStream<< "Tetgen/";

    Output::CreateFolder(nameFolderStream.str());

    nameFileStream<< nameFolderStream.str()<< nameFile;

    Output::CreateFolder(nameFolderStream.str());

    outputMeshPointer->firstnumber = 0;
    outputMeshPointer->save_nodes((char*)nameFileStream.str().c_str());
    outputMeshPointer->save_elements((char*)nameFileStream.str().c_str());
    outputMeshPointer->save_faces((char*)nameFileStream.str().c_str());
    outputMeshPointer->save_edges((char*)nameFileStream.str().c_str());

    return Output::Success;
  }
  // ***************************************************************************
  /*Output::ExitCodes TetgenVemInterface::TetgenToVemMesh(const tetgenio& tetgenMesh, VemMesh& mesh) const
    {
    /// Add Points
    mesh.InitializePoints(tetgenMesh.numberofpoints);
    for(unsigned int i = 0; i < tetgenMesh.numberofpoints; i++)
    {
    VemPoint& point = *(mesh.CreatePoint());
    point.SetCoordinates(tetgenMesh.pointlist[3*i],tetgenMesh.pointlist[3*i+1],tetgenMesh.pointlist[3*i+2]);
    point.SetMarker(tetgenMesh.pointmarkerlist[i]);
    mesh.AddPoint(&point);
    }

    /// Add Edges
    mesh.InitializeEdges(tetgenMesh.numberofedges);
    map< vector<unsigned int> , unsigned int > endpointsToEdge;
    vector<unsigned int> endpoints (2);
    for(unsigned int i = 0; i < tetgenMesh.numberofedges; i++)
    {
    VemEdge& edge = *(mesh.CreateEdge());
    endpoints[0] = tetgenMesh.edgelist[2*i];
    endpoints[1] = tetgenMesh.edgelist[2*i+1];
    edge.AddPoint(mesh.Point(endpoints[0]));
    edge.AddPoint(mesh.Point(endpoints[1]));
    sort(endpoints.begin(),endpoints.end());
    endpointsToEdge.insert(std::pair< vector<unsigned int>, unsigned int>(endpoints, edge.Id()));
    edge.SetMarker(tetgenMesh.edgemarkerlist[i]);
    mesh.AddEdge(&edge);
    }

    /// Add Faces
    map< vector<unsigned int> , unsigned int > verticesToFace;
    vector<unsigned int> vertices (3);
    mesh.InitializeFaces(tetgenMesh.numberoftrifaces);
    for(unsigned int i = 0; i < tetgenMesh.numberoftrifaces; i++)
    {
    VemFace& face = *(mesh.CreateFace());
    face.InitializePoints(3);
    /// Points
    for(unsigned int j = 0; j < 3; j++)
    {
    vertices[j] = tetgenMesh.trifacelist[3*i+j];
    face.AddPoint(mesh.Point(vertices[j]));
    }
    sort(vertices.begin(),vertices.end());
    verticesToFace.insert(std::pair< vector<unsigned int>, unsigned int >(vertices, face.Id()));
    /// Edges
    vector<unsigned int> edgeEndPoints (2);
    face.InitializeEdges(3);
    for(unsigned int j = 0; j < 3; j++)
    {
    edgeEndPoints[0] = face.Point(j)->Id();
    edgeEndPoints[1] = face.Point((j+1)%3)->Id();
    sort(edgeEndPoints.begin(),edgeEndPoints.end());
    if(endpointsToEdge.find(edgeEndPoints)==endpointsToEdge.end())
    {
    Output::PrintErrorMessage("%s: error retrieving an edge when adding edges to faces", true, __func__);
    exit(-1);
    }
    face.AddEdge(mesh.Edge(endpointsToEdge.find(edgeEndPoints)->second));
    }
    face.SetMarker(tetgenMesh.trifacemarkerlist[i]);
    face.InitializeCells(2);
    mesh.AddFace(&face);
    }

    /// Add Cells
    mesh.InitializeCells(tetgenMesh.numberoftetrahedra);
    map<vector<unsigned int> , unsigned int>::iterator verticesToFace_iterator;
    for(unsigned int i = 0; i < tetgenMesh.numberoftetrahedra; i++)
    {
    VemCell& cell = *(mesh.CreateCell3D());
    /// Points
    cell.InitializePoints(4);
    bool flag = false;
    cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i]));
    cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+1]));
    cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+2]));
    cell.AddPoint(mesh.Point(tetgenMesh.tetrahedronlist[tetgenMesh.numberofcorners*i+3]));
    /// Faces and Edges
    cell.InitializeFaces(4);
    cell.InitializeEdges(6);
    vector<unsigned int> faceVertices (3);
    for(unsigned int j = 0; j < 4; j++)
    {
    for(unsigned int k = 0; k < 3; k++)
    faceVertices[k] = cell.Point((j+k)%4)->Id();
    sort(faceVertices.begin(),faceVertices.end());
    verticesToFace_iterator = verticesToFace.find(faceVertices);
    if(verticesToFace_iterator != verticesToFace.end())
    {
    unsigned int faceId = verticesToFace_iterator->second;
    cell.AddFace(mesh.Face(faceId));
    mesh.AddFaceNeighbourCell( faceId, &cell );
    }
    else
    {
    Output::PrintErrorMessage("%s: error retrieving faces when adding faces to cells", true, __func__);
    exit(-1);
    }
    if(j < 3)
    {
    for(unsigned int k = 0; k < 3; k++)
    {
    bool flag = true;
    unsigned int pos = 0;
    while(pos < cell.NumberOfEdges() && flag)
    flag = cell.Face(j)->Edge(k) != cell.Edge(pos++);
    if(flag)
    cell.AddEdge(cell.Face(j)->Edge(k));
    }
    }
    }
    if(cell.NumberOfEdges() != 6)
    {
    Output::PrintErrorMessage("%s: error retrieving edges when adding edges to cells", true, __func__);
    exit(-1);
    }
    mesh.AddCell(&cell);
    }
    return Output::Success;
    }*/
  // ***************************************************************************
}
