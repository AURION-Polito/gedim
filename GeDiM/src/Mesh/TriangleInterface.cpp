#include <math.h>
#include <iomanip>
#include <sys/stat.h>

#include "TriangleInterface.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  TriangleInterface::TriangleInterface()
  {
  }
  TriangleInterface::~TriangleInterface()
  {
  }
  // ***************************************************************************
  void TriangleInterface::CreateTriangleInput(const MatrixXd& polygonVertices,
                                              struct triangulateio& triangleInput,
                                              const MatrixXd& constrainedPoints,
                                              const MatrixXi& constrainedSegments) const
  {
    const unsigned int& numVertices = polygonVertices.cols();
    unsigned int numberOfConstrainedPoints = constrainedPoints.cols();
    const unsigned int& numberOfConstrainedSegments = constrainedSegments.cols();
    const unsigned int& numberOfEdges = numVertices;

    triangleInput.pointlist = new double[2 * (numVertices + numberOfConstrainedPoints)];
    triangleInput.pointattributelist = NULL;
    triangleInput.pointmarkerlist = new int[numVertices + numberOfConstrainedPoints];
    triangleInput.numberofpoints = numVertices + numberOfConstrainedPoints;
    triangleInput.numberofpointattributes = 0;
    triangleInput.numberofsegments = numberOfEdges + numberOfConstrainedSegments;
    triangleInput.trianglelist = NULL;
    triangleInput.triangleattributelist = NULL;
    triangleInput.trianglearealist = NULL;
    triangleInput.neighborlist = NULL;
    triangleInput.numberoftriangles = 0;
    triangleInput.numberofcorners = 0;
    triangleInput.numberoftriangleattributes = 0;
    triangleInput.segmentlist = new int[2 * (numberOfEdges + numberOfConstrainedSegments)];
    triangleInput.segmentmarkerlist = new int[numberOfEdges + numberOfConstrainedSegments];
    triangleInput.holelist = NULL;
    triangleInput.numberofholes = 0;
    triangleInput.regionlist = NULL;
    triangleInput.numberofregions = 0;
    triangleInput.edgelist = NULL;
    triangleInput.edgemarkerlist = NULL;
    triangleInput.normlist = NULL;
    triangleInput.numberofedges = 0;

    double* point_list = triangleInput.pointlist;
    int* point_markerlist = triangleInput.pointmarkerlist;
    int* segment_list = triangleInput.segmentlist;
    int* segment_markerlist = triangleInput.segmentmarkerlist;

    for (unsigned int j = 0; j < numVertices; j++)
    {
      point_list[2 * j] = polygonVertices.col(j).x();
      point_list[2 * j + 1] = polygonVertices.col(j).y();

      point_markerlist[j] = j + 1;
    }

    for (unsigned int e = 0; e < numberOfEdges; e++)
    {
      segment_list[2 * e] = e;
      segment_list[2 * e + 1] = (e + 1) % numberOfEdges;

      segment_markerlist[e]= numVertices + numberOfConstrainedPoints + e + 1;
    }

    if(numberOfConstrainedPoints > 0)
    {
      for (unsigned int j = 0; j < numberOfConstrainedPoints; j++)
      {
        point_list[2 * (numVertices + j)] = constrainedPoints.col(j).x();
        point_list[2 * (numVertices + j) + 1] = constrainedPoints.col(j).y();

        point_markerlist[numVertices + j] = numVertices + j + 1;
      }

      for (unsigned int e = 0; e < numberOfConstrainedSegments; e++)
      {
        segment_list[2 * (numberOfEdges + e)] = constrainedSegments(0, e);
        segment_list[2 * (numberOfEdges + e) + 1] = constrainedSegments(1, e);

        segment_markerlist[numberOfEdges + e] = numVertices + numberOfConstrainedPoints + numberOfEdges + e + 1;
      }
    }
  }
  // ***************************************************************************
  void TriangleInterface::CreateTriangleOutput(const double& maxTriangleArea,
                                               struct triangulateio& triangleInput,
                                               struct triangulateio& triangleOutput,
                                               const string& triangleOptions) const
  {
    ostringstream options;
    options.precision(std::numeric_limits<double>::digits10);
    options<< triangleOptions;
    options<< maxTriangleArea;

    size_t sizeOptions = options.str().size();
    char* optionPointer = new char[sizeOptions + 1];
    options.str().copy(optionPointer, sizeOptions);
    optionPointer[sizeOptions] = '\0';

    triangulate(optionPointer,
                &triangleInput,
                &triangleOutput,
                (struct triangulateio*) NULL);

    delete[] optionPointer;
  }
  // ***************************************************************************
  void TriangleInterface::DeleteTriangleStructure(struct triangulateio& triangleInput,
                                                  struct triangulateio& triangleOutput) const
  {
    delete[] triangleInput.pointlist;
    triangleInput.pointlist = NULL;
    delete[] triangleInput.pointmarkerlist;
    triangleInput.pointmarkerlist = NULL;
    delete[] triangleInput.segmentlist;
    triangleInput.segmentlist = NULL;
    delete[] triangleInput.segmentmarkerlist;
    triangleInput.segmentmarkerlist = NULL;

    free(triangleOutput.pointlist);
    free(triangleOutput.pointattributelist);
    free(triangleOutput.pointmarkerlist);
    free(triangleOutput.trianglelist);
    free(triangleOutput.triangleattributelist);
    free(triangleOutput.trianglearealist);
    free(triangleOutput.neighborlist);
    free(triangleOutput.segmentlist);
    free(triangleOutput.segmentmarkerlist);
    free(triangleOutput.edgelist);
    free(triangleOutput.edgemarkerlist);
  }
  // ***************************************************************************
  void TriangleInterface::CreateMesh(const Eigen::MatrixXd& polygonVertices,
                                     const double& maxTriangleArea,
                                     IMeshDAO& mesh,
                                     const string& triangleOptions) const
  {
    struct triangulateio* triangleInput = new triangulateio();
    struct triangulateio* triangleOutput = new triangulateio();

    CreateTriangleInput(polygonVertices,
                        *triangleInput);
    CreateTriangleOutput(maxTriangleArea,
                         *triangleInput,
                         *triangleOutput,
                         triangleOptions);

    /// <li>	Fill mesh structures
    unsigned int numberOfCell2Ds = triangleOutput->numberoftriangles;
    unsigned int numberOfCell1Ds = triangleOutput->numberofedges;
    unsigned int numberOfCell0Ds = triangleOutput->numberofpoints;

    mesh.InitializeDimension(2);
    mesh.Cell2DsInitialize(numberOfCell2Ds);
    mesh.Cell1DsInitialize(numberOfCell1Ds);
    mesh.Cell0DsInitialize(numberOfCell0Ds);

    /// <li> Set Cell0Ds
    for (unsigned int p = 0; p < numberOfCell0Ds; p++)
    {
      mesh.Cell0DInsertCoordinates(p,
                                   Vector3d(triangleOutput->pointlist[2 * p],
                                   triangleOutput->pointlist[2 * p + 1],
          0.0));
      mesh.Cell0DSetState(p, true);
      mesh.Cell0DSetMarker(p,
                           triangleOutput->pointmarkerlist[p]);
    }

    /// <li> Set Cell1Ds
    struct Edge
    {
        unsigned int Cell1DIndex;
        int Cell0DEnd;
    };

    std::vector<std::list<Edge>> cell0DsCell1Ds(mesh.Cell0DTotalNumber());

    for (unsigned int ed = 0; ed < numberOfCell1Ds; ed++)
    {
      mesh.Cell1DInsertExtremes(ed,
                                triangleOutput->edgelist[2 * ed],
          triangleOutput->edgelist[2 * ed + 1]);
      mesh.Cell1DSetState(ed, true);
      mesh.Cell1DSetMarker(ed,
                           triangleOutput->edgemarkerlist[ed]);

      cell0DsCell1Ds[triangleOutput->edgelist[2 * ed]].push_back({ ed, triangleOutput->edgelist[2 * ed + 1] });
    }

    /// <li> Set Cell2Ds
    mesh.Cell2DsInitializeVertices(3);
    mesh.Cell2DsInitializeEdges(3);
    for (unsigned int c = 0; c < numberOfCell2Ds; c++)
    {
      vector<unsigned int> cell2DVertices(3);
      vector<unsigned int> cell2DEdges(3);

      cell2DVertices[0] = triangleOutput->trianglelist[3 * c];
      cell2DVertices[1] = triangleOutput->trianglelist[3 * c + 1];
      cell2DVertices[2] = triangleOutput->trianglelist[3 * c + 2];

      for (unsigned int e = 0; e < 3; e++)
      {
        const int vertexId = cell2DVertices[e];
        const int nextVertexId = cell2DVertices[(e + 1) % 3];
        unsigned int edgeId = numberOfCell1Ds;

        if (edgeId == numberOfCell1Ds)
        {
          const std::list<Edge>& originCell1Ds = cell0DsCell1Ds.at(vertexId);
          std::list<Edge>::const_iterator findOriginEdge = std::find_if(originCell1Ds.begin(),
                                                                        originCell1Ds.end(),
                                                                        [&] (const Edge& edge)
          { return edge.Cell0DEnd == nextVertexId; });

          if (findOriginEdge != originCell1Ds.end())
            edgeId = findOriginEdge->Cell1DIndex;
        }

        if (edgeId == numberOfCell1Ds)
        {
          const std::list<Edge>& endCell1Ds = cell0DsCell1Ds.at(nextVertexId);
          std::list<Edge>::const_iterator findEndEdge = std::find_if(endCell1Ds.begin(),
                                                                     endCell1Ds.end(),
                                                                     [&] (const Edge& edge)
          { return edge.Cell0DEnd == vertexId; });

          if (findEndEdge != endCell1Ds.end())
            edgeId = findEndEdge->Cell1DIndex;
        }

        Gedim::Output::Assert(edgeId != numberOfCell1Ds);
        cell2DEdges[e] = edgeId;
      }

      mesh.Cell2DInsertVertices(c,
                                cell2DVertices);
      mesh.Cell2DInsertEdges(c,
                             cell2DEdges);
      mesh.Cell2DSetState(c, true);
      mesh.Cell2DSetMarker(c,
                           0);
    }

    DeleteTriangleStructure(*triangleInput,
                            *triangleOutput);
    delete triangleInput;
    delete triangleOutput;
  }
  // ***************************************************************************
  void TriangleInterface::ExportTriangleMesh(const struct triangulateio& triangleInput,
                                             const struct triangulateio& triangleOutput,
                                             const string& nameFolder,
                                             const string& nameFile) const
  {
    ostringstream nameFolderStream, nameFileStream;

    nameFolderStream<< nameFolder<< "/";
    nameFolderStream<< "Triangle/";

    Output::CreateFolder(nameFolderStream.str());

    nameFileStream<< nameFolderStream.str()<< nameFile;

    ofstream file;
    char basefilename[50];
    char filename[50];
    strcpy(basefilename, nameFileStream.str().c_str());

    file<< setprecision(16);

    strcpy (filename, basefilename);
    strcat (filename, ".poly");
    file.open(filename);
    file<< triangleInput.numberofpoints<< " "<< "2"<< " "<< "0"<< " "<< "1"<< endl;
    for(int i = 0; i < triangleInput.numberofpoints; i++)
      file<< i + 1<< " "
          << triangleInput.pointlist[2 * i]<< " "
          << triangleInput.pointlist[2 * i + 1]<< " "
          << triangleInput.pointmarkerlist[i]<< endl;

    file<< triangleInput.numberofsegments<< " "<< " 0 "<< endl;

    for(int i = 0; i < triangleInput.numberofsegments; i++)
      file<< i + 1<< " "
          << triangleInput.segmentlist[2 * i] + 1<< " "
          << triangleInput.segmentlist[2 * i + 1] + 1<< " "
          << triangleInput.segmentmarkerlist[i]<< endl;

    file<< "0";
    file.close();

    strcpy (filename, basefilename);
    strcat (filename, ".1.node");
    file.open(filename);
    file<< triangleOutput.numberofpoints<< " "<< "2"<< " "<< "0"<< " "<< "1"<< endl;
    for(int i = 0; i < triangleOutput.numberofpoints; i++)
      file<< i + 1<< " "<< triangleOutput.pointlist[2 * i]<< " "
          << triangleOutput.pointlist[2 * i + 1]<< " "
          << triangleOutput.pointmarkerlist[i]<< endl;
    file.close();

    strcpy (filename, basefilename);
    strcat (filename, ".1.neigh");
    file.open(filename);
    file<< triangleOutput.numberoftriangles<< " "<< "3"<< endl;
    for(int i = 0; i < triangleOutput.numberoftriangles; i++)
      file<< i + 1<< " "
          << ((triangleOutput.neighborlist[3 * i] != -1) ? (triangleOutput.neighborlist[3 * i] + 1) : -1)<< " "
        << ((triangleOutput.neighborlist[3 * i + 1] != -1) ? (triangleOutput.neighborlist[3 * i + 1] + 1) : -1)<< " "
      << ((triangleOutput.neighborlist[3 * i + 2] != -1) ? (triangleOutput.neighborlist[3 * i + 2] + 1) : -1)<< endl;
    file.close();

    strcpy (filename,basefilename);
    strcat (filename,".1.edge");
    file.open(filename);
    file<< triangleOutput.numberofedges<< " "<< "1"<< endl;
    for(int i = 0; i < triangleOutput.numberofedges; i++)
      file<< i + 1<< " "<< triangleOutput.edgelist[2 * i] + 1
          << " "<< triangleOutput.edgelist[2 * i + 1] + 1
          << " "<< triangleOutput.edgemarkerlist[i]<< endl;
    file.close();

    strcpy (filename, basefilename);
    strcat (filename,".1.ele");
    file.open(filename);
    file<< triangleOutput.numberoftriangles<< " "<< "3"<< " "<< "0"<< endl;
    for(int i = 0; i < triangleOutput.numberoftriangles; i++)
      file<< i + 1<< " "
          << triangleOutput.trianglelist[3 * i] + 1<< " "
          << triangleOutput.trianglelist[3 * i + 1] + 1<< " "
          << triangleOutput.trianglelist[3 * i + 2] + 1<< endl;
    file.close();

    strcpy (filename, basefilename);
    strcat (filename, ".1.poly");
    file.open(filename);
    file<< " 0 "<< "2 "<< "0 "<< "1"<< endl;
    file<< triangleOutput.numberofsegments<< " 1"<< endl;
    for(int i = 0; i < triangleOutput.numberofsegments; i++)
    {
      file<< i + 1<< " "
          << triangleOutput.segmentlist[2 * i] + 1<< " "
          << triangleOutput.segmentlist[2 * i + 1] + 1<< " "
          << triangleOutput.segmentmarkerlist[i]<< endl;
    }

    file<< "0";
    file.close();
  }
  // ***************************************************************************
}
