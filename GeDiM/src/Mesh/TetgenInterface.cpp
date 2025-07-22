// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "TetgenInterface.hpp"
#include "CommonUtilities.hpp"
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
void TetgenInterface::CreateDelaunay(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, IMeshDAO &mesh) const
{
    Gedim::Utilities::Unused(points);
    Gedim::Utilities::Unused(points_marker);
    Gedim::Utilities::Unused(mesh);

#if ENABLE_TETGEN == 1
    tetgenio *tetgenInput = new tetgenio();
    tetgenio *tetgenOutput = new tetgenio();

    CreateDelaunayInput(points, points_marker, *tetgenInput);
    CreateTetgenOutput(*tetgenInput, *tetgenOutput, "Qfezn");

    ConvertTetgenOutputToMeshDAO(*tetgenOutput, mesh);

    DeleteTetgenStructure(*tetgenInput, *tetgenOutput);
    delete tetgenInput;
    delete tetgenOutput;
#endif
}
// ***************************************************************************
void TetgenInterface::CreateMesh(const Eigen::MatrixXd &polyhedronVertices,
                                 const Eigen::MatrixXi &polyhedronEdges,
                                 const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                 const double &maxTetrahedronVolume,
                                 IMeshDAO &mesh,
                                 const std::string &tetgenOptions) const
{
    Gedim::Utilities::Unused(polyhedronVertices);
    Gedim::Utilities::Unused(polyhedronEdges);
    Gedim::Utilities::Unused(polyhedronFaces);
    Gedim::Utilities::Unused(maxTetrahedronVolume);
    Gedim::Utilities::Unused(mesh);
    Gedim::Utilities::Unused(tetgenOptions);

#if ENABLE_TETGEN == 1
    tetgenio *tetgenInput = new tetgenio();
    tetgenio *tetgenOutput = new tetgenio();

    CreateTetgenInput(polyhedronVertices, polyhedronEdges, polyhedronFaces, *tetgenInput);
    CreateTetgenOutput(maxTetrahedronVolume, *tetgenInput, *tetgenOutput, tetgenOptions);

    ConvertTetgenOutputToMeshDAO(*tetgenOutput, mesh);

    DeleteTetgenStructure(*tetgenInput, *tetgenOutput);
    delete tetgenInput;
    delete tetgenOutput;
#endif
}
// ***************************************************************************
void TetgenInterface::CreateMesh(const Eigen::MatrixXd &points,
                                 const std::vector<std::vector<unsigned int>> &facets,
                                 const double &maxTetrahedronVolume,
                                 IMeshDAO &mesh,
                                 const std::string &tetgenOptions) const
{
    Gedim::Utilities::Unused(points);
    Gedim::Utilities::Unused(facets);
    Gedim::Utilities::Unused(maxTetrahedronVolume);
    Gedim::Utilities::Unused(mesh);
    Gedim::Utilities::Unused(tetgenOptions);

#if ENABLE_TETGEN == 1
    tetgenio *tetgenInput = new tetgenio();
    tetgenio *tetgenOutput = new tetgenio();

    CreateTetgenInput(points, facets, {}, *tetgenInput);
    CreateTetgenOutput(maxTetrahedronVolume, *tetgenInput, *tetgenOutput, tetgenOptions);

    ConvertTetgenOutputToMeshDAO(*tetgenOutput, mesh);

    DeleteTetgenStructure(*tetgenInput, *tetgenOutput);
    delete tetgenInput;
    delete tetgenOutput;
#endif
}
// ***************************************************************************
void TetgenInterface::CreateMesh(const Eigen::MatrixXd &points,
                const std::vector<std::vector<unsigned int>> &facets,
                const std::vector<TetgenInterface::Region> &regions,
                IMeshDAO &mesh,
                const std::string &tetgenOptions) const
{
  Gedim::Utilities::Unused(points);
  Gedim::Utilities::Unused(facets);
  Gedim::Utilities::Unused(regions);
  Gedim::Utilities::Unused(mesh);
  Gedim::Utilities::Unused(tetgenOptions);

#if ENABLE_TETGEN == 1
  tetgenio *tetgenInput = new tetgenio();
  tetgenio *tetgenOutput = new tetgenio();

  CreateTetgenInput(points, facets, regions, *tetgenInput);
  CreateTetgenOutput(0.0, *tetgenInput, *tetgenOutput, tetgenOptions);

  ConvertTetgenOutputToMeshDAO(*tetgenOutput, mesh);

  DeleteTetgenStructure(*tetgenInput, *tetgenOutput);
  delete tetgenInput;
  delete tetgenOutput;
#endif
}
// ***************************************************************************
#if ENABLE_TETGEN == 1
void TetgenInterface::DeleteTetgenStructure(tetgenio &tetgenInput, tetgenio &) const
{
    delete[] tetgenInput.pointlist;
    tetgenInput.pointlist = NULL;
    delete[] tetgenInput.pointmarkerlist;
    tetgenInput.pointmarkerlist = NULL;

    delete[] tetgenInput.edgelist;
    tetgenInput.edgelist = NULL;
    delete[] tetgenInput.edgemarkerlist;
    tetgenInput.edgemarkerlist = NULL;

    for (int f = 0; f < tetgenInput.numberoffacets; f++)
    {
        tetgenio::facet *tetgenFace = &tetgenInput.facetlist[f];

        tetgenio::polygon *tetgenPolygon = &tetgenFace->polygonlist[0];
        delete[] tetgenPolygon->vertexlist;
        tetgenPolygon->vertexlist = NULL;

        delete[] tetgenFace->polygonlist;
        tetgenFace->polygonlist = NULL;
    }

    delete[] tetgenInput.facetlist;
    tetgenInput.facetlist = NULL;
}
// ***************************************************************************
void TetgenInterface::CreateTetgenInput(const Eigen::MatrixXd &polyhedronVertices,
                                        const Eigen::MatrixXi &polyhedronEdges,
                                        const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                        tetgenio &tetgenInput,
                                        const MatrixXd &constrainedPoints,
                                        const std::vector<Eigen::VectorXi> &constrainedFaces) const
{
    const unsigned int &numberOfVertices = polyhedronVertices.cols();
    const unsigned int numberOfConstrainedPoints = constrainedPoints.cols();

    const unsigned int &numberOfEdges = polyhedronEdges.cols();
    const unsigned int &numberOfFaces = polyhedronFaces.size();
    const unsigned int numberOfConstrainedFaces = constrainedFaces.size();

    Output::Assert(numberOfVertices > 0 && numberOfEdges > 0 && numberOfFaces > 0);

    tetgenInput.firstnumber = 0;
    tetgenInput.numberofpoints = numberOfVertices + numberOfConstrainedPoints;
    tetgenInput.pointlist = new REAL[(numberOfVertices + numberOfConstrainedPoints) * 3];
    tetgenInput.pointmarkerlist = new int[numberOfVertices + numberOfConstrainedPoints];

    tetgenInput.numberofedges = numberOfEdges;
    tetgenInput.edgelist = new int[(numberOfEdges) * 2];
    tetgenInput.edgemarkerlist = new int[(numberOfEdges)];

    tetgenInput.numberoffacets = numberOfFaces + numberOfConstrainedFaces;
    tetgenInput.facetlist = new tetgenio::facet[numberOfFaces + numberOfConstrainedFaces];
    tetgenInput.facetmarkerlist = new int[numberOfFaces + numberOfConstrainedFaces];

    double *point_list = tetgenInput.pointlist;
    int *point_markerlist = tetgenInput.pointmarkerlist;

    int *edge_list = tetgenInput.edgelist;
    int *edge_markerlist = tetgenInput.edgemarkerlist;

    tetgenio::facet *face_list = tetgenInput.facetlist;
    int *face_markerlist = tetgenInput.facetmarkerlist;

    for (unsigned int v = 0; v < numberOfVertices; v++)
    {
        const Eigen::Vector3d &point = polyhedronVertices.col(v);
        point_list[3 * v] = point(0);
        point_list[3 * v + 1] = point(1);
        point_list[3 * v + 2] = point(2);

        point_markerlist[v] = v + 1;
    }

    if (numberOfConstrainedPoints > 0)
    {
        unsigned int localOffset = numberOfVertices;
        for (unsigned int j = 0; j < constrainedPoints.cols(); j++)
        {
            point_list[3 * (localOffset + j)] = constrainedPoints(0, j);
            point_list[3 * (localOffset + j) + 1] = constrainedPoints(1, j);
            point_list[3 * (localOffset + j) + 2] = constrainedPoints(2, j);

            point_markerlist[(localOffset + j)] = numberOfVertices + j + 1;
        }
    }

    for (unsigned int e = 0; e < numberOfEdges; e++)
    {
        const unsigned int &vId1 = polyhedronEdges(0, e);
        const unsigned int &vId2 = polyhedronEdges(1, e);

        Output::Assert(vId1 < numberOfVertices && vId2 < numberOfVertices);

        edge_list[2 * e] = vId1;
        edge_list[2 * e + 1] = vId2;

        edge_markerlist[e] = numberOfVertices + numberOfConstrainedPoints + e + 1;
    }

    for (unsigned int f = 0; f < numberOfFaces; f++)
    {
        tetgenio::facet *tetgenFace = &face_list[f];

        tetgenFace->numberofpolygons = 1;
        tetgenFace->polygonlist = new tetgenio::polygon[1];
        tetgenFace->numberofholes = 0;
        tetgenFace->holelist = NULL;

        tetgenio::polygon *tetgenPolygon = &tetgenFace->polygonlist[0];

        const MatrixXi &face = polyhedronFaces[f];
        const size_t numberFacePoints = face.cols();
        tetgenPolygon->numberofvertices = numberFacePoints;
        tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

        for (unsigned int v = 0; v < numberFacePoints; v++)
        {
            const unsigned int &vId = face(0, v);
            Output::Assert(vId < numberOfVertices);

            tetgenPolygon->vertexlist[v] = vId;
        }

        face_markerlist[f] = numberOfVertices + numberOfConstrainedPoints + numberOfEdges + f + 1;
    }

    if (constrainedFaces.size() > 0)
    {
        unsigned int localOffset = numberOfFaces;
        unsigned int localOffsetPoints = numberOfVertices;
        for (unsigned int numFac = 0; numFac < constrainedFaces.size(); numFac++)
        {
            const VectorXi &localIds = constrainedFaces[numFac];
            tetgenio::facet *tetgenFace = &face_list[localOffset + numFac];

            tetgenFace->numberofpolygons = 1;
            tetgenFace->polygonlist = new tetgenio::polygon[1];
            tetgenFace->numberofholes = 0;
            tetgenFace->holelist = NULL;

            tetgenio::polygon *tetgenPolygon = &tetgenFace->polygonlist[0];

            const size_t numberFacePoints = localIds.size();
            tetgenPolygon->numberofvertices = numberFacePoints;
            tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

            for (unsigned int v = 0; v < numberFacePoints; v++)
            {
                unsigned int vId = localIds(v) + localOffsetPoints;
                tetgenPolygon->vertexlist[v] = vId;
            }
            face_markerlist[localOffset + numFac] =
                numberOfVertices + numberOfConstrainedPoints + numberOfEdges + numberOfFaces + numFac + 1;
        }
    }
}
// ***************************************************************************
void TetgenInterface::CreateTetgenInput(const Eigen::MatrixXd &points,
                                        const std::vector<std::vector<unsigned int>> &facets,
                                        const std::vector<TetgenInterface::Region> &regions,
                                        tetgenio &tetgenInput) const
{
    const unsigned int num_points = points.cols();

    Output::Assert(num_points > 0);

    tetgenInput.firstnumber = 0;
    tetgenInput.numberofpoints = num_points;
    tetgenInput.pointlist = new REAL[num_points * 3];
    tetgenInput.pointmarkerlist = new int[num_points];

    const unsigned int num_facets = facets.size();

    tetgenInput.numberoffacets = num_facets;
    tetgenInput.facetlist = new tetgenio::facet[num_facets];
    tetgenInput.facetmarkerlist = new int[num_facets];

    double *point_list = tetgenInput.pointlist;
    int *point_markerlist = tetgenInput.pointmarkerlist;

    tetgenio::facet *face_list = tetgenInput.facetlist;
    int *face_markerlist = tetgenInput.facetmarkerlist;

    for (unsigned int v = 0; v < num_points; v++)
    {
        const Eigen::Vector3d &point = points.col(v);
        point_list[3 * v] = point(0);
        point_list[3 * v + 1] = point(1);
        point_list[3 * v + 2] = point(2);

        point_markerlist[v] = v + 1;
    }

    for (unsigned int f = 0; f < num_facets; f++)
    {
        tetgenio::facet *tetgenFace = &face_list[f];

        tetgenFace->numberofpolygons = 1;
        tetgenFace->polygonlist = new tetgenio::polygon[1];
        tetgenFace->numberofholes = 0;
        tetgenFace->holelist = NULL;

        tetgenio::polygon *tetgenPolygon = &tetgenFace->polygonlist[0];

        const auto &face = facets[f];
        const size_t numberFacePoints = face.size();
        tetgenPolygon->numberofvertices = numberFacePoints;
        tetgenPolygon->vertexlist = new int[tetgenPolygon->numberofvertices];

        for (unsigned int v = 0; v < numberFacePoints; v++)
        {
            const unsigned int &vId = face[v];
            Output::Assert(vId < num_points);

            tetgenPolygon->vertexlist[v] = vId;
        }

        face_markerlist[f] = num_points + f + 1;
    }

    if (regions.size() > 0)
    {
      const unsigned int num_regions = regions.size();

      tetgenInput.numberofregions = num_regions;
      tetgenInput.regionlist = new REAL[num_regions * 5];

      for (unsigned int r = 0; r < num_regions; ++r)
      {
        const auto& region = regions[r];

        tetgenInput.regionlist[0 + 5 * r] = region.centroid.x();
        tetgenInput.regionlist[1 + 5 * r] = region.centroid.y();
        tetgenInput.regionlist[2 + 5 * r] = region.centroid.z();
        tetgenInput.regionlist[3 + 5 * r] = region.id;
        tetgenInput.regionlist[4 + 5 * r] = region.max_volume;
      }
    }
}
// ***************************************************************************
void TetgenInterface::CreateDelaunayInput(const Eigen::MatrixXd &points, const std::vector<unsigned int> &points_marker, tetgenio &tetgenInput) const
{
    const unsigned int &number_of_points = points.cols();

    Output::Assert(number_of_points > 0);

    tetgenInput.firstnumber = 0;
    tetgenInput.numberofpoints = number_of_points;
    tetgenInput.pointlist = new REAL[(number_of_points) * 3];
    tetgenInput.pointmarkerlist = new int[number_of_points];

    double *point_list = tetgenInput.pointlist;
    int *point_markerlist = tetgenInput.pointmarkerlist;

    for (unsigned int v = 0; v < number_of_points; v++)
    {
        const Eigen::Vector3d &point = points.col(v);
        point_list[3 * v] = point(0);
        point_list[3 * v + 1] = point(1);
        point_list[3 * v + 2] = point(2);

        point_markerlist[v] = points_marker.at(v);
    }
}
// ***************************************************************************
void TetgenInterface::CreateTetgenOutput(const double &maxTetrahedronArea,
                                         tetgenio &tetgenInput,
                                         tetgenio &tetgenOutput,
                                         const std::string &tetgenOptions) const
{
    ostringstream options;
    options.precision(16);
    options << tetgenOptions;

    if (maxTetrahedronArea > 0.0)
      options << maxTetrahedronArea;

    CreateTetgenOutput(tetgenInput, tetgenOutput, options.str());
}
// ***************************************************************************
void TetgenInterface::CreateTetgenOutput(tetgenio &tetgenInput, tetgenio &tetgenOutput, const std::string &tetgenOptions) const
{
    tetgenbehavior b;

    size_t sizeOptions = tetgenOptions.size();
    char *optionPointer = new char[sizeOptions + 1];
    tetgenOptions.copy(optionPointer, sizeOptions);
    optionPointer[sizeOptions] = '\0';

    b.parse_commandline(optionPointer);

    tetrahedralize(&b, &tetgenInput, &tetgenOutput);

    delete[] optionPointer;
}
// ***************************************************************************
void TetgenInterface::ConvertTetgenOutputToMeshDAO(const tetgenio &tetgenOutput, IMeshDAO &mesh) const
{
    /// <li>	Fill mesh structures
    const unsigned int numberOfCellsMesh = tetgenOutput.numberoftetrahedra;
    const unsigned int numberOfFacesMesh = tetgenOutput.numberoftrifaces;
    const unsigned int numberOfEgdesMesh = tetgenOutput.numberofedges;
    const unsigned int numberOfPointsMesh = tetgenOutput.numberofpoints;
    const unsigned int num_tetra_attribute = tetgenOutput.numberoftetrahedronattributes;


    mesh.InitializeDimension(3);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEgdesMesh);
    mesh.Cell2DsInitialize(numberOfFacesMesh);
    mesh.Cell3DsInitialize(numberOfCellsMesh);

    /// <li> Set Cell0Ds
    for (unsigned int p = 0; p < numberOfPointsMesh; p++)
    {
        mesh.Cell0DSetState(p, true);
        mesh.Cell0DInsertCoordinates(
            p,
            Vector3d(tetgenOutput.pointlist[3 * p], tetgenOutput.pointlist[3 * p + 1], tetgenOutput.pointlist[3 * p + 2]));
        mesh.Cell0DSetMarker(p, tetgenOutput.pointmarkerlist[p]);
    }

    /// <li> Set Cell1Ds
    struct Edge
    {
        unsigned int Cell1DIndex;
        int Cell0DEnd;
    };

    std::vector<std::list<Edge>> cell0DsCell1Ds(mesh.Cell0DTotalNumber());
    for (unsigned int e = 0; e < numberOfEgdesMesh; e++)
    {
        mesh.Cell1DInsertExtremes(e, tetgenOutput.edgelist[2 * e + 0], tetgenOutput.edgelist[2 * e + 1]);
        mesh.Cell1DSetState(e, true);
        mesh.Cell1DSetMarker(e, tetgenOutput.edgemarkerlist[e]);
        cell0DsCell1Ds[tetgenOutput.edgelist[2 * e + 0]].push_back({e, tetgenOutput.edgelist[2 * e + 1]});
    }

    /// <li> Set Faces
    vector<unsigned int> vertices(3);
    vector<unsigned int> edgeEndPoints(2);

    mesh.Cell2DsInitializeVertices(3);
    mesh.Cell2DsInitializeEdges(3);

    unsigned long estimatedValue = 2 * numberOfPointsMesh * numberOfPointsMesh + 2 * numberOfPointsMesh + 1;
    SparseMatrix<int> connectivityPointsFaces(numberOfPointsMesh, estimatedValue);
    vector<Triplet<int, unsigned long>> tripletsFaces;
    tripletsFaces.reserve(numberOfFacesMesh);
    for (unsigned int f = 0; f < numberOfFacesMesh; f++)
    {
        const unsigned int numCell2DVertices = 3;

        for (unsigned int v = 0; v < numCell2DVertices; v++)
        {
            const int vertexId = tetgenOutput.trifacelist[3 * f + v];
            const int nextVertexId = tetgenOutput.trifacelist[3 * f + (v + 1) % 3];
            unsigned int edgeId = numberOfEgdesMesh;

            if (edgeId == numberOfEgdesMesh)
            {
                const std::list<Edge> &originCell1Ds = cell0DsCell1Ds.at(vertexId);
                std::list<Edge>::const_iterator findOriginEdge =
                    std::find_if(originCell1Ds.begin(), originCell1Ds.end(), [&](const Edge &edge) {
                        return edge.Cell0DEnd == nextVertexId;
                    });

                if (findOriginEdge != originCell1Ds.end())
                    edgeId = findOriginEdge->Cell1DIndex;
            }

            if (edgeId == numberOfEgdesMesh)
            {
                const std::list<Edge> &endCell1Ds = cell0DsCell1Ds.at(nextVertexId);
                std::list<Edge>::const_iterator findEndEdge =
                    std::find_if(endCell1Ds.begin(), endCell1Ds.end(), [&](const Edge &edge) {
                        return edge.Cell0DEnd == vertexId;
                    });

                if (findEndEdge != endCell1Ds.end())
                    edgeId = findEndEdge->Cell1DIndex;
            }

            Gedim::Output::Assert(edgeId != numberOfEgdesMesh);

            vertices[v] = vertexId;

            mesh.Cell2DInsertVertex(f, v, vertexId);
            mesh.Cell2DInsertEdge(f, v, edgeId);
            if (mesh.Cell1DMarker(edgeId) == 0)
                mesh.Cell1DSetMarker(edgeId, tetgenOutput.trifacemarkerlist[f]);
        }

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, tetgenOutput.trifacemarkerlist[f]);

        sort(vertices.begin(), vertices.end());
        unsigned long indexI = vertices[0];
        unsigned long indexJK = (vertices[1] + vertices[2]) * (vertices[1] + vertices[2] + 1) * 0.5 + vertices[2] + 1;
        tripletsFaces.push_back(Triplet<int, unsigned long>(indexI, indexJK, f + 1));
    }

    connectivityPointsFaces.setFromTriplets(tripletsFaces.begin(), tripletsFaces.end());
    connectivityPointsFaces.makeCompressed();

    /// <li> Set Cells
    mesh.Cell3DsInitializeVertices(std::vector<unsigned int>(numberOfCellsMesh, 4));
    mesh.Cell3DsInitializeEdges(std::vector<unsigned int>(numberOfCellsMesh, 6));
    mesh.Cell3DsInitializeFaces(std::vector<unsigned int>(numberOfCellsMesh, 4));

    if (num_tetra_attribute > 0)
    {
      mesh.Cell3DInitializeDoubleProperties(num_tetra_attribute);
      for (unsigned int c_p = 0; c_p < num_tetra_attribute; ++c_p)
        mesh.Cell3DAddDoubleProperty("Attribute_" + std::to_string(c_p));
    }

    vector<unsigned int> faceVertices(3);
    for (unsigned int c = 0; c < numberOfCellsMesh; c++)
    {
        const unsigned int numVertices = 4;
        const unsigned int numFaces = 4;

        for (unsigned int v = 0; v < numVertices; v++)
            mesh.Cell3DInsertVertex(c, v, tetgenOutput.tetrahedronlist[tetgenOutput.numberofcorners * c + v]);

        unordered_set<unsigned int> cell3DEdges; // find cell edges and faces
        for (unsigned int j = 0; j < numFaces; j++)
        {
            for (unsigned int k = 0; k < 3; k++)
            {
                const unsigned int faceVertexId = (mesh.Cell3DVertex(c, (j + k) % 4));
                faceVertices[k] = faceVertexId;
            }

            sort(faceVertices.begin(), faceVertices.end());
            const unsigned int indexI = faceVertices[0];
            const unsigned int indexJK =
                (faceVertices[1] + faceVertices[2]) * (faceVertices[1] + faceVertices[2] + 1) * 0.5 + faceVertices[2] + 1;

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

        if (num_tetra_attribute > 0)
        {
          for (unsigned int c_p = 0; c_p < num_tetra_attribute; ++c_p)
          {
            mesh.Cell3DInitializeDoublePropertyValues(c, c_p, 1);
            mesh.Cell3DInsertDoublePropertyValue(c, c_p, 0, tetgenOutput.tetrahedronattributelist[c + c_p * numberOfCellsMesh]);
          }
        }
    }
}
// ***************************************************************************
void TetgenInterface::ExportTetgenOutput(const string &nameFolder, const string &nameFile, tetgenio &tetgenOutput) const
{

    ostringstream nameFolderStream, nameFileStream;

    nameFolderStream << nameFolder << "/";
    nameFolderStream << "Tetgen/";

    Output::CreateFolder(nameFolderStream.str());

    nameFileStream << nameFolderStream.str() << nameFile;

    Output::CreateFolder(nameFolderStream.str());

    tetgenOutput.firstnumber = 0;
    tetgenOutput.save_nodes((char *)nameFileStream.str().c_str());
    tetgenOutput.save_elements((char *)nameFileStream.str().c_str());
    tetgenOutput.save_faces((char *)nameFileStream.str().c_str());
    tetgenOutput.save_edges((char *)nameFileStream.str().c_str());
}
#endif
// ***************************************************************************
} // namespace Gedim
