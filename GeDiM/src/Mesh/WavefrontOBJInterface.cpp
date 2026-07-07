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

#include "WavefrontOBJInterface.hpp"
#include "FileTextReader.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
// ***************************************************************************
WavefrontOBJInterface::WavefrontOBJInterface()
{
}
WavefrontOBJInterface::~WavefrontOBJInterface()
{
}
// ***************************************************************************
WavefrontOBJInterface::OBJMesh WavefrontOBJInterface::StringsToOBJMesh(const std::vector<std::string> &fileLines) const
{
    OBJMesh mesh;

    std::list<Eigen::Vector3d> vertices;
    std::list<Eigen::VectorXi> faces;

    for (const std::string &line : fileLines)
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);

        std::string tag;
        iss >> tag;

        if (tag == "v")
        {
            double x, y, z;
            iss >> x >> y >> z;

            vertices.push_back(Eigen::Vector3d(x, y, z));
        }
        else if (tag == "f")
        {
            if (line.find('/') != std::string::npos)
                throw std::runtime_error("obj format v/vt v//vn v/vt/vn  not supported");

            std::list<int> faceIndices;
            std::string token;

            while (iss >> token)
            {
                // Handle:
                // v
                // v/vt (NOT supported)
                // v//vn (NOT supported)
                // v/vt/vn (NOT supported)
                int idx = std::stoi(token);

                // Negative OBJ indexing
                if (idx < 0)
                    idx = static_cast<int>(vertices.size()) + idx + 1;

                --idx;

                faceIndices.push_back(idx);
            }

            faces.push_back(Eigen::VectorXi(faceIndices.size()));
            auto &face = faces.back();

            unsigned int fv = 0;
            for (const auto f_v : faceIndices)
                face[fv++] = f_v;
        }
    }

    mesh.NumCell0Ds = static_cast<unsigned int>(vertices.size());
    mesh.NumCell2Ds = static_cast<unsigned int>(faces.size());
    mesh.Cell0Ds.resize(3, mesh.NumCell0Ds);

    unsigned int v = 0;
    for (const auto &vertex : vertices)
        mesh.Cell0Ds.col(v++) = vertex;

    mesh.Cell2Ds = std::vector<Eigen::VectorXi>(faces.begin(), faces.end());

    return mesh;
}
// ***************************************************************************
std::vector<string> WavefrontOBJInterface::OBJMeshToStrings(const OBJMesh &mesh) const
{
    list<string> lines;

    for (unsigned int v = 0; v < mesh.NumCell0Ds; v++)
    {
        ostringstream stream;
        stream.precision(16);
        stream << "v" << " ";
        stream << scientific << mesh.Cell0Ds(0, v) << " ";
        stream << scientific << mesh.Cell0Ds(1, v) << " ";
        stream << scientific << mesh.Cell0Ds(2, v);

        lines.push_back(stream.str());
    }

    for (unsigned int f = 0; f < mesh.NumCell2Ds; f++)
    {
        ostringstream stream;
        stream << "f" << " ";
        for (unsigned int fv = 0; fv < mesh.Cell2Ds[f].size(); fv++)
            stream << " " << mesh.Cell2Ds[f][fv] + 1;

        lines.push_back(stream.str());
    }

    return std::vector<string>(lines.begin(), lines.end());
}
// ***************************************************************************
WavefrontOBJInterface::OBJMesh WavefrontOBJInterface::MeshDAOToOBJMesh(const IMeshDAO &originalMesh) const
{
    WavefrontOBJInterface::OBJMesh convertedMesh;

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
void WavefrontOBJInterface::OBJMeshToMeshDAO(const OBJMesh &originalMesh, const MeshUtilities &meshUtilities, IMeshDAO &convertedMesh) const
{
    convertedMesh.InitializeDimension(2);

    // Create Cell0Ds
    convertedMesh.Cell0DsInitialize(originalMesh.NumCell0Ds);
    convertedMesh.Cell0DsInsertCoordinates(originalMesh.Cell0Ds);
    for (unsigned int v = 0; v < originalMesh.NumCell0Ds; v++)
        convertedMesh.Cell0DSetState(v, true);

    // Compute Cell1Ds
    const MeshUtilities::ComputeMesh2DCell1DsResult cell1Ds =
        meshUtilities.ComputeMesh2DCell1Ds(originalMesh.Cell0Ds, originalMesh.Cell2Ds);

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
            const Eigen::MatrixXi &OBJ_face = cell1Ds.Cell2Ds[f];
            facesNumberVertices[f] = OBJ_face.cols();
        }

        convertedMesh.Cell2DsInitialize(originalMesh.NumCell2Ds);
        convertedMesh.Cell2DsInitializeVertices(facesNumberVertices);
        convertedMesh.Cell2DsInitializeEdges(facesNumberVertices);
        for (unsigned int f = 0; f < originalMesh.NumCell2Ds; f++)
        {
            const Eigen::MatrixXi &OBJ_face = cell1Ds.Cell2Ds[f];
            const unsigned int numFaceVertices = OBJ_face.cols();

            for (unsigned int fv = 0; fv < numFaceVertices; fv++)
            {
                convertedMesh.Cell2DInsertVertex(f, fv, OBJ_face(0, fv));
                convertedMesh.Cell2DInsertEdge(f, fv, OBJ_face(1, fv));
            }

            convertedMesh.Cell2DSetState(f, true);
        }
    }
}
// ***************************************************************************
std::vector<std::string> WavefrontOBJInterface::ImportMeshFromFile(const std::string &OBJFilePath,
                                                                   const MeshUtilities &meshUtilities,
                                                                   IMeshDAO &mesh) const
{
    vector<string> fileLines;

    {
        FileReader fileReader(OBJFilePath);

        if (!fileReader.Open())
            throw runtime_error("File " + OBJFilePath + " not found");

        fileReader.GetAllLines(fileLines);
        fileReader.Close();
    }

    const OBJMesh meshImported = StringsToOBJMesh(fileLines);
    OBJMeshToMeshDAO(meshImported, meshUtilities, mesh);

    return fileLines;
}
// ***************************************************************************
void WavefrontOBJInterface::ExportMeshToFile(const IMeshDAO &mesh, const std::string &OBJFilePath) const
{
    const OBJMesh meshToExport = MeshDAOToOBJMesh(mesh);
    const vector<string> fileLines = OBJMeshToStrings(meshToExport);

    ofstream exportFile(OBJFilePath);
    if (!exportFile.is_open())
        throw runtime_error("Unable to export file " + OBJFilePath);

    for (const string &line : fileLines)
        exportFile << line << endl;
    exportFile.close();
}
// ***************************************************************************
} // namespace Gedim
