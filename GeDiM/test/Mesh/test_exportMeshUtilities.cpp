#include "test_exportMeshUtilities.hpp"
#include "ImportExportUtilities.hpp"
#include <fstream>
#include <sys/stat.h>

namespace GedimUnitTesting
{
// ***************************************************************************
ExportMeshData ExportMeshUtilities::ImportMesh2DFromText(const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    ExportMeshData mesh;

    std::ifstream inFile;
    inFile.open(file_path.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Import file '" + file_path + "' not opened");

    inFile >> mesh.Cell0Ds;
    inFile >> mesh.Cell1Ds;
    inFile >> mesh.Cell2Ds;
    inFile >> mesh.Cell0DsMarker;
    inFile >> mesh.Cell1DsMarker;
    inFile >> mesh.Cell2DsMarker;

    inFile.close();

    return mesh;
}
// ***************************************************************************
void ExportMeshUtilities::ExportMesh2DToText(const ExportMeshData &mesh_data, const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    std::ofstream file(file_path);

    if (!file.is_open())
        throw std::runtime_error("Export file '" + file_path + "' not opened");

    file.precision(16);

    file << std::scientific << mesh_data.Cell0Ds << std::endl;
    file << std::scientific << mesh_data.Cell1Ds << std::endl;
    file << std::scientific << mesh_data.Cell2Ds << std::endl;
    file << std::scientific << mesh_data.Cell0DsMarker << std::endl;
    file << std::scientific << mesh_data.Cell1DsMarker << std::endl;
    file << std::scientific << mesh_data.Cell2DsMarker << std::endl;

    file.close();
}
// ***************************************************************************
ExportMeshGeometricData2D ExportMeshUtilities::ImportMeshGeometricData2DFromText(const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    ExportMeshGeometricData2D mesh_geometric_data;

    std::ifstream inFile;
    inFile.open(file_path.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Import file '" + file_path + "' not opened");

    inFile >> mesh_geometric_data.PolygonsVertices;
    inFile >> mesh_geometric_data.PolygonsTriangulations;
    inFile >> mesh_geometric_data.PolygonsArea;
    inFile >> mesh_geometric_data.PolygonsCentroid;
    inFile >> mesh_geometric_data.PolygonsDiameter;
    inFile >> mesh_geometric_data.PolygonsEdgesDirection;
    inFile >> mesh_geometric_data.PolygonsEdgesLength;
    inFile >> mesh_geometric_data.PolygonsEdgesTangent;
    inFile >> mesh_geometric_data.PolygonsEdgesTangentNormalized;
    inFile >> mesh_geometric_data.PolygonsEdgesCentroid;
    inFile >> mesh_geometric_data.PolygonsEdgesNormal;

    return mesh_geometric_data;
}
// ***************************************************************************
void ExportMeshUtilities::ExportMeshGeometricData2DToText(const ExportMeshGeometricData2D &mesh_geometric_data,
                                                          const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    std::ofstream file(file_path);

    if (!file.is_open())
        throw std::runtime_error("Export file '" + file_path + "' not opened");

    file.precision(16);

    file << std::scientific << mesh_geometric_data.PolygonsVertices << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsTriangulations << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsArea << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsCentroid << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsDiameter << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesDirection << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesLength << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesTangent << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesTangentNormalized << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesCentroid << std::endl;
    file << std::scientific << mesh_geometric_data.PolygonsEdgesNormal << std::endl;

    file.close();
}
// ***************************************************************************
ExportMeshData ExportMeshUtilities::ImportMesh3DFromText(const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    ExportMeshData mesh;

    std::ifstream inFile;
    inFile.open(file_path.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Import file '" + file_path + "' not opened");

    inFile >> mesh.Cell0Ds;
    inFile >> mesh.Cell1Ds;
    inFile >> mesh.Cell2Ds;
    inFile >> mesh.Cell3DsVertices;
    inFile >> mesh.Cell3DsEdges;
    inFile >> mesh.Cell3DsFaces;
    inFile >> mesh.Cell0DsMarker;
    inFile >> mesh.Cell1DsMarker;
    inFile >> mesh.Cell2DsMarker;
    inFile >> mesh.Cell3DsMarker;

    inFile.close();

    return mesh;
}
// ***************************************************************************
void ExportMeshUtilities::ExportMesh3DToText(const ExportMeshData &mesh_data, const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    std::ofstream file(file_path);

    if (!file.is_open())
        throw std::runtime_error("Export file '" + file_path + "' not opened");

    file.precision(16);

    file << std::scientific << mesh_data.Cell0Ds << std::endl;
    file << std::scientific << mesh_data.Cell1Ds << std::endl;
    file << std::scientific << mesh_data.Cell2Ds << std::endl;
    file << std::scientific << mesh_data.Cell3DsVertices << std::endl;
    file << std::scientific << mesh_data.Cell3DsEdges << std::endl;
    file << std::scientific << mesh_data.Cell3DsFaces << std::endl;
    file << std::scientific << mesh_data.Cell0DsMarker << std::endl;
    file << std::scientific << mesh_data.Cell1DsMarker << std::endl;
    file << std::scientific << mesh_data.Cell2DsMarker << std::endl;
    file << std::scientific << mesh_data.Cell3DsMarker << std::endl;

    file.close();
}
// ***************************************************************************
ExportMeshGeometricData3D ExportMeshUtilities::ImportMeshGeometricData3DFromText(const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    ExportMeshGeometricData3D mesh_geometric_data;

    std::ifstream inFile;
    inFile.open(file_path.c_str());

    if (!inFile.is_open())
        throw std::runtime_error("Import file '" + file_path + "' not opened");

    inFile >> mesh_geometric_data.PolyhedronsVertices;
    inFile >> mesh_geometric_data.PolyhedronsEdges;
    inFile >> mesh_geometric_data.PolyhedronsFaces;
    inFile >> mesh_geometric_data.PolyhedronsVolume;
    inFile >> mesh_geometric_data.PolyhedronsDiameter;
    inFile >> mesh_geometric_data.PolyhedronsCentroid;
    inFile >> mesh_geometric_data.PolyhedronsTetrahedronsVertices;
    inFile >> mesh_geometric_data.PolyhedronsFacesTranslation;
    inFile >> mesh_geometric_data.PolyhedronsFacesRotationMatrix;
    inFile >> mesh_geometric_data.PolyhedronsFacesNormal;
    inFile >> mesh_geometric_data.PolyhedronsFacesNormalDirection;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdgesDirection;
    inFile >> mesh_geometric_data.PolyhedronsFaces2DVertices;
    inFile >> mesh_geometric_data.PolyhedronsFacesTriangulations2DVertices;
    inFile >> mesh_geometric_data.PolyhedronsFacesArea;
    inFile >> mesh_geometric_data.PolyhedronsFaces2DCentroid;
    inFile >> mesh_geometric_data.PolyhedronsFacesDiameter;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdgesLength;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdges2DCentroid;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdges2DTangent;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdges2DTangentNormalized;
    inFile >> mesh_geometric_data.PolyhedronsFacesEdges2DNormal;

    return mesh_geometric_data;
}
// ***************************************************************************
void ExportMeshUtilities::ExportMeshGeometricData3DToText(const ExportMeshGeometricData3D &mesh_geometric_data,
                                                          const std::string &file_path)
{
    using namespace Gedim_ImportExport_Utilities;

    std::ofstream file(file_path);

    if (!file.is_open())
        throw std::runtime_error("Export file '" + file_path + "' not opened");

    file.precision(16);

    file << std::scientific << mesh_geometric_data.PolyhedronsVertices << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsEdges << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFaces << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsVolume << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsDiameter << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsCentroid << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsTetrahedronsVertices << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesTranslation << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesRotationMatrix << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesNormal << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesNormalDirection << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdgesDirection << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFaces2DVertices << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesTriangulations2DVertices << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesArea << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFaces2DCentroid << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesDiameter << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdgesLength << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdges2DCentroid << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdges2DTangent << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdges2DTangentNormalized << std::endl;
    file << std::scientific << mesh_geometric_data.PolyhedronsFacesEdges2DNormal << std::endl;

    file.close();
}
// ***************************************************************************
} // namespace GedimUnitTesting
