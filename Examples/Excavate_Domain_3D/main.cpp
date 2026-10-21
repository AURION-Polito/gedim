#include "CreateMeshWithHoles.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_GaussLobatto1D.hpp"
#include "TetgenInterface.hpp"
#include "UCDUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VoroInterface.hpp"
#include "VtkMeshExporter.hpp"
#include "program_configuration.hpp"
#include <numeric>

using namespace std;
using namespace Eigen;

void ExportDomain3DToUCD(const string &file_path, const std::vector<Gedim::GeometryUtilities::Polyhedron> &domains)
{
    for (unsigned int i = 0; i < domains.size(); i++)
    {
        Gedim::UCDUtilities ucdUtilities;

        std::vector<std::vector<unsigned int>> polyhedra_vertices(1);
        polyhedra_vertices[0].resize(domains[i].Vertices.cols());
        for (unsigned int i = 0; i < domains[i].Vertices.cols(); i++)
            polyhedra_vertices[0][i] = i;

        Eigen::VectorXi materials(1);
        materials << i;

        ucdUtilities.ExportPolyhedra(file_path + "/Domain_" + std::to_string(i) + ".inp", domains[i].Vertices, polyhedra_vertices, {}, {}, materials);

        Gedim::Output::PrintGenericMessage(Gedim::Output::MagentaColor + "Export 3D domain to paraview. File path: " + file_path +
                                               "/Domain_" + std::to_string(i) + ".inp" + Gedim::Output::EndColor,
                                           true);
    }
}

void AddEdgeInternalPoints(const Gedim::GeometryUtilities &geometryUtilities,
                           const unsigned int &num_internal_points,
                           const InternalPointsType &internal_points_type,
                           Gedim::IMeshDAO &mesh)
{
    if (num_internal_points == 0)
        return;

    std::map<unsigned int, std::vector<double>> hanging_per_edges;

    switch (internal_points_type)
    {
    case InternalPointsType::Equispaced: {
        std::vector<double> reference_points = geometryUtilities.EquispaceCoordinates(num_internal_points, 0.0, 1.0, false);

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
            hanging_per_edges.insert({e, reference_points});
    }
    break;
    case InternalPointsType::GaussLobatto: {
        const Gedim::Quadrature::QuadratureData reference_data =
            Gedim::Quadrature::Quadrature_GaussLobatto1D::FillPointsAndWeights(2 * num_internal_points + 1);

        std::vector<double> reference_points(num_internal_points);
        for (unsigned int i = 0; i < num_internal_points; i++)
            reference_points[i] = reference_data.Points(0, i + 1);

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
            hanging_per_edges.insert({e, reference_points});
    }
    break;
    default:
        throw std::runtime_error("not valid type: (1) equispaced - (2) lobatto");
    }

    Gedim::MeshUtilities::AddHangingNodesToMesh(hanging_per_edges, mesh);
}

void CreateMesh3D(const Gedim::GeometryUtilities &geometryUtilities,
                  const Program_configuration &config,
                  const IDomain &domain_data,
                  Gedim::IMeshDAO &mesh)
{

    switch (config.MeshGenerator())
    {
    case MeshGenerator_Types_3D::Tetrahedral: {
#if ENABLE_TRIANGLE == 0
        throw std::runtime_error("Triangle library not active");
#endif
        Gedim::External::MeshUtilities::CreateTetrahedralMesh(domain_data.domains[0].Vertices,
                                                              domain_data.domains[0].Edges,
                                                              domain_data.domains[0].Faces,
                                                              config.MeshMaxVolume() * domain_data.volume,
                                                              mesh);
    }
    break;
    case MeshGenerator_Types_3D::Minimal: {
        std::vector<unsigned int> marker_vertices(domain_data.domains[0].Vertices.cols());
        std::iota(marker_vertices.begin(), marker_vertices.end(), 1);
        std::vector<unsigned int> marker_edges(domain_data.domains[0].Edges.cols());
        std::iota(marker_edges.begin(), marker_edges.end(), domain_data.domains[0].Vertices.cols() + 1);
        std::vector<unsigned int> marker_faces(domain_data.domains[0].Faces.size());
        std::iota(marker_faces.begin(),
                  marker_faces.end(),
                  domain_data.domains[0].Vertices.cols() + domain_data.domains[0].Edges.cols() + 1);
        Gedim::MeshUtilities::Mesh3DFromPolyhedron(domain_data.domains[0].Vertices,
                                                   domain_data.domains[0].Edges,
                                                   domain_data.domains[0].Faces,
                                                   marker_vertices,
                                                   marker_edges,
                                                   marker_faces,
                                                   mesh);
    }
    break;
    case MeshGenerator_Types_3D::Polyhedral: {
#if ENABLE_VORO == 0
        throw std::runtime_error("Voro library not active");
#endif
        const unsigned num_cells = static_cast<unsigned int>(std::max(1.0, 1.0 / config.MeshMaxVolume()));

        Gedim::External::MeshUtilities::CreatePolyhedralMesh(geometryUtilities,
                                                             domain_data.domains[0].Vertices,
                                                             domain_data.domains[0].Edges,
                                                             domain_data.domains[0].Faces,
                                                             num_cells,
                                                             10,
                                                             mesh,
                                                             10);
    }
    break;
    case MeshGenerator_Types_3D::Cubic: {
        const Eigen::Vector3d domain_origin = domain_data.domains[0].Vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = domain_data.domains[0].Vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_width_tangent = domain_data.domains[0].Vertices.col(4) - domain_origin;
        const Eigen::Vector3d domain_heigth_tangent = domain_data.domains[0].Vertices.col(3) - domain_origin;

        Gedim::MeshUtilities::CreateParallelepipedMesh(
            domain_origin,
            domain_base_tangent,
            domain_heigth_tangent,
            domain_width_tangent,
            geometryUtilities.EquispaceCoordinates(config.NumCellsBase() + 1, 0.0, 1.0, true),
            geometryUtilities.EquispaceCoordinates(config.NumCellsHeight() + 1, 0.0, 1.0, true),
            geometryUtilities.EquispaceCoordinates(config.NumCellsWidth() + 1, 0.0, 1.0, true),
            mesh);
    }
    break;
    default:
        throw std::runtime_error("not valid mesh generator type");
    }

#ifdef EXPORT_TEST_VTU
    {
        std::string mesh_3D_export_folder = "./TEST";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(full_mesh, mesh_3D_export_folder, "Full_Mesh");
    }
#endif
}

void ExcavateMesh3D(const Gedim::GeometryUtilities &geometryUtilities,
                    const Program_configuration &config,
                    const IDomain &domain_data,
                    Gedim::IMeshDAO &mesh)
{

    if (domain_data.domains.size() == 1)
    {
        CreateMesh3D(geometryUtilities, config, domain_data, mesh);

        return;
    }

    Gedim::MeshMatrices full_mesh_data;
    Gedim::MeshMatricesDAO full_mesh(full_mesh_data);

    CreateMesh3D(geometryUtilities, config, domain_data, full_mesh);

#ifdef EXPORT_TEST_VTU
    {
        std::string mesh_3D_export_folder = "./TEST";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(full_mesh, mesh_3D_export_folder, "Full_Mesh");
    }
#endif

    Gedim::External::MeshUtilities::CreateMeshWithHoles(geometryUtilities, domain_data.domains, full_mesh_data, full_mesh, mesh);
}

void CreateTetrahedralMesh3D(const IDomain &domain_data, const double &max_volume, Gedim::IMeshDAO &mesh)
{
#if ENABLE_TRIANGLE == 0
    throw std::runtime_error("Triangle library not active");
#endif

    if (domain_data.domains.size() == 1)
    {
        Gedim::External::MeshUtilities::CreateTetrahedralMesh(domain_data.domains[0].Vertices,
                                                              domain_data.domains[0].Edges,
                                                              domain_data.domains[0].Faces,
                                                              max_volume,
                                                              mesh);
    }
    else
    {
        unsigned int num_vertices = 0;
        unsigned int num_faces = 0;
        for (unsigned int d = 0; d < domain_data.domains.size(); d++)
        {
            num_vertices += domain_data.domains[d].Vertices.cols();
            num_faces += domain_data.domains[d].Faces.size();
        }

        Eigen::MatrixXd points = Eigen::MatrixXd::Zero(3, num_vertices);
        std::vector<std::vector<unsigned int>> facets(num_faces);
        unsigned int offset_vertices = 0;
        unsigned int offset_faces = 0;
        for (unsigned int d = 0; d < domain_data.domains.size(); d++)
        {
            unsigned int num_local_vertices = domain_data.domains[d].Vertices.cols();
            points.middleCols(offset_vertices, num_local_vertices) = domain_data.domains[d].Vertices;

            for (unsigned int f = 0; f < domain_data.domains[d].Faces.size(); f++)
            {
                facets[offset_faces].resize(domain_data.domains[d].Faces[f].cols());
                for (unsigned int v = 0; v < domain_data.domains[d].Faces[f].cols(); v++)
                    facets[offset_faces][v] = domain_data.domains[d].Faces[f](0, v) + offset_vertices;

                offset_faces++;
            }

            offset_vertices += num_local_vertices;
        }

        std::vector<Eigen::Vector3d> holes_points(domain_data.domains.size() - 1);
        for (unsigned int d = 0; d < domain_data.domains.size() - 1; d++)
            holes_points[d] = domain_data.internal_points[1 + d];

        std::vector<Gedim::TetgenInterface::Region> regions(1);
        regions[0].centroid = domain_data.internal_points[0];
        regions[0].id = 1;
        regions[0].max_volume = max_volume * domain_data.volume;

        Gedim::TetgenInterface tetgen_interface;
        tetgen_interface.CreateMesh(points, facets, regions, holes_points, mesh, "QpqfezAa");

        for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
        {
            if (mesh.Cell0DMarker(v) == std::numeric_limits<unsigned int>::max())
                mesh.Cell0DSetMarker(v, num_faces + num_vertices + 1);
        }

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
        {
            if (mesh.Cell1DMarker(e) == std::numeric_limits<unsigned int>::max())
                mesh.Cell1DSetMarker(e, num_faces + num_vertices + 1);
        }
    }
}

int main(int argc, char **argv)
{
    // This is a work in progress...

    Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini", false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    const std::string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const std::string exportCsvFolder = exportFolder + "/Csv";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const std::string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const std::string exportUCDFolder = exportFolder + "/UCD";
    Gedim::Output::CreateFolder(exportUCDFolder);
    const std::string exportMeditFolder = exportFolder + "/Medit";
    Gedim::Output::CreateFolder(exportMeditFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    const Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    IDomain domain_data;
    switch (config.TestDomain())
    {
    case TestDomainTypes::SquareDomainWithOneHole:
        domain_data = SquareDomainWithOneHole();
        break;
    }

    if (config.ExportFormat()[1])
    {
        for (unsigned int d = 0; d < domain_data.domains.size(); d++)
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolyhedron(domain_data.domains[d].Vertices,
                                       domain_data.domains[d].Edges,
                                       domain_data.domains[d].Faces);
            vtkUtilities.Export(exportVtuFolder + "/Domain" + std::to_string(d) + ".vtu");
            Gedim::Output::PrintGenericMessage(Gedim::Output::MagentaColor + "Domain " + std::to_string(d) +
                                                   " is exported in: " + exportVtuFolder + "/Domain.vtu" + Gedim::Output::EndColor,
                                               true);
        }
    }

    if (config.ExportFormat()[2])
    {
        ExportDomain3DToUCD(exportUCDFolder, domain_data.domains);
    }

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    switch (config.MeshGenerator())
    {
    case MeshGenerator_Types_3D::PureTetrahedral: {
        CreateTetrahedralMesh3D(domain_data, config.MeshMaxVolume(), mesh);
    }
    break;
    case MeshGenerator_Types_3D::Tetrahedral:
    case MeshGenerator_Types_3D::Minimal:
    case MeshGenerator_Types_3D::Polyhedral:
        throw std::runtime_error("Case not supported yet");
    case MeshGenerator_Types_3D::Cubic: {
        if (config.NumCellsBase() % 4 != 0 || config.NumCellsWidth() % 4 != 0 || config.NumCellsHeight() % 4 != 0)
            throw std::runtime_error("Case not supported yet");

        ExcavateMesh3D(geometryUtilities, config, domain_data, mesh);
    }
    break;
    default:
        throw std::runtime_error("not valid mesh generator type");
    }

    AddEdgeInternalPoints(geometryUtilities, config.MethodOrder() - 1, config.InternalPoints(), mesh);

    if (config.ExportFormat()[0])
    {
        Gedim::MeshUtilities meshUtilities;
        const Gedim::MeshFromCsvUtilities utilities;
        Gedim::MeshDAOExporterToCsv exporter(utilities);
        Gedim::MeshFromCsvUtilities::Configuration configuration;
        configuration.Folder = exportCsvFolder;
        exporter.Export(configuration, mesh);

        Gedim::Output::PrintGenericMessage(
            Gedim::Output::MagentaColor + "Export 3D mesh to Csv. File path: " + exportCsvFolder + Gedim::Output::EndColor,
            true);
    }

    if (config.ExportFormat()[1])
    {
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
        Gedim::Output::PrintGenericMessage(
            Gedim::Output::MagentaColor + "Export 3D mesh to VTU. File path:: " + exportVtuFolder + Gedim::Output::EndColor,
            true);
    }

    if (config.ExportFormat()[2])
    {

        Gedim::MeshUtilities::ExportMeshToUCD(mesh, exportUCDFolder, "Domain_Mesh", false);

        Gedim::Output::PrintGenericMessage(
            Gedim::Output::MagentaColor + "Export 3D mesh to UCD. File path: " + exportUCDFolder + Gedim::Output::EndColor,
            true);
    }

    if (config.ExportFormat()[3])
    {
        Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
        Gedim::MeshUtilities::ExportMeshToMEDIT(mesh, exportMeditFolder, "Domain_Mesh", false, config.MethodOrder());

        Gedim::Output::PrintGenericMessage(
            Gedim::Output::MagentaColor + "Export 3D mesh to MEDIT. File path: " + exportMeditFolder + Gedim::Output::EndColor,
            true);
    }

    return 0;
}
