#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "MeshFromCsvUtilities.hpp"
#include "MeshUtilities.hpp"
#include "VtkMeshExporter.hpp"
#include "program_configuration.hpp"
#include <numeric>

int main(int argc, char **argv)
{
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

  Gedim::MeshMatrices mesh_data;
  Gedim::MeshMatricesDAO mesh(mesh_data);

  switch (config.MeshImportType())
  {
    case Mesh_Import_Types::CsvImporter:
    {
      Gedim::MeshFromCsvUtilities importerUtilities;
      Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
      meshImporterConfiguration.Folder = config.MeshImportFilePath();
      meshImporterConfiguration.Separator = ';';
      Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
      importer.Import(meshImporterConfiguration, mesh);
    }
      break;
    case Mesh_Import_Types::OVMImporter:
    {
      std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
      Gedim::MeshUtilities::ImportOpenVolumeMesh(config.MeshImportFilePath(), mesh, meshCell3DsFacesOrientation);
    }
      break;
    default:
      throw std::runtime_error("MeshImportType " + std::to_string((unsigned int)config.MeshImportType()) + " not supported");
  }


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
    Gedim::MeshUtilities::ExportMeshToMEDIT(mesh, exportMeditFolder, "Domain_Mesh", false, 0);

    Gedim::Output::PrintGenericMessage(
          Gedim::Output::MagentaColor + "Export 3D mesh to MEDIT. File path: " + exportMeditFolder + Gedim::Output::EndColor,
          true);
  }

  return 0;
}
