#ifndef __program_configuration_hpp
#define __program_configuration_hpp

#include "Configurations.hpp"

enum struct Mesh_Import_Types
{
    CsvImporter = 0, ///< import mesh from CSV format
    OVMImporter = 1, ///< import mesh from OVM format
};

struct Program_configuration final
{
    Program_configuration()
    {
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Export", "Folder where to export data (Default: ./Export)");
        Gedim::Configurations::AddProperty("ExportFormat",
                                           std::vector<unsigned int>({1, 0, 0, 0}),
                                           "A boolean vector of kind of desired export type for solution and mesh "
                                           "[Csv, Vtu, UCD, Medit] (Default: [1,0,0,0])");

        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshImportType",
                                           static_cast<unsigned int>(Mesh_Import_Types::CsvImporter),
                                           "Mesh import type, 0 - CsvImporter; 1 - OVMImporter (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
  }


    inline std::string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("ExportFolder");
    }
    inline std::vector<unsigned int> ExportFormat() const
    {
        return Gedim::Configurations::GetPropertyValue<std::vector<unsigned int>>("ExportFormat");
    }
    inline Mesh_Import_Types MeshImportType() const
    {
        return (Mesh_Import_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshImportType");
    }
    inline std::string MeshImportFilePath() const
    {
      return Gedim::Configurations::GetPropertyValue<std::string>("MeshImportFilePath");
    }
};

#endif // __program_configuration_hpp
