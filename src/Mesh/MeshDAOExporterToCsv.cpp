#include "MeshDAOExporterToCsv.hpp"
#include <iostream>
#include <fstream>

namespace Gedim
{
  // ***************************************************************************
  MeshDAOExporterToCsv::MeshDAOExporterToCsv(const MeshFromCsvUtilities& utilities) :
    utilities(utilities)
  {
  }
  MeshDAOExporterToCsv::~MeshDAOExporterToCsv()
  {
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::Export(const MeshFromCsvUtilities::Configuration& configuration,
                                    const IMeshDAO& mesh) const
  {
    utilities.ExportCell0Ds(configuration.Folder + "/" +
                            configuration.FileCell0DsName + "." +
                            configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell1Ds(configuration.Folder + "/" +
                            configuration.FileCell1DsName + "." +
                            configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell2Ds(configuration.Folder + "/" +
                            configuration.FileCell2DsName + "." +
                            configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell3Ds(configuration.Folder + "/" +
                            configuration.FileCell3DsName + "." +
                            configuration.FileExtension,
                            configuration.Separator,
                            mesh);

    utilities.ExportCell0DProperties(configuration.Folder,
                                     configuration.FileCell0DPropertiesName,
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);

    utilities.ExportCell1DProperties(configuration.Folder,
                                     configuration.FileCell1DPropertiesName,
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);

    utilities.ExportCell2DProperties(configuration.Folder,
                                     configuration.FileCell2DPropertiesName,
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);

    utilities.ExportCell3DProperties(configuration.Folder,
                                     configuration.FileCell3DPropertiesName,
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);

    utilities.ExportCell0DNeighbours(configuration.Folder + "/" +
                                     configuration.FileCell0DNeighboursName + "." +
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell1DNeighbours(configuration.Folder + "/" +
                                     configuration.FileCell1DNeighboursName + "." +
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell2DNeighbours(configuration.Folder + "/" +
                                     configuration.FileCell2DNeighboursName + "." +
                                     configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell2DSubDivisions(configuration.Folder + "/" +
                                       configuration.FileCell2DSubDivisionsName + "." +
                                       configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);
  }
  // ***************************************************************************
}
