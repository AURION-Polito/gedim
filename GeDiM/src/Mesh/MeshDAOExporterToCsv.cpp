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

#include "MeshDAOExporterToCsv.hpp"
#include <fstream>
#include <iostream>

namespace Gedim
{
// ***************************************************************************
MeshDAOExporterToCsv::MeshDAOExporterToCsv(const Gedim::MeshFromCsvUtilities &utilities) : utilities(utilities)
{
}
MeshDAOExporterToCsv::~MeshDAOExporterToCsv()
{
}
// ***************************************************************************
void MeshDAOExporterToCsv::Export(const Gedim::MeshFromCsvUtilities::Configuration &configuration, const Gedim::IMeshDAO &mesh) const
{
    utilities.ExportCell0Ds(configuration.Folder + "/" + configuration.FileCell0DsName + "." + configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell1Ds(configuration.Folder + "/" + configuration.FileCell1DsName + "." + configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell2Ds(configuration.Folder + "/" + configuration.FileCell2DsName + "." + configuration.FileExtension,
                            configuration.Separator,
                            mesh);
    utilities.ExportCell3Ds(configuration.Folder + "/" + configuration.FileCell3DsName + "." + configuration.FileExtension,
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

    utilities.ExportCell0DNeighbours(configuration.Folder + "/" + configuration.FileCell0DNeighboursName + "." +
                                         configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell1DNeighbours(configuration.Folder + "/" + configuration.FileCell1DNeighboursName + "." +
                                         configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell2DNeighbours(configuration.Folder + "/" + configuration.FileCell2DNeighboursName + "." +
                                         configuration.FileExtension,
                                     configuration.Separator,
                                     mesh);
    utilities.ExportCell2DSubDivisions(configuration.Folder + "/" + configuration.FileCell2DSubDivisionsName + "." +
                                           configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);

    utilities.ExportCell0DUpdatedCells(configuration.Folder + "/" + configuration.FileCell0DUpdatedCellsName + "." +
                                           configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);
    utilities.ExportCell1DUpdatedCells(configuration.Folder + "/" + configuration.FileCell1DUpdatedCellsName + "." +
                                           configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);
    utilities.ExportCell2DUpdatedCells(configuration.Folder + "/" + configuration.FileCell2DUpdatedCellsName + "." +
                                           configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);
    utilities.ExportCell3DUpdatedCells(configuration.Folder + "/" + configuration.FileCell3DUpdatedCellsName + "." +
                                           configuration.FileExtension,
                                       configuration.Separator,
                                       mesh);
}
// ***************************************************************************
} // namespace Gedim
