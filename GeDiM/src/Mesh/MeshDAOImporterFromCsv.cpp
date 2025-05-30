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

#include "MeshDAOImporterFromCsv.hpp"
#include "FileTextReader.hpp"
#include <fstream>
#include <iostream>

namespace Gedim
{
// ***************************************************************************
MeshDAOImporterFromCsv::MeshDAOImporterFromCsv(const MeshFromCsvUtilities &utilities) : utilities(utilities)
{
}
MeshDAOImporterFromCsv::~MeshDAOImporterFromCsv()
{
}
// ***************************************************************************
void MeshDAOImporterFromCsv::Import(const MeshFromCsvUtilities::Configuration &configuration, IMeshDAO &mesh)
{
    Gedim::FileReader csvCell0DsFile(configuration.Folder + "/" + configuration.FileCell0DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DsFile(configuration.Folder + "/" + configuration.FileCell1DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DsFile(configuration.Folder + "/" + configuration.FileCell2DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell3DsFile(configuration.Folder + "/" + configuration.FileCell3DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DPropertiesFile(configuration.Folder + "/" + configuration.FileCell0DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DPropertiesFile(configuration.Folder + "/" + configuration.FileCell1DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DPropertiesFile(configuration.Folder + "/" + configuration.FileCell2DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell3DPropertiesFile(configuration.Folder + "/" + configuration.FileCell3DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DNeighboursFile(configuration.Folder + "/" + configuration.FileCell0DNeighboursName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DNeighboursFile(configuration.Folder + "/" + configuration.FileCell1DNeighboursName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DNeighboursFile(configuration.Folder + "/" + configuration.FileCell2DNeighboursName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DSubDivisionsFile(configuration.Folder + "/" + configuration.FileCell2DSubDivisionsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell0DUpdatedCellsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell1DUpdatedCellsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell2DUpdatedCellsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell3DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell3DUpdatedCellsName +
                                                "." + configuration.FileExtension);

    Gedim::Profiler::StartTime("ReadFiles");
    std::vector<MeshFromCsvUtilities::Cell0D> cell0Ds = utilities.ImportCell0Ds(csvCell0DsFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell1D> cell1Ds = utilities.ImportCell1Ds(csvCell1DsFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell2D> cell2Ds = utilities.ImportCell2Ds(csvCell2DsFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell3D> cell3Ds = utilities.ImportCell3Ds(csvCell3DsFile, configuration.Separator);
    Gedim::Profiler::StopTime("ReadFiles");

    Gedim::Profiler::StartTime("ConvertFiles");
    Gedim::Profiler::StartTime("Cell0Ds");
    utilities.ConvertCell0Ds(cell0Ds, mesh);
    Gedim::Profiler::StopTime("Cell0Ds");
    Gedim::Profiler::StartTime("Cell1Ds");
    utilities.ConvertCell1Ds(cell1Ds, mesh);
    Gedim::Profiler::StopTime("Cell1Ds");
    Gedim::Profiler::StartTime("Cell2Ds");
    utilities.ConvertCell2Ds(cell2Ds, mesh);
    Gedim::Profiler::StopTime("Cell2Ds");
    Gedim::Profiler::StartTime("Cell3Ds");
    utilities.ConvertCell3Ds(cell3Ds, mesh);
    Gedim::Profiler::StopTime("Cell3Ds");
    Gedim::Profiler::StopTime("ConvertFiles");

    unsigned int meshDimension = 0;
    if (mesh.Cell0DTotalNumber() > 0 && mesh.Cell1DTotalNumber() > 0 && mesh.Cell2DTotalNumber() == 0 &&
        mesh.Cell3DTotalNumber() == 0)
        meshDimension = 1;
    else if (mesh.Cell0DTotalNumber() > 0 && mesh.Cell1DTotalNumber() > 0 && mesh.Cell2DTotalNumber() > 0 &&
             mesh.Cell3DTotalNumber() == 0)
        meshDimension = 2;
    else if (mesh.Cell0DTotalNumber() > 0 && mesh.Cell1DTotalNumber() > 0 && mesh.Cell2DTotalNumber() > 0 &&
             mesh.Cell3DTotalNumber() > 0)
        meshDimension = 3;
    else
        throw std::runtime_error("Dimension of imported mesh not recognized");

    mesh.InitializeDimension(meshDimension);

    Gedim::Profiler::StartTime("ReadProperties");
    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell0DPropertiesFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell1DPropertiesFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell2DPropertiesFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell3DPropertiesFile, configuration.Separator);
    Gedim::Profiler::StopTime("ReadProperties");

    Gedim::Profiler::StartTime("ConvertProperties");
    utilities.ConvertCell0DDoubleProperties(cell0DDoubleProperties, mesh);
    utilities.ConvertCell1DDoubleProperties(cell1DDoubleProperties, mesh);
    utilities.ConvertCell2DDoubleProperties(cell2DDoubleProperties, mesh);
    utilities.ConvertCell3DDoubleProperties(cell3DDoubleProperties, mesh);
    Gedim::Profiler::StopTime("ConvertProperties");

    Gedim::Profiler::StartTime("ReadNeighbours");
    std::vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours =
        utilities.ImportCell0DNeighbours(csvCell0DNeighboursFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours =
        utilities.ImportCell1DNeighbours(csvCell1DNeighboursFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours =
        utilities.ImportCell2DNeighbours(csvCell2DNeighboursFile, configuration.Separator);

    Gedim::Profiler::StopTime("ReadNeighbours");
    Gedim::Profiler::StartTime("ConvertNeighbours");
    utilities.ConvertCell0DNeighbours(cell0DNeighbours, mesh);
    utilities.ConvertCell1DNeighbours(cell1DNeighbours, mesh);
    utilities.ConvertCell2DNeighbours(cell2DNeighbours, mesh);
    Gedim::Profiler::StopTime("ConvertNeighbours");

    Gedim::Profiler::StartTime("ReadSubdivision");
    std::vector<MeshFromCsvUtilities::Cell2DSubDivision> cell2DSubDivisions =
        utilities.ImportCell2DSubDivision(csvCell2DSubDivisionsFile, configuration.Separator);
    Gedim::Profiler::StopTime("ReadSubdivision");
    Gedim::Profiler::StartTime("ConvertSubdivision");
    utilities.ConvertCell2DSubDivisions(cell2DSubDivisions, mesh);
    Gedim::Profiler::StopTime("ConvertSubdivision");

    Gedim::Profiler::StartTime("ReadUpdated");
    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell0DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell0DUpdatedCellsFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell1DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell1DUpdatedCellsFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell2DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell2DUpdatedCellsFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell3DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell3DUpdatedCellsFile, configuration.Separator);
    Gedim::Profiler::StopTime("ReadUpdated");

    Gedim::Profiler::StartTime("ConvertUpdated");
    utilities.ConvertCell0DUpdatedCells(cell0DUpdatedCells, mesh);
    utilities.ConvertCell1DUpdatedCells(cell1DUpdatedCells, mesh);
    utilities.ConvertCell2DUpdatedCells(cell2DUpdatedCells, mesh);
    utilities.ConvertCell3DUpdatedCells(cell3DUpdatedCells, mesh);
    Gedim::Profiler::StopTime("ConvertUpdated");

    mesh.Compress();
}
// ***************************************************************************
void MeshDAOImporterFromCsv::ImportMesh2D(const MeshFromCsvUtilities::Configuration &configuration, IMeshDAO &mesh)
{
    Gedim::FileReader csvCell0DsFile(configuration.Folder + "/" + configuration.FileCell0DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DsFile(configuration.Folder + "/" + configuration.FileCell1DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DsFile(configuration.Folder + "/" + configuration.FileCell2DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DPropertiesFile(configuration.Folder + "/" + configuration.FileCell0DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DPropertiesFile(configuration.Folder + "/" + configuration.FileCell1DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DPropertiesFile(configuration.Folder + "/" + configuration.FileCell2DPropertiesName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DNeighboursFile(configuration.Folder + "/" + configuration.FileCell0DNeighboursName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DNeighboursFile(configuration.Folder + "/" + configuration.FileCell1DNeighboursName +
                                              "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DSubDivisionsFile(configuration.Folder + "/" + configuration.FileCell2DSubDivisionsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell0DUpdatedCellsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell1DUpdatedCellsName +
                                                "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DUpdatedCellsFile(configuration.Folder + "/" + configuration.FileCell2DUpdatedCellsName +
                                                "." + configuration.FileExtension);

    std::vector<MeshFromCsvUtilities::Cell0D> cell0Ds = utilities.ImportCell0Ds(csvCell0DsFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell1D> cell1Ds = utilities.ImportCell1Ds(csvCell1DsFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell2D> cell2Ds = utilities.ImportCell2Ds(csvCell2DsFile, configuration.Separator);

    utilities.ConvertMesh2D(cell0Ds, cell1Ds, cell2Ds, mesh);

    Output::Assert(mesh.Cell0DTotalNumber() > 0 && mesh.Cell1DTotalNumber() > 0 && mesh.Cell2DTotalNumber() > 0 &&
                   mesh.Cell3DTotalNumber() == 0);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell0DPropertiesFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell1DPropertiesFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties =
        utilities.ImportCellDoubleProperties(csvCell2DPropertiesFile, configuration.Separator);

    utilities.ConvertCell0DDoubleProperties(cell0DDoubleProperties, mesh);
    utilities.ConvertCell1DDoubleProperties(cell1DDoubleProperties, mesh);
    utilities.ConvertCell2DDoubleProperties(cell2DDoubleProperties, mesh);

    std::vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours =
        utilities.ImportCell0DNeighbours(csvCell0DNeighboursFile, configuration.Separator);
    std::vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours =
        utilities.ImportCell1DNeighbours(csvCell1DNeighboursFile, configuration.Separator);

    utilities.ConvertCell0DNeighbours(cell0DNeighbours, mesh);
    utilities.ConvertCell1DNeighbours(cell1DNeighbours, mesh);

    std::vector<MeshFromCsvUtilities::Cell2DSubDivision> cell2DSubDivisions =
        utilities.ImportCell2DSubDivision(csvCell2DSubDivisionsFile, configuration.Separator);
    utilities.ConvertCell2DSubDivisions(cell2DSubDivisions, mesh);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell0DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell0DUpdatedCellsFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell1DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell1DUpdatedCellsFile, configuration.Separator);

    std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell2DUpdatedCells =
        utilities.ImportCellUpdatedCells(csvCell2DUpdatedCellsFile, configuration.Separator);

    utilities.ConvertCell0DUpdatedCells(cell0DUpdatedCells, mesh);
    utilities.ConvertCell1DUpdatedCells(cell1DUpdatedCells, mesh);
    utilities.ConvertCell2DUpdatedCells(cell2DUpdatedCells, mesh);

    mesh.Compress();
}
// ***************************************************************************
} // namespace Gedim
