#include "MeshDAOImporterFromCsv.hpp"
#include <iostream>
#include <fstream>
#include "FileTextReader.hpp"

namespace Gedim
{
  // ***************************************************************************
  MeshDAOImporterFromCsv::MeshDAOImporterFromCsv(const MeshFromCsvUtilities& utilities) :
    utilities(utilities)
  {
  }
  MeshDAOImporterFromCsv::~MeshDAOImporterFromCsv()
  {
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::Import(MeshFromCsvUtilities::Configuration& configuration,
                                      IMeshDAO& mesh)
  {
    Gedim::FileReader csvCell0DsFile(configuration.Folder + "/" + configuration.FileCell0DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DsFile(configuration.Folder + "/" + configuration.FileCell1DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DsFile(configuration.Folder + "/" + configuration.FileCell2DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell3DsFile(configuration.Folder + "/" + configuration.FileCell3DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DPropertiesFile(configuration.Folder + "/" + configuration.FileCell0DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DPropertiesFile(configuration.Folder + "/" + configuration.FileCell1DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DPropertiesFile(configuration.Folder + "/" + configuration.FileCell2DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell3DPropertiesFile(configuration.Folder + "/" + configuration.FileCell3DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DNeighboursFile(configuration.Folder + "/" + configuration.FileCell0DNeighboursName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DNeighboursFile(configuration.Folder + "/" + configuration.FileCell1DNeighboursName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DNeighboursFile(configuration.Folder + "/" + configuration.FileCell2DNeighboursName + "." + configuration.FileExtension);

    vector<MeshFromCsvUtilities::Cell0D> cell0Ds = utilities.ImportCell0Ds(csvCell0DsFile,
                                                                           configuration.Separator);
    vector<MeshFromCsvUtilities::Cell1D> cell1Ds = utilities.ImportCell1Ds(csvCell1DsFile,
                                                                           configuration.Separator);
    vector<MeshFromCsvUtilities::Cell2D> cell2Ds = utilities.ImportCell2Ds(csvCell2DsFile,
                                                                           configuration.Separator);
    vector<MeshFromCsvUtilities::Cell3D> cell3Ds = utilities.ImportCell3Ds(csvCell3DsFile,
                                                                           configuration.Separator);


    utilities.ConvertCell0Ds(cell0Ds,
                             mesh);
    utilities.ConvertCell1Ds(cell1Ds,
                             mesh);
    utilities.ConvertCell2Ds(cell2Ds,
                             mesh);
    utilities.ConvertCell3Ds(cell3Ds,
                             mesh);

    unsigned int meshDimension = 0;
    if (mesh.Cell0DTotalNumber() > 0 &&
        mesh.Cell1DTotalNumber() > 0 &&
        mesh.Cell2DTotalNumber() == 0 &&
        mesh.Cell3DTotalNumber() == 0)
      meshDimension = 1;
    else if (mesh.Cell0DTotalNumber() > 0 &&
             mesh.Cell1DTotalNumber() > 0 &&
             mesh.Cell2DTotalNumber() > 0 &&
             mesh.Cell3DTotalNumber() == 0)
      meshDimension = 2;
    else if (mesh.Cell0DTotalNumber() > 0 &&
             mesh.Cell1DTotalNumber() > 0 &&
             mesh.Cell2DTotalNumber() > 0 &&
             mesh.Cell3DTotalNumber() > 0)
      meshDimension = 3;
    else
      throw runtime_error("Dimension of imported mesh not recognized");

    mesh.InitializeDimension(meshDimension);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell0DPropertiesFile,
                                                                                                                   configuration.Separator);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell1DPropertiesFile,
                                                                                                                   configuration.Separator);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell2DPropertiesFile,
                                                                                                                   configuration.Separator);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell3DPropertiesFile,
                                                                                                                   configuration.Separator);

    utilities.ConvertCell0DDoubleProperties(cell0DDoubleProperties,
                                            mesh);
    utilities.ConvertCell1DDoubleProperties(cell1DDoubleProperties,
                                            mesh);
    utilities.ConvertCell2DDoubleProperties(cell2DDoubleProperties,
                                            mesh);
    utilities.ConvertCell3DDoubleProperties(cell3DDoubleProperties,
                                            mesh);

    vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours = utilities.ImportCell0DNeighbours(csvCell0DNeighboursFile,
                                                                                                       configuration.Separator);
    vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours = utilities.ImportCell1DNeighbours(csvCell1DNeighboursFile,
                                                                                                       configuration.Separator);
    vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours = utilities.ImportCell2DNeighbours(csvCell2DNeighboursFile,
                                                                                                       configuration.Separator);

    utilities.ConvertCell0DNeighbours(cell0DNeighbours,
                                      mesh);
    utilities.ConvertCell1DNeighbours(cell1DNeighbours,
                                      mesh);
    utilities.ConvertCell2DNeighbours(cell2DNeighbours,
                                      mesh);

    mesh.Compress();
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportMesh2D(MeshFromCsvUtilities::Configuration& configuration,
                                            IMeshDAO& mesh)
  {
    Gedim::FileReader csvCell0DsFile(configuration.Folder + "/" + configuration.FileCell0DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DsFile(configuration.Folder + "/" + configuration.FileCell1DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DsFile(configuration.Folder + "/" + configuration.FileCell2DsName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DPropertiesFile(configuration.Folder + "/" + configuration.FileCell0DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DPropertiesFile(configuration.Folder + "/" + configuration.FileCell1DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell2DPropertiesFile(configuration.Folder + "/" + configuration.FileCell2DPropertiesName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell0DNeighboursFile(configuration.Folder + "/" + configuration.FileCell0DNeighboursName + "." + configuration.FileExtension);
    Gedim::FileReader csvCell1DNeighboursFile(configuration.Folder + "/" + configuration.FileCell1DNeighboursName + "." + configuration.FileExtension);

    vector<MeshFromCsvUtilities::Cell0D> cell0Ds = utilities.ImportCell0Ds(csvCell0DsFile,
                                                                           configuration.Separator);
    vector<MeshFromCsvUtilities::Cell1D> cell1Ds = utilities.ImportCell1Ds(csvCell1DsFile,
                                                                           configuration.Separator);
    vector<MeshFromCsvUtilities::Cell2D> cell2Ds = utilities.ImportCell2Ds(csvCell2DsFile,
                                                                           configuration.Separator);

    utilities.ConvertMesh2D(cell0Ds,
                            cell1Ds,
                            cell2Ds,
                            mesh);

    Output::Assert(mesh.Cell0DTotalNumber() > 0 &&
                   mesh.Cell1DTotalNumber() > 0 &&
                   mesh.Cell2DTotalNumber() > 0 &&
                   mesh.Cell3DTotalNumber() == 0);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell0DPropertiesFile,
                                                                                                                   configuration.Separator);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell1DPropertiesFile,
                                                                                                                   configuration.Separator);

    vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties = utilities.ImportCellDoubleProperties(csvCell2DPropertiesFile,
                                                                                                                   configuration.Separator);

    utilities.ConvertCell0DDoubleProperties(cell0DDoubleProperties,
                                            mesh);
    utilities.ConvertCell1DDoubleProperties(cell1DDoubleProperties,
                                            mesh);
    utilities.ConvertCell2DDoubleProperties(cell2DDoubleProperties,
                                            mesh);

    vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours = utilities.ImportCell0DNeighbours(csvCell0DNeighboursFile,
                                                                                                       configuration.Separator);
    vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours = utilities.ImportCell1DNeighbours(csvCell1DNeighboursFile,
                                                                                                       configuration.Separator);

    utilities.ConvertCell0DNeighbours(cell0DNeighbours,
                                      mesh);
    utilities.ConvertCell1DNeighbours(cell1DNeighbours,
                                      mesh);

    mesh.Compress();
  }
  // ***************************************************************************
}
