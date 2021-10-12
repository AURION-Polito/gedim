#include "MeshDAOImporterFromCsv.hpp"
#include <iostream>
#include <fstream>
#include "FileTextReader.hpp"

namespace Gedim
{
  // ***************************************************************************
  MeshDAOImporterFromCsv::MeshDAOImporterFromCsv()
  {
  }
  MeshDAOImporterFromCsv::~MeshDAOImporterFromCsv()
  {
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::Import(MeshImporterFromCsvUtilities::Configuration& configuration,
                                      MeshImporterFromCsvUtilities& importerUtilities,
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

    vector<MeshImporterFromCsvUtilities::Cell0D> cell0Ds = importerUtilities.ImportCell0Ds(csvCell0DsFile,
                                                                                           configuration.Separator);
    vector<MeshImporterFromCsvUtilities::Cell1D> cell1Ds = importerUtilities.ImportCell1Ds(csvCell1DsFile,
                                                                                           configuration.Separator);
    vector<MeshImporterFromCsvUtilities::Cell2D> cell2Ds = importerUtilities.ImportCell2Ds(csvCell2DsFile,
                                                                                           configuration.Separator);
    vector<MeshImporterFromCsvUtilities::Cell3D> cell3Ds = importerUtilities.ImportCell3Ds(csvCell3DsFile,
                                                                                           configuration.Separator);


    importerUtilities.ConvertCell0Ds(cell0Ds,
                                     mesh);
    importerUtilities.ConvertCell1Ds(cell1Ds,
                                     mesh);
    importerUtilities.ConvertCell2Ds(cell2Ds,
                                     mesh);
    importerUtilities.ConvertCell3Ds(cell3Ds,
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

    vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties = importerUtilities.ImportCellDoubleProperties(csvCell0DPropertiesFile,
                                                                                                                                   configuration.Separator);

    vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties = importerUtilities.ImportCellDoubleProperties(csvCell1DPropertiesFile,
                                                                                                                                   configuration.Separator);

    vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties = importerUtilities.ImportCellDoubleProperties(csvCell2DPropertiesFile,
                                                                                                                                   configuration.Separator);

    vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties = importerUtilities.ImportCellDoubleProperties(csvCell3DPropertiesFile,
                                                                                                                                   configuration.Separator);

    importerUtilities.ConvertCell0DDoubleProperties(cell0DDoubleProperties,
                                                    mesh);
    importerUtilities.ConvertCell1DDoubleProperties(cell1DDoubleProperties,
                                                    mesh);
    importerUtilities.ConvertCell2DDoubleProperties(cell2DDoubleProperties,
                                                    mesh);
    importerUtilities.ConvertCell3DDoubleProperties(cell3DDoubleProperties,
                                                    mesh);

    vector<MeshImporterFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours = importerUtilities.ImportCell0DNeighbours(csvCell0DNeighboursFile,
                                                                                                                       configuration.Separator);
    vector<MeshImporterFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours = importerUtilities.ImportCell1DNeighbours(csvCell1DNeighboursFile,
                                                                                                                       configuration.Separator);
    vector<MeshImporterFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours = importerUtilities.ImportCell2DNeighbours(csvCell2DNeighboursFile,
                                                                                                                       configuration.Separator);

    importerUtilities.ConvertCell0DNeighbours(cell0DNeighbours,
                                              mesh);
    importerUtilities.ConvertCell1DNeighbours(cell1DNeighbours,
                                              mesh);
    importerUtilities.ConvertCell2DNeighbours(cell2DNeighbours,
                                              mesh);

    mesh.Compress();
  }
  // ***************************************************************************
}
