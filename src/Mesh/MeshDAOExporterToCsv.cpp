#include "MeshDAOExporterToCsv.hpp"
#include <iostream>
#include <fstream>

namespace Gedim
{
  // ***************************************************************************
  MeshDAOExporterToCsv::MeshDAOExporterToCsv()
  {
  }
  MeshDAOExporterToCsv::~MeshDAOExporterToCsv()
  {
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell0Ds(const string& filePath,
                                           const char& separator,
                                           const IMeshDAO& mesh) const
  {
    /// Export Cell0Ds
    ofstream fileCell0Ds;

    fileCell0Ds.open(filePath);
    fileCell0Ds.precision(16);

    if (fileCell0Ds.fail())
      throw runtime_error("Error on mesh cell0Ds file");

    fileCell0Ds<< "Id"<< separator;
    fileCell0Ds<< "Marker"<< separator;
    fileCell0Ds<< "Active"<< separator;
    fileCell0Ds<< "X"<< separator;
    fileCell0Ds<< "Y"<< separator;
    fileCell0Ds<< "Z"<< endl;
    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      fileCell0Ds<< scientific<< mesh.Cell0DId(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DMarker(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DIsActive(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateX(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateY(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateZ(v)<< endl;
    }

    fileCell0Ds.close();
  }

  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell1Ds(const string& filePath,
                                           const char& separator,
                                           const IMeshDAO& mesh) const
  {
    /// Export Cell1Ds
    ofstream fileCell1Ds;

    fileCell1Ds.open(filePath);
    fileCell1Ds.precision(16);

    if (fileCell1Ds.fail())
      throw runtime_error("Error on mesh cell1Ds file");

    fileCell1Ds<< "Id"<< separator;
    fileCell1Ds<< "Marker"<< separator;
    fileCell1Ds<< "Active"<< separator;
    fileCell1Ds<< "Origin"<< separator;
    fileCell1Ds<< "End"<< endl;
    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      fileCell1Ds<< scientific<< mesh.Cell1DId(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DMarker(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DIsActive(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DOrigin(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DEnd(e)<< endl;
    }

    fileCell1Ds.close();
  }

  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell2Ds(const string& filePath,
                                           const char& separator,
                                           const IMeshDAO& mesh) const
  {
    /// Export Cell2Ds
    ofstream fileCell2Ds;

    fileCell2Ds.open(filePath);
    fileCell2Ds.precision(16);

    if (fileCell2Ds.fail())
      throw runtime_error("Error on mesh cell2Ds file");

    fileCell2Ds<< "Id"<< separator;
    fileCell2Ds<< "Marker"<< separator;
    fileCell2Ds<< "Active"<< separator;
    fileCell2Ds<< "NumVertices"<< separator;
    fileCell2Ds<< "Vertices"<< separator;
    fileCell2Ds<< "NumEdges"<< separator;
    fileCell2Ds<< "Edges"<< endl;
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
    {
      fileCell2Ds<< scientific<< mesh.Cell2DId(f)<< separator;
      fileCell2Ds<< scientific<< mesh.Cell2DMarker(f)<< separator;
      fileCell2Ds<< scientific<< mesh.Cell2DIsActive(f)<< separator;

      fileCell2Ds<< scientific<< mesh.Cell2DNumberVertices(f);
      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(f); v++)
        fileCell2Ds<< scientific<< separator<< mesh.Cell2DVertex(f, v);
      fileCell2Ds<< separator;

      fileCell2Ds<< scientific<< mesh.Cell2DNumberEdges(f);
      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(f); e++)
        fileCell2Ds<< scientific<< separator<< mesh.Cell2DEdge(f, e);
      fileCell2Ds<< endl;
    }

    fileCell2Ds.close();
  }

  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell3Ds(const string& filePath,
                                           const char& separator,
                                           const IMeshDAO& mesh) const
  {
    /// Export Cell3Ds
    ofstream fileCell3Ds;

    fileCell3Ds.open(filePath);
    fileCell3Ds.precision(16);

    if (fileCell3Ds.fail())
      throw runtime_error("Error on mesh cell3Ds file");

    fileCell3Ds<< "Id"<< separator;
    fileCell3Ds<< "Marker"<< separator;
    fileCell3Ds<< "Active"<< separator;
    fileCell3Ds<< "NumVertices"<< separator;
    fileCell3Ds<< "Vertices"<< separator;
    fileCell3Ds<< "NumEdges"<< separator;
    fileCell3Ds<< "Edges"<< separator;
    fileCell3Ds<< "NumFaces"<< separator;
    fileCell3Ds<< "Faces"<< endl;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
      fileCell3Ds<< scientific<< mesh.Cell3DId(c)<< separator;
      fileCell3Ds<< scientific<< mesh.Cell3DMarker(c)<< separator;
      fileCell3Ds<< scientific<< mesh.Cell3DIsActive(c)<< separator;

      fileCell3Ds<< scientific<< mesh.Cell3DNumberVertices(c);
      for (unsigned int v = 0; v < mesh.Cell3DNumberVertices(c); v++)
        fileCell3Ds<< scientific<< separator<< mesh.Cell3DVertex(c, v);
      fileCell3Ds<< separator;

      fileCell3Ds<< scientific<< mesh.Cell3DNumberEdges(c);
      for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(c); e++)
        fileCell3Ds<< scientific<< separator<< mesh.Cell3DEdge(c, e);
      fileCell3Ds<< separator;

      fileCell3Ds<< scientific<< mesh.Cell3DNumberFaces(c);
      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(c); f++)
        fileCell3Ds<< scientific<< separator<< mesh.Cell3DFace(c, f);
      fileCell3Ds<< endl;
    }

    fileCell3Ds.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell0DProperties(const string& exportFolder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell0D Properties
    ofstream fileCell0DProperties;

    fileCell0DProperties.open(exportFolder + "/" + propertyFileName + "." + propertyFileExtension);
    fileCell0DProperties.precision(16);

    if (fileCell0DProperties.fail())
      throw runtime_error("Error on mesh cell0DProperties file");

    fileCell0DProperties<< "Id"<< separator;
    fileCell0DProperties<< "FilePath"<< endl;
    for (unsigned int p = 0; p < mesh.Cell0DNumberDoubleProperties(); p++)
    {
      const string propertyId = mesh.Cell0DDoublePropertyId(p);
      string propertyFilePath = propertyFileName + "_" +
                                propertyId + "." +
                                propertyFileExtension;

      fileCell0DProperties<< scientific<< propertyId<< separator;
      fileCell0DProperties<< scientific<< propertyFilePath<< endl;

      ExportCell0DProperty(p,
                           exportFolder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell0DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell0DProperty(const unsigned int& propertyIndex,
                                                  const string& filePath,
                                                  const char& separator,
                                                  const IMeshDAO& mesh) const
  {
    /// Export Cell0D Properties
    ofstream fileCell0DProperties;

    fileCell0DProperties.open(filePath);
    fileCell0DProperties.precision(16);

    if (fileCell0DProperties.fail())
      throw runtime_error("Error on mesh cell0DProperties file");

    fileCell0DProperties<< "Id"<< separator;
    fileCell0DProperties<< "PropertySize"<< separator;
    fileCell0DProperties<< "PropertyValues"<< endl;
    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      fileCell0DProperties<< scientific<< mesh.Cell0DId(v)<< separator;

      fileCell0DProperties<< scientific<< mesh.Cell0DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell0DDoublePropertySize(v, propertyIndex); n++)
        fileCell0DProperties<< scientific<< separator<< mesh.Cell0DDoublePropertyValue(v, propertyIndex, n);
      fileCell0DProperties<< endl;
    }

    fileCell0DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell1DProperties(const string& exportFolder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell1D Properties
    ofstream fileCell1DProperties;

    fileCell1DProperties.open(exportFolder + "/" + propertyFileName + "." + propertyFileExtension);
    fileCell1DProperties.precision(16);

    if (fileCell1DProperties.fail())
      throw runtime_error("Error on mesh cell1DProperties file");

    fileCell1DProperties<< "Id"<< separator;
    fileCell1DProperties<< "FilePath"<< endl;
    for (unsigned int p = 0; p < mesh.Cell1DNumberDoubleProperties(); p++)
    {
      const string propertyId = mesh.Cell1DDoublePropertyId(p);
      string propertyFilePath = propertyFileName + "_" +
                                propertyId + "." +
                                propertyFileExtension;

      fileCell1DProperties<< scientific<< propertyId<< separator;
      fileCell1DProperties<< scientific<< propertyFilePath<< endl;

      ExportCell1DProperty(p,
                           exportFolder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell1DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell1DProperty(const unsigned int& propertyIndex,
                                                  const string& filePath,
                                                  const char& separator,
                                                  const IMeshDAO& mesh) const
  {
    /// Export Cell1D Properties
    ofstream fileCell1DProperties;

    fileCell1DProperties.open(filePath);
    fileCell1DProperties.precision(16);

    if (fileCell1DProperties.fail())
      throw runtime_error("Error on mesh cell1DProperties file");

    fileCell1DProperties<< "Id"<< separator;
    fileCell1DProperties<< "PropertySize"<< separator;
    fileCell1DProperties<< "PropertyValues"<< endl;
    for (unsigned int v = 0; v < mesh.Cell1DTotalNumber(); v++)
    {
      fileCell1DProperties<< scientific<< mesh.Cell1DId(v)<< separator;

      fileCell1DProperties<< scientific<< mesh.Cell1DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell1DDoublePropertySize(v, propertyIndex); n++)
        fileCell1DProperties<< scientific<< separator<< mesh.Cell1DDoublePropertyValue(v, propertyIndex, n);
      fileCell1DProperties<< endl;
    }

    fileCell1DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell2DProperties(const string& exportFolder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell2D Properties
    ofstream fileCell2DProperties;

    fileCell2DProperties.open(exportFolder + "/" + propertyFileName + "." + propertyFileExtension);
    fileCell2DProperties.precision(16);

    if (fileCell2DProperties.fail())
      throw runtime_error("Error on mesh cell2DProperties file");

    fileCell2DProperties<< "Id"<< separator;
    fileCell2DProperties<< "FilePath"<< endl;
    for (unsigned int p = 0; p < mesh.Cell2DNumberDoubleProperties(); p++)
    {
      const string propertyId = mesh.Cell2DDoublePropertyId(p);
      string propertyFilePath = propertyFileName + "_" +
                                propertyId + "." +
                                propertyFileExtension;

      fileCell2DProperties<< scientific<< propertyId<< separator;
      fileCell2DProperties<< scientific<< propertyFilePath<< endl;

      ExportCell2DProperty(p,
                           exportFolder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell2DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell2DProperty(const unsigned int& propertyIndex,
                                                  const string& filePath,
                                                  const char& separator,
                                                  const IMeshDAO& mesh) const
  {
    /// Export Cell2D Properties
    ofstream fileCell2DProperties;

    fileCell2DProperties.open(filePath);
    fileCell2DProperties.precision(16);

    if (fileCell2DProperties.fail())
      throw runtime_error("Error on mesh cell2DProperties file");

    fileCell2DProperties<< "Id"<< separator;
    fileCell2DProperties<< "PropertySize"<< separator;
    fileCell2DProperties<< "PropertyValues"<< endl;
    for (unsigned int v = 0; v < mesh.Cell2DTotalNumber(); v++)
    {
      fileCell2DProperties<< scientific<< mesh.Cell2DId(v)<< separator;

      fileCell2DProperties<< scientific<< mesh.Cell2DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell2DDoublePropertySize(v, propertyIndex); n++)
        fileCell2DProperties<< scientific<< separator<< mesh.Cell2DDoublePropertyValue(v, propertyIndex, n);
      fileCell2DProperties<< endl;
    }

    fileCell2DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell3DProperties(const string& exportFolder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell3D Properties
    ofstream fileCell3DProperties;

    fileCell3DProperties.open(exportFolder + "/" + propertyFileName + "." + propertyFileExtension);
    fileCell3DProperties.precision(16);

    if (fileCell3DProperties.fail())
      throw runtime_error("Error on mesh cell3DProperties file");

    fileCell3DProperties<< "Id"<< separator;
    fileCell3DProperties<< "FilePath"<< endl;
    for (unsigned int p = 0; p < mesh.Cell3DNumberDoubleProperties(); p++)
    {
      const string propertyId = mesh.Cell3DDoublePropertyId(p);
      string propertyFilePath = propertyFileName + "_" +
                                propertyId + "." +
                                propertyFileExtension;

      fileCell3DProperties<< scientific<< propertyId<< separator;
      fileCell3DProperties<< scientific<< propertyFilePath<< endl;

      ExportCell3DProperty(p,
                           exportFolder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell3DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell3DProperty(const unsigned int& propertyIndex,
                                                  const string& filePath,
                                                  const char& separator,
                                                  const IMeshDAO& mesh) const
  {
    /// Export Cell3D Properties
    ofstream fileCell3DProperties;

    fileCell3DProperties.open(filePath);
    fileCell3DProperties.precision(16);

    if (fileCell3DProperties.fail())
      throw runtime_error("Error on mesh cell3DProperties file");

    fileCell3DProperties<< "Id"<< separator;
    fileCell3DProperties<< "PropertySize"<< separator;
    fileCell3DProperties<< "PropertyValues"<< endl;
    for (unsigned int v = 0; v < mesh.Cell3DTotalNumber(); v++)
    {
      fileCell3DProperties<< scientific<< mesh.Cell3DId(v)<< separator;

      fileCell3DProperties<< scientific<< mesh.Cell3DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell3DDoublePropertySize(v, propertyIndex); n++)
        fileCell3DProperties<< scientific<< separator<< mesh.Cell3DDoublePropertyValue(v, propertyIndex, n);
      fileCell3DProperties<< endl;
    }

    fileCell3DProperties.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell0DNeighbours(const string& filePath,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell0D Neigbours
    ofstream fileCell0DNeighbours;

    fileCell0DNeighbours.open(filePath);
    fileCell0DNeighbours.precision(16);

    if (fileCell0DNeighbours.fail())
      throw runtime_error("Error on mesh cell0DNeighbours file");

    fileCell0DNeighbours<< "Id"<< separator;
    fileCell0DNeighbours<< "Num1DNeighbours"<< separator;
    fileCell0DNeighbours<< "1DNeighbours"<< separator;
    fileCell0DNeighbours<< "Num2DNeighbours"<< separator;
    fileCell0DNeighbours<< "2DNeighbours"<< separator;
    fileCell0DNeighbours<< "Num3DNeighbours"<< separator;
    fileCell0DNeighbours<< "3DNeighbours"<< endl;
    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      fileCell0DNeighbours<< scientific<< mesh.Cell0DId(v)<< separator;

      fileCell0DNeighbours<< scientific<< mesh.Cell0DNumberNeighbourCell1D(v);
      for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell1D(v); n++)
        fileCell0DNeighbours<< scientific<< separator<< mesh.Cell0DNeighbourCell1D(v, n);
      fileCell0DNeighbours<< separator;

      fileCell0DNeighbours<< scientific<< mesh.Cell0DNumberNeighbourCell2D(v);
      for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell2D(v); n++)
        fileCell0DNeighbours<< scientific<< separator<< mesh.Cell0DNeighbourCell2D(v, n);
      fileCell0DNeighbours<< separator;

      fileCell0DNeighbours<< scientific<< mesh.Cell0DNumberNeighbourCell3D(v);
      for (unsigned int n = 0; n < mesh.Cell0DNumberNeighbourCell3D(v); n++)
        fileCell0DNeighbours<< scientific<< separator<< mesh.Cell0DNeighbourCell3D(v, n);
      fileCell0DNeighbours<< endl;
    }

    fileCell0DNeighbours.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell1DNeighbours(const string& filePath,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell1D Neigbours
    ofstream fileCell1DNeighbours;

    fileCell1DNeighbours.open(filePath);
    fileCell1DNeighbours.precision(16);

    if (fileCell1DNeighbours.fail())
      throw runtime_error("Error on mesh cell1DNeighbours file");

    fileCell1DNeighbours<< "Id"<< separator;
    fileCell1DNeighbours<< "Num2DNeighbours"<< separator;
    fileCell1DNeighbours<< "2DNeighbours"<< separator;
    fileCell1DNeighbours<< "Num3DNeighbours"<< separator;
    fileCell1DNeighbours<< "3DNeighbours"<< endl;
    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      fileCell1DNeighbours<< scientific<< mesh.Cell1DId(e)<< separator;

      fileCell1DNeighbours<< scientific<< mesh.Cell1DNumberNeighbourCell2D(e);
      for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(e); n++)
        fileCell1DNeighbours<< scientific<< separator<< mesh.Cell1DNeighbourCell2D(e, n);
      fileCell1DNeighbours<< separator;

      fileCell1DNeighbours<< scientific<< mesh.Cell1DNumberNeighbourCell3D(e);
      for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell3D(e); n++)
        fileCell1DNeighbours<< scientific<< separator<< mesh.Cell1DNeighbourCell3D(e, n);
      fileCell1DNeighbours<< endl;
    }

    fileCell1DNeighbours.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::ExportCell2DNeighbours(const string& filePath,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell2D Neigbours
    ofstream fileCell2DNeighbours;

    fileCell2DNeighbours.open(filePath);
    fileCell2DNeighbours.precision(16);

    if (fileCell2DNeighbours.fail())
      throw runtime_error("Error on mesh cell1DNeighbours file");

    fileCell2DNeighbours<< "Id"<< separator;
    fileCell2DNeighbours<< "Num3DNeighbours"<< separator;
    fileCell2DNeighbours<< "3DNeighbours"<< endl;
    for (unsigned int f = 0; f < mesh.Cell3DTotalNumber(); f++)
    {
      fileCell2DNeighbours<< scientific<< mesh.Cell2DId(f)<< separator;

      fileCell2DNeighbours<< scientific<< mesh.Cell2DNumberNeighbourCell3D(f);
      for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(f); n++)
        fileCell2DNeighbours<< scientific<< separator<< mesh.Cell2DNeighbourCell3D(f, n);
      fileCell2DNeighbours<< endl;
    }

    fileCell2DNeighbours.close();
  }
  // ***************************************************************************
  void MeshDAOExporterToCsv::Export(const MeshDAOExporterToCsv::Configuration& configuration,
                                    const IMeshDAO& mesh) const
  {

    ExportCell0Ds(configuration.ExportFolder + "/" +
                  configuration.FileCell0DsName + "." +
                  configuration.FileExtension,
                  configuration.Separator,
                  mesh);
    ExportCell1Ds(configuration.ExportFolder + "/" +
                  configuration.FileCell1DsName + "." +
                  configuration.FileExtension,
                  configuration.Separator,
                  mesh);
    ExportCell2Ds(configuration.ExportFolder + "/" +
                  configuration.FileCell2DsName + "." +
                  configuration.FileExtension,
                  configuration.Separator,
                  mesh);
    ExportCell3Ds(configuration.ExportFolder + "/" +
                  configuration.FileCell3DsName + "." +
                  configuration.FileExtension,
                  configuration.Separator,
                  mesh);

    ExportCell0DProperties(configuration.ExportFolder,
                           configuration.FileCell0DPropertiesName,
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);

    ExportCell1DProperties(configuration.ExportFolder,
                           configuration.FileCell1DPropertiesName,
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);

    ExportCell2DProperties(configuration.ExportFolder,
                           configuration.FileCell2DPropertiesName,
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);

    ExportCell3DProperties(configuration.ExportFolder,
                           configuration.FileCell3DPropertiesName,
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);

    ExportCell0DNeighbours(configuration.ExportFolder + "/" +
                           configuration.FileCell0DNeighboursName + "." +
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);
    ExportCell1DNeighbours(configuration.ExportFolder + "/" +
                           configuration.FileCell1DNeighboursName + "." +
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);
    ExportCell2DNeighbours(configuration.ExportFolder + "/" +
                           configuration.FileCell2DNeighboursName + "." +
                           configuration.FileExtension,
                           configuration.Separator,
                           mesh);


  }
  // ***************************************************************************
}
