#include "MeshDAOImporter2DFromCsv.hpp"
#include <iostream>
#include <fstream>
#include "FileTextReader.hpp"

namespace Gedim
{
  // ***************************************************************************
  MeshDAOImporter2DFromCsv::MeshDAOImporter2DFromCsv()
  {
  }
  MeshDAOImporter2DFromCsv::~MeshDAOImporter2DFromCsv()
  {
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell0Ds(IFileReader& csvFileReader,
                                               const char& separator,
                                               vector<Cell0D>& cell0Ds) const
  {
    /// Import Cell0Ds
    {
      vector<string> cell0DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell0Ds file");

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0Ds = cell0DsLines.size() - 1;

      if (numCell0Ds > 0)
      {
        cell0Ds.resize(numCell0Ds);
        for (unsigned int v = 0; v < numCell0Ds; v++)
        {
          istringstream converter(cell0DsLines[v + 1]);

          Cell0D& cell0D = cell0Ds[v];

          char temp;
          converter >> cell0D.Id;
          if (separator != ' ')
            converter >> temp;
          converter >> cell0D.Marker;
          if (separator != ' ')
            converter >> temp;
          converter >> cell0D.Active;
          if (separator != ' ')
            converter >> temp;
          converter >> cell0D.X;
          if (separator != ' ')
            converter >> temp;
          converter >> cell0D.Y;
          if (separator != ' ')
            converter >> temp;
          converter >> cell0D.Z;
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell1Ds(IFileReader& csvFileReader,
                                               const char& separator,
                                               vector<Cell1D>& cell1Ds) const
  {
    /// Import Cell1Ds
    {
      vector<string> cell1DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell1Ds file");

      csvFileReader.GetAllLines(cell1DsLines);
      csvFileReader.Close();

      unsigned int numCell1Ds = cell1DsLines.size() - 1;

      if (numCell1Ds > 0)
      {
        cell1Ds.resize(numCell1Ds);
        for (unsigned int e = 0; e < numCell1Ds; e++)
        {
          istringstream converter(cell1DsLines[e + 1]);

          Cell1D& cell1D = cell1Ds[e];

          char temp;
          converter >> cell1D.Id;
          if (separator != ' ')
            converter >> temp;
          converter >> cell1D.Marker;
          if (separator != ' ')
            converter >> temp;
          converter >> cell1D.Active;
          if (separator != ' ')
            converter >> temp;
          converter >> cell1D.Origin;
          if (separator != ' ')
            converter >> temp;
          converter >> cell1D.End;
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell2Ds(IFileReader& csvFileReader,
                                               const char& separator,
                                               vector<Cell2D>& cell2Ds) const
  {
    /// Import Cell2Ds
    {
      vector<string> cell2DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell2Ds file");

      csvFileReader.GetAllLines(cell2DsLines);
      csvFileReader.Close();

      unsigned int numCell2Ds = cell2DsLines.size() - 1;

      if (numCell2Ds > 0)
      {
        cell2Ds.resize(numCell2Ds);
        for (unsigned int f = 0; f < numCell2Ds; f++)
        {
          istringstream converter(cell2DsLines[f + 1]);

          Cell2D& cell2D = cell2Ds[f];

          char temp;
          converter >> cell2D.Id;
          if (separator != ' ')
            converter >> temp;

          converter >> cell2D.Marker;
          if (separator != ' ')
            converter >> temp;

          converter >> cell2D.Active;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCellVertices;
          converter >> numCellVertices;
          if (separator != ' ')
            converter >> temp;

          cell2D.Vertices.resize(numCellVertices);
          for (unsigned int v = 0; v < numCellVertices; v++)
          {
            converter >> cell2D.Vertices[v];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCellEdges;
          converter >> numCellEdges;
          if (separator != ' ')
            converter >> temp;

          cell2D.Edges.resize(numCellEdges);
          for (unsigned int e = 0; e < numCellEdges; e++)
          {
            converter >> cell2D.Edges[e];
            if (separator != ' ')
              converter >> temp;
          }
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::CreateMesh2D(const vector<Cell0D>& cell0Ds,
                                              const vector<Cell1D>& cell1Ds,
                                              const vector<Cell2D>& cell2Ds,
                                              IMeshDAO& mesh)
  {
    unsigned int numCell0Ds = cell0Ds.size();
    unsigned int numCell1Ds = cell1Ds.size();
    unsigned int numCell2Ds = cell2Ds.size();
    Eigen::MatrixXd meshCell0Ds(3, numCell0Ds);
    Eigen::MatrixXi meshCell1Ds(2, numCell1Ds);
    vector<Eigen::MatrixXi> meshCell2Ds(numCell2Ds);

    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const Cell0D& cell0D = cell0Ds[v];
      meshCell0Ds(0, v) = cell0D.X;
      meshCell0Ds(1, v) = cell0D.Y;
      meshCell0Ds(2, v) = cell0D.Z;
    }

    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      const Cell1D& cell1D = cell1Ds[e];
      meshCell1Ds(0, e) = cell1D.Origin;
      meshCell1Ds(1, e) = cell1D.End;
    }

    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const Cell2D& cell2D = cell2Ds[f];
      Output::Assert(cell2D.Vertices.size() == cell2D.Edges.size());
      const unsigned int numVertices = cell2D.Vertices.size();
      Eigen::MatrixXi& polygon = meshCell2Ds[f];
      polygon.resize(2,
                     numVertices);
      for (unsigned int v = 0; v < numVertices; v++)
      {
        polygon(0, v) = cell2D.Vertices[v];
        polygon(1, v) = cell2D.Edges[v];
      }
    }

    mesh.FillMesh2D(meshCell0Ds,
                    meshCell1Ds,
                    meshCell2Ds);

    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const Cell0D& cell0D = cell0Ds[v];
      mesh.Cell0DSetId(v, cell0D.Id);
      mesh.Cell0DSetMarker(v, cell0D.Marker);
      mesh.Cell0DSetState(v, cell0D.Active);
    }

    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      const Cell1D& cell1D = cell1Ds[e];
      mesh.Cell1DSetId(e, cell1D.Id);
      mesh.Cell1DSetMarker(e, cell1D.Marker);
      mesh.Cell1DSetState(e, cell1D.Active);
    }

    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const Cell2D& cell2D = cell2Ds[f];
      mesh.Cell2DSetId(f, cell2D.Id);
      mesh.Cell2DSetMarker(f, cell2D.Marker);
      mesh.Cell2DSetState(f, cell2D.Active);
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell0DProperties(IFileReader& csvFileReader,
                                                        const char& separator,
                                                        IMeshDAO& mesh) const
  {
    /// Import Cell0DProperties
    {
      vector<string> cell0DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell0DProperties file");

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0DProperties = cell0DsLines.size() - 1;

      if (numCell0DProperties > 0)
      {
        CellProperty cellProperty;
        mesh.Cell0DInitializeDoubleProperties(numCell0DProperties);
        for (unsigned int p = 0; p < numCell0DProperties; p++)
        {
          istringstream converter(cell0DsLines[p + 1]);

          if (separator == ' ')
          {
            char temp;
            converter >> cellProperty.Id;
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.FilePath;
          }
          else
          {
            string tempStr;
            converter >> tempStr;
            vector<string> result = Output::StringSplit(tempStr, separator);
            Output::Assert(result.size() == 2);

            cellProperty.Id = result[0];
            cellProperty.FilePath = result[1];
          }

          unsigned int pIndex = mesh.Cell0DAddDoubleProperty(cellProperty.Id);
          FileReader propertyFileReader(cellProperty.FilePath);
          ImportCell0DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell0DProperty(const unsigned int& propertyIndex,
                                                      IFileReader& csvFileReader,
                                                      const char& separator,
                                                      IMeshDAO& mesh) const
  {
    /// Import Cell0DProperty
    {
      vector<string> cell0DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell0DProperty file");

      csvFileReader.GetAllLines(cell0DsLines);
      csvFileReader.Close();

      unsigned int numCell0DProperty = cell0DsLines.size() - 1;

      if (numCell0DProperty > 0)
      {
        CellProperty::Value cellProperty;
        mesh.Cell0DInitializeDoubleProperties(numCell0DProperty);
        for (unsigned int p = 0; p < numCell0DProperty; p++)
        {
          istringstream converter(cell0DsLines[p + 1]);

          char temp;
          converter >> cellProperty.CellId;
          if (separator != ' ')
            converter >> temp;
          unsigned int numValues;
          converter >> numValues;
          cellProperty.Values.resize(numValues);
          for (unsigned int v = 0; v < numValues; v++)
          {
            converter >> cellProperty.Values[v];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell0DInitializeDoublePropertyValues(cellProperty.CellId,
                                                    propertyIndex,
                                                    numValues);
          for (unsigned int v = 0; v < numValues; v++)
            mesh.Cell0DInsertDoublePropertyValue(cellProperty.CellId,
                                                 propertyIndex,
                                                 v,
                                                 cellProperty.Values[v]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell1DProperties(IFileReader& csvFileReader,
                                                        const char& separator,
                                                        IMeshDAO& mesh) const
  {
    /// Import Cell1DProperties
    {
      vector<string> cell1DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell1DProperties file");

      csvFileReader.GetAllLines(cell1DsLines);
      csvFileReader.Close();

      unsigned int numCell1DProperties = cell1DsLines.size() - 1;

      if (numCell1DProperties > 0)
      {
        CellProperty cellProperty;
        mesh.Cell1DInitializeDoubleProperties(numCell1DProperties);
        for (unsigned int p = 0; p < numCell1DProperties; p++)
        {
          istringstream converter(cell1DsLines[p + 1]);

          if (separator == ' ')
          {
            char temp;
            converter >> cellProperty.Id;
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.FilePath;
          }
          else
          {
            string tempStr;
            converter >> tempStr;
            vector<string> result = Output::StringSplit(tempStr, separator);
            Output::Assert(result.size() == 2);

            cellProperty.Id = result[0];
            cellProperty.FilePath = result[1];
          }

          unsigned int pIndex = mesh.Cell1DAddDoubleProperty(cellProperty.Id);
          FileReader propertyFileReader(cellProperty.FilePath);
          ImportCell1DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell1DProperty(const unsigned int& propertyIndex,
                                                      IFileReader& csvFileReader,
                                                      const char& separator,
                                                      IMeshDAO& mesh) const
  {
    /// Import Cell1DProperty
    {
      vector<string> cell1DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell1DProperty file");

      csvFileReader.GetAllLines(cell1DsLines);
      csvFileReader.Close();

      unsigned int numCell1DProperty = cell1DsLines.size() - 1;

      if (numCell1DProperty > 0)
      {
        CellProperty::Value cellProperty;
        mesh.Cell1DInitializeDoubleProperties(numCell1DProperty);
        for (unsigned int p = 0; p < numCell1DProperty; p++)
        {
          istringstream converter(cell1DsLines[p + 1]);

          char temp;
          converter >> cellProperty.CellId;
          if (separator != ' ')
            converter >> temp;
          unsigned int numValues;
          converter >> numValues;
          cellProperty.Values.resize(numValues);
          for (unsigned int v = 0; v < numValues; v++)
          {
            converter >> cellProperty.Values[v];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell1DInitializeDoublePropertyValues(cellProperty.CellId,
                                                    propertyIndex,
                                                    numValues);
          for (unsigned int v = 0; v < numValues; v++)
            mesh.Cell1DInsertDoublePropertyValue(cellProperty.CellId,
                                                 propertyIndex,
                                                 v,
                                                 cellProperty.Values[v]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell2DProperties(IFileReader& csvFileReader,
                                                        const char& separator,
                                                        IMeshDAO& mesh) const
  {
    /// Import Cell2DProperties
    {
      vector<string> cell2DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell2DProperties file");

      csvFileReader.GetAllLines(cell2DsLines);
      csvFileReader.Close();

      unsigned int numCell2DProperties = cell2DsLines.size() - 1;

      if (numCell2DProperties > 0)
      {
        CellProperty cellProperty;
        mesh.Cell2DInitializeDoubleProperties(numCell2DProperties);
        for (unsigned int p = 0; p < numCell2DProperties; p++)
        {
          istringstream converter(cell2DsLines[p + 1]);

          if (separator == ' ')
          {
            char temp;
            converter >> cellProperty.Id;
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.FilePath;
          }
          else
          {
            string tempStr;
            converter >> tempStr;
            vector<string> result = Output::StringSplit(tempStr, separator);
            Output::Assert(result.size() == 2);

            cellProperty.Id = result[0];
            cellProperty.FilePath = result[1];
          }

          unsigned int pIndex = mesh.Cell2DAddDoubleProperty(cellProperty.Id);
          FileReader propertyFileReader(cellProperty.FilePath);
          ImportCell2DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell2DProperty(const unsigned int& propertyIndex,
                                                      IFileReader& csvFileReader,
                                                      const char& separator,
                                                      IMeshDAO& mesh) const
  {
    /// Import Cell2DProperty
    {
      vector<string> cell2DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell2DProperty file");

      csvFileReader.GetAllLines(cell2DsLines);
      csvFileReader.Close();

      unsigned int numCell2DProperty = cell2DsLines.size() - 1;

      if (numCell2DProperty > 0)
      {
        CellProperty::Value cellProperty;
        mesh.Cell2DInitializeDoubleProperties(numCell2DProperty);
        for (unsigned int p = 0; p < numCell2DProperty; p++)
        {
          istringstream converter(cell2DsLines[p + 1]);

          char temp;
          converter >> cellProperty.CellId;
          if (separator != ' ')
            converter >> temp;
          unsigned int numValues;
          converter >> numValues;
          cellProperty.Values.resize(numValues);
          for (unsigned int v = 0; v < numValues; v++)
          {
            converter >> cellProperty.Values[v];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell2DInitializeDoublePropertyValues(cellProperty.CellId,
                                                    propertyIndex,
                                                    numValues);
          for (unsigned int v = 0; v < numValues; v++)
            mesh.Cell2DInsertDoublePropertyValue(cellProperty.CellId,
                                                 propertyIndex,
                                                 v,
                                                 cellProperty.Values[v]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell0DNeighbours(IFileReader& csvFileReader,
                                                        const char& separator,
                                                        IMeshDAO& mesh) const
  {
    /// Import Cell0DNeighbours
    {
      vector<string> cell0DNeighboursLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell0DNeighbours file");

      csvFileReader.GetAllLines(cell0DNeighboursLines);
      csvFileReader.Close();

      unsigned int numCell0DNeighbours = cell0DNeighboursLines.size() - 1;

      if (numCell0DNeighbours > 0)
      {
        Cell0D cell0D;
        for (unsigned int v = 0; v < numCell0DNeighbours; v++)
        {
          istringstream converter(cell0DNeighboursLines[v + 1]);

          unsigned int cell0Did;
          char temp;
          converter >> cell0Did;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCell1DNeighbours;
          converter >> numCell1DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell0D.Cell1DNeighbours.resize(numCell1DNeighbours);
          for (unsigned int n = 0; n < numCell1DNeighbours; n++)
          {
            converter >> cell0D.Cell1DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCell2DNeighbours;
          converter >> numCell2DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell0D.Cell2DNeighbours.resize(numCell2DNeighbours);
          for (unsigned int n = 0; n < numCell2DNeighbours; n++)
          {
            converter >> cell0D.Cell2DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCell3DNeighbours;
          converter >> numCell3DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell0D.Cell3DNeighbours.resize(numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            converter >> cell0D.Cell3DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell0DInitializeNeighbourCell1Ds(v, numCell1DNeighbours);
          for (unsigned int n = 0; n < numCell1DNeighbours; n++)
          {
            if (!mesh.Cell0DHasNeighbourCell1D(v, n))
              continue;

            mesh.Cell0DInsertNeighbourCell1D(v, n, cell0D.Cell1DNeighbours[n]);
          }

          mesh.Cell0DInitializeNeighbourCell2Ds(v, numCell2DNeighbours);
          for (unsigned int n = 0; n < numCell2DNeighbours; n++)
          {
            if (!mesh.Cell0DHasNeighbourCell2D(v, n))
              continue;

            mesh.Cell0DInsertNeighbourCell2D(v, n, cell0D.Cell2DNeighbours[n]);
          }

          mesh.Cell0DInitializeNeighbourCell3Ds(v, numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            if (!mesh.Cell0DHasNeighbourCell3D(v, n))
              continue;

            mesh.Cell0DInsertNeighbourCell3D(v, n, cell0D.Cell3DNeighbours[n]);
          }

        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::ImportCell1DNeighbours(IFileReader& csvFileReader,
                                                        const char& separator,
                                                        IMeshDAO& mesh) const
  {
    /// Import Cell1DNeighbours
    {
      vector<string> cell1DNeighboursLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell1DNeighbours file");

      csvFileReader.GetAllLines(cell1DNeighboursLines);
      csvFileReader.Close();

      unsigned int numCell1DNeighbours = cell1DNeighboursLines.size() - 1;

      if (numCell1DNeighbours > 0)
      {
        Cell1D cell1D;
        for (unsigned int e = 0; e < numCell1DNeighbours; e++)
        {
          istringstream converter(cell1DNeighboursLines[e + 1]);

          unsigned int cell1Did;
          char temp;
          converter >> cell1Did;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCell2DNeighbours;
          converter >> numCell2DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell1D.Cell2DNeighbours.resize(numCell2DNeighbours);
          for (unsigned int n = 0; n < numCell2DNeighbours; n++)
          {
            converter >> cell1D.Cell2DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCell3DNeighbours;
          converter >> numCell3DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell1D.Cell3DNeighbours.resize(numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            converter >> cell1D.Cell3DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell1DInitializeNeighbourCell2Ds(e, numCell2DNeighbours);
          for (unsigned int n = 0; n < numCell2DNeighbours; n++)
          {
            if (!mesh.Cell1DHasNeighbourCell2D(e, n))
              continue;

            mesh.Cell1DInsertNeighbourCell2D(e, n, cell1D.Cell2DNeighbours[n]);
          }

          mesh.Cell1DInitializeNeighbourCell3Ds(e, numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            if (!mesh.Cell1DHasNeighbourCell3D(e, n))
              continue;

            mesh.Cell1DInsertNeighbourCell3D(e, n, cell1D.Cell3DNeighbours[n]);
          }
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporter2DFromCsv::Import(Configuration& configuration,
                                        IMeshDAO& mesh)
  {
    vector<Cell0D> cell0Ds;
    vector<Cell1D> cell1Ds;
    vector<Cell2D> cell2Ds;

    ImportCell0Ds(configuration.CsvCell0DsFile,
                  configuration.Separator,
                  cell0Ds);
    ImportCell1Ds(configuration.CsvCell1DsFile,
                  configuration.Separator,
                  cell1Ds);
    ImportCell2Ds(configuration.CsvCell2DsFile,
                  configuration.Separator,
                  cell2Ds);

    CreateMesh2D(cell0Ds,
                 cell1Ds,
                 cell2Ds,
                 mesh);

    Output::Assert(mesh.Cell0DTotalNumber() > 0 &&
                   mesh.Cell1DTotalNumber() > 0 &&
                   mesh.Cell2DTotalNumber() > 0 &&
                   mesh.Cell3DTotalNumber() == 0);

    ImportCell0DProperties(configuration.CsvCell0DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell1DProperties(configuration.CsvCell1DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell2DProperties(configuration.CsvCell2DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell0DNeighbours(configuration.CsvCell0DNeighboursFile,
                           configuration.Separator,
                           mesh);
    ImportCell1DNeighbours(configuration.CsvCell1DNeighboursFile,
                           configuration.Separator,
                           mesh);

    mesh.Compress();
  }
  // ***************************************************************************
}
