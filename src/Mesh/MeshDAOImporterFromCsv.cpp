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
  void MeshDAOImporterFromCsv::ImportCell0Ds(IFileReader& csvFileReader,
                                             const char& separator,
                                             IMeshDAO& mesh) const
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
        Cell0D cell0D;
        mesh.Cell0DsInitialize(numCell0Ds);
        for (unsigned int v = 0; v < numCell0Ds; v++)
        {
          istringstream converter(cell0DsLines[v + 1]);

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

          mesh.Cell0DSetId(v, cell0D.Id);
          mesh.Cell0DSetMarker(v, cell0D.Marker);
          mesh.Cell0DSetState(v, cell0D.Active);
          mesh.Cell0DInsertCoordinates(v, Vector3d(cell0D.X, cell0D.Y, cell0D.Z));
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell1Ds(IFileReader& csvFileReader,
                                             const char& separator,
                                             IMeshDAO& mesh) const
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
        Cell1D cell1D;
        mesh.Cell1DsInitialize(numCell1Ds);
        for (unsigned int e = 0; e < numCell1Ds; e++)
        {
          istringstream converter(cell1DsLines[e + 1]);

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

          mesh.Cell1DSetId(e, cell1D.Id);
          mesh.Cell1DSetMarker(e, cell1D.Marker);
          mesh.Cell1DSetState(e, cell1D.Active);
          mesh.Cell1DInsertExtremes(e,
                                    cell1D.Origin,
                                    cell1D.End);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell2Ds(IFileReader& csvFileReader,
                                             const char& separator,
                                             IMeshDAO& mesh) const
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
        Cell2D cell2D;

        mesh.Cell2DsInitialize(numCell2Ds);
        for (unsigned int f = 0; f < numCell2Ds; f++)
        {
          istringstream converter(cell2DsLines[f + 1]);

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

          mesh.Cell2DSetId(f, cell2D.Id);
          mesh.Cell2DSetMarker(f, cell2D.Marker);
          mesh.Cell2DSetState(f, cell2D.Active);
          mesh.Cell2DInitializeVertices(f, numCellVertices);
          for (unsigned int v = 0; v < numCellVertices; v++)
            mesh.Cell2DInsertVertex(f, v, cell2D.Vertices[v]);

          mesh.Cell2DInitializeEdges(f, numCellEdges);
          for (unsigned int e = 0; e < numCellEdges; e++)
            mesh.Cell2DInsertEdge(f, e, cell2D.Edges[e]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell3Ds(IFileReader& csvFileReader,
                                             const char& separator,
                                             IMeshDAO& mesh) const
  {
    /// Import Cell3Ds
    {
      vector<string> cell3DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell3Ds file");

      csvFileReader.GetAllLines(cell3DsLines);
      csvFileReader.Close();

      unsigned int numCell3Ds = cell3DsLines.size() - 1;

      if (numCell3Ds > 0)
      {
        Cell3D cell3D;

        mesh.Cell3DsInitialize(numCell3Ds);
        for (unsigned int c = 0; c < numCell3Ds; c++)
        {
          istringstream converter(cell3DsLines[c + 1]);

          char temp;
          converter >> cell3D.Id;
          if (separator != ' ')
            converter >> temp;

          converter >> cell3D.Marker;
          if (separator != ' ')
            converter >> temp;

          converter >> cell3D.Active;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCellVertices;
          converter >> numCellVertices;
          if (separator != ' ')
            converter >> temp;

          cell3D.Vertices.resize(numCellVertices);
          for (unsigned int v = 0; v < numCellVertices; v++)
          {
            converter >> cell3D.Vertices[v];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCellEdges;
          converter >> numCellEdges;
          if (separator != ' ')
            converter >> temp;

          cell3D.Edges.resize(numCellEdges);
          for (unsigned int e = 0; e < numCellEdges; e++)
          {
            converter >> cell3D.Edges[e];
            if (separator != ' ')
              converter >> temp;
          }

          unsigned int numCellFaces;
          converter >> numCellFaces;
          if (separator != ' ')
            converter >> temp;

          cell3D.Faces.resize(numCellFaces);
          for (unsigned int f = 0; f < numCellFaces; f++)
          {
            converter >> cell3D.Faces[f];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell3DSetId(c, cell3D.Id);
          mesh.Cell3DSetMarker(c, cell3D.Marker);
          mesh.Cell3DSetState(c, cell3D.Active);
          mesh.Cell3DInitializeVertices(c, numCellVertices);
          for (unsigned int v = 0; v < numCellVertices; v++)
            mesh.Cell3DInsertVertex(c, v, cell3D.Vertices[v]);

          mesh.Cell3DInitializeEdges(c, numCellEdges);
          for (unsigned int e = 0; e < numCellEdges; e++)
            mesh.Cell3DInsertEdge(c, e, cell3D.Edges[e]);

          mesh.Cell3DInitializeFaces(c, numCellFaces);
          for (unsigned int f = 0; f < numCellFaces; f++)
            mesh.Cell3DInsertVertex(c, f, cell3D.Faces[f]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell0DProperties(IFileReader& csvFileReader,
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

          string fileFolder, fileName, fileExtension;
          Gedim::Output::GetFilePath(csvFileReader.Path(), fileFolder, fileName, fileExtension);
          FileReader propertyFileReader(fileFolder + cellProperty.FilePath);
          ImportCell0DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell0DProperty(const unsigned int& propertyIndex,
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
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.Values[v];
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
  void MeshDAOImporterFromCsv::ImportCell1DProperties(IFileReader& csvFileReader,
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

          string fileFolder, fileName, fileExtension;
          Gedim::Output::GetFilePath(csvFileReader.Path(), fileFolder, fileName, fileExtension);
          FileReader propertyFileReader(fileFolder + cellProperty.FilePath);
          ImportCell1DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell1DProperty(const unsigned int& propertyIndex,
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
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.Values[v];
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
  void MeshDAOImporterFromCsv::ImportCell2DProperties(IFileReader& csvFileReader,
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

          string fileFolder, fileName, fileExtension;
          Gedim::Output::GetFilePath(csvFileReader.Path(), fileFolder, fileName, fileExtension);
          FileReader propertyFileReader(fileFolder + cellProperty.FilePath);
          ImportCell2DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell2DProperty(const unsigned int& propertyIndex,
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
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.Values[v];
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
  void MeshDAOImporterFromCsv::ImportCell3DProperties(IFileReader& csvFileReader,
                                                      const char& separator,
                                                      IMeshDAO& mesh) const
  {
    /// Import Cell3DProperties
    {
      vector<string> cell3DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell3DProperties file");

      csvFileReader.GetAllLines(cell3DsLines);
      csvFileReader.Close();

      unsigned int numCell3DProperties = cell3DsLines.size() - 1;

      if (numCell3DProperties > 0)
      {
        CellProperty cellProperty;
        mesh.Cell3DInitializeDoubleProperties(numCell3DProperties);
        for (unsigned int p = 0; p < numCell3DProperties; p++)
        {
          istringstream converter(cell3DsLines[p + 1]);

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

          unsigned int pIndex = mesh.Cell3DAddDoubleProperty(cellProperty.Id);

          string fileFolder, fileName, fileExtension;
          Gedim::Output::GetFilePath(csvFileReader.Path(), fileFolder, fileName, fileExtension);
          FileReader propertyFileReader(fileFolder + cellProperty.FilePath);
          ImportCell3DProperty(pIndex,
                               propertyFileReader,
                               separator,
                               mesh);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell3DProperty(const unsigned int& propertyIndex,
                                                    IFileReader& csvFileReader,
                                                    const char& separator,
                                                    IMeshDAO& mesh) const
  {
    /// Import Cell3DProperty
    {
      vector<string> cell3DsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell3DProperty file");

      csvFileReader.GetAllLines(cell3DsLines);
      csvFileReader.Close();

      unsigned int numCell3DProperty = cell3DsLines.size() - 1;

      if (numCell3DProperty > 0)
      {
        CellProperty::Value cellProperty;
        mesh.Cell3DInitializeDoubleProperties(numCell3DProperty);
        for (unsigned int p = 0; p < numCell3DProperty; p++)
        {
          istringstream converter(cell3DsLines[p + 1]);

          char temp;
          converter >> cellProperty.CellId;
          if (separator != ' ')
            converter >> temp;
          unsigned int numValues;
          converter >> numValues;
          cellProperty.Values.resize(numValues);
          for (unsigned int v = 0; v < numValues; v++)
          {
            if (separator != ' ')
              converter >> temp;
            converter >> cellProperty.Values[v];
          }

          mesh.Cell3DInitializeDoublePropertyValues(cellProperty.CellId,
                                                    propertyIndex,
                                                    numValues);
          for (unsigned int v = 0; v < numValues; v++)
            mesh.Cell3DInsertDoublePropertyValue(cellProperty.CellId,
                                                 propertyIndex,
                                                 v,
                                                 cellProperty.Values[v]);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell0DNeighbours(IFileReader& csvFileReader,
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
            if (cell0D.Cell1DNeighbours[n] >= mesh.Cell1DTotalNumber())
              continue;

            mesh.Cell0DInsertNeighbourCell1D(v, n, cell0D.Cell1DNeighbours[n]);
          }

          mesh.Cell0DInitializeNeighbourCell2Ds(v, numCell2DNeighbours);
          for (unsigned int n = 0; n < numCell2DNeighbours; n++)
          {
            if (cell0D.Cell2DNeighbours[n] >= mesh.Cell2DTotalNumber())
              continue;

            mesh.Cell0DInsertNeighbourCell2D(v, n, cell0D.Cell2DNeighbours[n]);
          }

          mesh.Cell0DInitializeNeighbourCell3Ds(v, numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            if (cell0D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
              continue;

            mesh.Cell0DInsertNeighbourCell3D(v, n, cell0D.Cell3DNeighbours[n]);
          }

        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell1DNeighbours(IFileReader& csvFileReader,
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
            if (cell1D.Cell2DNeighbours[n] >= mesh.Cell2DTotalNumber())
              continue;

            mesh.Cell1DInsertNeighbourCell2D(e, n, cell1D.Cell2DNeighbours[n]);
          }

          mesh.Cell1DInitializeNeighbourCell3Ds(e, numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            if (cell1D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
              continue;

            mesh.Cell1DInsertNeighbourCell3D(e, n, cell1D.Cell3DNeighbours[n]);
          }
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::ImportCell2DNeighbours(IFileReader& csvFileReader,
                                                      const char& separator,
                                                      IMeshDAO& mesh) const
  {
    /// Import Cell2DNeighbours
    {
      vector<string> cell2DNeighboursLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cell2DNeighbours file");

      csvFileReader.GetAllLines(cell2DNeighboursLines);
      csvFileReader.Close();

      unsigned int numCell2DNeighbours = cell2DNeighboursLines.size() - 1;

      if (numCell2DNeighbours > 0)
      {
        Cell2D cell2D;
        for (unsigned int f = 0; f < numCell2DNeighbours; f++)
        {
          istringstream converter(cell2DNeighboursLines[f + 1]);

          unsigned int cell2Did;
          char temp;
          converter >> cell2Did;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCell3DNeighbours;
          converter >> numCell3DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell2D.Cell3DNeighbours.resize(numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            converter >> cell2D.Cell3DNeighbours[n];
            if (separator != ' ')
              converter >> temp;
          }

          mesh.Cell2DInitializeNeighbourCell3Ds(f, numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            if (cell2D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
              continue;

            mesh.Cell2DInsertNeighbourCell3D(f, n, cell2D.Cell3DNeighbours[n]);
          }
        }
      }
    }
  }
  // ***************************************************************************
  void MeshDAOImporterFromCsv::Import(Configuration& configuration,
                                      IMeshDAO& mesh)
  {
    ImportCell0Ds(configuration.CsvCell0DsFile,
                  configuration.Separator,
                  mesh);
    ImportCell1Ds(configuration.CsvCell1DsFile,
                  configuration.Separator,
                  mesh);
    ImportCell2Ds(configuration.CsvCell2DsFile,
                  configuration.Separator,
                  mesh);
    ImportCell3Ds(configuration.CsvCell3DsFile,
                  configuration.Separator,
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

    ImportCell0DProperties(configuration.CsvCell0DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell1DProperties(configuration.CsvCell1DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell2DProperties(configuration.CsvCell2DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell3DProperties(configuration.CsvCell3DPropertiesFile,
                           configuration.Separator,
                           mesh);

    ImportCell0DNeighbours(configuration.CsvCell0DNeighboursFile,
                           configuration.Separator,
                           mesh);
    ImportCell1DNeighbours(configuration.CsvCell1DNeighboursFile,
                           configuration.Separator,
                           mesh);
    ImportCell2DNeighbours(configuration.CsvCell2DNeighboursFile,
                           configuration.Separator,
                           mesh);

    mesh.Compress();
  }
  // ***************************************************************************
}
