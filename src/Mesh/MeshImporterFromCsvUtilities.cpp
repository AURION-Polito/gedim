#include "MeshImporterFromCsvUtilities.hpp"
#include <iostream>
#include <fstream>
#include "FileTextReader.hpp"

namespace Gedim
{
  // ***************************************************************************
  MeshImporterFromCsvUtilities::MeshImporterFromCsvUtilities()
  {
  }
  MeshImporterFromCsvUtilities::~MeshImporterFromCsvUtilities()
  {
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell0D> MeshImporterFromCsvUtilities::ImportCell0Ds(IFileReader& csvFileReader,
                                                                                           const char& separator) const
  {
    vector<Cell0D> cell0Ds;

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

    return cell0Ds;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell1D> MeshImporterFromCsvUtilities::ImportCell1Ds(IFileReader& csvFileReader,
                                                                                           const char& separator) const
  {
    vector<Cell1D> cell1Ds;

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

    return cell1Ds;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell2D> MeshImporterFromCsvUtilities::ImportCell2Ds(IFileReader& csvFileReader,
                                                                                           const char& separator) const
  {
    vector<Cell2D> cell2Ds;

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

    return cell2Ds;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell3D> MeshImporterFromCsvUtilities::ImportCell3Ds(IFileReader& csvFileReader,
                                                                                           const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::Cell3D> cell3Ds;

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
        cell3Ds.resize(numCell3Ds);
        for (unsigned int c = 0; c < numCell3Ds; c++)
        {
          Cell3D& cell3D = cell3Ds[c];

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
        }
      }
    }

    return cell3Ds;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::CellDoubleProperty> MeshImporterFromCsvUtilities::ImportCellDoubleProperties(IFileReader& csvFileReader,
                                                                                                                    const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cellProperties;

    /// Import CellProperties
    {
      vector<string> cellsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cellProperties file");

      csvFileReader.GetAllLines(cellsLines);
      csvFileReader.Close();

      unsigned int numCellProperties = cellsLines.size() - 1;

      if (numCellProperties > 0)
      {
        cellProperties.resize(numCellProperties);
        for (unsigned int p = 0; p < numCellProperties; p++)
        {
          CellDoubleProperty& cellProperty = cellProperties[p];

          istringstream converter(cellsLines[p + 1]);

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

          string fileFolder, fileName, fileExtension;
          Gedim::Output::GetFilePath(csvFileReader.Path(), fileFolder, fileName, fileExtension);
          FileReader propertyFileReader(fileFolder + cellProperty.FilePath);

          cellProperty.Values = ImportCellProperty(propertyFileReader,
                                                   separator);
        }
      }
    }

    return cellProperties;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::CellDoubleProperty::Value> MeshImporterFromCsvUtilities::ImportCellProperty(IFileReader& csvFileReader, const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::CellDoubleProperty::Value> cellPropertyValues;

    /// Import CellProperty
    {
      vector<string> cellsLines;

      if (!csvFileReader.Open())
        throw runtime_error("Error on mesh cellProperty file");

      csvFileReader.GetAllLines(cellsLines);
      csvFileReader.Close();

      unsigned int numCellProperty = cellsLines.size() - 1;

      if (numCellProperty > 0)
      {
        cellPropertyValues.resize(numCellProperty);

        for (unsigned int p = 0; p < numCellProperty; p++)
        {
          CellDoubleProperty::Value& cellProperty = cellPropertyValues[p];

          istringstream converter(cellsLines[p + 1]);

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
        }
      }
    }

    return cellPropertyValues;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell0DNeighbours> MeshImporterFromCsvUtilities::ImportCell0DNeighbours(IFileReader& csvFileReader,
                                                                                                              const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours;

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
        cell0DNeighbours.resize(numCell0DNeighbours);

        for (unsigned int v = 0; v < numCell0DNeighbours; v++)
        {
          MeshImporterFromCsvUtilities::Cell0DNeighbours& cell0D = cell0DNeighbours[v];

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
        }
      }
    }

    return cell0DNeighbours;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell1DNeighbours> MeshImporterFromCsvUtilities::ImportCell1DNeighbours(IFileReader& csvFileReader,
                                                                                                              const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours;

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
        cell1DNeighbours.resize(numCell1DNeighbours);

        Cell1D cell1D;
        for (unsigned int e = 0; e < numCell1DNeighbours; e++)
        {
          MeshImporterFromCsvUtilities::Cell1DNeighbours& cell1D = cell1DNeighbours[e];

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
        }
      }
    }

    return cell1DNeighbours;
  }
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell2DNeighbours> MeshImporterFromCsvUtilities::ImportCell2DNeighbours(IFileReader& csvFileReader,
                                                                                                              const char& separator) const
  {
    vector<MeshImporterFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours;

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
        cell2DNeighbours.resize(numCell2DNeighbours);

        for (unsigned int f = 0; f < numCell2DNeighbours; f++)
        {
          MeshImporterFromCsvUtilities::Cell2DNeighbours& cell2D = cell2DNeighbours[f];

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
        }
      }
    }

    return cell2DNeighbours;
  }
  // ***************************************************************************
}
