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
  void MeshImporterFromCsvUtilities::ConvertCell0Ds(const vector<MeshImporterFromCsvUtilities::Cell0D> cell0Ds,
                                                    IMeshDAO& mesh) const
  {
    const unsigned int numCell0Ds = cell0Ds.size();

    if (numCell0Ds == 0)
      return;

    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const Cell0D& cell0D = cell0Ds[v];

      mesh.Cell0DSetId(v, cell0D.Id);
      mesh.Cell0DSetMarker(v, cell0D.Marker);
      mesh.Cell0DSetState(v, cell0D.Active);
      mesh.Cell0DInsertCoordinates(v, Vector3d(cell0D.X, cell0D.Y, cell0D.Z));
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell1Ds(const vector<MeshImporterFromCsvUtilities::Cell1D> cell1Ds,
                                                    IMeshDAO& mesh) const
  {
    const unsigned int numCell1Ds = cell1Ds.size();

    if (numCell1Ds == 0)
      return;

    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      const Cell1D& cell1D = cell1Ds[e];

      mesh.Cell1DSetId(e, cell1D.Id);
      mesh.Cell1DSetMarker(e, cell1D.Marker);
      mesh.Cell1DSetState(e, cell1D.Active);
      mesh.Cell1DInsertExtremes(e,
                                cell1D.Origin,
                                cell1D.End);
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell2Ds(const vector<MeshImporterFromCsvUtilities::Cell2D> cell2Ds,
                                                    IMeshDAO& mesh) const
  {
    const unsigned int numCell2Ds = cell2Ds.size();

    if (numCell2Ds == 0)
      return;

    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const Cell2D& cell2D = cell2Ds[f];

      mesh.Cell2DSetId(f, cell2D.Id);
      mesh.Cell2DSetMarker(f, cell2D.Marker);
      mesh.Cell2DSetState(f, cell2D.Active);

      const unsigned int numCellVertices = cell2D.Vertices.size();
      const unsigned int numCellEdges = cell2D.Edges.size();
      Output::Assert(numCellVertices == numCellEdges);

      mesh.Cell2DInitializeVertices(f, numCellVertices);
      for (unsigned int v = 0; v < numCellVertices; v++)
        mesh.Cell2DInsertVertex(f, v, cell2D.Vertices[v]);

      mesh.Cell2DInitializeEdges(f, numCellEdges);
      for (unsigned int e = 0; e < numCellEdges; e++)
        mesh.Cell2DInsertEdge(f, e, cell2D.Edges[e]);
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell3Ds(const vector<MeshImporterFromCsvUtilities::Cell3D> cell3Ds,
                                                    IMeshDAO& mesh) const
  {
    const unsigned int numCell3Ds = cell3Ds.size();

    if (numCell3Ds == 0)
      return;

    mesh.Cell3DsInitialize(numCell3Ds);
    for (unsigned int c = 0; c < numCell3Ds; c++)
    {
      const Cell3D& cell3D = cell3Ds[c];

      mesh.Cell3DSetId(c, cell3D.Id);
      mesh.Cell3DSetMarker(c, cell3D.Marker);
      mesh.Cell3DSetState(c, cell3D.Active);

      const unsigned int numCellVertices = cell3D.Vertices.size();
      const unsigned int numCellEdges = cell3D.Edges.size();
      const unsigned int numCellFaces = cell3D.Faces.size();

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
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell0DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours,
                                                             IMeshDAO& mesh) const
  {
    const unsigned int numCell0Ds = cell0DNeighbours.size();

    if (numCell0Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const MeshImporterFromCsvUtilities::Cell0DNeighbours& cell0D = cell0DNeighbours[v];

      const unsigned int numCell1DNeighbours = cell0D.Cell1DNeighbours.size();

      mesh.Cell0DInitializeNeighbourCell1Ds(v, numCell1DNeighbours);
      for (unsigned int n = 0; n < numCell1DNeighbours; n++)
      {
        if (cell0D.Cell1DNeighbours[n] >= mesh.Cell1DTotalNumber())
          continue;

        mesh.Cell0DInsertNeighbourCell1D(v, n, cell0D.Cell1DNeighbours[n]);
      }

      const unsigned int numCell2DNeighbours = cell0D.Cell2DNeighbours.size();

      mesh.Cell0DInitializeNeighbourCell2Ds(v, numCell2DNeighbours);
      for (unsigned int n = 0; n < numCell2DNeighbours; n++)
      {
        if (cell0D.Cell2DNeighbours[n] >= mesh.Cell2DTotalNumber())
          continue;

        mesh.Cell0DInsertNeighbourCell2D(v, n, cell0D.Cell2DNeighbours[n]);
      }

      const unsigned int numCell3DNeighbours = cell0D.Cell3DNeighbours.size();

      mesh.Cell0DInitializeNeighbourCell3Ds(v, numCell3DNeighbours);
      for (unsigned int n = 0; n < numCell3DNeighbours; n++)
      {
        if (cell0D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
          continue;

        mesh.Cell0DInsertNeighbourCell3D(v, n, cell0D.Cell3DNeighbours[n]);
      }
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell1DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours,
                                                             IMeshDAO& mesh) const
  {
    const unsigned int numCell1Ds = cell1DNeighbours.size();

    if (numCell1Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell1Ds; v++)
    {
      const MeshImporterFromCsvUtilities::Cell1DNeighbours& cell1D = cell1DNeighbours[v];

      const unsigned int numCell2DNeighbours = cell1D.Cell2DNeighbours.size();

      mesh.Cell1DInitializeNeighbourCell2Ds(v, numCell2DNeighbours);
      for (unsigned int n = 0; n < numCell2DNeighbours; n++)
      {
        if (cell1D.Cell2DNeighbours[n] >= mesh.Cell2DTotalNumber())
          continue;

        mesh.Cell1DInsertNeighbourCell2D(v, n, cell1D.Cell2DNeighbours[n]);
      }

      const unsigned int numCell3DNeighbours = cell1D.Cell3DNeighbours.size();

      mesh.Cell1DInitializeNeighbourCell3Ds(v, numCell3DNeighbours);
      for (unsigned int n = 0; n < numCell3DNeighbours; n++)
      {
        if (cell1D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
          continue;

        mesh.Cell1DInsertNeighbourCell3D(v, n, cell1D.Cell3DNeighbours[n]);
      }
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell2DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours,
                                                             IMeshDAO& mesh) const
  {
    const unsigned int numCell2Ds = cell2DNeighbours.size();

    if (numCell2Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell2Ds; v++)
    {
      const MeshImporterFromCsvUtilities::Cell2DNeighbours& cell2D = cell2DNeighbours[v];

      const unsigned int numCell3DNeighbours = cell2D.Cell3DNeighbours.size();

      mesh.Cell2DInitializeNeighbourCell3Ds(v, numCell3DNeighbours);
      for (unsigned int n = 0; n < numCell3DNeighbours; n++)
      {
        if (cell2D.Cell3DNeighbours[n] >= mesh.Cell3DTotalNumber())
          continue;

        mesh.Cell2DInsertNeighbourCell3D(v, n, cell2D.Cell3DNeighbours[n]);
      }
    }
  }
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell0DDoubleProperties(const vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties,
                                                                   IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell0DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell0DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshImporterFromCsvUtilities::CellDoubleProperty& cellsProperty = cell0DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell0DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshImporterFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

        const unsigned int numValues = cellProperty.Values.size();

        if (numValues == 0)
          continue;

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
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell1DDoubleProperties(const vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties,
                                                                   IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell1DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell1DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshImporterFromCsvUtilities::CellDoubleProperty& cellsProperty = cell1DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell1DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshImporterFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

        const unsigned int numValues = cellProperty.Values.size();

        if (numValues == 0)
          continue;

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
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell2DDoubleProperties(const vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties,
                                                                   IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell2DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell2DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshImporterFromCsvUtilities::CellDoubleProperty& cellsProperty = cell2DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell2DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshImporterFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

        const unsigned int numValues = cellProperty.Values.size();

        if (numValues == 0)
          continue;

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
  // ***************************************************************************
  void MeshImporterFromCsvUtilities::ConvertCell3DDoubleProperties(const vector<MeshImporterFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties,
                                                                   IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell3DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell3DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshImporterFromCsvUtilities::CellDoubleProperty& cellsProperty = cell3DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell3DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshImporterFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

        const unsigned int numValues = cellProperty.Values.size();

        if (numValues == 0)
          continue;

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
  // ***************************************************************************
  vector<MeshImporterFromCsvUtilities::Cell0D> MeshImporterFromCsvUtilities::ImportCell0Ds(IFileReader& csvFileReader,
                                                                                           const char& separator) const
  {
    vector<Cell0D> cell0Ds;

    /// Import Cell0Ds
    {
      vector<string> cell0DsLines;

      if (!csvFileReader.Open())
        return cell0Ds;

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
        return cell1Ds;

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
        return cell2Ds;

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
        return cell3Ds;

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
        return cellProperties;

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
        return cell0DNeighbours;

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
        return cell1DNeighbours;

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
        return cell2DNeighbours;

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
