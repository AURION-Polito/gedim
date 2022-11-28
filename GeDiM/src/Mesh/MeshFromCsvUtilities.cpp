#include "MeshFromCsvUtilities.hpp"
#include <iostream>
#include <fstream>
#include "FileTextReader.hpp"
#include "MeshUtilities.hpp"

using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  MeshFromCsvUtilities::MeshFromCsvUtilities()
  {
  }
  MeshFromCsvUtilities::~MeshFromCsvUtilities()
  {
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertMesh2D(const vector<MeshFromCsvUtilities::Cell0D>& cell0Ds,
                                           const vector<MeshFromCsvUtilities::Cell1D>& cell1Ds,
                                           const vector<MeshFromCsvUtilities::Cell2D>& cell2Ds,
                                           IMeshDAO& mesh) const
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

    MeshUtilities meshUtilities;

    meshUtilities.FillMesh2D(meshCell0Ds,
                             meshCell1Ds,
                             meshCell2Ds,
                             mesh);

    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const Cell0D& cell0D = cell0Ds[v];
      mesh.Cell0DSetMarker(v, cell0D.Marker);
      mesh.Cell0DSetState(v, cell0D.Active);
    }

    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      const Cell1D& cell1D = cell1Ds[e];
      mesh.Cell1DSetMarker(e, cell1D.Marker);
      mesh.Cell1DSetState(e, cell1D.Active);
    }

    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const Cell2D& cell2D = cell2Ds[f];
      mesh.Cell2DSetMarker(f, cell2D.Marker);
      mesh.Cell2DSetState(f, cell2D.Active);
    }
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertCell0Ds(const vector<MeshFromCsvUtilities::Cell0D> cell0Ds,
                                            IMeshDAO& mesh) const
  {
    const unsigned int numCell0Ds = cell0Ds.size();

    if (numCell0Ds == 0)
      return;

    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const Cell0D& cell0D = cell0Ds[v];

      mesh.Cell0DSetMarker(v, cell0D.Marker);
      mesh.Cell0DSetState(v, cell0D.Active);
      mesh.Cell0DInsertCoordinates(v, Vector3d(cell0D.X, cell0D.Y, cell0D.Z));
    }
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertCell1Ds(const vector<MeshFromCsvUtilities::Cell1D> cell1Ds,
                                            IMeshDAO& mesh) const
  {
    const unsigned int numCell1Ds = cell1Ds.size();

    if (numCell1Ds == 0)
      return;

    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      const Cell1D& cell1D = cell1Ds[e];

      mesh.Cell1DSetMarker(e, cell1D.Marker);
      mesh.Cell1DSetState(e, cell1D.Active);
      mesh.Cell1DInsertExtremes(e,
                                cell1D.Origin,
                                cell1D.End);
    }
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertCell2Ds(const vector<MeshFromCsvUtilities::Cell2D> cell2Ds,
                                            IMeshDAO& mesh) const
  {
    const unsigned int numCell2Ds = cell2Ds.size();

    if (numCell2Ds == 0)
      return;

    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const Cell2D& cell2D = cell2Ds[f];

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
  void MeshFromCsvUtilities::ConvertCell3Ds(const vector<MeshFromCsvUtilities::Cell3D> cell3Ds,
                                            IMeshDAO& mesh) const
  {
    const unsigned int numCell3Ds = cell3Ds.size();

    if (numCell3Ds == 0)
      return;

    mesh.Cell3DsInitialize(numCell3Ds);
    for (unsigned int c = 0; c < numCell3Ds; c++)
    {
      const Cell3D& cell3D = cell3Ds[c];

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
        mesh.Cell3DInsertFace(c, f, cell3D.Faces[f]);
    }
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertCell0DNeighbours(const vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours,
                                                     IMeshDAO& mesh) const
  {
    const unsigned int numCell0Ds = cell0DNeighbours.size();

    if (numCell0Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      const MeshFromCsvUtilities::Cell0DNeighbours& cell0D = cell0DNeighbours[v];

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
  void MeshFromCsvUtilities::ConvertCell1DNeighbours(const vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours,
                                                     IMeshDAO& mesh) const
  {
    const unsigned int numCell1Ds = cell1DNeighbours.size();

    if (numCell1Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell1Ds; v++)
    {
      const MeshFromCsvUtilities::Cell1DNeighbours& cell1D = cell1DNeighbours[v];

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
  void MeshFromCsvUtilities::ConvertCell2DNeighbours(const vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours,
                                                     IMeshDAO& mesh) const
  {
    const unsigned int numCell2Ds = cell2DNeighbours.size();

    if (numCell2Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell2Ds; v++)
    {
      const MeshFromCsvUtilities::Cell2DNeighbours& cell2D = cell2DNeighbours[v];

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
  void MeshFromCsvUtilities::ConvertCell2DSubDivisions(const vector<Cell2DSubDivision> cell2DSubDivisions,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numCell2Ds = cell2DSubDivisions.size();

    if (numCell2Ds == 0)
      return;

    for (unsigned int v = 0; v < numCell2Ds; v++)
    {
      const MeshFromCsvUtilities::Cell2DSubDivision& cell2D = cell2DSubDivisions[v];

      const unsigned int numCell2DSubDivision = cell2D.SubDivision.size();

      mesh.Cell2DInitializeSubDivision(v, numCell2DSubDivision);
      for (unsigned int n = 0; n < numCell2DSubDivision; n++)
        mesh.Cell2DInsertSubDivision(v, n, cell2D.SubDivision[n]);
    }
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ConvertCell0DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties,
                                                           IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell0DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell0DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshFromCsvUtilities::CellDoubleProperty& cellsProperty = cell0DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell0DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

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
  void MeshFromCsvUtilities::ConvertCell1DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties,
                                                           IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell1DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell1DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshFromCsvUtilities::CellDoubleProperty& cellsProperty = cell1DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell1DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

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
  void MeshFromCsvUtilities::ConvertCell2DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties,
                                                           IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell2DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell2DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshFromCsvUtilities::CellDoubleProperty& cellsProperty = cell2DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell2DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

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
  void MeshFromCsvUtilities::ConvertCell3DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties,
                                                           IMeshDAO& mesh) const
  {
    const unsigned int numCellProperties = cell3DDoubleProperties.size();

    if (numCellProperties == 0)
      return;

    mesh.Cell3DInitializeDoubleProperties(numCellProperties);

    for (unsigned int p = 0; p < numCellProperties; p++)
    {
      const MeshFromCsvUtilities::CellDoubleProperty& cellsProperty = cell3DDoubleProperties[p];
      const unsigned int numCells = cellsProperty.Values.size();

      if (numCells == 0)
        continue;

      unsigned int propertyIndex = mesh.Cell3DAddDoubleProperty(cellsProperty.Id);

      for (unsigned int c = 0; c < numCells; c++)
      {
        const MeshFromCsvUtilities::CellDoubleProperty::Value& cellProperty = cellsProperty.Values[c];

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
  vector<MeshFromCsvUtilities::Cell0D> MeshFromCsvUtilities::ImportCell0Ds(IFileReader& csvFileReader,
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
  vector<MeshFromCsvUtilities::Cell1D> MeshFromCsvUtilities::ImportCell1Ds(IFileReader& csvFileReader,
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
  vector<MeshFromCsvUtilities::Cell2D> MeshFromCsvUtilities::ImportCell2Ds(IFileReader& csvFileReader,
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
  vector<MeshFromCsvUtilities::Cell3D> MeshFromCsvUtilities::ImportCell3Ds(IFileReader& csvFileReader,
                                                                           const char& separator) const
  {
    vector<MeshFromCsvUtilities::Cell3D> cell3Ds;

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
  vector<MeshFromCsvUtilities::CellDoubleProperty> MeshFromCsvUtilities::ImportCellDoubleProperties(IFileReader& csvFileReader,
                                                                                                    const char& separator) const
  {
    vector<MeshFromCsvUtilities::CellDoubleProperty> cellProperties;

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
  vector<MeshFromCsvUtilities::CellDoubleProperty::Value> MeshFromCsvUtilities::ImportCellProperty(IFileReader& csvFileReader, const char& separator) const
  {
    vector<MeshFromCsvUtilities::CellDoubleProperty::Value> cellPropertyValues;

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
  vector<MeshFromCsvUtilities::Cell0DNeighbours> MeshFromCsvUtilities::ImportCell0DNeighbours(IFileReader& csvFileReader,
                                                                                              const char& separator) const
  {
    vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours;

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
          MeshFromCsvUtilities::Cell0DNeighbours& cell0D = cell0DNeighbours[v];

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
  vector<MeshFromCsvUtilities::Cell1DNeighbours> MeshFromCsvUtilities::ImportCell1DNeighbours(IFileReader& csvFileReader,
                                                                                              const char& separator) const
  {
    vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours;

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

        for (unsigned int e = 0; e < numCell1DNeighbours; e++)
        {
          MeshFromCsvUtilities::Cell1DNeighbours& cell1D = cell1DNeighbours[e];

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
  vector<MeshFromCsvUtilities::Cell2DNeighbours> MeshFromCsvUtilities::ImportCell2DNeighbours(IFileReader& csvFileReader,
                                                                                              const char& separator) const
  {
    vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours;

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
          MeshFromCsvUtilities::Cell2DNeighbours& cell2D = cell2DNeighbours[f];

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
  vector<MeshFromCsvUtilities::Cell2DSubDivision> MeshFromCsvUtilities::ImportCell2DSubDivision(IFileReader& csvFileReader,
                                                                                                const char& separator) const
  {
    vector<MeshFromCsvUtilities::Cell2DSubDivision> cell2DSubDivision;

    /// Import Cell2DSubDivision
    {
      vector<string> cell2DSubDivisionLines;

      if (!csvFileReader.Open())
        return cell2DSubDivision;

      csvFileReader.GetAllLines(cell2DSubDivisionLines);
      csvFileReader.Close();

      unsigned int numCell2DSubDivision = cell2DSubDivisionLines.size() - 1;

      if (numCell2DSubDivision > 0)
      {
        cell2DSubDivision.resize(numCell2DSubDivision);

        for (unsigned int f = 0; f < numCell2DSubDivision; f++)
        {
          MeshFromCsvUtilities::Cell2DSubDivision& cell2D = cell2DSubDivision[f];

          istringstream converter(cell2DSubDivisionLines[f + 1]);

          unsigned int cell2Did;
          char temp;
          converter >> cell2Did;
          if (separator != ' ')
            converter >> temp;

          unsigned int numCell3DNeighbours;
          converter >> numCell3DNeighbours;
          if (separator != ' ')
            converter >> temp;

          cell2D.SubDivision.resize(numCell3DNeighbours);
          for (unsigned int n = 0; n < numCell3DNeighbours; n++)
          {
            converter >> cell2D.SubDivision[n];
            if (separator != ' ')
              converter >> temp;
          }
        }
      }
    }

    return cell2DSubDivision;
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell0Ds(const string& filePath,
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
      fileCell0Ds<< scientific<< v<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DMarker(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DIsActive(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateX(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateY(v)<< separator;
      fileCell0Ds<< scientific<< mesh.Cell0DCoordinateZ(v)<< endl;
    }

    fileCell0Ds.close();
  }

  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell1Ds(const string& filePath,
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
      fileCell1Ds<< scientific<< e<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DMarker(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DIsActive(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DOrigin(e)<< separator;
      fileCell1Ds<< scientific<< mesh.Cell1DEnd(e)<< endl;
    }

    fileCell1Ds.close();
  }

  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell2Ds(const string& filePath,
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
      fileCell2Ds<< scientific<< f<< separator;
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
  void MeshFromCsvUtilities::ExportCell3Ds(const string& filePath,
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
      fileCell3Ds<< scientific<< c<< separator;
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
  void MeshFromCsvUtilities::ExportCell0DProperties(const string& Folder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell0D Properties
    ofstream fileCell0DProperties;

    fileCell0DProperties.open(Folder + "/" + propertyFileName + "." + propertyFileExtension);
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
                           Folder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell0DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell0DProperty(const unsigned int& propertyIndex,
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
      fileCell0DProperties<< scientific<< v<< separator;

      fileCell0DProperties<< scientific<< mesh.Cell0DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell0DDoublePropertySize(v, propertyIndex); n++)
        fileCell0DProperties<< scientific<< separator<< mesh.Cell0DDoublePropertyValue(v, propertyIndex, n);
      fileCell0DProperties<< endl;
    }

    fileCell0DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell1DProperties(const string& Folder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell1D Properties
    ofstream fileCell1DProperties;

    fileCell1DProperties.open(Folder + "/" + propertyFileName + "." + propertyFileExtension);
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
                           Folder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell1DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell1DProperty(const unsigned int& propertyIndex,
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
      fileCell1DProperties<< scientific<< v<< separator;

      fileCell1DProperties<< scientific<< mesh.Cell1DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell1DDoublePropertySize(v, propertyIndex); n++)
        fileCell1DProperties<< scientific<< separator<< mesh.Cell1DDoublePropertyValue(v, propertyIndex, n);
      fileCell1DProperties<< endl;
    }

    fileCell1DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell2DProperties(const string& Folder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell2D Properties
    ofstream fileCell2DProperties;

    fileCell2DProperties.open(Folder + "/" + propertyFileName + "." + propertyFileExtension);
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
                           Folder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell2DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell2DProperty(const unsigned int& propertyIndex,
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
      fileCell2DProperties<< scientific<< v<< separator;

      fileCell2DProperties<< scientific<< mesh.Cell2DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell2DDoublePropertySize(v, propertyIndex); n++)
        fileCell2DProperties<< scientific<< separator<< mesh.Cell2DDoublePropertyValue(v, propertyIndex, n);
      fileCell2DProperties<< endl;
    }

    fileCell2DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell3DProperties(const string& Folder,
                                                    const string& propertyFileName,
                                                    const string& propertyFileExtension,
                                                    const char& separator,
                                                    const IMeshDAO& mesh) const
  {
    /// Export Cell3D Properties
    ofstream fileCell3DProperties;

    fileCell3DProperties.open(Folder + "/" + propertyFileName + "." + propertyFileExtension);
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
                           Folder + "/" + propertyFilePath,
                           separator,
                           mesh);
    }

    fileCell3DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell3DProperty(const unsigned int& propertyIndex,
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
      fileCell3DProperties<< scientific<< v<< separator;

      fileCell3DProperties<< scientific<< mesh.Cell3DDoublePropertySize(v, propertyIndex);
      for (unsigned int n = 0; n < mesh.Cell3DDoublePropertySize(v, propertyIndex); n++)
        fileCell3DProperties<< scientific<< separator<< mesh.Cell3DDoublePropertyValue(v, propertyIndex, n);
      fileCell3DProperties<< endl;
    }

    fileCell3DProperties.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell0DNeighbours(const string& filePath,
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
      fileCell0DNeighbours<< scientific<< v<< separator;

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
  void MeshFromCsvUtilities::ExportCell1DNeighbours(const string& filePath,
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
      fileCell1DNeighbours<< scientific<< e<< separator;

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
  void MeshFromCsvUtilities::ExportCell2DNeighbours(const string& filePath,
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
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
    {
      fileCell2DNeighbours<< scientific<< f<< separator;

      fileCell2DNeighbours<< scientific<< mesh.Cell2DNumberNeighbourCell3D(f);
      for (unsigned int n = 0; n < mesh.Cell2DNumberNeighbourCell3D(f); n++)
        fileCell2DNeighbours<< scientific<< separator<< mesh.Cell2DNeighbourCell3D(f, n);
      fileCell2DNeighbours<< endl;
    }

    fileCell2DNeighbours.close();
  }
  // ***************************************************************************
  void MeshFromCsvUtilities::ExportCell2DSubDivisions(const string& filePath,
                                                      const char& separator,
                                                      const IMeshDAO& mesh) const
  {
    /// Export Cell2D Neigbours
    ofstream fileCell2DSubDivisions;

    fileCell2DSubDivisions.open(filePath);
    fileCell2DSubDivisions.precision(16);

    if (fileCell2DSubDivisions.fail())
      throw runtime_error("Error on mesh cell2DSubDivisions file");

    fileCell2DSubDivisions<< "Id"<< separator;
    fileCell2DSubDivisions<< "NumSubDivision"<< separator;
    fileCell2DSubDivisions<< "SubDivisions"<< endl;
    for (unsigned int f = 0; f < mesh.Cell2DTotalNumber(); f++)
    {
      fileCell2DSubDivisions<< scientific<< f<< separator;

      fileCell2DSubDivisions<< scientific<< mesh.Cell2DNumberSubDivision(f);
      for (unsigned int n = 0; n < mesh.Cell2DNumberSubDivision(f); n++)
        fileCell2DSubDivisions<< scientific<< separator<< mesh.Cell2DSubDivisionCell0D(f, n);
      fileCell2DSubDivisions<< endl;
    }

    fileCell2DSubDivisions.close();
  }
  // ***************************************************************************
}
