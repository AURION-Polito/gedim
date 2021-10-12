#ifndef __MeshImporterFromCsvUtilities_H
#define __MeshImporterFromCsvUtilities_H

#include "IMeshDAO.hpp"
#include "IFileTextReader.hpp"

using namespace std;

namespace Gedim
{
  /// \brief MeshImporterFromCsvUtilities
  /// \note each file could be EmptyFileReader if not necessary
  /// \copyright See top level LICENSE file for details
  class MeshImporterFromCsvUtilities final {
    public:
      struct Configuration {
          string Folder = "./";
          string FileCell0DsName = "Cell0Ds";
          string FileCell1DsName = "Cell1Ds";
          string FileCell2DsName = "Cell2Ds";
          string FileCell3DsName = "Cell3Ds";
          string FileCell0DNeighboursName = "Cell0DNeighbours";
          string FileCell1DNeighboursName = "Cell1DNeighbours";
          string FileCell2DNeighboursName = "Cell2DNeighbours";
          string FileCell0DPropertiesName = "Cell0DProperties";
          string FileCell1DPropertiesName = "Cell1DProperties";
          string FileCell2DPropertiesName = "Cell2DProperties";
          string FileCell3DPropertiesName = "Cell3DProperties";
          char Separator = ';';
          string FileExtension = "csv";
      };

      struct CellDoubleProperty {
          struct Value {
              unsigned int CellId;
              vector<double> Values;
          };

          string Id;
          string FilePath;
          vector<Value> Values;
      };

      struct Cell0D {
          unsigned int Id;
          double X;
          double Y;
          double Z;
          unsigned int Marker;
          bool Active;
      };

      struct Cell0DNeighbours {
          unsigned int Id;
          vector<unsigned int> Cell1DNeighbours;
          vector<unsigned int> Cell2DNeighbours;
          vector<unsigned int> Cell3DNeighbours;
      };

      struct Cell1D {
          unsigned int Id;
          unsigned int Origin;
          unsigned int End;
          unsigned int Marker;
          bool Active;
      };

      struct Cell1DNeighbours {
          unsigned int Id;
          vector<unsigned int> Cell2DNeighbours;
          vector<unsigned int> Cell3DNeighbours;
      };

      struct Cell2D {
          unsigned int Id;
          vector<unsigned int> Vertices;
          vector<unsigned int> Edges;
          unsigned int Marker;
          bool Active;
      };

      struct Cell2DNeighbours {
          unsigned int Id;
          vector<unsigned int> Cell3DNeighbours;
      };

      struct Cell3D {
          unsigned int Id;
          vector<unsigned int> Vertices;
          vector<unsigned int> Edges;
          vector<unsigned int> Faces;
          unsigned int Marker;
          bool Active;
      };

    private:
      /// \brief Import Cell0DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<CellDoubleProperty::Value> ImportCellProperty(IFileReader& csvFileReader,
                                                           const char& separator) const;

    public:
      MeshImporterFromCsvUtilities();
      ~MeshImporterFromCsvUtilities();

      /// \brief Convert the imported Cell0Ds to mesh
      /// \param cell0Ds the container of cell0Ds
      /// \param mesh the mesh
      void ConvertCell0Ds(const vector<MeshImporterFromCsvUtilities::Cell0D> cell0Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell1Ds to mesh
      /// \param cell1Ds the container of cell1Ds
      /// \param mesh the mesh
      void ConvertCell1Ds(const vector<MeshImporterFromCsvUtilities::Cell1D> cell1Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell2Ds to mesh
      /// \param cell2Ds the container of cell2Ds
      /// \param mesh the mesh
      void ConvertCell2Ds(const vector<MeshImporterFromCsvUtilities::Cell2D> cell2Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell3Ds to mesh
      /// \param cell3Ds the container of cell3Ds
      /// \param mesh the mesh
      void ConvertCell3Ds(const vector<MeshImporterFromCsvUtilities::Cell3D> cell3Ds,
                          IMeshDAO& mesh) const;

      /// \brief Convert the imported Cell0D neighbours to mesh
      /// \param cell0DNeighbours the container of cell0D neighbours
      /// \param mesh the mesh
      void ConvertCell0DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours,
                                   IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell1D neighbours to mesh
      /// \param cell1DNeighbours the container of cell1D neighbours
      /// \param mesh the mesh
      void ConvertCell1DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours,
                                   IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell2D neighbours to mesh
      /// \param cell2DNeighbours the container of cell2D neighbours
      /// \param mesh the mesh
      void ConvertCell2DNeighbours(const vector<MeshImporterFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours,
                                   IMeshDAO& mesh) const;

      /// \brief Import Cell0Ds; format: Id, Marker, Active, X, Y, Z
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      vector<Cell0D> ImportCell0Ds(IFileReader& csvFileReader,
                                   const char& separator) const;
      /// \brief Import Cell1Ds; format: Id, Marker, Active, Origin, End
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell1D> ImportCell1Ds(IFileReader& csvFileReader,
                                   const char& separator) const;
      /// \brief Import Cell2Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell2D> ImportCell2Ds(IFileReader& csvFileReader,
                                   const char& separator) const;
      /// \brief Import Cell3Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges, NumFaces, Faces
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell3D> ImportCell3Ds(IFileReader& csvFileReader,
                                   const char& separator) const;

      /// \brief Import Cell0DNeighbours; format: Id, Num1DNeighbours, 1DNeighbours, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell0DNeighbours> ImportCell0DNeighbours(IFileReader& csvFileReader,
                                                      const char& separator) const;
      /// \brief Import Cell1DNeighbours; format: Id, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell1DNeighbours> ImportCell1DNeighbours(IFileReader& csvFileReader,
                                                      const char& separator) const;

      /// \brief Import Cell2DNeighbours; format: Id, Num3DNeighbours, 3DNeighbours
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell2DNeighbours> ImportCell2DNeighbours(IFileReader& csvFileReader,
                                                      const char& separator) const;

      /// \brief Import CellProperties; format: Id, FilePath
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<CellDoubleProperty> ImportCellDoubleProperties(IFileReader& csvFileReader,
                                                            const char& separator) const;
  };

}

#endif // __MeshImporterFromCsvUtilities_H
