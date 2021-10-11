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
