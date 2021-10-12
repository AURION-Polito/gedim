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
  class MeshFromCsvUtilities final
  {
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
          string FileCell2DSubDivisionsName = "Cell2DSubDivisions";
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

      struct Cell2DSubDivision {
          unsigned int Id;
          vector<unsigned int> SubDivision;
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
      MeshFromCsvUtilities();
      ~MeshFromCsvUtilities();

      /// \brief Convert a 2D Mesh
      /// \param cell0Ds the container of cell0Ds
      /// \param cell1Ds the container of cell1Ds
      /// \param cell2Ds the container of cell2Ds
      /// \param mesh the resulting mesh
      void ConvertMesh2D(const vector<MeshFromCsvUtilities::Cell0D>& cell0Ds,
                         const vector<MeshFromCsvUtilities::Cell1D>& cell1Ds,
                         const vector<MeshFromCsvUtilities::Cell2D>& cell2Ds,
                         IMeshDAO& mesh) const;

      /// \brief Convert the imported Cell0Ds to mesh
      /// \param cell0Ds the container of cell0Ds
      /// \param mesh the mesh
      void ConvertCell0Ds(const vector<MeshFromCsvUtilities::Cell0D> cell0Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell1Ds to mesh
      /// \param cell1Ds the container of cell1Ds
      /// \param mesh the mesh
      void ConvertCell1Ds(const vector<MeshFromCsvUtilities::Cell1D> cell1Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell2Ds to mesh
      /// \param cell2Ds the container of cell2Ds
      /// \param mesh the mesh
      void ConvertCell2Ds(const vector<MeshFromCsvUtilities::Cell2D> cell2Ds,
                          IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell3Ds to mesh
      /// \param cell3Ds the container of cell3Ds
      /// \param mesh the mesh
      void ConvertCell3Ds(const vector<MeshFromCsvUtilities::Cell3D> cell3Ds,
                          IMeshDAO& mesh) const;

      /// \brief Convert the imported Cell0D neighbours to mesh
      /// \param cell0DNeighbours the container of cell0D neighbours
      /// \param mesh the mesh
      void ConvertCell0DNeighbours(const vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours,
                                   IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell1D neighbours to mesh
      /// \param cell1DNeighbours the container of cell1D neighbours
      /// \param mesh the mesh
      void ConvertCell1DNeighbours(const vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours,
                                   IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell2D neighbours to mesh
      /// \param cell2DNeighbours the container of cell2D neighbours
      /// \param mesh the mesh
      void ConvertCell2DNeighbours(const vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours,
                                   IMeshDAO& mesh) const;

      /// \brief Convert the imported Cell2D subdivision to mesh
      /// \param cell2DSubDivisions the container of cell2D neighbours
      /// \param mesh the mesh
      void ConvertCell2DSubDivisions(const vector<MeshFromCsvUtilities::Cell2DSubDivision> cell2DSubDivisions,
                                     IMeshDAO& mesh) const;

      /// \brief Convert the imported Cell0D double properties to mesh
      /// \param cell0DDoubleProperties the container of cell0D double properties
      /// \param mesh the mesh
      void ConvertCell0DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties,
                                         IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell1D double properties to mesh
      /// \param cell1DDoubleProperties the container of cell1D double properties
      /// \param mesh the mesh
      void ConvertCell1DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties,
                                         IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell2D double properties to mesh
      /// \param cell2DDoubleProperties the container of cell2D double properties
      /// \param mesh the mesh
      void ConvertCell2DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties,
                                         IMeshDAO& mesh) const;
      /// \brief Convert the imported Cell3D double properties to mesh
      /// \param cell3DDoubleProperties the container of cell3D double properties
      /// \param mesh the mesh
      void ConvertCell3DDoubleProperties(const vector<MeshFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties,
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

      /// \brief Import Cell2DSubDivision; format: Id, NumSubDivision, SubDivisions
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<Cell2DSubDivision> ImportCell2DSubDivision(IFileReader& csvFileReader,
                                                        const char& separator) const;

      /// \brief Import CellProperties; format: Id, FilePath
      /// \param csvFileReader the file reader
      /// \param separator the file separator
      vector<CellDoubleProperty> ImportCellDoubleProperties(IFileReader& csvFileReader,
                                                            const char& separator) const;

      /// \brief Export Cell0Ds; format: Id, Marker, Active, X, Y, Z
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell0Ds(const string& filePath,
                         const char& separator,
                         const IMeshDAO& mesh) const;
      /// \brief Export Cell1Ds; format: Id, Marker, Active, Origin, End
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell1Ds(const string& filePath,
                         const char& separator,
                         const IMeshDAO& mesh) const;
      /// \brief Export Cell2Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell2Ds(const string& filePath,
                         const char& separator,
                         const IMeshDAO& mesh) const;
      /// \brief Export Cell3Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges, NumFaces, Faces
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell3Ds(const string& filePath,
                         const char& separator,
                         const IMeshDAO& mesh) const;

      /// \brief Export Cell0DProperties; format: Id, FilePath
      /// \param exportFolder the folder where to export the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell0DProperties(const string& exportFolder,
                                  const string& propertyFileName,
                                  const string& propertyFileExtension,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;

      /// \brief Export Cell0DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell0DProperty(const unsigned int& propertyIndex,
                                const string& filePath,
                                const char& separator,
                                const IMeshDAO& mesh) const;

      /// \brief Export Cell1DProperties; format: Id, FilePath
      /// \param exportFolder the folder where to export the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell1DProperties(const string& exportFolder,
                                  const string& propertyFileName,
                                  const string& propertyFileExtension,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;

      /// \brief Export Cell1DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell1DProperty(const unsigned int& propertyIndex,
                                const string& filePath,
                                const char& separator,
                                const IMeshDAO& mesh) const;

      /// \brief Export Cell2DProperties; format: Id, FilePath
      /// \param exportFolder the folder where to export the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell2DProperties(const string& exportFolder,
                                  const string& propertyFileName,
                                  const string& propertyFileExtension,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;

      /// \brief Export Cell2DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell2DProperty(const unsigned int& propertyIndex,
                                const string& filePath,
                                const char& separator,
                                const IMeshDAO& mesh) const;

      /// \brief Export Cell3DProperties; format: Id, FilePath
      /// \param exportFolder the folder where to export the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell3DProperties(const string& exportFolder,
                                  const string& propertyFileName,
                                  const string& propertyFileExtension,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;

      /// \brief Export Cell3DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell3DProperty(const unsigned int& propertyIndex,
                                const string& filePath,
                                const char& separator,
                                const IMeshDAO& mesh) const;

      /// \brief Export Cell0DNeighbours; format: Id, Num1DNeighbours, 1DNeighbours, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell0DNeighbours(const string& filePath,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;
      /// \brief Export Cell1DNeighbours; format: Id, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell1DNeighbours(const string& filePath,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;
      /// \brief Export Cell2DNeighbours; format: Id, Num3DNeighbours, 3DNeighbours
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell2DNeighbours(const string& filePath,
                                  const char& separator,
                                  const IMeshDAO& mesh) const;
      /// \brief Export Cell2DSubDivisions; format: Id, NumSubDivision, SubDivisions
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be exported
      void ExportCell2DSubDivisions(const string& filePath,
                                    const char& separator,
                                    const IMeshDAO& mesh) const;
  };

}

#endif // __MeshImporterFromCsvUtilities_H
