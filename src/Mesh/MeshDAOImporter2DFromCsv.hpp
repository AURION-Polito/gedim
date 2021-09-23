#ifndef __MeshDAOImporter2DFromCsv_H
#define __MeshDAOImporter2DFromCsv_H

#include "IMeshDAO.hpp"
#include "IFileTextReader.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  /// \brief MeshDAOImporter2DFromCsv
  /// \note each file could be EmptyFileReader if not necessary
  /// \copyright See top level LICENSE file for details
  class MeshDAOImporter2DFromCsv final {

    private:
      struct CellProperty {
          struct Value {
              unsigned int CellId;
              vector<double> Values;
          };

          string Id;
          string FilePath;
      };

      struct Cell0D {
          unsigned int Id;
          double X;
          double Y;
          double Z;
          unsigned int Marker;
          bool Active;
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
          vector<unsigned int> Cell2DNeighbours;
          vector<unsigned int> Cell3DNeighbours;
      };

      struct Cell2D {
          unsigned int Id;
          vector<unsigned int> Vertices;
          vector<unsigned int> Edges;
          unsigned int Marker;
          bool Active;
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

    public:
      class Configuration {
        public:
          IFileReader& CsvCell0DsFile;
          IFileReader& CsvCell1DsFile;
          IFileReader& CsvCell2DsFile;
          IFileReader& CsvCell0DPropertiesFile;
          IFileReader& CsvCell1DPropertiesFile;
          IFileReader& CsvCell2DPropertiesFile;
          IFileReader& CsvCell0DNeighboursFile;
          IFileReader& CsvCell1DNeighboursFile;
          char Separator = ';';

          Configuration(IFileReader& csvCell0DsFile,
                        IFileReader& csvCell1DsFile,
                        IFileReader& csvCell2DsFile,
                        IFileReader& csvCell0DPropertiesFile,
                        IFileReader& csvCell1DPropertiesFile,
                        IFileReader& csvCell2DPropertiesFile,
                        IFileReader& csvCell0DNeighboursFile,
                        IFileReader& csvCell1DNeighboursFile) :
            CsvCell0DsFile(csvCell0DsFile),
            CsvCell1DsFile(csvCell1DsFile),
            CsvCell2DsFile(csvCell2DsFile),
            CsvCell0DPropertiesFile(csvCell0DPropertiesFile),
            CsvCell1DPropertiesFile(csvCell1DPropertiesFile),
            CsvCell2DPropertiesFile(csvCell2DPropertiesFile),
            CsvCell0DNeighboursFile(csvCell0DNeighboursFile),
            CsvCell1DNeighboursFile(csvCell1DNeighboursFile)
          { }
      };

    public:
      MeshDAOImporter2DFromCsv();
      ~MeshDAOImporter2DFromCsv();

      /// \brief Import Cell0Ds; format: Id, Marker, Active, X, Y, Z
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell0Ds(IFileReader& csvFileReader,
                         const char& separator,
                         vector<Cell0D>& cell0Ds) const;
      /// \brief Import Cell1Ds; format: Id, Marker, Active, Origin, End
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell1Ds(IFileReader& csvFileReader,
                         const char& separator,
                         vector<Cell1D>& cell1Ds) const;
      /// \brief Import Cell2Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell2Ds(IFileReader& csvFileReader,
                         const char& separator,
                         vector<Cell2D>& cell2Ds) const;

      void CreateMesh2D(const vector<Cell0D>& cell0Ds,
                        const vector<Cell1D>& cell1Ds,
                        const vector<Cell2D>& cell2Ds,
                        IMeshDAO& mesh);

      /// \brief Import Cell0DProperties; format: Id, FilePath
      /// \param ImportFolder the folder where to Import the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell0DProperties(IFileReader& csvFileReader,
                                  const char& separator,
                                  IMeshDAO& mesh) const;

      /// \brief Import Cell0DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell0DProperty(const unsigned int& propertyIndex,
                                IFileReader& csvFileReader,
                                const char& separator,
                                IMeshDAO& mesh) const;

      /// \brief Import Cell1DProperties; format: Id, FilePath
      /// \param ImportFolder the folder where to Import the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell1DProperties(IFileReader& csvFileReader,
                                  const char& separator,
                                  IMeshDAO& mesh) const;

      /// \brief Import Cell1DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell1DProperty(const unsigned int& propertyIndex,
                                IFileReader& csvFileReader,
                                const char& separator,
                                IMeshDAO& mesh) const;

      /// \brief Import Cell2DProperties; format: Id, FilePath
      /// \param ImportFolder the folder where to Import the files
      /// \param propertyFileName the name of property file
      /// \param propertyFileExtension the extension of the files
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell2DProperties(IFileReader& csvFileReader,
                                  const char& separator,
                                  IMeshDAO& mesh) const;

      /// \brief Import Cell2DProperty identified by index; format: Id, PropertySize, PropertyValues
      /// \param propertyIndex the property index
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell2DProperty(const unsigned int& propertyIndex,
                                IFileReader& csvFileReader,
                                const char& separator,
                                IMeshDAO& mesh) const;

      /// \brief Import Cell0DNeighbours; format: Id, Num1DNeighbours, 1DNeighbours, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell0DNeighbours(IFileReader& csvFileReader,
                                  const char& separator,
                                  IMeshDAO& mesh) const;
      /// \brief Import Cell1DNeighbours; format: Id, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
      /// \param filePath the path of the file
      /// \param separator the file separator
      /// \param mesh the mesh to be Imported
      void ImportCell1DNeighbours(IFileReader& csvFileReader,
                                  const char& separator,
                                  IMeshDAO& mesh) const;

      void Import(Configuration& configuration,
                  IMeshDAO& mesh);
  };

}

#endif // __MeshDAOImporter2DFromCsv_H
