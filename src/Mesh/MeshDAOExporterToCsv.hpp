#ifndef __MeshDAOExporterToCsv_H
#define __MeshDAOExporterToCsv_H

#include "IMeshDAO.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  /// \brief MeshDAOExporterToCsv
  /// \copyright See top level LICENSE file for details.
  class MeshDAOExporterToCsv final {
    public:
      struct Configuration {
          string ExportFolder = "./";
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

    public:
      MeshDAOExporterToCsv();
      ~MeshDAOExporterToCsv();

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

      /// \brief Export the mesh in all parts
      /// \param configuration the configuration for export
      /// \param mesh the mesh to be exported
      void Export(const MeshDAOExporterToCsv::Configuration& configuration,
                  const IMeshDAO& mesh) const;
  };

}

#endif // __MeshDAOExporterToCsv_H
