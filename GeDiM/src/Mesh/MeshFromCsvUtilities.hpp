#ifndef __MeshImporterFromCsvUtilities_H
#define __MeshImporterFromCsvUtilities_H

#include "IFileTextReader.hpp"
#include "IMeshDAO.hpp"

namespace Gedim
{
/// \brief MeshImporterFromCsvUtilities
/// \note each file could be EmptyFileReader if not necessary
/// \copyright See top level LICENSE file for details
class MeshFromCsvUtilities final
{
  public:
    struct Configuration
    {
        std::string Folder = "./";
        std::string FileCell0DsName = "Cell0Ds";
        std::string FileCell1DsName = "Cell1Ds";
        std::string FileCell2DsName = "Cell2Ds";
        std::string FileCell3DsName = "Cell3Ds";
        std::string FileCell0DNeighboursName = "Cell0DNeighbours";
        std::string FileCell1DNeighboursName = "Cell1DNeighbours";
        std::string FileCell2DNeighboursName = "Cell2DNeighbours";
        std::string FileCell0DPropertiesName = "Cell0DProperties";
        std::string FileCell1DPropertiesName = "Cell1DProperties";
        std::string FileCell2DPropertiesName = "Cell2DProperties";
        std::string FileCell3DPropertiesName = "Cell3DProperties";
        std::string FileCell2DSubDivisionsName = "Cell2DSubDivisions";
        std::string FileCell0DUpdatedCellsName = "Cell0DUpdatedCells";
        std::string FileCell1DUpdatedCellsName = "Cell1DUpdatedCells";
        std::string FileCell2DUpdatedCellsName = "Cell2DUpdatedCells";
        std::string FileCell3DUpdatedCellsName = "Cell3DUpdatedCells";
        char Separator = ';';
        std::string FileExtension = "csv";
    };

    struct CellDoubleProperty
    {
        struct Value
        {
            unsigned int CellId;
            std::vector<double> Values;
        };

        std::string Id;
        std::string FilePath;
        std::vector<Value> Values;
    };

    struct Cell0D
    {
        unsigned int Id;
        double X;
        double Y;
        double Z;
        unsigned int Marker;
        bool Active;
    };

    struct Cell0DNeighbours
    {
        unsigned int Id;
        std::vector<unsigned int> Cell1DNeighbours;
        std::vector<unsigned int> Cell2DNeighbours;
        std::vector<unsigned int> Cell3DNeighbours;
    };

    struct CellUpdatedCells
    {
        unsigned int Id;
        std::vector<unsigned int> UpdatedCells;
    };

    struct Cell1D
    {
        unsigned int Id;
        unsigned int Origin;
        unsigned int End;
        unsigned int Marker;
        bool Active;
    };

    struct Cell1DNeighbours
    {
        unsigned int Id;
        std::vector<unsigned int> Cell2DNeighbours;
        std::vector<unsigned int> Cell3DNeighbours;
    };

    struct Cell2D
    {
        unsigned int Id;
        std::vector<unsigned int> Vertices;
        std::vector<unsigned int> Edges;
        unsigned int Marker;
        bool Active;
    };

    struct Cell2DNeighbours
    {
        unsigned int Id;
        std::vector<unsigned int> Cell3DNeighbours;
    };

    struct Cell2DSubDivision
    {
        unsigned int Id;
        std::vector<unsigned int> SubDivision;
    };

    struct Cell3D
    {
        unsigned int Id;
        std::vector<unsigned int> Vertices;
        std::vector<unsigned int> Edges;
        std::vector<unsigned int> Faces;
        unsigned int Marker;
        bool Active;
    };

  private:
    /// \brief Import Cell0DProperty identified by index; format: Id, PropertySize, PropertyValues
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<CellDoubleProperty::Value> ImportCellProperty(IFileReader &csvFileReader, const char &separator) const;

  public:
    MeshFromCsvUtilities();
    ~MeshFromCsvUtilities();

    /// \brief Convert a 2D Mesh
    /// \param cell0Ds the container of cell0Ds
    /// \param cell1Ds the container of cell1Ds
    /// \param cell2Ds the container of cell2Ds
    /// \param mesh the resulting mesh
    void ConvertMesh2D(const std::vector<MeshFromCsvUtilities::Cell0D> &cell0Ds,
                       const std::vector<MeshFromCsvUtilities::Cell1D> &cell1Ds,
                       const std::vector<MeshFromCsvUtilities::Cell2D> &cell2Ds,
                       IMeshDAO &mesh) const;

    /// \brief Convert the imported Cell0Ds to mesh
    /// \param cell0Ds the container of cell0Ds
    /// \param mesh the mesh
    void ConvertCell0Ds(const std::vector<MeshFromCsvUtilities::Cell0D> cell0Ds, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell1Ds to mesh
    /// \param cell1Ds the container of cell1Ds
    /// \param mesh the mesh
    void ConvertCell1Ds(const std::vector<MeshFromCsvUtilities::Cell1D> cell1Ds, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell2Ds to mesh
    /// \param cell2Ds the container of cell2Ds
    /// \param mesh the mesh
    void ConvertCell2Ds(const std::vector<MeshFromCsvUtilities::Cell2D> cell2Ds, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell3Ds to mesh
    /// \param cell3Ds the container of cell3Ds
    /// \param mesh the mesh
    void ConvertCell3Ds(const std::vector<MeshFromCsvUtilities::Cell3D> cell3Ds, IMeshDAO &mesh) const;

    /// \brief Convert the imported Cell0D neighbours to mesh
    /// \param cell0DNeighbours the container of cell0D neighbours
    /// \param mesh the mesh
    void ConvertCell0DNeighbours(const std::vector<MeshFromCsvUtilities::Cell0DNeighbours> cell0DNeighbours, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell1D neighbours to mesh
    /// \param cell1DNeighbours the container of cell1D neighbours
    /// \param mesh the mesh
    void ConvertCell1DNeighbours(const std::vector<MeshFromCsvUtilities::Cell1DNeighbours> cell1DNeighbours, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell2D neighbours to mesh
    /// \param cell2DNeighbours the container of cell2D neighbours
    /// \param mesh the mesh
    void ConvertCell2DNeighbours(const std::vector<MeshFromCsvUtilities::Cell2DNeighbours> cell2DNeighbours, IMeshDAO &mesh) const;

    /// \brief Convert the imported Cell2D subdivision to mesh
    /// \param cell2DSubDivisions the container of cell2D neighbours
    /// \param mesh the mesh
    void ConvertCell2DSubDivisions(const std::vector<MeshFromCsvUtilities::Cell2DSubDivision> cell2DSubDivisions, IMeshDAO &mesh) const;

    /// \brief Convert the imported Cell0D double properties to mesh
    /// \param cell0DDoubleProperties the container of cell0D double properties
    /// \param mesh the mesh
    void ConvertCell0DDoubleProperties(const std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell0DDoubleProperties,
                                       IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell1D double properties to mesh
    /// \param cell1DDoubleProperties the container of cell1D double properties
    /// \param mesh the mesh
    void ConvertCell1DDoubleProperties(const std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell1DDoubleProperties,
                                       IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell2D double properties to mesh
    /// \param cell2DDoubleProperties the container of cell2D double properties
    /// \param mesh the mesh
    void ConvertCell2DDoubleProperties(const std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell2DDoubleProperties,
                                       IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell3D double properties to mesh
    /// \param cell3DDoubleProperties the container of cell3D double properties
    /// \param mesh the mesh
    void ConvertCell3DDoubleProperties(const std::vector<MeshFromCsvUtilities::CellDoubleProperty> cell3DDoubleProperties,
                                       IMeshDAO &mesh) const;

    /// \brief Convert the imported Cell0D updated cells to mesh
    /// \param cell0DUpdatedCells the container of cell0D updated cells
    /// \param mesh the mesh
    void ConvertCell0DUpdatedCells(const std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell0DUpdatedCells, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell1D updated cells to mesh
    /// \param cell1DUpdatedCells the container of cell1D updated cells
    /// \param mesh the mesh
    void ConvertCell1DUpdatedCells(const std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell1DUpdatedCells, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell2D updated cells to mesh
    /// \param cell2DUpdatedCells the container of cell2D updated cells
    /// \param mesh the mesh
    void ConvertCell2DUpdatedCells(const std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell2DUpdatedCells, IMeshDAO &mesh) const;
    /// \brief Convert the imported Cell3D updated cells to mesh
    /// \param cell3DUpdatedCells the container of cell3D updated cells
    /// \param mesh the mesh
    void ConvertCell3DUpdatedCells(const std::vector<MeshFromCsvUtilities::CellUpdatedCells> cell3DUpdatedCells, IMeshDAO &mesh) const;

    /// \brief Import Cell0Ds; format: Id, Marker, Active, X, Y, Z
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    /// \param mesh the mesh to be Imported
    std::vector<Cell0D> ImportCell0Ds(IFileReader &csvFileReader, const char &separator) const;
    /// \brief Import Cell1Ds; format: Id, Marker, Active, Origin, End
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell1D> ImportCell1Ds(IFileReader &csvFileReader, const char &separator) const;
    /// \brief Import Cell2Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell2D> ImportCell2Ds(IFileReader &csvFileReader, const char &separator) const;
    /// \brief Import Cell3Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges, NumFaces, Faces
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell3D> ImportCell3Ds(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Import Cell0DNeighbours; format: Id, Num1DNeighbours, 1DNeighbours, Num2DNeighbours, 2DNeighbours,
    /// Num3DNeighbours, 3DNeighbours \param csvFileReader the file reader \param separator the file separator
    std::vector<Cell0DNeighbours> ImportCell0DNeighbours(IFileReader &csvFileReader, const char &separator) const;
    /// \brief Import Cell1DNeighbours; format: Id, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell1DNeighbours> ImportCell1DNeighbours(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Import Cell2DNeighbours; format: Id, Num3DNeighbours, 3DNeighbours
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell2DNeighbours> ImportCell2DNeighbours(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Import Cell2DSubDivision; format: Id, NumSubDivision, SubDivisions
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<Cell2DSubDivision> ImportCell2DSubDivision(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Import CellProperties; format: Id, FilePath
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<CellDoubleProperty> ImportCellDoubleProperties(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Import CellUpdatedCells; format: Id, NumUpdatedCells, UpdatedCells
    /// \param csvFileReader the file reader
    /// \param separator the file separator
    std::vector<CellUpdatedCells> ImportCellUpdatedCells(IFileReader &csvFileReader, const char &separator) const;

    /// \brief Export Cell0Ds; format: Id, Marker, Active, X, Y, Z
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell0Ds(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell1Ds; format: Id, Marker, Active, Origin, End
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell1Ds(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell2Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2Ds(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell3Ds; format: Id, Marker, Active, NumVertices, Vertices, NumEdges, Edges, NumFaces, Faces
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell3Ds(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell0DProperties; format: Id, FilePath
    /// \param exportFolder the folder where to export the files
    /// \param propertyFileName the name of property file
    /// \param propertyFileExtension the extension of the files
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell0DProperties(const std::string &exportFolder,
                                const std::string &propertyFileName,
                                const std::string &propertyFileExtension,
                                const char &separator,
                                const IMeshDAO &mesh) const;

    /// \brief Export Cell0DProperty identified by index; format: Id, PropertySize, PropertyValues
    /// \param propertyIndex the property index
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell0DProperty(const unsigned int &propertyIndex, const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell1DProperties; format: Id, FilePath
    /// \param exportFolder the folder where to export the files
    /// \param propertyFileName the name of property file
    /// \param propertyFileExtension the extension of the files
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell1DProperties(const std::string &exportFolder,
                                const std::string &propertyFileName,
                                const std::string &propertyFileExtension,
                                const char &separator,
                                const IMeshDAO &mesh) const;

    /// \brief Export Cell1DProperty identified by index; format: Id, PropertySize, PropertyValues
    /// \param propertyIndex the property index
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell1DProperty(const unsigned int &propertyIndex, const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell2DProperties; format: Id, FilePath
    /// \param exportFolder the folder where to export the files
    /// \param propertyFileName the name of property file
    /// \param propertyFileExtension the extension of the files
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2DProperties(const std::string &exportFolder,
                                const std::string &propertyFileName,
                                const std::string &propertyFileExtension,
                                const char &separator,
                                const IMeshDAO &mesh) const;

    /// \brief Export Cell2DProperty identified by index; format: Id, PropertySize, PropertyValues
    /// \param propertyIndex the property index
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2DProperty(const unsigned int &propertyIndex, const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell3DProperties; format: Id, FilePath
    /// \param exportFolder the folder where to export the files
    /// \param propertyFileName the name of property file
    /// \param propertyFileExtension the extension of the files
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell3DProperties(const std::string &exportFolder,
                                const std::string &propertyFileName,
                                const std::string &propertyFileExtension,
                                const char &separator,
                                const IMeshDAO &mesh) const;

    /// \brief Export Cell3DProperty identified by index; format: Id, PropertySize, PropertyValues
    /// \param propertyIndex the property index
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell3DProperty(const unsigned int &propertyIndex, const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell0DNeighbours; format: Id, Num1DNeighbours, 1DNeighbours, Num2DNeighbours, 2DNeighbours,
    /// Num3DNeighbours, 3DNeighbours \param filePath the path of the file \param separator the file separator \param
    /// mesh the mesh to be exported
    void ExportCell0DNeighbours(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell1DNeighbours; format: Id, Num2DNeighbours, 2DNeighbours, Num3DNeighbours, 3DNeighbours
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell1DNeighbours(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell2DNeighbours; format: Id, Num3DNeighbours, 3DNeighbours
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2DNeighbours(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell2DSubDivisions; format: Id, NumSubDivision, SubDivisions
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2DSubDivisions(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;

    /// \brief Export Cell0DUpdatedCells; format: Id, NumUpdatedCells, UpdatedCells
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell0DUpdatedCells(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell1DUpdatedCells; format: Id, NumUpdatedCells, UpdatedCells
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell1DUpdatedCells(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell2DUpdatedCells; format: Id, NumUpdatedCells, UpdatedCells
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell2DUpdatedCells(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
    /// \brief Export Cell3DUpdatedCells; format: Id, NumUpdatedCells, UpdatedCells
    /// \param filePath the path of the file
    /// \param separator the file separator
    /// \param mesh the mesh to be exported
    void ExportCell3DUpdatedCells(const std::string &filePath, const char &separator, const IMeshDAO &mesh) const;
};

} // namespace Gedim

#endif // __MeshImporterFromCsvUtilities_H
