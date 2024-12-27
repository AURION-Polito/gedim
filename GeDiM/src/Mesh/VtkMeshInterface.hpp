#ifndef __VtkMeshInterface_H
#define __VtkMeshInterface_H

#include "IMeshDAO.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
class VtkMeshInterface final
{
  public:
    struct VtkMesh final
    {
        struct Cell3D final
        {
            enum struct Types
            {
                Unknown = 0,
                Tetrahedron = 1,
                Hexahedron = 2
            };

            std::vector<unsigned int> Cell0D_Id;
            Types Type;
        };

        unsigned int NumCell0Ds;
        unsigned int NumCell3Ds;

        Eigen::MatrixXd Cell0Ds;
        std::vector<Cell3D> Cell3Ds;
    };

    struct VtkMesh3D final
    {
        std::array<std::vector<unsigned int>, 4> Markers;

        Eigen::MatrixXd Cell0Ds;
        Eigen::MatrixXi Cell1Ds;
        std::vector<Eigen::MatrixXi> Cell2Ds;
        std::vector<Gedim::MeshUtilities::Mesh3DPolyhedron> Cell3Ds;
    };

  private:
    struct Mesh_Face final
    {
        Eigen::MatrixXi Extremes;
        unsigned int Index;
        unsigned int Marker;
    };

    unsigned int InsertNewEdge(const unsigned int origin,
                               const unsigned int end,
                               std::map<std::pair<unsigned int, unsigned int>, unsigned int> &edges) const;
    unsigned int InsertNewFace(const std::vector<unsigned int> &face_vertices,
                               const std::vector<unsigned int> &face_edges,
                               std::map<std::pair<unsigned int, unsigned int>, Mesh_Face> &faces) const;
    VtkMesh3D ComputeVtkMesh3D(VtkMesh &originalMesh) const;

  public:
    VtkMeshInterface();
    ~VtkMeshInterface();

    VtkMesh3D ImportMesh3DFromFile(const std::string &vtkFilePath) const;
};
} // namespace Gedim

#endif // __OpenVolumeMeshInterface_H
