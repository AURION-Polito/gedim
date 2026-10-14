#ifndef __program_configuration_hpp
#define __program_configuration_hpp

#include "Configurations.hpp"
#include "GeometryUtilities.hpp"

enum struct InternalPointsType
{
    Equispaced = 1,
    GaussLobatto = 2
};

struct IDomain
{
    // The first represents the domain, the others represent holes.
    std::vector<Gedim::GeometryUtilities::Polyhedron> domains;
    std::vector<Eigen::Vector3d> internal_points; // point internal to the domain
    double volume;
};

enum struct TestDomainTypes
{
    SquareDomainWithOneHole = 1
};

struct SquareDomainWithOneHole : IDomain
{
    SquareDomainWithOneHole()
    {
        this->domains.resize(2);
        domains[0].Vertices = Eigen::MatrixXd::Zero(3, 8);
        domains[0].Vertices << -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0,
            -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0;

        domains[0].Edges = Eigen::MatrixXi::Zero(2, 12);
        domains[0].Edges << 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7;

        domains[0].Faces.resize(6, Eigen::MatrixXi::Zero(2, 4));

        domains[0].Faces[0] << 0, 1, 2, 3, 0, 1, 2, 3;
        domains[0].Faces[1] << 4, 5, 6, 7, 4, 5, 6, 7;
        domains[0].Faces[2] << 0, 3, 7, 4, 3, 11, 7, 8;
        domains[0].Faces[3] << 1, 2, 6, 5, 1, 10, 5, 9;
        domains[0].Faces[4] << 0, 1, 5, 4, 0, 9, 4, 8;
        domains[0].Faces[5] << 3, 2, 6, 7, 2, 10, 6, 11;

        domains[1].Vertices = Eigen::MatrixXd::Zero(3, 8);
        domains[1].Vertices << -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5,
            -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5;

        domains[1].Edges = Eigen::MatrixXi::Zero(2, 12);
        domains[1].Edges << 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7;

        domains[1].Faces.resize(6, Eigen::MatrixXi::Zero(2, 4));

        domains[1].Faces[0] << 0, 1, 2, 3, 0, 1, 2, 3;
        domains[1].Faces[1] << 4, 5, 6, 7, 4, 5, 6, 7;
        domains[1].Faces[2] << 0, 3, 7, 4, 3, 11, 7, 8;
        domains[1].Faces[3] << 1, 2, 6, 5, 1, 10, 5, 9;
        domains[1].Faces[4] << 0, 1, 5, 4, 0, 9, 4, 8;
        domains[1].Faces[5] << 3, 2, 6, 7, 2, 10, 6, 11;

        internal_points.resize(2);

        internal_points[0] << 0.0, 0.0, -0.7;
        internal_points[1] << 0.0, 0.0, 0.0;

        volume = 8.0;
    }
};

enum class MeshGenerator_Types_3D
{
    PureTetrahedral = 0, ///< generated tetrahedral mesh which remains tetrahedral after excavation
    Tetrahedral = 1,     ///< generated tetrahedral mesh
    Minimal = 2,         ///< generated minimal mesh
    Polyhedral = 3,      ///< generated voronoi polyhedral mesh
    Cubic = 4            ///< cubic mesh
};

struct Program_configuration final
{
    Program_configuration()
    {
        Gedim::Configurations::AddProperty("TestDomain",
                                           static_cast<unsigned int>(TestDomainTypes::SquareDomainWithOneHole),
                                           "Test Type: 1 - SquareDomainWithOneHole (Default: 1)");

        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Export", "Folder where to export data (Default: ./Export)");
        Gedim::Configurations::AddProperty("ExportFormat",
                                           std::vector<unsigned int>({1, 0, 0, 0}),
                                           "A boolean vector of kind of desired export type for solution and mesh "
                                           "[Csv, Vtu, UCD, Medit] (Default: [1,0,0,0])");

        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshGenerator",
                                           static_cast<unsigned int>(MeshGenerator_Types_3D::PureTetrahedral),
                                           "Mesh 3D gereator type, 0 - PureTetrahedral; 1 - Tetrahedral; 2 - "
                                           "Minimal; 3 - Polyhedral; 4 - Cubic (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxVolume", 0.1, "Mesh 3D maximum relative cell volume for tetrahedral and voronoi mesh (Default: 0.1)");
        Gedim::Configurations::AddProperty("NumCellsHeight",
                                           static_cast<unsigned int>(12),
                                           "Mesh 3D number of cells Height for parallelpipedal mesh (Default: 12)");
        Gedim::Configurations::AddProperty("NumCellsBase",
                                           static_cast<unsigned int>(12),
                                           "Mesh 3D number of cells Base for parallelpipedal mesh (Default: 12)");
        Gedim::Configurations::AddProperty("NumCellsWidth",
                                           static_cast<unsigned int>(12),
                                           "Mesh 3D number of cells Width for parallelpipedal mesh (Default: 12)");
        Gedim::Configurations::AddProperty("GeometricTolerance1D", 1.0e-12, "Geometric Tolerance 1D (Default: 1.0e-12)");
        Gedim::Configurations::AddProperty("GeometricTolerance2D", 1.0e-14, "Geometric Tolerance 2D (Default: 1.0e-14)");
        Gedim::Configurations::AddProperty("GeometricTolerance3D", 1.0e-15, "Geometric Tolerance 3D (Default: 1.0e-15)");

        // Add internal points to edge to simulate edge dofs
        Gedim::Configurations::AddProperty("MethodOrder", static_cast<unsigned int>(1), "Method order (Default: 1)");
        Gedim::Configurations::AddProperty("InternalPointsType",
                                           static_cast<unsigned int>(InternalPointsType::Equispaced),
                                           "Internal Points Type:; 1 - Equispaced; 2 - Lobatto (Default: 1)");
    }

    inline TestDomainTypes TestDomain() const
    {
        return (TestDomainTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("TestDomain");
    }

    inline InternalPointsType InternalPoints() const
    {
        return (InternalPointsType)Gedim::Configurations::GetPropertyValue<unsigned int>("InternalPointsType");
    }

    inline std::string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("ExportFolder");
    }

    inline MeshGenerator_Types_3D MeshGenerator() const
    {
        return (MeshGenerator_Types_3D)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator");
    }

    inline std::string MeshImportFilePath() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("MeshImportFilePath");
    }
    inline double MeshMaxVolume() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("MeshMaxVolume");
    }
    inline double GeometricTolerance1D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D");
    }
    inline double GeometricTolerance2D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance2D");
    }
    inline double GeometricTolerance3D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance3D");
    }

    inline unsigned int MethodOrder() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("MethodOrder");
    }

    inline unsigned int NumCellsHeight() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("NumCellsHeight");
    }

    inline unsigned int NumCellsBase() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("NumCellsBase");
    }

    inline unsigned int NumCellsWidth() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("NumCellsWidth");
    }

    inline std::vector<unsigned int> ExportFormat() const
    {
        return Gedim::Configurations::GetPropertyValue<std::vector<unsigned int>>("ExportFormat");
    }
};

#endif // __program_configuration_hpp
