#ifndef __SphereMeshUtilities_H
#define __SphereMeshUtilities_H

#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "PlatonicSolid.hpp"
#include <numbers>
#include <numeric>

namespace Gedim
{
/// \brief MeshUtilities
/// \copyright See top level LICENSE file for details.
///
/// https://danielsieger.com/blog/2021/03/27/generating-spheres.html
class SphereMeshUtilities final
{
    const GeometryUtilities &geometryUtilities;
    const MeshUtilities &meshUtilities;

  public:
    SphereMeshUtilities(const GeometryUtilities &geometryUtilities, const MeshUtilities &meshUtilities)
        : geometryUtilities(geometryUtilities), meshUtilities(meshUtilities)
    {
    }

    virtual ~SphereMeshUtilities(){};

    GeometryUtilities::Polyhedron uv_sphere(const unsigned int &meridians, const unsigned int &parallels) const
    {
        Output::Assert(meridians >= 4 && parallels >= 2);

        GeometryUtilities::Polyhedron polyhedron;
        polyhedron.Vertices = Eigen::MatrixXd::Zero(3, 2 + (parallels - 1) * meridians);

        unsigned int v = 0;
        polyhedron.Vertices.col(v++) << 0.0, 1.0, 0.0; // add top vertex// add top vertex

        // generate vertices per stack / slice
        for (int i = 0; i < parallels - 1; i++)
        {
            const double phi = std::numbers::pi * double(i + 1) / double(parallels);
            for (int j = 0; j < meridians; j++)
            {
                const double theta = 2.0 * std::numbers::pi * double(j) / double(meridians);
                polyhedron.Vertices.col(v++) << std::sin(phi) * std::cos(theta), std::cos(phi), std::sin(phi) * std::sin(theta);
            }
        }

        // add bottom vertex
        polyhedron.Vertices.col(v++) << 0.0, -1.0, 0.0;

        const unsigned int idLastVertices = v - 1;
        std::vector<Eigen::VectorXi> faces;

        // add top / bottom triangles
        for (unsigned int i = 0; i < meridians; ++i)
        {
            Eigen::VectorXi face_top(3);
            face_top << 0, i + 1, (i + 1) % meridians + 1;
            faces.push_back(face_top);

            Eigen::VectorXi face_bottom(3);
            face_bottom << idLastVertices, i + meridians * (parallels - 2) + 1,
                (i + 1) % meridians + meridians * (parallels - 2) + 1;
            faces.push_back(face_bottom);
        }

        // add quads per stack / slice
        for (unsigned int j = 0; j < parallels - 2; j++)
        {
            const unsigned int j0 = j * meridians + 1;
            const unsigned int j1 = (j + 1) * meridians + 1;
            for (int i = 0; i < meridians; i++)
            {
                const unsigned int i0 = j0 + i;
                const unsigned int i1 = j0 + (i + 1) % meridians;
                const unsigned int i2 = j1 + (i + 1) % meridians;
                const unsigned int i3 = j1 + i;

                Eigen::VectorXi face(4);
                face << i0, i1, i2, i3;
                faces.push_back(face);
            }
        }

        Gedim::MeshUtilities::ComputeMesh2DCell1DsResult result = meshUtilities.ComputeMesh2DCell1Ds(polyhedron.Vertices, faces);

        polyhedron.Edges = result.Cell1Ds;
        polyhedron.Faces = result.Cell2Ds;

        return polyhedron;
    }
};

} // namespace Gedim

#endif // __SphereMeshUtilities_H
