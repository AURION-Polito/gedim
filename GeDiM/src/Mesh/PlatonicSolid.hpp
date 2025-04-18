#ifndef __PlatonicSolid_H
#define __PlatonicSolid_H

#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include <numeric>

namespace Gedim
{
/// \brief MeshUtilities
/// \copyright See top level LICENSE file for details.
///
/// https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html
class PlatonicSolid final
{
    const GeometryUtilities &geometryUtilities;
    const MeshUtilities &meshUtilities;

  public:
    PlatonicSolid(const GeometryUtilities &geometryUtilities, const MeshUtilities &meshUtilities)
        : geometryUtilities(geometryUtilities), meshUtilities(meshUtilities)
    {
    }

    virtual ~PlatonicSolid() {};

    void project_to_unit_sphere(GeometryUtilities::Polyhedron &polyhedron) const
    {
        polyhedron.Vertices.colwise().normalize();
    }

    GeometryUtilities::Polyhedron dual_polyhedron(const GeometryUtilities::Polyhedron &polyhedron) const;

    GeometryUtilities::Polyhedron first_class_geodesic_polyhedron(const GeometryUtilities::Polyhedron &starting_polyhedron,
                                                                  const unsigned int &frequency) const;

    GeometryUtilities::Polyhedron goldberg_polyhedron(const GeometryUtilities::Polyhedron &geodesic_polyhedron) const
    {

        GeometryUtilities::Polyhedron polyhedron = dual_polyhedron(geodesic_polyhedron);
        project_to_unit_sphere(polyhedron);

        return polyhedron;
    }

    GeometryUtilities::Polyhedron tetrahedron() const
    {
        // choose coordinates on the unit sphere
        const double a = 1.0 / 3.0;
        const double b = sqrt(8.0 / 9.0);
        const double c = sqrt(2.0 / 9.0);
        const double d = sqrt(2.0 / 3.0);

        GeometryUtilities::Polyhedron polyhedron;

        // add the 4 vertices
        polyhedron.Vertices.resize(3, 4);
        polyhedron.Vertices.col(0) << 0.0, 0.0, 1.0;
        polyhedron.Vertices.col(1) << -c, d, -a;
        polyhedron.Vertices.col(2) << -c, -d, -a;
        polyhedron.Vertices.col(3) << b, 0.0, -a;

        // add the 6 edges
        polyhedron.Edges.resize(2, 6);
        polyhedron.Edges.col(0) << 0, 1; // 0
        polyhedron.Edges.col(1) << 1, 2; // 1
        polyhedron.Edges.col(2) << 0, 2; // 2
        polyhedron.Edges.col(3) << 2, 3; // 3
        polyhedron.Edges.col(4) << 0, 3; // 4
        polyhedron.Edges.col(5) << 1, 3; // 5

        // add the 4 faces
        polyhedron.Faces.resize(4);
        polyhedron.Faces[0].resize(2, 3);
        polyhedron.Faces[0] << 0, 1, 2, 0, 1, 2;

        polyhedron.Faces[1].resize(2, 3);
        polyhedron.Faces[1] << 0, 2, 3, 2, 3, 4;

        polyhedron.Faces[2].resize(2, 3);
        polyhedron.Faces[2] << 0, 3, 1, 4, 5, 0;

        polyhedron.Faces[3].resize(2, 3);
        polyhedron.Faces[3] << 3, 2, 1, 3, 1, 5;

        return polyhedron;
    }

    GeometryUtilities::Polyhedron hexahedron() const
    {

        // choose coordinates on the unit sphere
        const double a = 1.0 / sqrt(3.0);

        GeometryUtilities::Polyhedron polyhedron;

        // add the 8 vertices
        polyhedron.Vertices.resize(3, 8);
        polyhedron.Vertices.col(0) << -a, -a, -a;
        polyhedron.Vertices.col(1) << a, -a, -a;
        polyhedron.Vertices.col(2) << a, a, -a;
        polyhedron.Vertices.col(3) << -a, a, -a;
        polyhedron.Vertices.col(4) << -a, -a, a;
        polyhedron.Vertices.col(5) << a, -a, a;
        polyhedron.Vertices.col(6) << a, a, a;
        polyhedron.Vertices.col(7) << -a, a, a;

        // add the 12 edges
        polyhedron.Edges.resize(2, 12);
        polyhedron.Edges.col(0) << 0, 1;  // 0
        polyhedron.Edges.col(1) << 1, 2;  // 1
        polyhedron.Edges.col(2) << 2, 3;  // 2
        polyhedron.Edges.col(3) << 0, 3;  // 3
        polyhedron.Edges.col(4) << 2, 6;  // 4
        polyhedron.Edges.col(5) << 5, 6;  // 5
        polyhedron.Edges.col(6) << 1, 5;  // 6
        polyhedron.Edges.col(7) << 6, 7;  // 7
        polyhedron.Edges.col(8) << 4, 7;  // 8
        polyhedron.Edges.col(9) << 4, 5;  // 9
        polyhedron.Edges.col(10) << 0, 4; // 10
        polyhedron.Edges.col(11) << 3, 7; // 11

        // add the 6 faces
        polyhedron.Faces.resize(6);
        polyhedron.Faces[0].resize(2, 4);
        polyhedron.Faces[0] << 3, 2, 1, 0, 2, 1, 0, 3;

        polyhedron.Faces[1].resize(2, 4);
        polyhedron.Faces[1] << 2, 6, 5, 1, 4, 5, 6, 1;

        polyhedron.Faces[2].resize(2, 4);
        polyhedron.Faces[2] << 5, 6, 7, 4, 5, 7, 8, 9;

        polyhedron.Faces[3].resize(2, 4);
        polyhedron.Faces[3] << 0, 4, 7, 3, 10, 8, 11, 3;

        polyhedron.Faces[4].resize(2, 4);
        polyhedron.Faces[4] << 3, 7, 6, 2, 11, 7, 4, 2;

        polyhedron.Faces[5].resize(2, 4);
        polyhedron.Faces[5] << 1, 5, 4, 0, 6, 9, 10, 0;

        return polyhedron;
    }

    GeometryUtilities::Polyhedron octahedron() const
    {
        const GeometryUtilities::Polyhedron polyhedron = hexahedron();
        GeometryUtilities::Polyhedron dual = dual_polyhedron(polyhedron);
        project_to_unit_sphere(dual);
        return dual;
    }

    GeometryUtilities::Polyhedron icosahedron() const
    {
        GeometryUtilities::Polyhedron polyhedron;

        const double phi = (1.0 + sqrt(5.0)) * 0.5; // golden ratio
        const double a = 1.0;
        const double b = 1.0 / phi;

        // add the 12 vertices
        polyhedron.Vertices.resize(3, 12);
        polyhedron.Vertices.col(0) << 0.0, b, -a;
        polyhedron.Vertices.col(1) << b, a, 0.0;
        polyhedron.Vertices.col(2) << -b, a, 0.0;
        polyhedron.Vertices.col(3) << 0.0, b, a;
        polyhedron.Vertices.col(4) << 0.0, -b, a;
        polyhedron.Vertices.col(5) << -a, 0.0, b;
        polyhedron.Vertices.col(6) << 0.0, -b, -a;
        polyhedron.Vertices.col(7) << a, 0.0, -b;
        polyhedron.Vertices.col(8) << a, 0.0, b;
        polyhedron.Vertices.col(9) << -a, 0.0, -b;
        polyhedron.Vertices.col(10) << b, -a, 0.0;
        polyhedron.Vertices.col(11) << -b, -a, 0.0;

        // add the 12 edges
        polyhedron.Edges.resize(2, 30);
        polyhedron.Edges.col(0) << 1, 2;    // 0
        polyhedron.Edges.col(1) << 0, 1;    // 1
        polyhedron.Edges.col(2) << 0, 2;    // 2
        polyhedron.Edges.col(3) << 2, 3;    // 3
        polyhedron.Edges.col(4) << 1, 3;    // 4
        polyhedron.Edges.col(5) << 4, 5;    // 5
        polyhedron.Edges.col(6) << 3, 4;    // 6
        polyhedron.Edges.col(7) << 3, 5;    // 7
        polyhedron.Edges.col(8) << 4, 8;    // 8
        polyhedron.Edges.col(9) << 3, 8;    // 9
        polyhedron.Edges.col(10) << 6, 7;   // 10
        polyhedron.Edges.col(11) << 0, 6;   // 11
        polyhedron.Edges.col(12) << 0, 7;   // 12
        polyhedron.Edges.col(13) << 6, 9;   // 13
        polyhedron.Edges.col(14) << 0, 9;   // 14
        polyhedron.Edges.col(15) << 10, 11; // 15
        polyhedron.Edges.col(16) << 6, 11;  // 16
        polyhedron.Edges.col(17) << 6, 10;  // 17
        polyhedron.Edges.col(18) << 4, 10;  // 18
        polyhedron.Edges.col(19) << 4, 11;  // 19
        polyhedron.Edges.col(20) << 5, 9;   // 20
        polyhedron.Edges.col(21) << 2, 5;   // 21
        polyhedron.Edges.col(22) << 2, 9;   // 22
        polyhedron.Edges.col(23) << 9, 11;  // 23
        polyhedron.Edges.col(24) << 5, 11;  // 24
        polyhedron.Edges.col(25) << 7, 8;   // 25
        polyhedron.Edges.col(26) << 1, 7;   // 26
        polyhedron.Edges.col(27) << 1, 8;   // 27
        polyhedron.Edges.col(28) << 8, 10;  // 28
        polyhedron.Edges.col(29) << 7, 10;  // 29

        // add the 6 faces
        polyhedron.Faces.resize(20);
        polyhedron.Faces[0].resize(2, 3);
        polyhedron.Faces[0] << 2, 1, 0, 0, 1, 2;

        polyhedron.Faces[1].resize(2, 3);
        polyhedron.Faces[1] << 1, 2, 3, 0, 3, 4;

        polyhedron.Faces[2].resize(2, 3);
        polyhedron.Faces[2] << 5, 4, 3, 5, 6, 7;

        polyhedron.Faces[3].resize(2, 3);
        polyhedron.Faces[3] << 4, 8, 3, 8, 9, 6;

        polyhedron.Faces[4].resize(2, 3);
        polyhedron.Faces[4] << 7, 6, 0, 10, 11, 12;

        polyhedron.Faces[5].resize(2, 3);
        polyhedron.Faces[5] << 6, 9, 0, 13, 14, 11;

        polyhedron.Faces[6].resize(2, 3);
        polyhedron.Faces[6] << 11, 10, 4, 15, 18, 19;

        polyhedron.Faces[7].resize(2, 3);
        polyhedron.Faces[7] << 10, 11, 6, 15, 16, 17;

        polyhedron.Faces[8].resize(2, 3);
        polyhedron.Faces[8] << 9, 5, 2, 20, 21, 22;

        polyhedron.Faces[9].resize(2, 3);
        polyhedron.Faces[9] << 5, 9, 11, 20, 23, 24;

        polyhedron.Faces[10].resize(2, 3);
        polyhedron.Faces[10] << 8, 7, 1, 25, 26, 27;

        polyhedron.Faces[11].resize(2, 3);
        polyhedron.Faces[11] << 7, 8, 10, 25, 28, 29;

        polyhedron.Faces[12].resize(2, 3);
        polyhedron.Faces[12] << 2, 5, 3, 21, 7, 3;

        polyhedron.Faces[13].resize(2, 3);
        polyhedron.Faces[13] << 8, 1, 3, 27, 4, 9;

        polyhedron.Faces[14].resize(2, 3);
        polyhedron.Faces[14] << 9, 2, 0, 22, 2, 14;

        polyhedron.Faces[15].resize(2, 3);
        polyhedron.Faces[15] << 1, 7, 0, 26, 12, 1;

        polyhedron.Faces[16].resize(2, 3);
        polyhedron.Faces[16] << 11, 9, 6, 23, 13, 16;

        polyhedron.Faces[17].resize(2, 3);
        polyhedron.Faces[17] << 7, 10, 6, 29, 17, 10;

        polyhedron.Faces[18].resize(2, 3);
        polyhedron.Faces[18] << 5, 11, 4, 24, 19, 5;

        polyhedron.Faces[19].resize(2, 3);
        polyhedron.Faces[19] << 10, 8, 4, 28, 8, 18;

        project_to_unit_sphere(polyhedron);

        return polyhedron;
    }

    GeometryUtilities::Polyhedron dodecahedron() const
    {
        const GeometryUtilities::Polyhedron polyhedron = icosahedron();
        GeometryUtilities::Polyhedron dual = dual_polyhedron(polyhedron);
        project_to_unit_sphere(dual);
        return dual;
    }
};

} // namespace Gedim

#endif // __PlatonicSolid_H
