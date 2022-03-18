#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateTetrahedronWithOrigin(const Eigen::Vector3d& origin,
                                                                               const Eigen::Vector3d& lengthVector,
                                                                               const Eigen::Vector3d& heightVector,
                                                                               const Eigen::Vector3d& widthVector) const
  {
    GeometryUtilities::Polyhedron tetrahedron;

    // create vertices
    tetrahedron.Vertices.resize(3, 4);
    tetrahedron.Vertices.col(0) << origin;
    tetrahedron.Vertices.col(1) << origin + lengthVector;
    tetrahedron.Vertices.col(2) << origin + widthVector;
    tetrahedron.Vertices.col(3) << origin + heightVector;

    // create edges
    tetrahedron.Edges.resize(2, 6);
    tetrahedron.Edges.col(0) << 0, 1;
    tetrahedron.Edges.col(1) << 0, 2;
    tetrahedron.Edges.col(2) << 1, 2;
    tetrahedron.Edges.col(3) << 0, 3;
    tetrahedron.Edges.col(4) << 1, 3;
    tetrahedron.Edges.col(5) << 2, 3;

    // create faces
    tetrahedron.Faces.resize(4, MatrixXi::Zero(2, 3));
    tetrahedron.Faces[0].row(0)<< 0, 1, 2;
    tetrahedron.Faces[1].row(0)<< 0, 1, 3;
    tetrahedron.Faces[2].row(0)<< 0, 2, 3;
    tetrahedron.Faces[3].row(0)<< 1, 2, 3;

    tetrahedron.Faces[0].row(1)<< 0, 2, 1;
    tetrahedron.Faces[1].row(1)<< 0, 4, 3;
    tetrahedron.Faces[2].row(1)<< 1, 5, 3;
    tetrahedron.Faces[3].row(1)<< 2, 5, 4;

    return tetrahedron;
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateTetrahedronWithVertices(const Eigen::Vector3d& v1,
                                                                                 const Eigen::Vector3d& v2,
                                                                                 const Eigen::Vector3d& v3,
                                                                                 const Eigen::Vector3d& v4) const
  {
    GeometryUtilities::Polyhedron tetrahedron;

    // create vertices
    tetrahedron.Vertices.resize(3, 4);
    tetrahedron.Vertices.col(0) << v1;
    tetrahedron.Vertices.col(1) << v2;
    tetrahedron.Vertices.col(2) << v3;
    tetrahedron.Vertices.col(3) << v4;

    // create edges
    tetrahedron.Edges.resize(2, 6);
    tetrahedron.Edges.col(0)<< 0, 1;
    tetrahedron.Edges.col(1)<< 0, 2;
    tetrahedron.Edges.col(2)<< 1, 2;
    tetrahedron.Edges.col(3)<< 0, 3;
    tetrahedron.Edges.col(4)<< 1, 3;
    tetrahedron.Edges.col(5)<< 2, 3;

    // create faces
    tetrahedron.Faces.resize(4, MatrixXi::Zero(2, 3));
    tetrahedron.Faces[0].row(0)<< 0, 1, 2;
    tetrahedron.Faces[1].row(0)<< 0, 1, 3;
    tetrahedron.Faces[2].row(0)<< 0, 2, 3;
    tetrahedron.Faces[3].row(0)<< 1, 2, 3;

    tetrahedron.Faces[0].row(1)<< 0, 2, 1;
    tetrahedron.Faces[1].row(1)<< 0, 4, 3;
    tetrahedron.Faces[2].row(1)<< 1, 5, 3;
    tetrahedron.Faces[3].row(1)<< 2, 5, 4;

    return tetrahedron;
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron GeometryUtilities::CreateCubeWithOrigin(const Eigen::Vector3d& origin,
                                                                        const Eigen::Vector3d& lengthVector,
                                                                        const Eigen::Vector3d& heightVector,
                                                                        const Eigen::Vector3d& widthVector) const
  {
    using namespace Eigen;

    Gedim::GeometryUtilities::Polyhedron cube;

    // create vertices
    cube.Vertices.setZero(3, 8);
    cube.Vertices.col(0)<< origin;
    cube.Vertices.col(1)<< origin + lengthVector;
    cube.Vertices.col(2)<< origin + lengthVector + widthVector;
    cube.Vertices.col(3)<< origin + widthVector;
    cube.Vertices.col(4)<< origin + heightVector;
    cube.Vertices.col(5)<< origin + heightVector + lengthVector;
    cube.Vertices.col(6)<< origin + heightVector + lengthVector + widthVector;
    cube.Vertices.col(7)<< origin + heightVector + widthVector;

    // create edges
    cube.Edges.setZero(2, 12);
    cube.Edges.col(0)<< 0, 1;
    cube.Edges.col(1)<< 1, 2;
    cube.Edges.col(2)<< 2, 3;
    cube.Edges.col(3)<< 3, 0;
    cube.Edges.col(4)<< 4, 5;
    cube.Edges.col(5)<< 5, 6;
    cube.Edges.col(6)<< 6, 7;
    cube.Edges.col(7)<< 7, 4;
    cube.Edges.col(8)<< 0, 4;
    cube.Edges.col(9)<< 1, 5;
    cube.Edges.col(10)<< 2, 6;
    cube.Edges.col(11)<< 3, 7;

    // create faces
    cube.Faces.resize(6, MatrixXi::Zero(2, 4));
    cube.Faces[0].row(0)<< 0, 1, 2, 3;
    cube.Faces[0].row(1)<< 0, 1, 2, 3;

    cube.Faces[1].row(0)<< 4, 5, 6, 7;
    cube.Faces[1].row(1)<< 4, 5, 6, 7;

    cube.Faces[2].row(0)<< 0, 3, 7, 4;
    cube.Faces[2].row(1)<< 3, 11, 7, 8;

    cube.Faces[3].row(0)<< 1, 2, 6, 5;
    cube.Faces[3].row(1)<< 1, 10, 5, 9;

    cube.Faces[4].row(0)<< 0, 1, 5, 4;
    cube.Faces[4].row(1)<< 0, 9, 4, 8;

    cube.Faces[5].row(0)<< 3, 2, 6, 7;
    cube.Faces[5].row(1)<< 2, 10, 6, 11;

    return cube;
  }
  // ***************************************************************************
}
