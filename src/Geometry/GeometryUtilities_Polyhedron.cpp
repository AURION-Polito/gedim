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
  GeometryUtilities::Polyhedron GeometryUtilities::CreateParallelepipedWithOrigin(const Eigen::Vector3d& origin,
                                                                                  const Eigen::Vector3d& lengthVector,
                                                                                  const Eigen::Vector3d& heightVector,
                                                                                  const Eigen::Vector3d& widthVector) const
  {
    Gedim::GeometryUtilities::Polyhedron parallelpiped;

    // create vertices
    parallelpiped.Vertices.setZero(3, 8);
    parallelpiped.Vertices.col(0)<< origin;
    parallelpiped.Vertices.col(1)<< origin + lengthVector;
    parallelpiped.Vertices.col(2)<< origin + lengthVector + widthVector;
    parallelpiped.Vertices.col(3)<< origin + widthVector;
    parallelpiped.Vertices.col(4)<< origin + heightVector;
    parallelpiped.Vertices.col(5)<< origin + heightVector + lengthVector;
    parallelpiped.Vertices.col(6)<< origin + heightVector + lengthVector + widthVector;
    parallelpiped.Vertices.col(7)<< origin + heightVector + widthVector;

    // create edges
    parallelpiped.Edges.setZero(2, 12);
    parallelpiped.Edges.col(0)<< 0, 1;
    parallelpiped.Edges.col(1)<< 1, 2;
    parallelpiped.Edges.col(2)<< 2, 3;
    parallelpiped.Edges.col(3)<< 3, 0;
    parallelpiped.Edges.col(4)<< 4, 5;
    parallelpiped.Edges.col(5)<< 5, 6;
    parallelpiped.Edges.col(6)<< 6, 7;
    parallelpiped.Edges.col(7)<< 7, 4;
    parallelpiped.Edges.col(8)<< 0, 4;
    parallelpiped.Edges.col(9)<< 1, 5;
    parallelpiped.Edges.col(10)<< 2, 6;
    parallelpiped.Edges.col(11)<< 3, 7;

    // create faces
    parallelpiped.Faces.resize(6, MatrixXi::Zero(2, 4));
    parallelpiped.Faces[0].row(0)<< 0, 1, 2, 3;
    parallelpiped.Faces[0].row(1)<< 0, 1, 2, 3;

    parallelpiped.Faces[1].row(0)<< 4, 5, 6, 7;
    parallelpiped.Faces[1].row(1)<< 4, 5, 6, 7;

    parallelpiped.Faces[2].row(0)<< 0, 3, 7, 4;
    parallelpiped.Faces[2].row(1)<< 3, 11, 7, 8;

    parallelpiped.Faces[3].row(0)<< 1, 2, 6, 5;
    parallelpiped.Faces[3].row(1)<< 1, 10, 5, 9;

    parallelpiped.Faces[4].row(0)<< 0, 1, 5, 4;
    parallelpiped.Faces[4].row(1)<< 0, 9, 4, 8;

    parallelpiped.Faces[5].row(0)<< 3, 2, 6, 7;
    parallelpiped.Faces[5].row(1)<< 2, 10, 6, 11;

    return parallelpiped;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::PolyhedronFaceVertices(const Eigen::MatrixXd& polyhedronVertices,
                                                             const Eigen::MatrixXi& polyhedronEdges,
                                                             const vector<Eigen::MatrixXi> polyhedronFaces) const
  {
    vector<Eigen::MatrixXd> facesVertices(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceVerticesIndices = polyhedronFaces[f].row(0);
      Eigen::MatrixXd& faceVertices = facesVertices[f];
      faceVertices.setZero(3, faceVerticesIndices.size());
      for (unsigned int v = 0; v < faceVerticesIndices.size(); v++)
        faceVertices.col(v)<< polyhedronVertices.col(faceVerticesIndices[v]);
    }

    return facesVertices;
  }
  // ***************************************************************************
  vector<Matrix3d> GeometryUtilities::PolyhedronFaceRotationMatrices(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                     const vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                     const vector<Eigen::Vector3d>& polyhedronFaceTranslations) const
  {
    vector<Matrix3d> rotationMatrices;
    rotationMatrices.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
    {
      rotationMatrices.push_back(PolygonRotationMatrix(polyhedronFaceVertices[f],
                                                       polyhedronFaceNormals[f],
                                                       polyhedronFaceTranslations[f]));
    }

    return rotationMatrices;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceTranslations(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> translations;
    translations.reserve(polyhedronFaceVertices.size());
    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      translations.push_back(PolygonTranslation(polyhedronFaceVertices[f]));

    return translations;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceNormals(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                            const Eigen::Vector3d& pointInsidePolyhedron) const
  {
    vector<Vector3d> faceNormals;
    faceNormals.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
    {
      const Eigen::Vector3d normal = PolygonNormal(polyhedronFaceVertices[f]);
      const PointPlanePositionTypes pointFacePosition = PointPlanePosition(pointInsidePolyhedron,
                                                                           normal,
                                                                           polyhedronFaceVertices[f].col(0));
      Output::Assert(pointFacePosition == PointPlanePositionTypes::Negative ||
                     pointFacePosition == PointPlanePositionTypes::Positive);

      const double normalOutgoingDirection = (pointFacePosition == PointPlanePositionTypes::Positive) ? -1.0 :
                                                                                                        1.0;

      faceNormals.push_back(normalOutgoingDirection * normal);
    }

    return faceNormals;
  }
  // ***************************************************************************
}
