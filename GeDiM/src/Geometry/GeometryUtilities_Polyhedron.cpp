#include "IOUtilities.hpp"
#include "GeometryUtilities.hpp"
#include "MapTetrahedron.hpp"
#include "MapTriangle.hpp"
#include "VTKUtilities.hpp"
#include <queue>
#include "CommonUtilities.hpp"

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
  GeometryUtilities::Polyhedron GeometryUtilities::CreatePolyhedronWithExtrusion(const Eigen::MatrixXd& polygonVertices,
                                                                                 const std::vector<Eigen::Vector3d>& heightVectors) const
  {
    Output::Assert(static_cast<unsigned int>(polygonVertices.cols()) ==
                   static_cast<unsigned int>(heightVectors.size()) &&
                   polygonVertices.rows() == 3);

    Gedim::GeometryUtilities::Polyhedron polyhedron;

    const unsigned int numPolygonVertices = polygonVertices.cols();

    // create vertices
    const unsigned int numPolyhedronVertices = 2 * numPolygonVertices;
    polyhedron.Vertices.setZero(3, numPolyhedronVertices);
    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      polyhedron.Vertices.col(v)<< polygonVertices.col(v);
      polyhedron.Vertices.col(numPolygonVertices + v)<< polygonVertices.col(v) + heightVectors[v];
    }

    // create edges
    const unsigned int numPolyhedronEdges = 3 * numPolygonVertices;
    polyhedron.Edges.setZero(2, numPolyhedronEdges);
    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      polyhedron.Edges.col(v)<< v, (v + 1) % numPolygonVertices;
      polyhedron.Edges.col(numPolygonVertices + v)<< numPolygonVertices + v, numPolygonVertices + (v + 1) % numPolygonVertices;
      polyhedron.Edges.col(2 * numPolygonVertices + v)<< v, numPolygonVertices + v;
    }

    // create faces
    const unsigned int numPolyhedronFaces = 2 + numPolygonVertices;
    polyhedron.Faces.reserve(numPolyhedronFaces);

    polyhedron.Faces.push_back(MatrixXi::Zero(2, numPolygonVertices));
    polyhedron.Faces.push_back(MatrixXi::Zero(2, numPolygonVertices));
    for (unsigned int v = 0; v < numPolygonVertices; v++)
    {
      polyhedron.Faces[0](0, v) = v;
      polyhedron.Faces[0](1, v) = v;
      polyhedron.Faces[1](0, v) = numPolygonVertices + v;
      polyhedron.Faces[1](1, v) = numPolygonVertices + v;

      polyhedron.Faces.push_back(MatrixXi::Zero(2, 4));
      polyhedron.Faces[2 + v].row(0)<< v, (v + 1) % numPolygonVertices, numPolygonVertices + (v + 1) % numPolygonVertices, numPolygonVertices + v;
      polyhedron.Faces[2 + v].row(1)<< v, 2 * numPolygonVertices + (v + 1) % numPolygonVertices, numPolygonVertices + v, 2 * numPolygonVertices + v;
    }

    return polyhedron;

  }
  // ***************************************************************************
  double GeometryUtilities::PolyhedronVolumeByBoundaryIntegral(const std::vector<std::vector<Eigen::Matrix3d> >& polyhedronFaceRotatedTriangulationPoints,
                                                               const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                               const std::vector<bool>& polyhedronFaceNormalDirections,
                                                               const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                               const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                                               const MatrixXd& referenceQuadraturePoints,
                                                               const VectorXd& referenceQuadratureWeights) const
  {
    double volume = 0.0;
    const unsigned int numFaces = polyhedronFaceRotatedTriangulationPoints.size();

    const Gedim::MapTriangle mapTriangle;
    for (unsigned int f = 0; f < numFaces; f++)
    {
      const unsigned int numFaceTriangles = polyhedronFaceRotatedTriangulationPoints[f].size();

      const Eigen::Matrix3d& faceRotationMatrix = polyhedronFaceRotationMatrices[f];
      const Eigen::Vector3d& faceTranslation = polyhedronFaceTranslations[f];
      const Eigen::Vector3d& faceNormal = polyhedronFaceNormalDirections[f] ? polyhedronFaceNormals[f] :
                                                                              -1.0 * polyhedronFaceNormals[f];

      for (unsigned int t = 0; t < numFaceTriangles; t++)
      {
        const Eigen::Matrix3d& face2DTriangle = polyhedronFaceRotatedTriangulationPoints[f][t];

        const Gedim::MapTriangle::MapTriangleData mapData = mapTriangle.Compute(face2DTriangle);

        const MatrixXd faceWeightsTimesNormal = (mapTriangle.DetJ(mapData,
                                                                  referenceQuadratureWeights).array().abs() *
                                                 referenceQuadratureWeights.array()).matrix() *
                                                faceNormal.transpose();
        const Eigen::MatrixXd faceQuadraturePoints =
            faceRotationMatrix * mapData.B * referenceQuadraturePoints +
            faceRotationMatrix * mapData.b +
            faceTranslation;

        volume += 1.0 / 3.0 * (faceQuadraturePoints.row(0) * faceWeightsTimesNormal.col(0) +
                               faceQuadraturePoints.row(1) * faceWeightsTimesNormal.col(1) +
                               faceQuadraturePoints.row(2) * faceWeightsTimesNormal.col(2)).sum();
      }
    }

    return volume;
  }
  // ***************************************************************************
  double GeometryUtilities::PolyhedronVolumeByInternalIntegral(const std::vector<Eigen::MatrixXd>& polyhedronTetrahedronVertices,
                                                               const Eigen::VectorXd& referenceQuadratureWeights) const
  {
    const unsigned int numPolyhedronTetra = polyhedronTetrahedronVertices.size();

    const unsigned int numReferenceQuadraturePoints = referenceQuadratureWeights.size();
    const unsigned int numQuadraturePoints = numPolyhedronTetra *
                                             numReferenceQuadraturePoints;

    Eigen::VectorXd quadratureWeights = Eigen::VectorXd::Zero(numQuadraturePoints);

    const Gedim::MapTetrahedron mapTetra(*this);

    for (unsigned int t = 0; t < numPolyhedronTetra; t++)
    {
      const Eigen::MatrixXd& tetraVertices = polyhedronTetrahedronVertices[t];

      const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapTetra.Compute(tetraVertices);
      quadratureWeights.segment(numReferenceQuadraturePoints * t,
                                numReferenceQuadraturePoints) = referenceQuadratureWeights.array() *
                                                                mapTetra.DetJ(mapData,
                                                                              referenceQuadratureWeights).array().abs();

    }

    return quadratureWeights.sum();
  }
  // ***************************************************************************
  Eigen::Vector3d GeometryUtilities::PolyhedronCentroid(const std::vector<std::vector<Eigen::Matrix3d> >& polyhedronFaceRotatedTriangulationPoints,
                                                        const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                        const std::vector<bool>& polyhedronFaceNormalDirections,
                                                        const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                        const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices,
                                                        const double& polyhedronVolume) const
  {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero(3);
    const unsigned int numFaces = polyhedronFaceRotatedTriangulationPoints.size();

    Eigen::Matrix3d quadraturePoints;
    Eigen::Vector3d quadratureWeights;
    quadraturePoints.setZero();
    quadraturePoints(0,0) = 6.666666666666666297e-01;
    quadraturePoints(1,0) = 1.666666666666666574e-01;
    quadraturePoints(0,1) = 1.666666666666666574e-01;
    quadraturePoints(1,1) = 6.666666666666666297e-01;
    quadraturePoints(0,2) = 1.666666666666666574e-01;
    quadraturePoints(1,2) = 1.666666666666666574e-01;
    quadratureWeights[0] = 1.666666666666666574e-01;
    quadratureWeights[1] = 1.666666666666666574e-01;
    quadratureWeights[2] = 1.666666666666666574e-01;

    for (unsigned int f = 0; f < numFaces; f++)
    {
      const unsigned int numFaceTriangles = polyhedronFaceRotatedTriangulationPoints[f].size();
      const Eigen::Matrix3d& faceRotationMatrix = polyhedronFaceRotationMatrices[f];
      const Eigen::Vector3d& faceTranslation = polyhedronFaceTranslations[f];
      const Eigen::Vector3d& faceNormal = polyhedronFaceNormalDirections[f] ? polyhedronFaceNormals[f] :
                                                                              -1.0 * polyhedronFaceNormals[f];

      for (unsigned int t = 0; t < numFaceTriangles; t++)
      {
        const Eigen::Matrix3d& face2DTriangle = polyhedronFaceRotatedTriangulationPoints[f][t];

        Gedim::MapTriangle mapping;
        const Gedim::MapTriangle::MapTriangleData mapData = mapping.Compute(face2DTriangle);

        const Eigen::Matrix3d mappedQuadraturePoints = mapping.F(mapData,
                                                                 quadraturePoints);
        const Eigen::Matrix3d quadraturePoints3D = RotatePointsFrom2DTo3D(mappedQuadraturePoints,
                                                                          faceRotationMatrix,
                                                                          faceTranslation);

        centroid(0) += quadraturePoints3D.row(0) *
                       (faceNormal(0) * quadratureWeights * mapData.B.determinant()).asDiagonal() *
                       quadraturePoints3D.row(0).transpose();
        centroid(1) += quadraturePoints3D.row(1) *
                       (faceNormal(1) * quadratureWeights * mapData.B.determinant()).asDiagonal() *
                       quadraturePoints3D.row(1).transpose();
        centroid(2) += quadraturePoints3D.row(2) *
                       (faceNormal(2) * quadratureWeights * mapData.B.determinant()).asDiagonal() *
                       quadraturePoints3D.row(2).transpose();
      }
    }

    return 0.5 * centroid / polyhedronVolume;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolyhedronEdgesCentroid(const Eigen::MatrixXd& polyhedronVertices,
                                                      const Eigen::MatrixXi& polyhedronEdges) const
  {
    MatrixXd edgesCentroid(3, polyhedronEdges.cols());

    for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      edgesCentroid.col(e) = SegmentBarycenter(polyhedronVertices.col(polyhedronEdges(0, e)),
                                               polyhedronVertices.col(polyhedronEdges(1, e)));
    }

    return edgesCentroid;
  }
  // ***************************************************************************
  VectorXd GeometryUtilities::PolyhedronEdgesLength(const Eigen::MatrixXd& polyhedronVertices,
                                                    const Eigen::MatrixXi& polyhedronEdges) const
  {
    VectorXd edgesLength(polyhedronEdges.cols());

    for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      edgesLength[e] = SegmentLength(polyhedronVertices.col(polyhedronEdges(0, e)),
                                     polyhedronVertices.col(polyhedronEdges(1, e)));
    }

    return edgesLength;
  }
  // ***************************************************************************
  MatrixXd GeometryUtilities::PolyhedronEdgeTangents(const Eigen::MatrixXd& polyhedronVertices,
                                                     const Eigen::MatrixXi& polyhedronEdges) const
  {
    MatrixXd edgeTangents(3, polyhedronEdges.cols());

    for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      edgeTangents.col(e) = SegmentTangent(polyhedronVertices.col(polyhedronEdges(0, e)),
                                           polyhedronVertices.col(polyhedronEdges(1, e)));
    }

    return edgeTangents;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::PolyhedronFaceVertices(const Eigen::MatrixXd& polyhedronVertices,
                                                             const vector<Eigen::MatrixXi>& polyhedronFaces) const
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
  vector<vector<bool>> GeometryUtilities::PolyhedronFaceEdgeDirections(const Eigen::MatrixXd& ,
                                                                       const Eigen::MatrixXi& polyhedronEdges,
                                                                       const vector<Eigen::MatrixXi>& polyhedronFaces) const
  {
    vector<vector<bool>> facesEdgeDirections(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceEdgesIndices = polyhedronFaces[f].row(1);
      facesEdgeDirections[f].resize(faceEdgesIndices.size(), true);
      for (unsigned int e = 0; e < faceEdgesIndices.size(); e++)
        facesEdgeDirections[f][e] = (polyhedronEdges(0, faceEdgesIndices[e]) == polyhedronFaces[f](0, e));
    }

    return facesEdgeDirections;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::PolyhedronFaceEdgeTangents(const Eigen::MatrixXd& ,
                                                                 const Eigen::MatrixXi& ,
                                                                 const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                 const vector<vector<bool>>& polyhedronFaceEdgeDirections,
                                                                 const Eigen::MatrixXd& polyhedronEdgeTangents) const
  {
    vector<Eigen::MatrixXd> facesEdgeTangents(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      const Eigen::VectorXi& faceEdgesIndices = polyhedronFaces[f].row(1);
      Eigen::MatrixXd& faceEdgeTangents = facesEdgeTangents[f];
      faceEdgeTangents.setZero(3, faceEdgesIndices.size());
      for (unsigned int e = 0; e < faceEdgesIndices.size(); e++)
        faceEdgeTangents.col(e)<< (polyhedronFaceEdgeDirections[f][e] ? 1.0 : -1.0) *
                                  polyhedronEdgeTangents.col(faceEdgesIndices[e]);
    }

    return facesEdgeTangents;
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
  std::vector<MatrixXd> GeometryUtilities::PolyhedronFaceRotatedVertices(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                         const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                         const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const
  {
    const unsigned int numFaces = polyhedronFaceVertices.size();
    std::vector<MatrixXd> faceRotatedVertices(numFaces);

    for (unsigned int f = 0; f < numFaces; f++)
    {
      faceRotatedVertices[f] = RotatePointsFrom3DTo2D(polyhedronFaceVertices[f],
                                                      polyhedronFaceRotationMatrices[f].transpose(),
                                                      polyhedronFaceTranslations[f]);
    }

    return faceRotatedVertices;
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int>> GeometryUtilities::PolyhedronFacesUnalignedVertices(const std::vector<Eigen::MatrixXd>& polyhedronFacesRotatedVertices) const
  {
    const unsigned int numFaces = polyhedronFacesRotatedVertices.size();
    std::vector<std::vector<unsigned int>> faceUnalignedVertices(numFaces);

    for (unsigned int f = 0; f < numFaces; f++)
      faceUnalignedVertices[f] = UnalignedPoints(polyhedronFacesRotatedVertices[f]);

    return faceUnalignedVertices;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceNormals(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> faceNormals;
    faceNormals.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      faceNormals.push_back(PolygonNormal(polyhedronFaceVertices[f]));

    return faceNormals;
  }
  // ***************************************************************************
  std::vector<std::array<Eigen::Vector3d, 2>> GeometryUtilities::PolyhedronFaceTangents(const std::vector<Eigen::MatrixXd>& polyhedronFacesVertices,
                                                                                        const std::vector<Eigen::Vector3d>& polyhedronFacesNormal,
                                                                                        const std::vector<bool>& polyhedronFacesNormalDirection) const
  {
    vector<std::array<Eigen::Vector3d, 2>> facesTangents;
    facesTangents.reserve(polyhedronFacesVertices.size());

    for (unsigned int f = 0; f < polyhedronFacesVertices.size(); f++)
    {
      const double normal_direction = polyhedronFacesNormalDirection.at(f) ?
                                        +1.0 : -1.0;
      facesTangents.push_back(PolygonTangents(polyhedronFacesVertices.at(f),
                                              normal_direction *
                                              polyhedronFacesNormal.at(f)));
    }

    return facesTangents;
  }
  // ***************************************************************************
  vector<Vector3d> GeometryUtilities::PolyhedronFaceBarycenter(const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<Vector3d> faceBarycenters;
    faceBarycenters.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
      faceBarycenters.push_back(PolygonBarycenter(polyhedronFaceVertices[f]));

    return faceBarycenters;
  }
  // ***************************************************************************
  Matrix3d GeometryUtilities::PolyhedronInertia(const Eigen::Vector3d& polyhedronCentroid,
                                                const std::vector<Eigen::MatrixXd>& polyhedronTetrahedraPoints) const
  {
    // Create reference quadrature points
    Eigen::MatrixXd referencePoints = Eigen::MatrixXd::Zero(3, 4);
    referencePoints(0,0) = 3.6531451881463450e-01;
    referencePoints(0,1) = 4.5746158708559548e-01;
    referencePoints(0,2) = 3.7551502872933986e-04;
    referencePoints(0,3) = 1.2366680032845828e-01;
    referencePoints(1,0) = 1.8002969351036546e-01;
    referencePoints(1,1) = 1.5593312049918590e-01;
    referencePoints(1,2) = 2.1607642918484790e-01;
    referencePoints(1,3) = 8.2157254096761967e-01;
    referencePoints(2,0) = 6.9232355736274656e-03;
    referencePoints(2,1) = 3.8176535606934675e-01;
    referencePoints(2,2) = 4.3070170707783589e-01;
    referencePoints(2,3) = 3.9933048641498409e-02;
    Eigen::VectorXd referenceWeights(4);
    referenceWeights[0] = 5.0086823222829355e-02;
    referenceWeights[1] = 4.6462929447761252e-02;
    referenceWeights[2] = 5.3182322583579078e-02;
    referenceWeights[3] = 1.6934591412496775e-02;

    // Create polygon quadrature points
    const unsigned int numPolyhedronTetra = polyhedronTetrahedraPoints.size();

    const unsigned int numTetraQuadraturePoints = referencePoints.cols();
    const unsigned int numPolyhedronQuadraturePoints = numPolyhedronTetra *
                                                       numTetraQuadraturePoints;

    Eigen::MatrixXd polyhedronQuadraturePoints = Eigen::MatrixXd::Zero(3, numPolyhedronQuadraturePoints);
    Eigen::VectorXd polyhedronQuadratureWeights = Eigen::VectorXd::Zero(numPolyhedronQuadraturePoints);

    Gedim::MapTetrahedron mapTetra(*this);

    for (unsigned int t = 0; t < numPolyhedronTetra; t++)
    {
      const Eigen::MatrixXd& tetraVertices = polyhedronTetrahedraPoints[t];

      Gedim::MapTetrahedron::MapTetrahedronData mapData = mapTetra.Compute(tetraVertices);
      polyhedronQuadraturePoints.block(0,
                                       numTetraQuadraturePoints * t,
                                       3,
                                       numTetraQuadraturePoints) = mapTetra.F(mapData,
                                                                              referencePoints);
      polyhedronQuadratureWeights.segment(numTetraQuadraturePoints * t,
                                          numTetraQuadraturePoints) = referenceWeights *
                                                                      abs(mapData.DetQ);

    }

    // Compute Inertia Matrix
    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();

    inertia(0, 0) = ((
                       (polyhedronQuadraturePoints.row(1).array() - polyhedronCentroid.y()).square() +
                       (polyhedronQuadraturePoints.row(2).array() - polyhedronCentroid.z()).square()
                       ) *
                     polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(1, 1) = ((
                       (polyhedronQuadraturePoints.row(0).array() - polyhedronCentroid.x()).square() +
                       (polyhedronQuadraturePoints.row(2).array() - polyhedronCentroid.z()).square()
                       ) *
                     polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(2, 2) = ((
                       (polyhedronQuadraturePoints.row(1).array() - polyhedronCentroid.y()).square() +
                       (polyhedronQuadraturePoints.row(0).array() - polyhedronCentroid.x()).square()
                       ) *
                     polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(0, 1) = - ((polyhedronQuadraturePoints.row(0).array() - polyhedronCentroid.x()) *
                       (polyhedronQuadraturePoints.row(1).array() - polyhedronCentroid.y()) *
                       polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(0, 2) = - ((polyhedronQuadraturePoints.row(0).array() - polyhedronCentroid.x()) *
                       (polyhedronQuadraturePoints.row(2).array() - polyhedronCentroid.z()) *
                       polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(1, 2) = - ((polyhedronQuadraturePoints.row(1).array() - polyhedronCentroid.y()) *
                       (polyhedronQuadraturePoints.row(2).array() - polyhedronCentroid.z()) *
                       polyhedronQuadratureWeights.transpose().array()).sum();
    inertia(1, 0) = inertia(0, 1);
    inertia(2, 0) = inertia(0, 2);
    inertia(2, 1) = inertia(1, 2);

    return inertia;
  }
  // ***************************************************************************
  bool GeometryUtilities::PolyhedronIsConvex(const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                             const std::vector<MatrixXd>& polyhedronFaceRotatedVertices,
                                             const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                             const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                             const std::vector<bool>& polyhedronFaceNormalDirections,
                                             const Eigen::Vector3d& polyhedronInternalPoint) const
  {
    /// check convexity with the intersections between the faces and
    /// a line starting from the polyhedron internal point
    /// to each face internal point.
    const unsigned int numPolyhedronFaces = polyhedronFaceVertices.size();

    for (unsigned int f1 = 0; f1 < numPolyhedronFaces; f1++)
    {
      // Check face convexity
      const vector<unsigned int> faceConvexHull = ConvexHull(polyhedronFaceRotatedVertices[f1]);
      if (!PolygonIsConvex(polyhedronFaceRotatedVertices[f1],
                           ExtractPoints(polyhedronFaceRotatedVertices[f1],
                                         faceConvexHull)))
      {
        return false;
      }

      for (unsigned int f2 = 0; f2 < numPolyhedronFaces; f2++)
      {
        if (f2 == f1)
          continue;

        // Check line intersection with other faces
        const Eigen::Vector3d faceNormal = polyhedronFaceNormalDirections[f2] ? polyhedronFaceNormals[f2] :
                                                                                -1.0 * polyhedronFaceNormals[f2];
        const IntersectionSegmentPlaneResult intersection = IntersectionSegmentPlane(polyhedronInternalPoint,
                                                                                     polyhedronFaceInternalPoints[f1],
                                                                                     faceNormal,
                                                                                     polyhedronFaceVertices[f2].col(0));

        switch (intersection.Type)
        {
          case IntersectionSegmentPlaneResult::Types::NoIntersection:
            continue;
          case IntersectionSegmentPlaneResult::Types::SingleIntersection:
          {
            if (intersection.SingleIntersection.Type == PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
              continue;

            switch (intersection.SingleIntersection.Type)
            {
              case PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
              case PointSegmentPositionTypes::OnSegmentLineAfterEnd:
              case PointSegmentPositionTypes::OnSegmentOrigin:
              case PointSegmentPositionTypes::OnSegmentEnd:
                continue;
              case PointSegmentPositionTypes::InsideSegment:
                return false;
              default:
                throw runtime_error("intersection.SingleIntersection.Type not expected");
            }
          }
            break;
          default:
            throw runtime_error("Intersection not expected");
        }
      }
    }

    return true;
  }
  // ***************************************************************************
  vector<bool> GeometryUtilities::PolyhedronFaceNormalDirections(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                 const Eigen::Vector3d& pointInsidePolyhedron,
                                                                 const vector<Vector3d>& polyhedronFaceNormals) const
  {
    vector<bool> faceDirections;
    faceDirections.reserve(polyhedronFaceVertices.size());

    for (unsigned int f = 0; f < polyhedronFaceVertices.size(); f++)
    {
      const Eigen::Vector3d& normal = polyhedronFaceNormals[f];
      const PointPlanePositionTypes pointFacePosition = PointPlanePosition(
                                                          PointPlaneDistance(pointInsidePolyhedron,
                                                                             normal,
                                                                             polyhedronFaceVertices[f].col(0)));
      Output::Assert(pointFacePosition == PointPlanePositionTypes::Negative ||
                     pointFacePosition == PointPlanePositionTypes::Positive);

      faceDirections.push_back(!(pointFacePosition == PointPlanePositionTypes::Positive));
    }

    return faceDirections;
  }
  // ***************************************************************************
  std::vector<bool> GeometryUtilities::PolyhedronFaceNormalDirections(const Eigen::MatrixXd& polyhedronVertices,
                                                                      const Eigen::MatrixXi& polyhedronEdges,
                                                                      const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                      const std::vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                      const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                      const std::vector<Eigen::MatrixXd>& polyhedronFaceRotatedVertices,
                                                                      const std::vector<Eigen::Vector3d>& polyhedronFaceNormals,
                                                                      const std::vector<Eigen::Vector3d>& polyhedronFaceTranslations,
                                                                      const std::vector<Eigen::Matrix3d>& polyhedronFaceRotationMatrices) const
  {
    std::vector<bool> vertices_intersection(polyhedronVertices.cols(), false);
    std::vector<bool> edges_intersection(polyhedronEdges.cols(), false);
    vector<bool> faceDirections(polyhedronFaceVertices.size(), true);

    for (unsigned int f1 = 0; f1 < polyhedronFaceVertices.size(); f1++)
    {
      const Eigen::Vector3d& normal = polyhedronFaceNormals[f1];
      const Eigen::Vector3d segmentOrigin = polyhedronFaceInternalPoints[f1];
      const Eigen::Vector3d segmentEnd = segmentOrigin + normal;
      list<double> curvilinearCoordinates;

      for (unsigned int f2 = 0; f2 < polyhedronFaceVertices.size(); f2++)
      {
        if (f2 == f1)
          continue;

        const IntersectionSegmentPlaneResult intersections = IntersectionSegmentPlane(segmentOrigin,
                                                                                      segmentEnd,
                                                                                      polyhedronFaceNormals[f2],
                                                                                      polyhedronFaceVertices[f2].col(0));

        switch (intersections.Type)
        {
          case IntersectionSegmentPlaneResult::Types::NoIntersection:
            continue;
          case IntersectionSegmentPlaneResult::Types::MultipleIntersections:
            continue;
          case IntersectionSegmentPlaneResult::Types::SingleIntersection:
          {
            if (intersections.SingleIntersection.Type == PointSegmentPositionTypes::OnSegmentLineBeforeOrigin)
              continue;

            switch (intersections.SingleIntersection.Type)
            {
              case PointSegmentPositionTypes::OnSegmentLineBeforeOrigin:
              case PointSegmentPositionTypes::OnSegmentOrigin:
                continue;
              case PointSegmentPositionTypes::OnSegmentLineAfterEnd:
              case PointSegmentPositionTypes::OnSegmentEnd:
              case PointSegmentPositionTypes::InsideSegment:
              {
                bool intersectionAlreadyConsidered = false;
                for (const double& curvilinearCoordinate : curvilinearCoordinates)
                {
                  if (AreValuesEqual(curvilinearCoordinate,
                                     intersections.SingleIntersection.CurvilinearCoordinate,
                                     Tolerance1D()))
                  {
                    intersectionAlreadyConsidered = true;
                    break;
                  }
                }

                if (intersectionAlreadyConsidered)
                  continue;

                const Eigen::Vector3d pointIntersection = segmentOrigin +
                                                          intersections.SingleIntersection.CurvilinearCoordinate * normal;

                const Eigen::Vector3d pointIntersection2D = RotatePointsFrom3DTo2D(pointIntersection,
                                                                                   polyhedronFaceRotationMatrices[f2].transpose(),
                                                                                   polyhedronFaceTranslations[f2]);
                const PointPolygonPositionResult& intersectionPosition = PointPolygonPosition_RayCasting(pointIntersection2D,
                                                                                                         polyhedronFaceRotatedVertices[f2]);

                switch (intersectionPosition.Type)
                {
                  case PointPolygonPositionResult::Types::Outside:
                    continue;
                  case PointPolygonPositionResult::Types::BorderVertex:
                  {
                    const unsigned int faceVertexIndex = polyhedronFaces[f2](0, intersectionPosition.BorderIndex);
                    if (!vertices_intersection[faceVertexIndex])
                    {
                      curvilinearCoordinates.push_back(intersections.SingleIntersection.CurvilinearCoordinate);
                      vertices_intersection[faceVertexIndex] = true;
                    }
                    continue;
                  }
                  case PointPolygonPositionResult::Types::BorderEdge:
                  {
                    const unsigned int faceEdgeIndex = polyhedronFaces[f2](1, intersectionPosition.BorderIndex);
                    if (!edges_intersection[faceEdgeIndex])
                    {
                      curvilinearCoordinates.push_back(intersections.SingleIntersection.CurvilinearCoordinate);
                      edges_intersection[faceEdgeIndex] = true;
                    }
                    continue;
                  }
                  case PointPolygonPositionResult::Types::Inside:
                    curvilinearCoordinates.push_back(intersections.SingleIntersection.CurvilinearCoordinate);
                    continue;
                  default:
                    throw runtime_error("intersectionPosition.Type not expected");
                }
              }
                break;
              default:
                throw runtime_error("intersection.SingleIntersection.Type not expected");
            }
          }
            break;
          default:
            throw runtime_error("Intersection not expected");
        }
      }

      faceDirections[f1] = (curvilinearCoordinates.size() % 2 == 0);
    }

    return faceDirections;
  }
  // ***************************************************************************
  vector<vector<unsigned int> > GeometryUtilities::PolyhedronFaceTriangulations(const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                const vector<vector<unsigned int>>& localFaceTriangulations) const
  {
    vector<vector<unsigned int>> facesTriangulation(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      vector<unsigned int>& faceTriangulation = facesTriangulation[f];

      vector<unsigned int> faceLocalTriangulation = localFaceTriangulations[f];
      faceTriangulation.resize(faceLocalTriangulation.size());

      for (unsigned int v = 0; v < faceLocalTriangulation.size(); v++)
        faceTriangulation[v] = polyhedronFaces[f](0, faceLocalTriangulation[v]);
    }

    return facesTriangulation;
  }
  // ***************************************************************************
  vector<vector<unsigned int>> GeometryUtilities::PolyhedronFaceTriangulationsByFirstVertex(const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                            const vector<Eigen::MatrixXd>& polyhedronFaceVertices) const
  {
    vector<vector<unsigned int>> polyhedronFacesTriangulations(polyhedronFaces.size());
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
      polyhedronFacesTriangulations[f] = PolygonTriangulationByFirstVertex(polyhedronFaceVertices[f]);

    return polyhedronFacesTriangulations;
  }
  // ***************************************************************************
  std::vector<std::vector<unsigned int> > GeometryUtilities::PolyhedronFaceTriangulationsByEarClipping(const unsigned int numPolyhedronFaces,
                                                                                                       const std::vector<Eigen::MatrixXd>& polyhedronFaces2DVertices) const
  {
    vector<vector<unsigned int>> polyhedronFacesTriangulations(numPolyhedronFaces);
    for (unsigned int f = 0; f < numPolyhedronFaces; f++)
      polyhedronFacesTriangulations[f] = PolygonTriangulationByEarClipping(polyhedronFaces2DVertices[f]);

    return polyhedronFacesTriangulations;
  }
  // ***************************************************************************
  std::vector<std::vector<Matrix3d>> GeometryUtilities::PolyhedronFaceExtractTriangulationPoints(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                                 const std::vector<std::vector<unsigned int> >& polyhedronFaceTriangulations) const
  {
    const unsigned int numFaces = polyhedronFaceTriangulations.size();
    vector<std::vector<Matrix3d>> faceTriangulationsPoints(numFaces);

    for (unsigned int f = 0; f < numFaces; f++)
      faceTriangulationsPoints[f] = ExtractTriangulationPoints(polyhedronFaceVertices[f],
                                                               polyhedronFaceTriangulations[f]);

    return faceTriangulationsPoints;
  }
  // ***************************************************************************
  std::vector<std::vector<Matrix3d> > GeometryUtilities::PolyhedronFaceTriangulationPointsByInternalPoint(const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                                          const std::vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                                                          const std::vector<std::vector<unsigned int> >& polyhedronFaceTriangulations) const
  {
    const unsigned int numFaces = polyhedronFaceTriangulations.size();
    vector<std::vector<Matrix3d>> faceTriangulationsPoints(numFaces);

    for (unsigned int f = 0; f < numFaces; f++)
      faceTriangulationsPoints[f] = ExtractTriangulationPointsByInternalPoint(polyhedronFaceVertices[f],
                                                                              polyhedronFaceInternalPoints[f],
                                                                              polyhedronFaceTriangulations[f]);

    return faceTriangulationsPoints;
  }
  // ***************************************************************************
  vector<vector<unsigned int>> GeometryUtilities::PolyhedronFaceTriangulationsByInternalPoint(const Eigen::MatrixXd& ,
                                                                                              const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                              const vector<Eigen::MatrixXd>& polyhedronFaceVertices,
                                                                                              const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints) const
  {
    vector<vector<unsigned int>> polyhedronFacesTriangulation(polyhedronFaces.size());

    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      polyhedronFacesTriangulation[f] = PolygonTriangulationByInternalPoint(polyhedronFaceVertices[f],
                                                                            polyhedronFaceInternalPoints[f]);
    }

    return polyhedronFacesTriangulation;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                                     const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                     const vector<vector<unsigned int>>& polyhedronFaceTriangulations,
                                                                                     const Eigen::Vector3d& ) const
  {
    list<unsigned int> tetrahedronList;

    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    for (unsigned int f = 0; f < polyhedronFaceTriangulations.size(); f++)
    {
      const unsigned int numFaceTriangulation = polyhedronFaceTriangulations[f].size() / 3;
      for (unsigned int ft = 0; ft < numFaceTriangulation; ft++)
      {
        tetrahedronList.push_back(polyhedronFaces[f](0, polyhedronFaceTriangulations[f][3 * ft]));
            tetrahedronList.push_back(polyhedronFaces[f](0, polyhedronFaceTriangulations[f][3 * ft + 1]));
            tetrahedronList.push_back(polyhedronFaces[f](0, polyhedronFaceTriangulations[f][3 * ft + 2]));
            tetrahedronList.push_back(numPolyhedronVertices);
      }
    }

    Output::Assert(tetrahedronList.size() % 4 == 0);

    return vector<unsigned int>(tetrahedronList.begin(), tetrahedronList.end());
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronTetrahedronsByFaceTriangulations(const Eigen::MatrixXd& polyhedronVertices,
                                                                                     const vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                                     const vector<vector<unsigned int> >& polyhedronFaceTriangulations,
                                                                                     const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                                                     const Eigen::Vector3d& ) const
  {
    list<unsigned int> tetrahedronList;

    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    for (unsigned int f = 0; f < polyhedronFaceTriangulations.size(); f++)
    {
      const unsigned int numFaceTriangulation = polyhedronFaceTriangulations[f].size() / 3;
      for (unsigned int ft = 0; ft < numFaceTriangulation; ft++)
      {
        tetrahedronList.push_back(numPolyhedronVertices + f);
        tetrahedronList.push_back(polyhedronFaces[f](0, polyhedronFaceTriangulations[f][3 * ft + 1]));
            tetrahedronList.push_back(polyhedronFaces[f](0, polyhedronFaceTriangulations[f][3 * ft + 2]));
            tetrahedronList.push_back(numPolyhedronVertices + polyhedronFaceInternalPoints.size());
      }
    }

    Output::Assert(tetrahedronList.size() % 4 == 0);

    return vector<unsigned int>(tetrahedronList.begin(), tetrahedronList.end());
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                               const Eigen::Vector3d& polyhedronInternalPoint,
                                                               const vector<unsigned int>& pointTetrahedrons) const
  {
    const unsigned int numTetra = pointTetrahedrons.size() / 4;
    vector<MatrixXd> tetra(numTetra, Eigen::MatrixXd::Zero(3, 4));

    for (unsigned int t = 0; t < numTetra; t++)
    {
      Eigen::MatrixXd& tetraVertices = tetra[t];
      tetraVertices.col(0)<< polyhedronVertices.col(pointTetrahedrons[4 * t]);
      tetraVertices.col(1)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 1]);
      tetraVertices.col(2)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 2]);
      tetraVertices.col(3)<< polyhedronInternalPoint;
    }

    return tetra;
  }
  // ***************************************************************************
  vector<MatrixXd> GeometryUtilities::ExtractTetrahedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                               const Eigen::Vector3d& polyhedronInternalPoint,
                                                               const vector<Eigen::Vector3d>& polyhedronFaceInternalPoints,
                                                               const vector<unsigned int>& pointTetrahedrons) const
  {
    const unsigned int numTetra = pointTetrahedrons.size() / 4;
    vector<MatrixXd> tetra(numTetra, Eigen::MatrixXd::Zero(3, 4));

    const unsigned int numPolyhedronVertices = polyhedronVertices.cols();
    for (unsigned int t = 0; t < numTetra; t++)
    {
      Eigen::MatrixXd& tetraVertices = tetra[t];

      tetraVertices.col(0)<< polyhedronFaceInternalPoints[pointTetrahedrons[4 * t] - numPolyhedronVertices];
      tetraVertices.col(1)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 1]);
      tetraVertices.col(2)<< polyhedronVertices.col(pointTetrahedrons[4 * t + 2]);
      tetraVertices.col(3)<< polyhedronInternalPoint;
    }

    return tetra;
  }
  // ***************************************************************************
  vector<unsigned int> GeometryUtilities::PolyhedronCoordinateSystem(const Eigen::MatrixXd& ,
                                                                     const Eigen::MatrixXi& polyhedronEdges)
  {
    vector<unsigned int> coordinateSystem(4);

    coordinateSystem[0] = 0;
    unsigned int vertexIndex = 1;
    for(unsigned int e = 0; e < polyhedronEdges.cols(); e++)
    {
      const unsigned int& originId = polyhedronEdges(0, e);
      const unsigned int& endId = polyhedronEdges(1, e);

      bool first = (originId == 0);
      bool second = (endId == 0);
      if (first || second)
        coordinateSystem[vertexIndex++] = first ? endId : originId;

      if (vertexIndex == 4)
        break;
    }

    Output::Assert(vertexIndex == 4);

    return coordinateSystem;
  }
  // ***************************************************************************
  void GeometryUtilities::ExportPolyhedronToVTU(const Eigen::MatrixXd& polyhedronVertices,
                                                const Eigen::MatrixXi& polyhedronEdges,
                                                const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                const std::string& exportFolder) const
  {
    vector<vector<unsigned int>> facesVertices(polyhedronFaces.size());
    for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
    {
      facesVertices[f].resize(polyhedronFaces[f].cols());
      for (unsigned int v = 0; v < polyhedronFaces[f].cols(); v++)
        facesVertices[f][v] = polyhedronFaces[f](0, v);
    }

    {
      // export vertices
      VTKUtilities exporter;
      vector<double> verticesIndex(polyhedronVertices.cols());
      for (unsigned int v = 0; v < polyhedronVertices.cols(); v++)
        verticesIndex[v] = v;

      exporter.AddPoints(polyhedronVertices,
                         {
                           {
                             "Index",
                             Gedim::VTPProperty::Formats::Cells,
                             static_cast<unsigned int>(verticesIndex.size()),
                             verticesIndex.data()
                           }
                         });
      exporter.Export(exportFolder +
                      "/Polyhedron_Vertices.vtu");
    }

    {
      // export edges
      VTKUtilities exporter;
      vector<double> edgesIndex(polyhedronEdges.cols());
      for (unsigned int e = 0; e < polyhedronEdges.cols(); e++)
        edgesIndex[e] = e;

      exporter.AddSegments(polyhedronVertices,
                           polyhedronEdges,
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(edgesIndex.size()),
                               edgesIndex.data()
                             }
                           });
      exporter.Export(exportFolder +
                      "/Polyhedron_Edges.vtu");
    }

    {
      // export faces
      VTKUtilities exporter;
      vector<double> facesIndex(polyhedronFaces.size());
      for (unsigned int f = 0; f < polyhedronFaces.size(); f++)
        facesIndex[f] = f;

      exporter.AddPolygons(polyhedronVertices,
                           facesVertices,
                           {
                             {
                               "Index",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(facesIndex.size()),
                               facesIndex.data()
                             }
                           });
      exporter.Export(exportFolder +
                      "/Polyhedron_Faces.vtu");
    }

    {
      // export polyhedron
      VTKUtilities exporter;
      exporter.AddPolyhedron(polyhedronVertices,
                             facesVertices);
      exporter.Export(exportFolder +
                      "/Polyhedron.vtu");
    }
  }
  // ***************************************************************************
  void GeometryUtilities::ExportPolyhedronToVTU(const unsigned int& index,
                                                const Eigen::MatrixXd& polyhedronVertices,
                                                const Eigen::MatrixXi& polyhedronEdges,
                                                const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                const std::vector<Eigen::MatrixXd>& polyhedronTetra,
                                                const double& polyhedronVolume,
                                                const Eigen::Vector3d& polyhedronCentroid,
                                                const std::vector<Eigen::MatrixXd>& polyhedronFaces3DVertices,
                                                const std::vector<double>& polyhedronFacesArea,
                                                const std::vector<Eigen::Vector3d>& polyhedronFaces2DCentroid,
                                                const std::vector<Eigen::Vector3d>& polyhedronFacesTranslation,
                                                const std::vector<Eigen::Matrix3d>& polyhedronFacesRotationMatrix,
                                                const std::vector<std::vector<Eigen::Matrix3d>>& polyhedronFaces3DTriangles,
                                                const std::vector<Eigen::Vector3d>& polyhedronFaces3DInternalPoint,
                                                const std::vector<Eigen::Vector3d>& polyhedronFaces3DNormal,
                                                const std::vector<bool>& polyhedronFaces3DNormalDirection,
                                                const std::string& exportFolder) const
  {
    {
      Gedim::VTKUtilities exporter;

      vector<double> id(1, index);
      vector<double> volume(1, polyhedronVolume);

      // Export Polyhedron
      exporter.AddPolyhedron(polyhedronVertices,
                             polyhedronEdges,
                             polyhedronFaces,
                             {
                               {
                                 "Id",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(id.size()),
                                 id.data()
                               },
                               {
                                 "Volume",
                                 Gedim::VTPProperty::Formats::Cells,
                                 static_cast<unsigned int>(volume.size()),
                                 volume.data()
                               }
                             });

      exporter.Export(exportFolder + "/" +
                      "Polyhedron.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      // Export Polyhedron tetra
      for (unsigned int t = 0; t < polyhedronTetra.size(); t++)
      {
        vector<double> id(1, t);

        const GeometryUtilities::Polyhedron polyhedron = CreateTetrahedronWithVertices(polyhedronTetra[t].col(0),
                                                                                       polyhedronTetra[t].col(1),
                                                                                       polyhedronTetra[t].col(2),
                                                                                       polyhedronTetra[t].col(3));

        exporter.AddPolyhedron(polyhedron.Vertices,
                               polyhedron.Edges,
                               polyhedron.Faces,
                               {
                                 {
                                   "Id",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(id.size()),
                                   id.data()
                                 }
                               });
      }

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_Tetra.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      for (unsigned int f = 0; f < polyhedronFaces3DVertices.size(); f++)
      {
        vector<double> id(1, f);
        vector<double> area(1, polyhedronFacesArea[f]);

        // Export faces
        exporter.AddPolygon(polyhedronFaces3DVertices[f],
                            {
                              {
                                "Id",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(id.size()),
                                id.data()
                              },
                              {
                                "Area",
                                Gedim::VTPProperty::Formats::Cells,
                                static_cast<unsigned int>(area.size()),
                                area.data()
                              }
                            });
      }

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_Faces.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      // Export faces triangles
      unsigned int numFaceTriangles = 0;
      for (unsigned int f = 0; f < polyhedronFaces3DTriangles.size(); f++)
      {
        vector<double> face(1, f);

        for (unsigned int t = 0; t < polyhedronFaces3DTriangles[f].size(); t++)
        {
          vector<double> id(1, numFaceTriangles++);
          exporter.AddPolygon(polyhedronFaces3DTriangles[f][t],
                              {
                                {
                                  "Face",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(face.size()),
                                  face.data()
                                },
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                }
                              });
        }
      }

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_FacesTriangles.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      // Export faces normal
      for (unsigned int f = 0; f < polyhedronFaces3DNormal.size(); f++)
      {
        vector<double> face(1, f);
        vector<double> normalDirection(1, polyhedronFaces3DNormalDirection[f] ? 1.0 : -1.0);

        exporter.AddSegment(polyhedronFaces3DInternalPoint[f],
                            polyhedronFaces3DInternalPoint[f] +
                            normalDirection[0] * polyhedronFaces3DNormal[f],
        {
          {
            "Face",
            Gedim::VTPProperty::Formats::Cells,
                static_cast<unsigned int>(face.size()),
                face.data()
          },
          {
            "NormalDirection",
            Gedim::VTPProperty::Formats::Cells,
                static_cast<unsigned int>(normalDirection.size()),
                normalDirection.data()
          }
        });
      }

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_FacesNormal.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      // Export Polyhedron centroid
      exporter.AddPoint(polyhedronCentroid);

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_Centroid.vtu");
    }

    {
      Gedim::VTKUtilities exporter;

      // Export faces centroids
      for (unsigned int f = 0; f < polyhedronFaces2DCentroid.size(); f++)
      {
        vector<double> face(1, f);

        const Eigen::Vector3d rotatedCentroid = RotatePointsFrom2DTo3D(polyhedronFaces2DCentroid[f],
                                                                       polyhedronFacesRotationMatrix[f],
                                                                       polyhedronFacesTranslation[f]);
        exporter.AddPoint(rotatedCentroid,
                          {
                            {
                              "Face",
                              Gedim::VTPProperty::Formats::Cells,
                              static_cast<unsigned int>(face.size()),
                              face.data()
                            }
                          });
      }

      exporter.Export(exportFolder + "/" +
                      "Polyhedron_FacesCentroid.vtu");
    }
  }
  // ***************************************************************************
  std::vector<unsigned int> GeometryUtilities::UnalignedPolyhedronPoints(const Eigen::MatrixXd& polyhedronVertices,
                                                                         const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                                         const std::vector<Eigen::Vector3d>& polyhedronFacesTranslation,
                                                                         const std::vector<Eigen::Matrix3d>& polyhedronFacesRotationMatrix,
                                                                         const std::vector<std::vector<unsigned int>>& polyhedronUnaligedFaces,
                                                                         const std::vector<std::vector<unsigned int>>& polyhedronFacesUnalignedVertices) const
  {
    std::set<unsigned int> unalignedVertices;

    for (const std::vector<unsigned int>& faces : polyhedronUnaligedFaces)
    {
      const unsigned int numAlignedFaces = faces.size();
      Gedim::Output::Assert(numAlignedFaces > 0);

      if (numAlignedFaces > 1)
      {
        std::vector<unsigned int> numFacesUnalignedVertices(numAlignedFaces + 1, 0);
        for (unsigned int f = 0; f < numAlignedFaces; f++)
          numFacesUnalignedVertices[f + 1] += numFacesUnalignedVertices[f] + polyhedronFacesUnalignedVertices[faces[f]].size();
        Eigen::MatrixXd unifiedUnalignedFacesVertices(3, numFacesUnalignedVertices[numAlignedFaces]);

        for (unsigned int f = 0; f < numAlignedFaces; f++)
        {
          const unsigned int polyhedronFace = faces[f];
          for (unsigned int v = 0; v < polyhedronFacesUnalignedVertices[polyhedronFace].size(); v++)
          {
            const unsigned int vertexIndex = polyhedronFaces[polyhedronFace](0, polyhedronFacesUnalignedVertices[polyhedronFace][v]);
            unifiedUnalignedFacesVertices.col(numFacesUnalignedVertices[f] + v) =
                polyhedronVertices.col(vertexIndex);
          }
        }

        const Eigen::MatrixXd unifiedUnalignedFacesVertices2D = RotatePointsFrom3DTo2D(unifiedUnalignedFacesVertices,
                                                                                       polyhedronFacesRotationMatrix[faces[0]].transpose(),
            polyhedronFacesTranslation[faces[0]]);

        const std::vector<unsigned int> convexHull = ConvexHull(unifiedUnalignedFacesVertices2D,
                                                                false);

        for (const unsigned int convexHullVertex : convexHull)
        {
          unsigned int f = 0;
          while (f < numAlignedFaces &&
                 convexHullVertex >= numFacesUnalignedVertices[f + 1])
            f++;

          Gedim::Output::Assert(convexHullVertex >= numFacesUnalignedVertices[f]);

          const unsigned int polyhedronFace = faces[f];
          const unsigned int unifyVertexIndex = convexHullVertex - numFacesUnalignedVertices[f];
          const unsigned int faceUnalignedVertexIndex = polyhedronFacesUnalignedVertices[polyhedronFace][unifyVertexIndex];
          const unsigned int polyhedronVertex = polyhedronFaces[polyhedronFace](0, faceUnalignedVertexIndex);

          if (unalignedVertices.find(polyhedronVertex) ==
              unalignedVertices.end())
            unalignedVertices.insert(polyhedronVertex);
        }
      }
      else
      {
        const unsigned int f = faces[0];
        for (const unsigned int unalignedFaceLocalPoint : polyhedronFacesUnalignedVertices[f])
        {
          const unsigned int unalignedFacePoint = polyhedronFaces[f](0, unalignedFaceLocalPoint);

          if (unalignedVertices.find(unalignedFacePoint) ==
              unalignedVertices.end())
            unalignedVertices.insert(unalignedFacePoint);
        }
      }
    }

    return std::vector<unsigned int>(unalignedVertices.begin(),
                                     unalignedVertices.end());
  }
  // ***************************************************************************
  GeometryUtilities::AlignedPolyhedronEdgesResult GeometryUtilities::AlignedPolyhedronEdges(const Eigen::MatrixXd& polyhedronVertices,
                                                                                            const std::vector<std::vector<unsigned int>>& verticesAdjacency,
                                                                                            const std::vector<std::vector<unsigned int>>& edgesAdjacency,
                                                                                            const std::vector<std::unordered_map<unsigned int, unsigned int>>& adjacencyVerticesMap,
                                                                                            const Eigen::MatrixXd& polyhedronEdgeTangents,
                                                                                            const Eigen::VectorXd& polyhedronEdgeSquaredLenghts) const
  {
    std::vector<bool> markedVertices(polyhedronVertices.cols(), false);
    std::vector<std::vector<bool>> markedSubEdges(verticesAdjacency.size());
    for (unsigned int v = 0; v < verticesAdjacency.size(); v++)
      markedSubEdges[v].resize(verticesAdjacency[v].size(), false);

    struct AlignedEdge
    {
        std::list<unsigned int> Vertices;
        Eigen::Vector3d LineOrigin;
        Eigen::Vector3d LineTangent;
        double LineTangentSquaredLength;
    };

    std::list<AlignedEdge> alignedEdges;

    std::queue<unsigned int> queue;

    unsigned int visitedVertex = 0;

    queue.push(visitedVertex);

    while (!queue.empty())
    {
      visitedVertex = queue.front();
      queue.pop();

      if (markedVertices[visitedVertex])
        continue;

      markedVertices[visitedVertex] = true;

      for (unsigned int av = 0; av < verticesAdjacency[visitedVertex].size(); av++)
      {
        if (markedSubEdges[visitedVertex][av])
          continue;

        const unsigned int adjacentVertex = verticesAdjacency[visitedVertex][av];
        const unsigned int adjacentEdge = edgesAdjacency[visitedVertex][av];
        const unsigned int visitedVertexAdjacentIndex = adjacencyVerticesMap[adjacentVertex].find(visitedVertex)->second;

        markedSubEdges[visitedVertex][av] = true;
        markedSubEdges[adjacentVertex][visitedVertexAdjacentIndex] = true;

        queue.push(adjacentVertex);

        AlignedEdge* alignedEdgeFound = nullptr;

        Eigen::MatrixXd edgeCoordinates(3, 2);
        edgeCoordinates.col(0)<< polyhedronVertices.col(visitedVertex);
        edgeCoordinates.col(1)<< polyhedronVertices.col(adjacentVertex);
        for (AlignedEdge& alignedEdge : alignedEdges)
        {
          const std::vector<bool> pointsAreOnLine = PointsAreOnLine(edgeCoordinates,
                                                                    alignedEdge.LineOrigin,
                                                                    alignedEdge.LineTangent);

          if (!pointsAreOnLine[0] ||
              !pointsAreOnLine[1])
            continue;

          alignedEdgeFound = &alignedEdge;
          break;
        }

        if (alignedEdgeFound == nullptr)
        {
          alignedEdges.push_back(AlignedEdge());

          AlignedEdge& alignedEdge = alignedEdges.back();
          alignedEdge.LineOrigin = polyhedronVertices.col(visitedVertex);
          alignedEdge.LineTangent = polyhedronEdgeTangents.col(adjacentEdge);
          alignedEdge.LineTangentSquaredLength = polyhedronEdgeSquaredLenghts[adjacentEdge];
          alignedEdge.Vertices.push_back(visitedVertex);
          alignedEdge.Vertices.push_back(adjacentVertex);
        }
        else
        {
          AlignedEdge& alignedEdge = *alignedEdgeFound;
          if (find(alignedEdge.Vertices.begin(),
                   alignedEdge.Vertices.end(),
                   visitedVertex) ==
              alignedEdge.Vertices.end())
            alignedEdge.Vertices.push_back(visitedVertex);

          if (find(alignedEdge.Vertices.begin(),
                   alignedEdge.Vertices.end(),
                   adjacentVertex) ==
              alignedEdge.Vertices.end())
            alignedEdge.Vertices.push_back(adjacentVertex);
        }
      }
    }

    GeometryUtilities::AlignedPolyhedronEdgesResult result;

    result.AlignedEdgesVertices.resize(alignedEdges.size());
    result.AlignedEdgesEdges.resize(alignedEdges.size());

    unsigned int alignedEdgeIndex = 0;
    for (const AlignedEdge& alignedEdge : alignedEdges)
    {
      std::vector<unsigned int> alignedVertices(alignedEdge.Vertices.begin(),
                                                alignedEdge.Vertices.end());
      std::vector<double> curvilinearCoordinates(alignedEdge.Vertices.size());

      for (unsigned int v = 0; v < alignedVertices.size(); v++)
        curvilinearCoordinates[v] = PointLineCurvilinearCoordinate(polyhedronVertices.col(alignedVertices.at(v)),
                                                                   alignedEdge.LineOrigin,
                                                                   alignedEdge.LineTangent,
                                                                   alignedEdge.LineTangentSquaredLength);

      const std::vector<unsigned int> orderedVerticesIndex = Gedim::Utilities::SortArrayIndices(curvilinearCoordinates);

      result.AlignedEdgesVertices[alignedEdgeIndex].resize(alignedVertices.size());
      result.AlignedEdgesEdges[alignedEdgeIndex].resize(alignedVertices.size() - 1);

      for (unsigned int ov = 0; ov < orderedVerticesIndex.size() - 1; ov++)
      {
        result.AlignedEdgesVertices[alignedEdgeIndex][ov] = alignedVertices.at(orderedVerticesIndex.at(ov));

        const unsigned int edgeOrigin = alignedVertices.at(orderedVerticesIndex.at(ov));
        const unsigned int edgeEnd = alignedVertices.at(orderedVerticesIndex.at(ov + 1));

        const unsigned int edgePosition = adjacencyVerticesMap[edgeOrigin].find(edgeEnd)->second;
        const unsigned int edgeIndex = edgesAdjacency[edgeOrigin][edgePosition];
        result.AlignedEdgesEdges[alignedEdgeIndex][ov] = edgeIndex;
      }
      result.AlignedEdgesVertices[alignedEdgeIndex][orderedVerticesIndex.size() - 1] = alignedVertices.at(orderedVerticesIndex.at(orderedVerticesIndex.size() - 1));


      alignedEdgeIndex++;
    }

    return result;
  }
  // ***************************************************************************
}
