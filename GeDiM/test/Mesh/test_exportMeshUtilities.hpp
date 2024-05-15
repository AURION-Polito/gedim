#ifndef __test_exportMeshUtilities_H
#define __test_exportMeshUtilities_H

#include "Eigen/Eigen"

namespace GedimUnitTesting
{
  struct ExportMeshData final
  {
      Eigen::MatrixXd Cell0Ds;
      Eigen::MatrixXi Cell1Ds;
      std::vector<Eigen::MatrixXi> Cell2Ds;
      std::vector<std::vector<unsigned int>> Cell3DsVertices;
      std::vector<std::vector<unsigned int>> Cell3DsEdges;
      std::vector<std::vector<unsigned int>> Cell3DsFaces;
      std::vector<unsigned int> Cell0DsMarker;
      std::vector<unsigned int> Cell1DsMarker;
      std::vector<unsigned int> Cell2DsMarker;
      std::vector<unsigned int> Cell3DsMarker;
  };

  struct ExportMeshGeometricData3D final
  {
      std::vector<Eigen::MatrixXd> PolyhedronsVertices;
      std::vector<Eigen::MatrixXi> PolyhedronsEdges;
      std::vector<std::vector<Eigen::MatrixXi>> PolyhedronsFaces;
      std::vector<double> PolyhedronsVolume;
      std::vector<double> PolyhedronsDiameter;
      std::vector<Eigen::Vector3d> PolyhedronsCentroid;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsTetrahedronsVertices;
      std::vector<std::vector<Eigen::Vector3d>> PolyhedronsFacesTranslation;
      std::vector<std::vector<Eigen::Matrix3d>> PolyhedronsFacesRotationMatrix;
      std::vector<Eigen::MatrixXd> PolyhedronsFacesNormal;
      std::vector<std::vector<bool>> PolyhedronsFacesNormalDirection;
      std::vector<std::vector<std::vector<bool>>> PolyhedronsFacesEdgesDirection;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFaces2DVertices;
      std::vector<std::vector<std::vector<Eigen::Matrix3d>>> PolyhedronsFacesTriangulations2DVertices;
      std::vector<std::vector<double>> PolyhedronsFacesArea;
      std::vector<std::vector<Eigen::Vector3d>> PolyhedronsFaces2DCentroid;
      std::vector<Eigen::VectorXd> PolyhedronsFacesDiameter;
      std::vector<std::vector<Eigen::VectorXd>> PolyhedronsFacesEdgesLength;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DCentroid;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DTangent;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DTangentNormalized;
      std::vector<std::vector<Eigen::MatrixXd>> PolyhedronsFacesEdges2DNormal;
  };

  struct ExportMeshGeometricData2D final
  {
      std::vector<Eigen::MatrixXd> PolygonsVertices;
      std::vector<std::vector<Eigen::Matrix3d>> PolygonsTriangulations;
      std::vector<double> PolygonsArea;
      std::vector<Eigen::Vector3d> PolygonsCentroid;
      std::vector<double> PolygonsDiameter;
      std::vector<std::vector<bool>> PolygonsEdgesDirection;
      std::vector<Eigen::VectorXd> PolygonsEdgesLength;
      std::vector<Eigen::MatrixXd> PolygonsEdgesTangent;
      std::vector<Eigen::MatrixXd> PolygonsEdgesTangentNormalized;
      std::vector<Eigen::MatrixXd> PolygonsEdgesCentroid;
      std::vector<Eigen::MatrixXd> PolygonsEdgesNormal;
  };

  class ExportMeshUtilities
  {
    public:
    public:
      ExportMeshUtilities(ExportMeshUtilities& other) = delete;
      void operator=(const ExportMeshUtilities&) = delete;

    public:
      ExportMeshUtilities() {}

      static ExportMeshData ImportMesh2DFromText(const std::string& file_path);
      static void ExportMesh2DToText(const ExportMeshData& mesh_data,
                                     const std::string& file_path);
      static ExportMeshGeometricData2D ImportMeshGeometricData2DFromText(const std::string& file_path);
      static void ExportMeshGeometricData2DToText(const ExportMeshGeometricData2D& mesh_geometric_data,
                                                  const std::string& file_path);

      static ExportMeshData ImportMesh3DFromText(const std::string& file_path);
      static void ExportMesh3DToText(const ExportMeshData& mesh_data,
                                     const std::string& file_path);
      static ExportMeshGeometricData3D ImportMeshGeometricData3DFromText(const std::string& file_path);
      static void ExportMeshGeometricData3DToText(const ExportMeshGeometricData3D& mesh_geometric_data,
                                                  const std::string& file_path);
  };

  template <class T>
  std::ostream& operator<<(std::ostream& out,
                           const std::vector<T>& elements)
  {
    out<< elements.size()<< ",";
    out<< "{";
    unsigned int i = 0;
    for (const auto& element : elements)
    {
      out<< (i != 0 ? "," : "")<< element;
      i++;
    }
    out<< "}";

    return out;
  }

  template <class T, int Rows, int Cols>
  std::ostream& operator<<(std::ostream& out,
                           const Eigen::Matrix<T, Rows, Cols>& matrix)
  {
    out<< matrix.rows()<< ","<< matrix.cols()<< ",";
    out<< "{";
    for (unsigned int c = 0; c < matrix.cols(); c++)
    {
      out<< (c != 0 ? ",{" : "{");
      for (unsigned int r = 0; r < matrix.rows(); r++)
        out<< (r != 0 ? ", " : "")<< matrix(r, c);
      out<< "}";
    }
    out<< "}";

    return out;
  }

  template <class T>
  std::istream& operator>>(std::istream& in,
                           std::vector<T>& elements)
  {
    char separator;
    unsigned int size = 0;
    in >> size;
    elements.resize(size);
    in >> separator;

    in >> separator;
    for (unsigned int v = 0; v < size; v++)
    {
      T element;
      in >> element;
      in >> separator;
      elements[v] = element;
    }

    return in;
  }

  template <class T, int Rows, int Cols>
  std::istream& operator>>(std::istream& in,
                           Eigen::Matrix<T, Rows, Cols>& matrix)
  {
    char separator;
    unsigned int rows = 0, cols = 0;
    in >> rows;
    in >> separator;
    in >> cols;
    matrix.resize(rows, cols);
    in >> separator;

    in >> separator;
    for (unsigned int c = 0; c < matrix.cols(); c++)
    {
      in >> separator;
      for (unsigned int r = 0; r < matrix.rows(); r++)
      {
        in >> matrix(r, c);
        in >> separator;
      }
      in >> separator;
    }

    return in;
  }
}

#endif // __IOUtilities_H
