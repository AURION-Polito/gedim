#ifndef __VTPExtension_H
#define __VTPExtension_H

#include "Macro.hpp"

#include <vector>
#include <map>

#include "Eigen/Eigen"

#if ENABLE_VTK == 1
#include "vtkCellType.h"
#include "vtkSmartPointer.h"
#include "vtkWeakPointer.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkCellTypes.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkAppendFilter.h"
#endif // ENABLE_VTK


namespace Gedim
{
  struct VTPPoint final
  {
      const Eigen::Vector3d& Point;

      VTPPoint(const Eigen::Vector3d& point) :
        Point(point)
      { }
  };

  struct VTPSegment final
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::VectorXi& Edge;

      VTPSegment(const Eigen::MatrixXd& vertices,
                 const Eigen::VectorXi& edge) :
        Vertices(vertices),
        Edge(edge)
      { }
  };

  struct VTPPolygon final
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::MatrixXi& Edges;
      const Eigen::MatrixXi& Face;

      VTPPolygon(const Eigen::MatrixXd& vertices,
                 const Eigen::MatrixXi& edges,
                 const Eigen::MatrixXi& face) :
        Vertices(vertices),
        Edges(edges),
        Face(face)
      { }
  };

  struct VTPPolyhedron final
  {
      const Eigen::MatrixXd& Vertices; /// size 3xnumVertices
      const Eigen::MatrixXi& Edges; /// size 2xnumEdges
      const std::vector<Eigen::MatrixXi>& Faces; /// size numFacesx2xnumFaceVertices

      VTPPolyhedron(const Eigen::MatrixXd& vertices,
                    const Eigen::MatrixXi& edges,
                    const std::vector<Eigen::MatrixXi>& faces) :
        Vertices(vertices),
        Edges(edges),
        Faces(faces)
      { }
  };

  struct VTPProperty final
  {
      enum Formats
      {
        Points = 0,
        Cells = 1
      };

      std::string Label;
      Formats Format;

      unsigned int Size;
      const double* Data;
  };

  class IGeometryToPolyData
  {
    public:
      virtual ~IGeometryToPolyData() { }

#if ENABLE_VTK == 1
      virtual void Convert(vtkNew<vtkAppendFilter>& exportData) const = 0;
#endif
  };

  class VTKUtilities final
  {
    public:
      enum ExportFormats
      {
        Binary = 0,
        Ascii = 1,
        Appended = 2
      };

    private:
      std::list<IGeometryToPolyData*> geometries;

#if ENABLE_VTK == 1
      vtkNew<vtkAppendFilter> exportData;
#endif

      void AddToExportData(IGeometryToPolyData& polyData);

    public:
      VTKUtilities();
      ~VTKUtilities();

      void AddPoint(const Eigen::Vector3d& point,
                    const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::Vector3d& origin,
                      const Eigen::Vector3d& end,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const Eigen::VectorXi& edge,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const Eigen::MatrixXi& edges,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolyhedron(const Eigen::MatrixXd& vertices,
                         const Eigen::MatrixXi& edges,
                         const std::vector<Eigen::MatrixXi>& faces,
                         const std::vector<VTPProperty>& properties = {});

      void Export(const std::string& filePath,
                  const ExportFormats& format = ExportFormats::Binary) const;
  };

  template <typename T>
  class GeometryToPolyData final : public IGeometryToPolyData
  {
    protected:
      const T& geometry;
      const std::vector<VTPProperty>& properties;

#if ENABLE_VTK == 1
      void AddPoint(const Eigen::Vector3d& point,
                    vtkNew<vtkPoints>& points) const;
      void AddVertex(const unsigned int& pointId,
                     vtkNew<vtkCellArray>& vertices) const;
      void AddLine(const unsigned int& originId,
                   const unsigned int& endId,
                   vtkNew<vtkCellArray>& lines) const;
      void AddPolygon(const Eigen::VectorXi& faceVerticesIds,
                      vtkNew<vtkCellArray>& faces) const;
      void AppendSolution(vtkDataSet* polyData) const;
#endif

    public:
      GeometryToPolyData(const T& geometry,
                         const std::vector<VTPProperty>& properties);
      ~GeometryToPolyData();

#if ENABLE_VTK == 1
      void Convert(vtkNew<vtkAppendFilter>& exportData) const;
#endif
  };
}

#endif // __VTPExtension_H
