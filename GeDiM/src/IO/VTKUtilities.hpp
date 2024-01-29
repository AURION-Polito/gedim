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

  struct VTPPoints final
  {
      const Eigen::MatrixXd& Points;

      VTPPoints(const Eigen::MatrixXd& points) :
        Points(points)
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

  struct VTPSegments final
  {
      const Eigen::MatrixXd& Vertices;
      const Eigen::MatrixXi& Edges;

      VTPSegments(const Eigen::MatrixXd& vertices,
                  const Eigen::MatrixXi& edges) :
        Vertices(vertices),
        Edges(edges)
      { }
  };

  struct VTPPolygon final
  {
      const Eigen::MatrixXd& Vertices;
      const std::vector<unsigned int>& Polygon;

      VTPPolygon(const Eigen::MatrixXd& vertices,
                 const std::vector<unsigned int>& polygon) :
        Vertices(vertices),
        Polygon(polygon)
      { }
  };

  struct VTPPolygons final
  {
      const Eigen::MatrixXd& Vertices;
      const std::vector<std::vector<unsigned int>>& Polygons;

      VTPPolygons(const Eigen::MatrixXd& vertices,
                  const std::vector<std::vector<unsigned int>>& polygons) :
        Vertices(vertices),
        Polygons(polygons)
      { }
  };

  struct VTPPolyhedron final
  {
      const Eigen::MatrixXd& Vertices; /// size 3xnumVertices
      const std::vector<std::vector<unsigned int>>& PolyhedronFaces; /// size numFaces x numFaceVertices

      VTPPolyhedron(const Eigen::MatrixXd& vertices,
                    const std::vector<std::vector<unsigned int>>& faces) :
        Vertices(vertices),
        PolyhedronFaces(faces)
      { }
  };

  struct VTPPolyhedrons final
  {
      const Eigen::MatrixXd& Vertices; /// size 3xnumVertices
      const std::vector<std::vector<std::vector<unsigned int>>>& PolyhedronsFaces; /// size numPolyhedrons x numPolyhedronFaces x numPolyhedronFaceVertices

      VTPPolyhedrons(const Eigen::MatrixXd& vertices,
                     const std::vector<std::vector<std::vector<unsigned int>>>& polyhedronsFaces) :
        Vertices(vertices),
        PolyhedronsFaces(polyhedronsFaces)
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
      void AddPoints(const Eigen::MatrixXd& points,
                     const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::Vector3d& origin,
                      const Eigen::Vector3d& end,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegment(const Eigen::MatrixXd& vertices,
                      const Eigen::VectorXi& edge,
                      const std::vector<VTPProperty>& properties = {});
      void AddSegments(const Eigen::MatrixXd& vertices,
                       const Eigen::MatrixXi& edges,
                       const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygon(const Eigen::MatrixXd& vertices,
                      const std::vector<unsigned int>& polygon,
                      const std::vector<VTPProperty>& properties = {});
      void AddPolygons(const Eigen::MatrixXd& vertices,
                       const std::vector<Eigen::MatrixXi>& polygons,
                       const std::vector<VTPProperty>& properties = {});
      void AddPolygons(const Eigen::MatrixXd& vertices,
                       const std::vector<std::vector<unsigned int>>& polygons,
                       const std::vector<VTPProperty>& properties = {});
      void AddPolyhedron(const Eigen::MatrixXd& vertices,
                         const Eigen::MatrixXi& edges,
                         const std::vector<Eigen::MatrixXi>& faces,
                         const std::vector<VTPProperty>& properties = {});
      void AddPolyhedron(const Eigen::MatrixXd& vertices,
                         const std::vector<std::vector<unsigned int>>& polyhedronFaces,
                         const std::vector<VTPProperty>& properties = {});
      void AddPolyhedrons(const Eigen::MatrixXd& vertices,
                          const std::vector<std::vector<std::vector<unsigned int>>>& polyhedronsFaces,
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
      void AddPoints(const Eigen::MatrixXd& points,
                     vtkNew<vtkPoints>& vtkPoints) const;
      void AddVertex(const unsigned int& pointId,
                     vtkNew<vtkCellArray>& vertices) const;
      void AddVertices(const std::vector<unsigned int>& pointIds,
                       vtkNew<vtkCellArray>& vertices) const;
      void AddLine(const unsigned int& originId,
                   const unsigned int& endId,
                   vtkNew<vtkCellArray>& lines) const;
      void AddLines(const Eigen::MatrixXi& edges,
                    vtkNew<vtkCellArray>& lines) const;
      void AddPolygon(const std::vector<unsigned int>& faceVerticesIds,
                      vtkNew<vtkCellArray>& faces) const;
      void AddPolygons(const std::vector<std::vector<unsigned int>>& facesVerticesIds,
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
