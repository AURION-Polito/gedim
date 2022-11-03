#include "VTKUtilities.hpp"

#include "IOUtilities.hpp"
#include "CommonUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  VTKUtilities::VTKUtilities()
  {
  }
  VTKUtilities::~VTKUtilities()
  {
  }
  // ***************************************************************************
  void VTKUtilities::AddToExportData(IGeometryToPolyData& polyData,
                                     const vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    if (properties.size() > 0)
    {
      polyData.InitializeSolutions(properties.size());
      unsigned int p = 0;
      for (const VTPProperty& property : properties)
      {
        polyData.SetSolutionOptions(p,
                                    property.Label,
                                    property.Format);
        polyData.AddSolution(p,
                             property.Size,
                             property.Data);
        p++;
      }
    }

    polyData.Convert(exportData);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPoint(const Eigen::Vector3d& point,
                              const vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    const VTPPoint vtpPoint(point);
    GeometryToPolyData<VTPPoint> polyData(vtpPoint);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddSegment(const Eigen::Vector3d& origin,
                                const Eigen::Vector3d& end,
                                const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    const Eigen::MatrixXd vertices =  (Eigen::MatrixXd(3, 2)<< origin, end).finished();
    const Eigen::VectorXi edge =  (Eigen::VectorXi(2)<< 0, 1).finished();
    const VTPSegment vtpSegment(vertices,
                                edge);
    GeometryToPolyData<VTPSegment> polyData(vtpSegment);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddSegment(const Eigen::MatrixXd& vertices,
                                const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    const Eigen::VectorXi edge =  (Eigen::VectorXi(2)<< 0, 1).finished();
    const VTPSegment vtpSegment(vertices,
                                edge);
    GeometryToPolyData<VTPSegment> polyData(vtpSegment);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPolygon(const Eigen::MatrixXd& vertices,
                                const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    Eigen::MatrixXi edges(2, vertices.cols());
    for (unsigned int v = 0; v < vertices.cols(); v++)
      edges.col(v)<< v, (v + 1) % vertices.cols();

    Eigen::MatrixXi face(2, vertices.cols());
    face.row(0)<< Eigen::VectorXi::LinSpaced(vertices.cols(),
                                             0.0,
                                             vertices.cols() - 1).transpose();
    face.row(1)<< Eigen::VectorXi::LinSpaced(vertices.cols(),
                                             0.0,
                                             vertices.cols() - 1).transpose();

    const VTPPolygon vtpPolygon(vertices,
                                edges,
                                face);
    GeometryToPolyData<VTPPolygon> polyData(vtpPolygon);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::AddPolyhedron(const Eigen::MatrixXd& vertices,
                                   const Eigen::MatrixXi& edges,
                                   const std::vector<Eigen::MatrixXi>& faces,
                                   const std::vector<VTPProperty>& properties)
  {
#if ENABLE_VTK == 1
    const VTPPolyhedron vtpPolyhedron(vertices,
                                      edges,
                                      faces);
    GeometryToPolyData<VTPPolyhedron> polyData(vtpPolyhedron);
    AddToExportData(polyData,
                    properties);
#endif
  }
  // ***************************************************************************
  void VTKUtilities::Export(const std::string& filePath,
                            const ExportFormats& format) const
  {
#if ENABLE_VTK == 1
    exportData->Update();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filePath.c_str());
    writer->SetInputData(exportData->GetOutput());

    switch (format)
    {
      case VTKUtilities::Binary:
        writer->SetDataModeToBinary();
        break;
      case VTKUtilities::Ascii:
        writer->SetDataModeToAscii();
        break;
      case VTKUtilities::Appended:
        writer->SetDataModeToAppended();
        break;
      default:
        throw runtime_error("Export Format not supported");
    }

    writer->Write();
#endif
  }
  // ***************************************************************************
  template class GeometryToPolyData<VTPPoint>;
  template class GeometryToPolyData<VTPSegment>;
  template class GeometryToPolyData<VTPPolygon>;
  template class GeometryToPolyData<VTPPolyhedron>;
  // ***************************************************************************
  template <typename T>
  GeometryToPolyData<T>::GeometryToPolyData(const T& geometry) :
    geometry(geometry)
  {
    numberSolutions = 0;
  }
  template <typename T>
  GeometryToPolyData<T>::~GeometryToPolyData()
  {
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::InitializeSolutions(const unsigned int& _numberSolutions)
  {
    Output::Assert(_numberSolutions > 0);
    numberSolutions = _numberSolutions;

    solutionLabels.resize(numberSolutions, "Solution");
    solutionFormats.resize(numberSolutions, VTPProperty::Formats::Points);
    solutionSizes.resize(numberSolutions, 0);
    solutions.resize(numberSolutions, NULL);
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::SetSolutionOptions(const unsigned int& solutionPosition,
                                                 const string& solutionLabel,
                                                 const VTPProperty::Formats& format)
  {
    Output::Assert(solutionPosition < numberSolutions);

    solutionLabels[solutionPosition] = solutionLabel;
    solutionFormats[solutionPosition] = format;
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddSolution(const unsigned int& solutionPosition,
                                          const unsigned int& solutionSize,
                                          const double* solution)
  {
    Output::Assert(solutionPosition < solutionLabels.size());

    if (solutions[solutionPosition] != NULL)
      throw runtime_error("Solution " + to_string(solutionPosition) + " already inserted");

    solutions[solutionPosition] = solution;
    solutionSizes[solutionPosition] = solutionSize;
  }
  // ***************************************************************************
#if ENABLE_VTK == 1
  template <typename T>
  void GeometryToPolyData<T>::AddPoint(const Eigen::Vector3d& point,
                                       vtkNew<vtkPoints>& points) const
  {
    points->InsertNextPoint(point(0),
                            point(1),
                            point(2));
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddVertex(const unsigned int& pointId,
                                        vtkNew<vtkCellArray>& vertices) const
  {
    vertices->InsertNextCell(1);
    vertices->InsertCellPoint(pointId);
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddLine(const unsigned int& originId,
                                      const unsigned int& endId,
                                      vtkNew<vtkCellArray>& lines) const
  {
    lines->InsertNextCell(2);
    lines->InsertCellPoint(originId);
    lines->InsertCellPoint(endId);
  }
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddPolygon(const Eigen::VectorXi& faceVerticesIds,
                                         vtkNew<vtkCellArray>& faces) const
  {
    faces->InsertNextCell(faceVerticesIds.size());
    for (unsigned int fp = 0 ; fp < faceVerticesIds.size(); fp++)
      faces->InsertCellPoint(faceVerticesIds[fp]);
  }
#endif
  // ***************************************************************************
  template <typename T>
  void GeometryToPolyData<T>::AddSolution(const unsigned int& solutionPosition,
                                          vtkNew<vtkDoubleArray>& vtkSolution) const
  {
    if (solutionPosition >= numberSolutions)
      throw runtime_error("");

    const unsigned int& solutionSize = solutionSizes[solutionPosition];
    const double* solution = solutions[solutionPosition];

    if (solutionSize == 0 || solution == NULL)
      throw runtime_error("");

#if ENABLE_VTK == 1
    vtkSolution->SetNumberOfValues(solutionSize);
    for (unsigned int p = 0; p < solutionSize; p++)
      vtkSolution->SetValue(p, solution[p]);
#endif // ENABLE_VTK
  }
  // ***************************************************************************
#if ENABLE_VTK == 1
  template<typename T>
  void GeometryToPolyData<T>::AppendSolution(vtkDataSet* polyData) const
  {
    if (numberSolutions > 0)
    {
      for (unsigned int s = 0; s < numberSolutions; s++)
      {
        vtkNew<vtkDoubleArray> vtkSolution;

        const string& solutionLabel = solutionLabels[s];
        VTPProperty::Formats solutionFormat = solutionFormats[s];

        vtkSolution->SetName(solutionLabel.c_str());

        switch (solutionFormat)
        {
          case VTPProperty::Formats::Points:
            polyData->GetPointData()->AddArray(vtkSolution);

            if (s == 0)
              polyData->GetPointData()->SetActiveScalars(solutionLabel.c_str());
            break;
          case VTPProperty::Formats::Cells:
            polyData->GetCellData()->AddArray(vtkSolution);

            if (s == 0)
              polyData->GetCellData()->SetActiveScalars(solutionLabel.c_str());
            break;
          default:
            throw runtime_error("Solution Format not supported");
        }

        AddSolution(s, vtkSolution);
      }
    }
  }
#endif // ENABLE_VTK
  // ***************************************************************************
#if ENABLE_VTK == 1
  template <>
  void GeometryToPolyData<VTPPoint>::Convert(vtkNew<vtkAppendFilter>& exportData) const
  {
    vtkNew<vtkPolyData> polyData;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> vertices;

    const unsigned int numberTotalPoints = 1;
    points->Allocate(numberTotalPoints);
    vertices->Allocate(numberTotalPoints);

    polyData->SetPoints(points);
    polyData->SetVerts(vertices);

    AddPoint(geometry.Point,
             points);
    AddVertex(0,
              vertices);

    AppendSolution(polyData);

    exportData->AddInputData(polyData);
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPSegment>::Convert(vtkNew<vtkAppendFilter>& exportData) const
  {
    vtkNew<vtkPolyData> polyData;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> lines;

    const unsigned int numberTotalPoints = 2;
    points->Allocate(numberTotalPoints);

    polyData->SetPoints(points);
    polyData->SetLines(lines);

    AddPoint(geometry.Vertices.col(0),
             points);
    AddPoint(geometry.Vertices.col(1),
             points);

    AddLine(geometry.Edge(0),
            geometry.Edge(1),
            lines);

    AppendSolution(polyData);

    exportData->AddInputData(polyData);
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolygon>::Convert(vtkNew<vtkAppendFilter>& exportData) const
  {
    vtkNew<vtkPolyData> polyData;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> faces;

    const unsigned int numberTotalPoints = geometry.Vertices.cols();
    unsigned int numberTotalFaces = 1;
    points->Allocate(numberTotalPoints);
    faces->Allocate(numberTotalFaces);

    polyData->SetPoints(points);
    polyData->SetPolys(faces);

    for (unsigned int i = 0 ; i < geometry.Vertices.cols(); i++)
    {
      AddPoint(geometry.Vertices.col(i),
               points);
    }

    AddPolygon(geometry.Face.row(0).transpose(),
               faces);

    AppendSolution(polyData);

    exportData->AddInputData(polyData);
  }
  // ***************************************************************************
  template <>
  void GeometryToPolyData<VTPPolyhedron>::Convert(vtkNew<vtkAppendFilter>& exportData) const
  {
    vtkNew<vtkUnstructuredGrid> polyData;

    vtkNew<vtkPoints> points;
    vtkNew<vtkIdList> faces;

    const unsigned int numPoints = geometry.Vertices.cols();
    unsigned int numFaces = geometry.Faces.size();
    points->Allocate(numPoints);
    faces->Allocate(numFaces);

    vector<vtkIdType> pointIds(numPoints);
    for (unsigned int i = 0 ; i < numPoints; i++)
    {
      AddPoint(geometry.Vertices.col(i),
               points);

      pointIds[i] = static_cast<vtkIdType>(i);
    }

    for (unsigned int f = 0; f < numFaces; f++)
    {
      const unsigned int numFaceVertices = geometry.Faces[f].cols();
      faces->InsertNextId(numFaceVertices);
      for (unsigned int v = 0; v < numFaceVertices; v++)
        faces->InsertNextId(geometry.Faces[f](0, v));
    }

    polyData->SetPoints(points);
    polyData->InsertNextCell(VTK_POLYHEDRON, numPoints, pointIds.data(), numFaces, faces->GetPointer(0));

    AppendSolution(polyData);

    exportData->AddInputData(polyData);
  }
#endif // ENABLE_VTK
  // ***************************************************************************
}
