#include "ConformMeshUtilities.hpp"

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  ConformMeshUtilities::ConformMeshUtilities(const GeometryUtilities& geometryUtilities,
                                             const MeshUtilities& meshUtilities) :
    geometryUtilities(geometryUtilities),
    meshUtilities(meshUtilities)
  {
  }
  ConformMeshUtilities::~ConformMeshUtilities()
  {
  }
  // ***************************************************************************
  void ConformMeshUtilities::ComputeDomainConformedMesh(const std::vector<std::list<double>>& segmentsAdditionalPoints,
                                                        const vector<Eigen::MatrixXd>& segmentsVertices,
                                                        const vector<Eigen::Vector3d>& segmentsTangent,
                                                        const vector<Eigen::Vector3d>& segmentsBarycenter,
                                                        const vector<double>& segmentsLength,
                                                        const vector<double>& segmentsSquaredLength,
                                                        IMeshDAO& domainMesh,
                                                        vector<IntersectorMesh2DSegment::IntersectionMesh>& segmentsIntersectionMesh,
                                                        std::vector<std::vector<double> >& segmentsCurvilinearCoordinatesMesh,
                                                        vector<UnionMeshSegment::UnionMesh>& segmentsUnionMesh,
                                                        vector<ConformerMeshSegment::ConformMesh>& segmentsConformMesh,
                                                        vector<ConformerMeshPolygon::ConformMesh>& segmentsConformMeshInfo,
                                                        const ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types& conformDomainMeshType,
                                                        const ComputeDomainConformedMeshOptions& options) const
  {
    const unsigned int numberOfInterfaces = segmentsVertices.size();

    for (unsigned int i = 0; i < numberOfInterfaces; i++)
    {
      if (options.PrintStatus)
      {
        Gedim::Output::PrintGenericMessage("Compute Interface Mesh %d / %d...", true,
                                           i,
                                           numberOfInterfaces - 1);
      }

      // Intersect interface with domain mesh
      const Eigen::MatrixXd& interface2D = segmentsVertices.at(i);
      const Eigen::Vector3d& interface2DTangent = segmentsTangent.at(i);
      const Eigen::Vector3d& interface2DBarycenter = segmentsBarycenter.at(i);
      const double& interface2DLength = segmentsLength.at(i);

      Gedim::IntersectorMesh2DSegment intersectorMesh(domainMesh,
                                                      geometryUtilities);

      intersectorMesh.CreateIntersectionMesh(interface2D.col(0),
                                             interface2D.col(1),
                                             interface2DTangent,
                                             interface2DBarycenter,
                                             interface2DLength,
                                             segmentsIntersectionMesh[i]);

      Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(segmentsIntersectionMesh[i],
                                                                segmentsCurvilinearCoordinatesMesh[i]);

      for (const Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshSegment& segment : segmentsIntersectionMesh[i].Segments)
        Gedim::Output::Assert(geometryUtilities.IsValue1DPositive(abs(segment.Points[1] - segment.Points[0])));

      if (options.PrintStatus)
      {
        Gedim::Output::PrintStatusProgram(" %d / %d: Intersect Interface Mesh",
                                          i,
                                          numberOfInterfaces - 1);
      }

      // Unify interface mesh
      Gedim::UnionMeshSegment unionMesher(geometryUtilities);
      unionMesher.CreateUnionMesh(segmentsCurvilinearCoordinatesMesh.at(i),
                                  segmentsCurvilinearCoordinatesMesh.at(i),
                                  segmentsUnionMesh.at(i));

      if (options.PrintStatus)
      {
        Gedim::Output::PrintStatusProgram(" %d / %d: Unify Interface Mesh",
                                          i,
                                          numberOfInterfaces - 1);
      }

      // Conform interface and domain mesh
      Gedim::ConformerMeshSegment conformMeshInterface(geometryUtilities);

      // Add mesh 1D additional points
      for (const double& additionalPoint : segmentsAdditionalPoints[i])
        conformMeshInterface.InsertExternalPoint(interface2D.col(0),
                                                 interface2D.col(1),
                                                 domainMesh,
                                                 additionalPoint,
                                                 segmentsConformMesh[i]);

      // Conform mesh 1D
      segmentsConformMesh[i].Segments.clear();
      conformMeshInterface.CreateConformMesh(segmentsIntersectionMesh[i],
                                             segmentsUnionMesh[i],
                                             0,
                                             segmentsConformMesh[i]);

      for (const Gedim::ConformerMeshSegment::ConformMesh::ConformMeshSegment& segment : segmentsConformMesh[i].Segments)
        Gedim::Output::Assert(geometryUtilities.IsValue1DPositive(abs(segment.Points[1] - segment.Points[0])));

      // Conform mesh 2D
      Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration conformMeshDomainConfiguration;
      conformMeshDomainConfiguration.Type = static_cast<Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types>(conformDomainMeshType);
      Gedim::ConformerMeshPolygon conformerMeshDomain(geometryUtilities,
                                                      conformMeshDomainConfiguration);
      conformerMeshDomain.CreateConformMesh(interface2D.col(0),
                                            interface2D.col(1),
                                            interface2DTangent,
                                            segmentsConformMesh[i],
                                            domainMesh,
                                            segmentsConformMeshInfo[i]);

      if (!options.VtkExportFolder.empty())
      {
        meshUtilities.ExportMeshToVTU(domainMesh,
                                      options.VtkExportFolder,
                                      "ConformedMesh");
      }

      // Update mesh 1D with updated mesh 2D
      for (unsigned int di = 0; di < i; di++)
      {
        conformMeshInterface.UpdateWithUpdatedMesh2D(domainMesh,
                                                     segmentsConformMesh[di]);
      }

      // Clean mesh 2D
      Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
      meshUtilities.ExtractActiveMesh(domainMesh,
                                      extractionData);

      if (!options.VtkExportFolder.empty())
      {
        meshUtilities.ExportMeshToVTU(domainMesh,
                                      options.VtkExportFolder,
                                      "ConformedMesh");
      }

      // Update mesh 1D with cleaned mesh 2D
      for (unsigned int di = 0; di < i + 1; di++)
      {
        conformMeshInterface.UpdateWithActiveMesh2D(extractionData,
                                                    segmentsConformMesh[di]);

        conformMeshInterface.AddMissingMesh2DCell0Ds(segmentsVertices[di].col(0),
                                                     segmentsTangent[di],
                                                     segmentsSquaredLength[di],
                                                     domainMesh,
                                                     segmentsConformMesh[di]);
      }

      if (options.PrintStatus)
      {
        Gedim::Output::PrintStatusProgram("Compute Interface Mesh %d / %d",
                                          i,
                                          numberOfInterfaces - 1);
      }
    }
  }
  // ***************************************************************************
}
