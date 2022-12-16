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
  void ConformMeshUtilities::ComputeConformedMeshWithSegments(const std::vector<std::list<double>>& segmentsAdditionalPoints,
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
      const Eigen::MatrixXd& segmentVertices = segmentsVertices.at(i);
      const Eigen::Vector3d& segmentTangent = segmentsTangent.at(i);
      const Eigen::Vector3d& segmentBarycenter = segmentsBarycenter.at(i);
      const double& segmentLength = segmentsLength.at(i);

      Gedim::IntersectorMesh2DSegment intersectorMesh(domainMesh,
                                                      geometryUtilities);

      intersectorMesh.CreateIntersectionMesh(segmentVertices.col(0),
                                             segmentVertices.col(1),
                                             segmentTangent,
                                             segmentBarycenter,
                                             segmentLength,
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
        conformMeshInterface.InsertExternalPoint(segmentVertices.col(0),
                                                 segmentVertices.col(1),
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
      conformerMeshDomain.CreateConformMesh(segmentVertices.col(0),
                                            segmentVertices.col(1),
                                            segmentTangent,
                                            segmentsConformMesh[i],
                                            domainMesh);

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
  void ConformMeshUtilities::AddConformedMeshProperties(IMeshDAO& networkMesh) const
  {
    networkMesh.Cell0DInitializeDoubleProperties(1);
    networkMesh.Cell0DAddDoubleProperty("marked");
    networkMesh.Cell1DInitializeDoubleProperties(2);
    networkMesh.Cell1DAddDoubleProperty("marked");
    networkMesh.Cell1DAddDoubleProperty("interface");

    for (unsigned int c = 0; c < networkMesh.Cell0DTotalNumber(); c++)
    {
      for (unsigned int p = 0; p < networkMesh.Cell0DNumberDoubleProperties(); p++)
      {
        networkMesh.Cell0DInitializeDoublePropertyValues(c, p, 1);
        networkMesh.Cell0DInsertDoublePropertyValue(c, p, 0, 0.0);
      }
    }
    for (unsigned int e = 0; e < networkMesh.Cell1DTotalNumber(); e++)
    {
      for (unsigned int p = 0; p < networkMesh.Cell1DNumberDoubleProperties(); p++)
      {
        networkMesh.Cell1DInitializeDoublePropertyValues(e, p, 1);
        networkMesh.Cell1DInsertDoublePropertyValue(e, p, 0, 0.0);
      }
    }
    for (unsigned int f = 0; f < networkMesh.Cell2DTotalNumber(); f++)
    {
      for (unsigned int p = 0; p < networkMesh.Cell2DNumberDoubleProperties(); p++)
      {
        networkMesh.Cell2DInitializeDoublePropertyValues(f, p, 1);
        networkMesh.Cell2DInsertDoublePropertyValue(f, p, 0, 0.0);
      }
    }

    //    if (networkMeshInformation.InterfacesCell0DsToNetworkCell0Ds.size() > 0 &&
    //        networkMeshInformation.InterfacesCell1DsToNetworkCell1Ds.size() > 0)
    //    {
    //      const unsigned int numInterfaces = networkMeshInformation.InterfacesCell0DsToNetworkCell0Ds.size();

    //      for (unsigned int i = 0; i < numInterfaces; i++)
    //      {
    //        const unsigned int interfaceId = i;

    //        for (map<unsigned int, unsigned int>::const_iterator cell1D_it = networkMeshInformation.InterfacesCell1DsToNetworkCell1Ds[i].begin();
    //             cell1D_it != networkMeshInformation.InterfacesCell1DsToNetworkCell1Ds[i].end();
    //             cell1D_it++)
    //        {
    //          const unsigned int networkCell1DId = cell1D_it->second;

    //          networkMesh.Cell1DInsertDoublePropertyValue(networkCell1DId,
    //                                                      0,
    //                                                      0,
    //                                                      1.0);
    //          networkMesh.Cell1DInsertDoublePropertyValue(networkCell1DId,
    //                                                      1,
    //                                                      0,
    //                                                      interfaceId + 1);

    //        }

    //        for (map<double, unsigned int>::const_iterator cell0D_it = networkMeshInformation.InterfacesCell0DsToNetworkCell0Ds[i].begin();
    //             cell0D_it != networkMeshInformation.InterfacesCell0DsToNetworkCell0Ds[i].end();
    //             cell0D_it++)
    //        {
    //          const unsigned int networkCell0DId = cell0D_it->second;

    //          networkMesh.Cell0DInsertDoublePropertyValue(networkCell0DId,
    //                                                      0,
    //                                                      0,
    //                                                      1.0);
    //        }
    //      }
    //    }
  }
  // ***************************************************************************
}
