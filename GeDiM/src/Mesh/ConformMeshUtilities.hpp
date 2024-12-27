#ifndef __ConformMeshUtilities_H
#define __ConformMeshUtilities_H

#include "IMeshDAO.hpp"
#include <list>

#include "ConformerMeshPolygon.hpp"
#include "MeshUtilities.hpp"

namespace Gedim
{
class ConformMeshUtilities final
{
  public:
    struct ComputeDomainConformedMeshOptions final
    {
        bool PrintStatus = false;
        std::string VtkExportFolder = "";
    };

  private:
    const GeometryUtilities &geometryUtilities;
    const MeshUtilities &meshUtilities;

  public:
    ConformMeshUtilities(const GeometryUtilities &geometryUtilities, const MeshUtilities &meshUtilities);
    ~ConformMeshUtilities();

    void ComputeConformedMeshWithSegments(const std::vector<std::list<double>> &segmentsAdditionalPoints,
                                          const std::vector<Eigen::MatrixXd> &segmentsVertices,
                                          const std::vector<Eigen::Vector3d> &segmentsTangent,
                                          const std::vector<Eigen::Vector3d> &segmentsBarycenter,
                                          const std::vector<double> &segmentsLength,
                                          const std::vector<double> &segmentsSquaredLength,
                                          IMeshDAO &domainMesh,
                                          std::vector<IntersectorMesh2DSegment::IntersectionMesh> &segmentsIntersectionMesh,
                                          std::vector<std::vector<double>> &segmentsCurvilinearCoordinatesMesh,
                                          std::vector<UnionMeshSegment::UnionMesh> &segmentsUnionMesh,
                                          std::vector<ConformerMeshSegment::ConformMesh> &segmentsConformMesh,
                                          const ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types &conformDomainMeshType,
                                          const ComputeDomainConformedMeshOptions &options) const;

    void AddConformedMeshProperties(const std::vector<ConformerMeshSegment::ConformMesh> &segmentsConformMesh,
                                    IMeshDAO &conformedMesh) const;
};
} // namespace Gedim

#endif
