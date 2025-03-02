#include "ConformMeshUtilities.hpp"

using namespace std;

namespace Gedim
{
// ***************************************************************************
ConformMeshUtilities::ConformMeshUtilities(const GeometryUtilities &geometryUtilities, const MeshUtilities &meshUtilities)
    : geometryUtilities(geometryUtilities), meshUtilities(meshUtilities)
{
}
ConformMeshUtilities::~ConformMeshUtilities()
{
}
// ***************************************************************************
void ConformMeshUtilities::ComputeConformedMeshWithSegments(const std::vector<std::list<double>> &segmentsAdditionalPoints,
                                                            const vector<Eigen::MatrixXd> &segmentsVertices,
                                                            const vector<Eigen::Vector3d> &segmentsTangent,
                                                            const vector<Eigen::Vector3d> &segmentsBarycenter,
                                                            const vector<double> &segmentsLength,
                                                            const vector<double> &segmentsSquaredLength,
                                                            IMeshDAO &domainMesh,
                                                            vector<IntersectorMesh2DSegment::IntersectionMesh> &segmentsIntersectionMesh,
                                                            std::vector<std::vector<double>> &segmentsCurvilinearCoordinatesMesh,
                                                            vector<UnionMeshSegment::UnionMesh> &segmentsUnionMesh,
                                                            vector<ConformerMeshSegment::ConformMesh> &segmentsConformMesh,
                                                            const ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types &conformDomainMeshType,
                                                            const ComputeDomainConformedMeshOptions &options) const
{
    const unsigned int numberOfInterfaces = segmentsVertices.size();

    for (unsigned int i = 0; i < numberOfInterfaces; i++)
    {
        if (options.PrintStatus)
        {
            Gedim::Output::PrintGenericMessage("Compute Interface Mesh %d / %d...", true, i, numberOfInterfaces - 1);
        }

        // Intersect interface with domain mesh
        const Eigen::MatrixXd &segmentVertices = segmentsVertices.at(i);
        const Eigen::Vector3d &segmentTangent = segmentsTangent.at(i);
        const Eigen::Vector3d &segmentBarycenter = segmentsBarycenter.at(i);
        const double &segmentLength = segmentsLength.at(i);

        Gedim::IntersectorMesh2DSegment intersectorMesh(domainMesh, geometryUtilities);

        intersectorMesh.CreateIntersectionMesh(segmentVertices.col(0),
                                               segmentVertices.col(1),
                                               segmentTangent,
                                               segmentBarycenter,
                                               segmentLength,
                                               segmentsIntersectionMesh[i]);

        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(segmentsIntersectionMesh[i],
                                                                  segmentsCurvilinearCoordinatesMesh[i]);

        for (const Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshSegment &segment :
             segmentsIntersectionMesh[i].Segments)
            Gedim::Output::Assert(geometryUtilities.IsValuePositive(std::abs(segment.Points[1] - segment.Points[0]),
                                                                    geometryUtilities.Tolerance1D()));

        if (options.PrintStatus)
        {
            Gedim::Output::PrintStatusProgram(" %d / %d: Intersect Interface Mesh", i, numberOfInterfaces - 1);
        }

        // Unify interface mesh
        Gedim::UnionMeshSegment unionMesher(geometryUtilities);
        unionMesher.CreateUnionMesh(segmentsCurvilinearCoordinatesMesh.at(i),
                                    segmentsCurvilinearCoordinatesMesh.at(i),
                                    segmentsUnionMesh.at(i));

        if (options.PrintStatus)
        {
            Gedim::Output::PrintStatusProgram(" %d / %d: Unify Interface Mesh", i, numberOfInterfaces - 1);
        }

        // Conform interface and domain mesh
        Gedim::ConformerMeshSegment conformMeshInterface(geometryUtilities);

        // Add mesh 1D additional points
        for (const double &additionalPoint : segmentsAdditionalPoints[i])
            conformMeshInterface.InsertExternalPoint(segmentVertices.col(0),
                                                     segmentVertices.col(1),
                                                     domainMesh,
                                                     additionalPoint,
                                                     segmentsConformMesh[i]);

        // Conform mesh 1D
        segmentsConformMesh[i].Segments.clear();
        conformMeshInterface.CreateConformMesh(segmentsIntersectionMesh[i], segmentsUnionMesh[i], 0, segmentsConformMesh[i]);

        for (const Gedim::ConformerMeshSegment::ConformMesh::ConformMeshSegment &segment : segmentsConformMesh[i].Segments)
            Gedim::Output::Assert(geometryUtilities.IsValuePositive(std::abs(segment.Points[1] - segment.Points[0]),
                                                                    geometryUtilities.Tolerance1D()));

        // Conform mesh 2D
        Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration conformMeshDomainConfiguration;
        conformMeshDomainConfiguration.Type =
            static_cast<Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types>(conformDomainMeshType);
        Gedim::ConformerMeshPolygon conformerMeshDomain(geometryUtilities, conformMeshDomainConfiguration);
        conformerMeshDomain.CreateConformMesh(segmentVertices.col(0), segmentVertices.col(1), segmentTangent, segmentsConformMesh[i], domainMesh);

        if (!options.VtkExportFolder.empty())
        {
            meshUtilities.ExportMeshToVTU(domainMesh, options.VtkExportFolder, "ConformedMesh");
        }

        // Update mesh 1D with updated mesh 2D
        for (unsigned int di = 0; di < i; di++)
        {
            conformMeshInterface.UpdateWithUpdatedMesh2D(domainMesh, segmentsConformMesh[di]);
        }

        // Clean mesh 2D
        Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
        meshUtilities.ExtractActiveMesh(domainMesh, extractionData);

        if (!options.VtkExportFolder.empty())
        {
            meshUtilities.ExportMeshToVTU(domainMesh, options.VtkExportFolder, "ConformedMesh");
        }

        // Update mesh 1D with cleaned mesh 2D
        for (unsigned int di = 0; di < i + 1; di++)
        {
            conformMeshInterface.UpdateWithActiveMesh2D(extractionData, segmentsConformMesh[di]);

            conformMeshInterface.AddMissingMesh2DCell0Ds(segmentsVertices[di].col(0),
                                                         segmentsTangent[di],
                                                         segmentsSquaredLength[di],
                                                         domainMesh,
                                                         segmentsConformMesh[di]);
        }

        if (options.PrintStatus)
        {
            Gedim::Output::PrintStatusProgram("Compute Interface Mesh %d / %d", i, numberOfInterfaces - 1);
        }
    }
}
// ***************************************************************************
void ConformMeshUtilities::AddConformedMeshProperties(const std::vector<ConformerMeshSegment::ConformMesh> &segmentsConformMesh,
                                                      IMeshDAO &conformedMesh) const
{
    conformedMesh.Cell0DInitializeDoubleProperties(1);
    conformedMesh.Cell0DAddDoubleProperty("marked");
    conformedMesh.Cell1DInitializeDoubleProperties(2);
    conformedMesh.Cell1DAddDoubleProperty("marked");
    conformedMesh.Cell1DAddDoubleProperty("segment");

    for (unsigned int p = 0; p < conformedMesh.Cell0DNumberDoubleProperties(); p++)
        conformedMesh.Cell0DsInitializeDoublePropertyValues(p, std::vector<unsigned int>(conformedMesh.Cell0DTotalNumber(), 1));

    for (unsigned int c = 0; c < conformedMesh.Cell0DTotalNumber(); c++)
    {
        for (unsigned int p = 0; p < conformedMesh.Cell0DNumberDoubleProperties(); p++)
            conformedMesh.Cell0DInsertDoublePropertyValue(c, p, 0, 0.0);
    }

    for (unsigned int p = 0; p < conformedMesh.Cell1DNumberDoubleProperties(); p++)
        conformedMesh.Cell1DsInitializeDoublePropertyValues(p, std::vector<unsigned int>(conformedMesh.Cell1DTotalNumber(), 1));

    for (unsigned int e = 0; e < conformedMesh.Cell1DTotalNumber(); e++)
    {
        for (unsigned int p = 0; p < conformedMesh.Cell1DNumberDoubleProperties(); p++)
            conformedMesh.Cell1DInsertDoublePropertyValue(e, p, 0, 0.0);
    }

    for (unsigned int s = 0; s < segmentsConformMesh.size(); s++)
    {
        const ConformerMeshSegment::ConformMesh &segmentConformMesh = segmentsConformMesh[s];

        for (const auto &point : segmentConformMesh.Points)
        {
            Output::Assert(point.second.Vertex2DIds.size() == 1);
            conformedMesh.Cell0DInsertDoublePropertyValue(*point.second.Vertex2DIds.begin(), 0, 0, 1.0);
        }

        for (const ConformerMeshSegment::ConformMesh::ConformMeshSegment &segment : segmentConformMesh.Segments)
        {
            Output::Assert(segment.Edge2DIds.size() == 1);
            conformedMesh.Cell1DInsertDoublePropertyValue(*segment.Edge2DIds.begin(), 0, 0, 1.0);
            conformedMesh.Cell1DInsertDoublePropertyValue(*segment.Edge2DIds.begin(), 1, 0, s + 1);
        }
    }
}
// ***************************************************************************
} // namespace Gedim
