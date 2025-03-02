#include "ConformerMeshSegment.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
// ***************************************************************************

ConformerMeshSegment::ConformerMeshSegment(const GeometryUtilities &geometryUtilities)
    : _geometryUtilities(geometryUtilities)
{
}

ConformerMeshSegment::~ConformerMeshSegment()
{
}

// ***************************************************************************
void ConformerMeshSegment::ToCurvilinearCoordinates(const ConformerMeshSegment::ConformMesh &conformMesh, vector<double> &curvilinearCoordinates)
{
    curvilinearCoordinates.reserve(conformMesh.Points.size());
    for (std::map<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint>::const_iterator it = conformMesh.Points.begin();
         it != conformMesh.Points.end();
         it++)
    {
        curvilinearCoordinates.push_back(it->first);
    }
}
// ***************************************************************************
string ConformerMeshSegment::ToString(const ConformerMeshSegment::ConformMesh &conformMesh)
{
    ostringstream converter;

    converter.precision(16);
    converter << "Points:" << endl;
    for_each(conformMesh.Points.begin(),
             conformMesh.Points.end(),
             [&converter](const std::pair<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint> &p) {
                 converter << scientific << "{ Key: " << p.first;
                 converter << "; Value:";
                 converter << " Type: " << (unsigned int)p.second.Type;
                 converter << " V: " << p.second.Vertex2DIds;
                 converter << " E: " << p.second.Edge2DIds;
                 converter << " C: " << p.second.Cell2DIds;
                 converter << " }\n";
             });
    converter << "Segments:" << endl;
    for_each(conformMesh.Segments.begin(),
             conformMesh.Segments.end(),
             [&converter](const ConformerMeshSegment::ConformMesh::ConformMeshSegment &p) {
                 converter << scientific << "{ Key: " << p.Points[0] << "->" << p.Points[1]
                           << "; Value: E: " << p.Edge2DIds << " C: " << p.Cell2DIds << " }\n";
             });

    return converter.str();
}
// ***************************************************************************
ConformerMeshSegment::ConformMesh::ConformMeshPoint &ConformerMeshSegment::InsertNewIntersection(const double &curvilinearCoordinate,
                                                                                                 ConformerMeshSegment::ConformMesh &result,
                                                                                                 bool &found) const
{
    double foundCoordinate = -1.0;
    for (std::map<double, ConformMesh::ConformMeshPoint>::const_iterator it = result.Points.begin(); it != result.Points.end(); it++)
    {
        if (!_geometryUtilities.IsValuePositive(std::abs(it->first - curvilinearCoordinate), _geometryUtilities.Tolerance1D()))
        {
            foundCoordinate = it->first;
            break;
        }
    }

    if (foundCoordinate != -1.0)
    {
        found = true;
        return result.Points[foundCoordinate];
    }

    result.Points.insert(pair<double, ConformMesh::ConformMeshPoint>(curvilinearCoordinate, ConformMesh::ConformMeshPoint()));
    found = false;
    return result.Points[curvilinearCoordinate];
}
// ***************************************************************************
void ConformerMeshSegment::CreateConformPoints(const Gedim::IntersectorMesh2DSegment::IntersectionMesh &meshIntersection,
                                               const Gedim::UnionMeshSegment::UnionMesh &meshUnion,
                                               const unsigned int &meshIntersectionPosition,
                                               ConformerMeshSegment::ConformMesh &result)
{
    map<double, Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>::const_iterator itIntersectionPoint =
        meshIntersection.Points.begin();
    unsigned int intersectionSegmentIndex = 0;
    map<double, UnionMeshSegment::UnionMesh::UnionMeshPoint>::const_iterator itUnionPoint = meshUnion.Points.begin();

    UnionMeshSegment::UnionMesh::UnionMeshPoint::Types otherIntersectionMeshUnionType =
        (meshIntersectionPosition == 0) ? UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Second
                                        : UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::First;
    for (unsigned int p = 0; p < meshUnion.Points.size(); p++)
    {
        const double &unionCurvilinearCoordinate = itUnionPoint->first;

        const UnionMeshSegment::UnionMesh::UnionMeshPoint &unionPoint = itUnionPoint->second;

        bool found;
        ConformerMeshSegment::ConformMesh::ConformMeshPoint &conformPoint =
            InsertNewIntersection(unionCurvilinearCoordinate, result, found);
        // check if union point is in intersection mesh
        if (unionPoint.Type != otherIntersectionMeshUnionType)
        {
            // union point belongs to mine intersection mesh
            // inherit properties from intersection mesh point
            const IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint &intersectionPoint =
                itIntersectionPoint->second;
            conformPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;

            conformPoint.Vertex2DIds = intersectionPoint.Vertex2DIds; // copy assignment
            conformPoint.Edge2DIds = intersectionPoint.Edge2DIds;     // copy assignment
            conformPoint.Cell2DIds = intersectionPoint.Cell2DIds;     // copy assignment

            itIntersectionPoint++;
            intersectionSegmentIndex++;
        }
        else
        {
            Output::Assert(intersectionSegmentIndex > 0 && intersectionSegmentIndex <= meshIntersection.Segments.size());

            // union point belongs to the other intersection mesh
            // add the new point to conform mesh and inherit properties from intersection mesh segment
            const IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshSegment &intersectionSegment =
                meshIntersection.Segments[intersectionSegmentIndex - 1];
            conformPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::Inherited;

            conformPoint.Edge2DIds = intersectionSegment.Edge2DIds; // copy assignment
            conformPoint.Cell2DIds = intersectionSegment.Cell2DIds; // copy assignment
        }

        itUnionPoint++;
    }
}
// ***************************************************************************
void ConformerMeshSegment::CreateConformSegments(ConformerMeshSegment::ConformMesh &result)
{
    result.Segments.resize(result.Points.size() - 1);
    map<double, ConformMesh::ConformMeshPoint>::const_iterator itPoint = result.Points.begin();
    map<double, ConformMesh::ConformMeshPoint>::const_iterator itPointNext = result.Points.begin();
    itPointNext++;
    for (unsigned int p = 0; p < result.Segments.size(); p++)
    {
        const double &curvilinearCoordinatePoint = itPoint->first;
        const double &curvilinearCoordinatePointNext = itPointNext->first;
        const ConformMesh::ConformMeshPoint &intersectionPoint = itPoint->second;
        const ConformMesh::ConformMeshPoint &intersectionPointNext = itPointNext->second;

        // fill origin and end of segment
        ConformMesh::ConformMeshSegment &meshSegment = result.Segments[p];
        meshSegment.Points.resize(2);
        meshSegment.Points[0] = curvilinearCoordinatePoint;
        meshSegment.Points[1] = curvilinearCoordinatePointNext;

        // fill the mesh 2D edges
        set_intersection(intersectionPoint.Edge2DIds.begin(),
                         intersectionPoint.Edge2DIds.end(),
                         intersectionPointNext.Edge2DIds.begin(),
                         intersectionPointNext.Edge2DIds.end(),
                         std::inserter(meshSegment.Edge2DIds, meshSegment.Edge2DIds.begin()));
        // fill the mesh 2D cells
        set_intersection(intersectionPoint.Cell2DIds.begin(),
                         intersectionPoint.Cell2DIds.end(),
                         intersectionPointNext.Cell2DIds.begin(),
                         intersectionPointNext.Cell2DIds.end(),
                         std::inserter(meshSegment.Cell2DIds, meshSegment.Cell2DIds.begin()));
        itPoint++;
        itPointNext++;
    }
}
// ***************************************************************************
void ConformerMeshSegment::Serialize(ostream &os, const ConformMesh &mesh)

{
    char separator = ' ';
    os.precision(16);

    os << scientific << mesh.Points.size() << separator;
    for_each(mesh.Points.begin(),
             mesh.Points.end(),
             [&os, &separator](const std::pair<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint> &p) {
                 os << scientific << p.first << separator;
                 os << scientific << (unsigned int)p.second.Type << separator;

                 os << scientific << p.second.Vertex2DIds.size() << separator;
                 for_each(p.second.Vertex2DIds.begin(), p.second.Vertex2DIds.end(), [&os, &separator](const unsigned int &id) {
                     os << scientific << id << separator;
                 });

                 os << scientific << p.second.Edge2DIds.size() << separator;
                 for_each(p.second.Edge2DIds.begin(), p.second.Edge2DIds.end(), [&os, &separator](const unsigned int &id) {
                     os << scientific << id << separator;
                 });

                 os << scientific << p.second.Cell2DIds.size() << separator;
                 for_each(p.second.Cell2DIds.begin(), p.second.Cell2DIds.end(), [&os, &separator](const unsigned int &id) {
                     os << scientific << id << separator;
                 });
             });

    os << scientific << mesh.Segments.size() << separator;
    for_each(mesh.Segments.begin(), mesh.Segments.end(), [&os, &separator](const ConformerMeshSegment::ConformMesh::ConformMeshSegment &s) {
        os << scientific << s.Points[0] << separator;
        os << scientific << s.Points[1] << separator;

        os << scientific << s.Edge2DIds.size() << separator;
        for_each(s.Edge2DIds.begin(), s.Edge2DIds.end(), [&os, &separator](const unsigned int &id) {
            os << scientific << id << separator;
        });

        os << scientific << s.Cell2DIds.size() << separator;
        for_each(s.Cell2DIds.begin(), s.Cell2DIds.end(), [&os, &separator](const unsigned int &id) {
            os << scientific << id << separator;
        });
    });
}
// ***************************************************************************
void ConformerMeshSegment::Deserialize(istream &is, ConformMesh &mesh)
{
    unsigned int numMeshPoints = 0;
    is >> numMeshPoints;

    for (unsigned int p = 0; p < numMeshPoints; p++)
    {
        double pointCurvilinearCoordinate;
        is >> pointCurvilinearCoordinate;

        mesh.Points.insert(pair<double, ConformMesh::ConformMeshPoint>(pointCurvilinearCoordinate, ConformMesh::ConformMeshPoint()));
        ConformMesh::ConformMeshPoint &point = mesh.Points.at(pointCurvilinearCoordinate);

        unsigned int pointType;
        is >> pointType;

        switch (pointType)
        {
        case 1:
            point.Type = ConformMesh::ConformMeshPoint::Types::Original;
            break;
        case 2:
            point.Type = ConformMesh::ConformMeshPoint::Types::Inherited;
            break;
        case 3:
            point.Type = ConformMesh::ConformMeshPoint::Types::External;
            break;
        default:
            throw runtime_error("Unknown point type");
        }

        unsigned int pointVerticesSize = 0;
        is >> pointVerticesSize;
        for (unsigned int pv = 0; pv < pointVerticesSize; pv++)
        {
            unsigned int pointVertex;
            is >> pointVertex;
            point.Vertex2DIds.insert(pointVertex);
        }

        unsigned int pointEdgesSize = 0;
        is >> pointEdgesSize;
        for (unsigned int pe = 0; pe < pointEdgesSize; pe++)
        {
            unsigned int pointEdge;
            is >> pointEdge;
            point.Edge2DIds.insert(pointEdge);
        }

        unsigned int pointCell2DSize = 0;
        is >> pointCell2DSize;
        for (unsigned int pc = 0; pc < pointCell2DSize; pc++)
        {
            unsigned int pointCell2D;
            is >> pointCell2D;
            point.Cell2DIds.insert(pointCell2D);
        }
    }

    unsigned int numMeshSegments = 0;
    is >> numMeshSegments;

    mesh.Segments.resize(numMeshSegments);
    for (unsigned int s = 0; s < numMeshSegments; s++)
    {
        ConformMesh::ConformMeshSegment &segment = mesh.Segments.at(s);
        segment.Points.resize(2);
        is >> segment.Points[0];
        is >> segment.Points[1];

        unsigned int segmentEdgesSize = 0;
        is >> segmentEdgesSize;
        for (unsigned int se = 0; se < segmentEdgesSize; se++)
        {
            unsigned int segmentEdge;
            is >> segmentEdge;
            segment.Edge2DIds.insert(segmentEdge);
        }

        unsigned int segmentCell2DSize = 0;
        is >> segmentCell2DSize;
        for (unsigned int sc = 0; sc < segmentCell2DSize; sc++)
        {
            unsigned int segmentCell2D;
            is >> segmentCell2D;
            segment.Cell2DIds.insert(segmentCell2D);
        }
    }
}
// ***************************************************************************
void ConformerMeshSegment::CreateConformMesh(const Gedim::IntersectorMesh2DSegment::IntersectionMesh &meshIntersection,
                                             const UnionMeshSegment::UnionMesh &meshUnion,
                                             const unsigned int &meshIntersectionPosition,
                                             ConformMesh &result)
{
    CreateConformPoints(meshIntersection, meshUnion, meshIntersectionPosition, result);
    CreateConformSegments(result);
}
// ***************************************************************************
void ConformerMeshSegment::InsertExternalPoint(const Vector3d &segmentOrigin,
                                               const Vector3d &segmentEnd,
                                               const Gedim::IMeshDAO &mesh2D,
                                               const double &curvilinearCoordinate,
                                               ConformerMeshSegment::ConformMesh &result)
{
    bool found;
    ConformerMeshSegment::ConformMesh::ConformMeshPoint &newPoint = InsertNewIntersection(curvilinearCoordinate, result, found);
    // Point already added
    if (found)
        return;

    newPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::External;
    const Vector3d newPoint2D = segmentOrigin + curvilinearCoordinate * _geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);

    // look for 2D cell which contains the new point
    for (unsigned int c = 0; c < mesh2D.Cell2DTotalNumber(); c++)
    {
        if (!mesh2D.Cell2DIsActive(c))
            continue;

        const MatrixXd cell2DVertices = mesh2D.Cell2DVerticesCoordinates(c);

        const Eigen::MatrixXd cell2D_bounding_box = _geometryUtilities.PointsBoundingBox(cell2DVertices);

        if (!_geometryUtilities.IsPointInBoundingBox(newPoint2D, cell2D_bounding_box))
            continue;

        GeometryUtilities::PointPolygonPositionResult pointPolygonPositionResult =
            _geometryUtilities.PointPolygonPosition(newPoint2D, cell2DVertices);

        if (pointPolygonPositionResult.Type == GeometryUtilities::PointPolygonPositionResult::Types::Outside)
            continue;

        newPoint.Cell2DIds.insert(c);
        if (pointPolygonPositionResult.Type == GeometryUtilities::PointPolygonPositionResult::Types::Inside)
            break;

        else if (pointPolygonPositionResult.Type == GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge)
        {
            const unsigned int edgeToInsert = mesh2D.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
            if (newPoint.Edge2DIds.find(edgeToInsert) == newPoint.Edge2DIds.end())
                newPoint.Edge2DIds.insert(edgeToInsert);
        }
        else if (pointPolygonPositionResult.Type == GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex)
        {
            const unsigned int vertexToInsert = mesh2D.Cell2DVertex(c, pointPolygonPositionResult.BorderIndex);
            int previousVertex = pointPolygonPositionResult.BorderIndex - 1;
            const unsigned int edgeToInsert = mesh2D.Cell2DEdge(c, pointPolygonPositionResult.BorderIndex);
            const unsigned int previousEdgeToInsert =
                mesh2D.Cell2DEdge(c, (previousVertex < 0) ? mesh2D.Cell2DNumberVertices(c) - 1 : previousVertex);
            if (newPoint.Vertex2DIds.find(vertexToInsert) == newPoint.Vertex2DIds.end())
                newPoint.Vertex2DIds.insert(vertexToInsert);
            if (newPoint.Edge2DIds.find(edgeToInsert) == newPoint.Edge2DIds.end())
                newPoint.Vertex2DIds.insert(edgeToInsert);
            if (newPoint.Edge2DIds.find(previousEdgeToInsert) == newPoint.Edge2DIds.end())
                newPoint.Vertex2DIds.insert(previousEdgeToInsert);
        }
    }

    // recreate conform mesh segments
    result.Segments.clear();
    CreateConformSegments(result);
}
// ***************************************************************************
void ConformerMeshSegment::UpdateWithUpdatedMesh2D(const IMeshDAO &mesh2D, ConformMesh &conformedMesh) const
{
    for (map<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint>::iterator it = conformedMesh.Points.begin();
         it != conformedMesh.Points.end();
         it++)
    {
        ConformerMeshSegment::ConformMesh::ConformMeshPoint &point = it->second;

        {
            list<unsigned int> newCell0DIndices;
            for (const unsigned int &v : point.Vertex2DIds)
            {
                if (!mesh2D.Cell0DHasUpdatedCell0Ds(v))
                    continue;

                list<unsigned int> updatedCell0Ds;
                mesh2D.Cell0DUpdatedCell0Ds(v, updatedCell0Ds);
                for (const unsigned int &nv : updatedCell0Ds)
                    newCell0DIndices.push_back(nv);
            }
            for (const unsigned int &nv : newCell0DIndices)
                point.Vertex2DIds.insert(nv);
        }

        Output::Assert(point.Vertex2DIds.size() == 1);
        const unsigned int mesh2DCell0DIndex = *point.Vertex2DIds.rbegin();

        {
            list<unsigned int> newCell1DIndices;
            for (const unsigned int &e : point.Edge2DIds)
            {
                if (!mesh2D.Cell1DHasUpdatedCell1Ds(e))
                    continue;

                list<unsigned int> updatedCell1Ds;
                mesh2D.Cell1DUpdatedCell1Ds(e, updatedCell1Ds);
                for (const unsigned int &ne : updatedCell1Ds)
                {
                    bool hasCell0D = false;
                    for (unsigned int nv = 0; nv < 2; nv++)
                    {
                        hasCell0D = (mesh2D.Cell1DVertex(ne, nv) == mesh2DCell0DIndex);

                        if (hasCell0D)
                            break;
                    }

                    if (!hasCell0D)
                        continue;

                    newCell1DIndices.push_back(ne);
                }
            }

            for (const unsigned int &ne : newCell1DIndices)
                point.Edge2DIds.insert(ne);
        }

        {
            list<unsigned int> newCell2DIndices;
            for (const unsigned int &c : point.Cell2DIds)
            {
                if (!mesh2D.Cell2DHasUpdatedCell2Ds(c))
                    continue;

                list<unsigned int> updatedCell2Ds;
                mesh2D.Cell2DUpdatedCell2Ds(c, updatedCell2Ds);
                for (const unsigned int &nc : updatedCell2Ds)
                {
                    bool hasCell0D = false;
                    for (unsigned int nv = 0; nv < mesh2D.Cell2DNumberVertices(nc); nv++)
                    {
                        hasCell0D = (mesh2D.Cell2DVertex(nc, nv) == mesh2DCell0DIndex);

                        if (hasCell0D)
                            break;
                    }

                    if (!hasCell0D)
                        continue;

                    newCell2DIndices.push_back(nc);
                }
            }

            for (const unsigned int &nc : newCell2DIndices)
                point.Cell2DIds.insert(nc);
        }
    }

    for (ConformerMeshSegment::ConformMesh::ConformMeshSegment &segment : conformedMesh.Segments)
    {
        {
            list<unsigned int> newCell1DIndices;
            for (const unsigned int &e : segment.Edge2DIds)
            {
                if (!mesh2D.Cell1DHasUpdatedCell1Ds(e))
                    continue;

                list<unsigned int> updatedCell1Ds;
                mesh2D.Cell1DUpdatedCell1Ds(e, updatedCell1Ds);
                for (const unsigned int &ne : updatedCell1Ds)
                    newCell1DIndices.push_back(ne);
            }

            for (const unsigned int &ne : newCell1DIndices)
                segment.Edge2DIds.insert(ne);
        }

        {
            list<unsigned int> newCell2DIndices;
            for (const unsigned int &c : segment.Cell2DIds)
            {
                if (!mesh2D.Cell2DHasUpdatedCell2Ds(c))
                    continue;

                list<unsigned int> updatedCell2Ds;
                mesh2D.Cell2DUpdatedCell2Ds(c, updatedCell2Ds);
                for (const unsigned int &nc : updatedCell2Ds)
                    newCell2DIndices.push_back(nc);
            }

            for (const unsigned int &nc : newCell2DIndices)
                segment.Cell2DIds.insert(nc);
        }
    }
}
// ***************************************************************************
void ConformerMeshSegment::UpdateWithActiveMesh2D(const MeshUtilities::ExtractActiveMeshData &activeMesh2DData,
                                                  ConformMesh &conformedMesh) const
{
    // update points
    for (map<double, ConformerMeshSegment::ConformMesh::ConformMeshPoint>::iterator it = conformedMesh.Points.begin();
         it != conformedMesh.Points.end();
         it++)
    {
        ConformerMeshSegment::ConformMesh::ConformMeshPoint &point = it->second;

        set<unsigned int> newVertex2DIds;
        for (const auto &vertexId : point.Vertex2DIds)
        {
            if (activeMesh2DData.OldCell0DToNewCell0D.find(vertexId) != activeMesh2DData.OldCell0DToNewCell0D.end())
                newVertex2DIds.insert(activeMesh2DData.OldCell0DToNewCell0D.at(vertexId));
        }
        point.Vertex2DIds.clear();
        point.Vertex2DIds = newVertex2DIds;

        set<unsigned int> newEdge2DIds;
        for (const auto &edgeId : point.Edge2DIds)
        {
            if (activeMesh2DData.OldCell1DToNewCell1D.find(edgeId) != activeMesh2DData.OldCell1DToNewCell1D.end())
                newEdge2DIds.insert(activeMesh2DData.OldCell1DToNewCell1D.at(edgeId));
        }
        point.Edge2DIds.clear();
        point.Edge2DIds = newEdge2DIds;

        set<unsigned int> newCell2DIds;
        for (const auto &cell2DIndex : point.Cell2DIds)
        {
            if (activeMesh2DData.OldCell2DToNewCell2D.find(cell2DIndex) != activeMesh2DData.OldCell2DToNewCell2D.end())
                newCell2DIds.insert(activeMesh2DData.OldCell2DToNewCell2D.at(cell2DIndex));
        }
        point.Cell2DIds.clear();
        point.Cell2DIds = newCell2DIds;
    }

    // update segments
    for (ConformerMeshSegment::ConformMesh::ConformMeshSegment &segment : conformedMesh.Segments)
    {
        set<unsigned int> newEdge2DIds;
        for (const auto &edgeId : segment.Edge2DIds)
        {
            if (activeMesh2DData.OldCell1DToNewCell1D.find(edgeId) != activeMesh2DData.OldCell1DToNewCell1D.end())
                newEdge2DIds.insert(activeMesh2DData.OldCell1DToNewCell1D.at(edgeId));
        }
        segment.Edge2DIds.clear();
        segment.Edge2DIds = newEdge2DIds;

        set<unsigned int> newCell2DIds;
        for (const auto &cell2DIndex : segment.Cell2DIds)
        {
            if (activeMesh2DData.OldCell2DToNewCell2D.find(cell2DIndex) != activeMesh2DData.OldCell2DToNewCell2D.end())
                newCell2DIds.insert(activeMesh2DData.OldCell2DToNewCell2D.at(cell2DIndex));
        }
        segment.Cell2DIds.clear();
        segment.Cell2DIds = newCell2DIds;
    }
}
// ***************************************************************************
vector<double> ConformerMeshSegment::AddMissingMesh2DCell0Ds(const Vector3d &segmentOrigin,
                                                             const Vector3d &segmentTangent,
                                                             const double &segmentSquaredLength,
                                                             const Gedim::IMeshDAO &mesh2D,
                                                             ConformMesh &conformedMesh) const
{
    list<double> newCoordinatesCreated;

    for (ConformerMeshSegment::ConformMesh::ConformMeshSegment &segment : conformedMesh.Segments)
    {
        if (segment.Edge2DIds.size() == 1)
            continue;

        Output::Assert(segment.Edge2DIds.size() == 2);
        // add missing Cell0D point
        const unsigned int cell1DOneIndex = *segment.Edge2DIds.begin();
        const unsigned int cell1DTwoIndex = *segment.Edge2DIds.rbegin();
        const unsigned int meshCell0DIndex = mesh2D.Cell1DOrigin(cell1DOneIndex) == mesh2D.Cell1DOrigin(cell1DTwoIndex) ||
                                                     mesh2D.Cell1DOrigin(cell1DOneIndex) == mesh2D.Cell1DEnd(cell1DTwoIndex)
                                                 ? mesh2D.Cell1DOrigin(cell1DOneIndex)
                                                 : mesh2D.Cell1DEnd(cell1DOneIndex);
        const Eigen::Vector3d cell0DCoordinate = mesh2D.Cell0DCoordinates(meshCell0DIndex);

        const double curvilinearCoordinate =
            _geometryUtilities.PointLineCurvilinearCoordinate(cell0DCoordinate, segmentOrigin, segmentTangent, segmentSquaredLength);

        bool found;
        ConformerMeshSegment::ConformMesh::ConformMeshPoint &conformPoint =
            InsertNewIntersection(curvilinearCoordinate, conformedMesh, found);

        Output::Assert(!found);

        conformPoint.Type = ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original;
        conformPoint.Vertex2DIds.insert(meshCell0DIndex);
        conformPoint.Edge2DIds.insert(cell1DOneIndex);
        conformPoint.Edge2DIds.insert(cell1DTwoIndex);

        if (mesh2D.Cell1DHasNeighbourCell2D(cell1DOneIndex, 0))
            conformPoint.Cell2DIds.insert(mesh2D.Cell1DNeighbourCell2D(cell1DOneIndex, 0));
        if (mesh2D.Cell1DHasNeighbourCell2D(cell1DOneIndex, 1))
            conformPoint.Cell2DIds.insert(mesh2D.Cell1DNeighbourCell2D(cell1DOneIndex, 1));
        if (mesh2D.Cell1DHasNeighbourCell2D(cell1DTwoIndex, 0))
            conformPoint.Cell2DIds.insert(mesh2D.Cell1DNeighbourCell2D(cell1DTwoIndex, 0));
        if (mesh2D.Cell1DHasNeighbourCell2D(cell1DTwoIndex, 1))
            conformPoint.Cell2DIds.insert(mesh2D.Cell1DNeighbourCell2D(cell1DTwoIndex, 1));

        newCoordinatesCreated.push_back(curvilinearCoordinate);
    }

    if (newCoordinatesCreated.size() > 0)
    {
        // recreate conform mesh segments
        conformedMesh.Segments.clear();
        CreateConformSegments(conformedMesh);
    }

    return vector<double>(newCoordinatesCreated.begin(), newCoordinatesCreated.end());
}
// ***************************************************************************
} // namespace Gedim
