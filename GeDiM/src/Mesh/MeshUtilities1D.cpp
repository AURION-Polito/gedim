#include "MeshUtilities.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::FillMesh1D(const GeometryUtilities&,
                                 const Vector3d& segmentOrigin,
                                 const Vector3d& segmentTangent,
                                 const vector<double>& coordinates,
                                 IMeshDAO& mesh) const
  {
    if (coordinates.size() == 0)
      return;

    mesh.InitializeDimension(1);

    const unsigned int& numCell0Ds = coordinates.size();
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int c = 0; c < numCell0Ds; c++)
    {
      mesh.Cell0DSetState(c, true);
      mesh.Cell0DInsertCoordinates(c, segmentOrigin + coordinates[c] * segmentTangent);
    }
    mesh.Cell0DSetMarker(0, 1);
    mesh.Cell0DSetMarker(numCell0Ds - 1, 2);

    const unsigned int numCell1Ds = numCell0Ds - 1;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                e,
                                e + 1);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, 0);
    }
  }
  // ***************************************************************************
  void MeshUtilities::Mesh1DFromSegment(const GeometryUtilities& geometryUtilities,
                                        const Eigen::MatrixXd& segmentVertices,
                                        const vector<unsigned int> vertexMarkers,
                                        IMeshDAO& mesh) const
  {
    FillMesh1D(geometryUtilities,
               segmentVertices.col(0),
               segmentVertices.col(1),
               { 0.0, 1.0 },
               mesh);

    mesh.Cell0DSetMarker(0, vertexMarkers[0]);
    mesh.Cell0DSetMarker(mesh.Cell0DTotalNumber() - 1, vertexMarkers[1]);
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData1D MeshUtilities::FillMesh1DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData1D result;

    result.Cell1DsVertices.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsTangents.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsSquaredLengths.resize(convexMesh.Cell1DTotalNumber());
    result.Cell1DsCentroids.resize(convexMesh.Cell1DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell1DTotalNumber(); c++)
    {
      // Extract original cell1D geometric information

      result.Cell1DsVertices[c] = convexMesh.Cell1DCoordinates(c);
      result.Cell1DsTangents[c] = geometryUtilities.SegmentTangent(result.Cell1DsVertices[c].col(0),
                                                                   result.Cell1DsVertices[c].col(1));
      result.Cell1DsLengths[c] = geometryUtilities.SegmentLength(result.Cell1DsVertices[c].col(0),
                                                                 result.Cell1DsVertices[c].col(1));
      result.Cell1DsSquaredLengths[c] = result.Cell1DsLengths[c] * result.Cell1DsLengths[c];
      result.Cell1DsCentroids[c] = geometryUtilities.SimplexBarycenter(result.Cell1DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  std::vector<unsigned int> MeshUtilities::SplitCell1D(const unsigned int& cell1DIndex,
                                                       const Eigen::MatrixXi subCell1Ds,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numSubCells = subCell1Ds.cols();
    unsigned int newCell1DsStartingIndex = mesh.Cell1DAppend(numSubCells);

    mesh.Cell1DSetState(cell1DIndex, false);

    vector<unsigned int> newCell1DsIndex(numSubCells);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell1DsIndex[c] = newCell1DsStartingIndex + c;

      const unsigned int& newCell1DIndex = newCell1DsIndex[c];

      mesh.Cell1DInsertExtremes(newCell1DIndex,
                                subCell1Ds(0, c),
                                subCell1Ds(1, c));

      mesh.Cell1DSetMarker(newCell1DIndex, mesh.Cell1DMarker(cell1DIndex));
      mesh.Cell1DSetState(newCell1DIndex, true);

      mesh.Cell1DInsertUpdatedCell1D(cell1DIndex,
                                     newCell1DIndex);

      const unsigned int numCell1DNumberNeighbourCell2D = mesh.Cell1DNumberNeighbourCell2D(cell1DIndex);
      if (numCell1DNumberNeighbourCell2D == 0)
        continue;

      mesh.Cell1DInitializeNeighbourCell2Ds(newCell1DIndex,
                                            numCell1DNumberNeighbourCell2D);

      for (unsigned int n = 0; n < numCell1DNumberNeighbourCell2D; n++)
      {
        if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
          continue;

        mesh.Cell1DInsertNeighbourCell2D(newCell1DIndex,
                                         n,
                                         mesh.Cell1DNeighbourCell2D(cell1DIndex,
                                                                    n));
      }
    }

    return newCell1DsIndex;
  }
  // ***************************************************************************
}
