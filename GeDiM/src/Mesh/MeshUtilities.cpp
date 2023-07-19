#include "MeshUtilities.hpp"

#include "TriangleInterface.hpp"
#include "TetgenInterface.hpp"
#include "VTKUtilities.hpp"
#include "MapTetrahedron.hpp"
#include "OpenVolumeMeshInterface.hpp"

#include "numeric"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  MeshUtilities::MeshUtilities()
  {
  }
  MeshUtilities::~MeshUtilities()
  {
  }
  // ***************************************************************************
  void MeshUtilities::ExtractActiveMesh(IMeshDAO& mesh,
                                        ExtractActiveMeshData& extractionData) const
  {
    // remove inactive Cell0Ds
    unsigned int numNewCell0Ds = 0;
    list<unsigned int> cell0DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
      if (!mesh.Cell0DIsActive(c))
      {
        cell0DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell0DToOldCell0D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell0Ds, c));
      extractionData.OldCell0DToNewCell0D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell0Ds));
      numNewCell0Ds++;
    }

    unsigned int removedCell0Ds = 0;
    for (const unsigned int& c : cell0DIdToRemove)
    {
      mesh.Cell0DRemove(c - removedCell0Ds);
      removedCell0Ds++;
    }

    // remove inactive Cell1D
    unsigned int numNewCell1Ds = 0;
    list<unsigned int> cell1DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
      if (!mesh.Cell1DIsActive(c))
      {
        cell1DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell1DToOldCell1D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell1Ds, c));
      extractionData.OldCell1DToNewCell1D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell1Ds));
      numNewCell1Ds++;
    }

    unsigned int removedCell1Ds = 0;
    for (const unsigned int& c : cell1DIdToRemove)
    {
      mesh.Cell1DRemove(c - removedCell1Ds);
      removedCell1Ds++;
    }

    // remove inactive Cell2Ds
    unsigned int numNewCell2Ds = 0;
    list<unsigned int> cell2DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
      {
        cell2DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell2DToOldCell2D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell2Ds, c));
      extractionData.OldCell2DToNewCell2D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell2Ds));
      numNewCell2Ds++;
    }

    unsigned int removedCell2Ds = 0;
    for (const unsigned int& c : cell2DIdToRemove)
    {
      mesh.Cell2DRemove(c - removedCell2Ds);
      removedCell2Ds++;
    }

    // remove inactive Cell3Ds
    unsigned int numNewCell3Ds = 0;
    list<unsigned int> cell3DIdToRemove;
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
      if (!mesh.Cell3DIsActive(c))
      {
        cell3DIdToRemove.push_back(c);
        continue;
      }

      extractionData.NewCell3DToOldCell3D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell3Ds, c));
      extractionData.OldCell3DToNewCell3D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell3Ds));
      numNewCell3Ds++;
    }

    unsigned int removedCell3Ds = 0;
    for (const unsigned int& c : cell3DIdToRemove)
    {
      mesh.Cell3DRemove(c - removedCell3Ds);
      removedCell3Ds++;
    }

    mesh.Compress();
  }
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
  void MeshUtilities::FillMesh2D(const MatrixXd& cell0Ds,
                                 const MatrixXi& cell1Ds,
                                 const vector<MatrixXi>& cell2Ds,
                                 IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    // Create Cell0Ds
    Output::Assert(cell0Ds.rows() == 3);
    const unsigned int& numCell0Ds = cell0Ds.cols();
    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell0DsInsertCoordinates(cell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
      mesh.Cell0DSetState(v, true);

    // Create Cell1Ds
    Output::Assert(cell1Ds.rows() == 2);
    unsigned int numCell1Ds = cell1Ds.cols();
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell1DsInsertExtremes(cell1Ds);
    for (int e = 0; e < cell1Ds.cols(); e++)
      mesh.Cell1DSetState(e, true);

    // Create Cell2Ds
    const unsigned int& numCell2Ds = cell2Ds.size();
    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const MatrixXi& polygon = cell2Ds[f];
      Output::Assert(polygon.rows() == 2);
      const unsigned int& numVertices = polygon.cols();

      mesh.Cell2DInitializeVertices(f, numVertices);
      mesh.Cell2DInitializeEdges(f, numVertices);

      for (unsigned int v = 0; v < numVertices; v++)
        mesh.Cell2DInsertVertex(f, v, polygon(0, v));
      for (unsigned int e = 0; e < numVertices; e++)
        mesh.Cell2DInsertEdge(f, e, polygon(1, e));

      mesh.Cell2DSetState(f, true);
    }
  }
  // ***************************************************************************
  MeshUtilities::ComputeMesh2DCell1DsResult MeshUtilities::ComputeMesh2DCell1Ds(const Eigen::MatrixXd& cell0Ds,
                                                                                const std::vector<Eigen::VectorXi>& cell2Ds) const
  {
    ComputeMesh2DCell1DsResult result;

    const unsigned int numVertices = cell0Ds.cols();
    const unsigned int numCell2Ds = cell2Ds.size();

    Eigen::SparseMatrix<unsigned int> edges;
    edges.resize(numVertices,
                 numVertices);

    std::list<Eigen::Triplet<unsigned int>> triplets;
    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];
        triplets.push_back(Eigen::Triplet<unsigned int>(origin, end, 1));
        triplets.push_back(Eigen::Triplet<unsigned int>(end, origin, 1));
      }
    }

    edges.setFromTriplets(triplets.begin(), triplets.end());
    edges.makeCompressed();

    unsigned int numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          it.valueRef() = 1 + numEdges++;
      }
    }

    result.Cell1Ds.resize(2, numEdges);

    numEdges = 0;
    for (int k = 0; k < edges.outerSize(); k++)
    {
      for (SparseMatrix<unsigned int>::InnerIterator it(edges, k); it; ++it)
      {
        if (it.row() < it.col())
          result.Cell1Ds.col(numEdges++)<< it.row(), it.col();
      }
    }

    result.Cell2Ds.resize(numCell2Ds);

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
      Eigen::MatrixXi& cell2D = result.Cell2Ds.at(c);
      const Eigen::VectorXi& cell2DVertices = cell2Ds.at(c);
      const unsigned int& numCell2DVertices = cell2DVertices.size();

      cell2D.resize(2, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
      {
        const unsigned int origin = cell2DVertices[v];
        const unsigned int end = cell2DVertices[(v + 1) % numCell2DVertices];

        cell2D(0, v) = origin;
        cell2D(1, v) = (origin < end) ? edges.coeff(origin, end) - 1 : edges.coeff(end, origin) - 1;
      }
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh2D(const CheckMesh2DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 2);

    // check Cell0D are 2D
    if (configuration.Cell0D_CheckCoordinates2D)
      Output::Assert(geometryUtilities.PointsAre2D(convexMesh.Cell0DsCoordinates()));

    // check Cell0D duplications
    if (configuration.Cell0D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell0DTotalNumber(); p1++)
      {
        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell0DTotalNumber(); p2++)
        {
          Output::Assert(!geometryUtilities.PointsAreCoincident(convexMesh.Cell0DCoordinates(p1),
                                                                convexMesh.Cell0DCoordinates(p2)));
        }
      }
    }

    if (configuration.Cell1D_CheckDuplications)
    {
      for (unsigned int e1 = 0; e1 < convexMesh.Cell1DTotalNumber(); e1++)
      {
        Output::Assert(convexMesh.Cell1DExists(convexMesh.Cell1DOrigin(e1),
                                               convexMesh.Cell1DEnd(e1)));
        Output::Assert(!convexMesh.Cell1DExists(convexMesh.Cell1DEnd(e1),
                                                convexMesh.Cell1DOrigin(e1)));

        for (unsigned int e2 = e1 + 1; e2 < convexMesh.Cell1DTotalNumber(); e2++)
        {
          Output::Assert(!(convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DEnd(e2)));
          Output::Assert(!(convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DEnd(e2)));
        }
      }
    }

    if (configuration.Cell1D_CheckNeighbours)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(convexMesh.Cell1DNumberNeighbourCell2D(e) == 2);

        if (convexMesh.Cell1DHasNeighbourCell2D(e, 0))
        {
          const unsigned int cell2DRight = convexMesh.Cell1DNeighbourCell2D(e, 0);
          const vector<unsigned int> cell2DEdges = convexMesh.Cell2DEdges(cell2DRight);

          // check edge orientation
          vector<unsigned int>::const_iterator it = std::find(cell2DEdges.begin(), cell2DEdges.end(), e);
          Output::Assert(it != cell2DEdges.end());

          const unsigned int cell2DEdgeIndex = std::distance(cell2DEdges.begin(), it);
          const unsigned int edgeOrigin = convexMesh.Cell2DVertex(cell2DRight,
                                                                  (cell2DEdgeIndex + 1) % cell2DEdges.size());
          const unsigned int edgeEnd = convexMesh.Cell2DVertex(cell2DRight,
                                                               cell2DEdgeIndex);

          Output::Assert(convexMesh.Cell1DExists(edgeOrigin,
                                                 edgeEnd) &&
                         convexMesh.Cell1DByExtremes(edgeOrigin,
                                                     edgeEnd) == e);
        }

        if (convexMesh.Cell1DHasNeighbourCell2D(e, 1))
        {
          const unsigned int cell2DLeft = convexMesh.Cell1DNeighbourCell2D(e, 1);
          const vector<unsigned int> cell2DEdges = convexMesh.Cell2DEdges(cell2DLeft);

          // check edge orientation
          vector<unsigned int>::const_iterator it = std::find(cell2DEdges.begin(), cell2DEdges.end(), e);
          Output::Assert(it != cell2DEdges.end());

          const unsigned int cell2DEdgeIndex = std::distance(cell2DEdges.begin(), it);
          const unsigned int edgeOrigin = convexMesh.Cell2DVertex(cell2DLeft,
                                                                  cell2DEdgeIndex);
          const unsigned int edgeEnd = convexMesh.Cell2DVertex(cell2DLeft,
                                                               (cell2DEdgeIndex + 1) % cell2DEdges.size());

          Output::Assert(convexMesh.Cell1DExists(edgeOrigin,
                                                 edgeEnd) &&
                         convexMesh.Cell1DByExtremes(edgeOrigin,
                                                     edgeEnd) == e);
        }
      }
    }

    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(geometryUtilities.IsValue1DPositive(
                         geometryUtilities.SegmentLength(convexMesh.Cell1DOriginCoordinates(e),
                                                         convexMesh.Cell1DEndCoordinates(e))));
      }
    }

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        for (unsigned int v = 0; v < convexMesh.Cell2DNumberVertices(p); v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % convexMesh.Cell2DNumberVertices(p));
          Output::Assert(convexMesh.Cell1DExists(eO, eE) || convexMesh.Cell1DExists(eE, eO));
          const unsigned int edgeFromVertices = convexMesh.Cell1DExists(eO, eE) ? convexMesh.Cell1DByExtremes(eO, eE) :
                                                                                  convexMesh.Cell1DByExtremes(eE, eO);
          Output::Assert(convexMesh.Cell2DEdge(p, v) == edgeFromVertices);
        }
      }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell2DTotalNumber(); p1++)
      {
        vector<unsigned int> cell2D1Vertices = convexMesh.Cell2DVertices(p1);
        sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
        vector<unsigned int> cell2D1Edges = convexMesh.Cell2DEdges(p1);
        sort(cell2D1Edges.begin(), cell2D1Edges.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell2DTotalNumber(); p2++)
        {
          vector<unsigned int> cell2D2Vertices = convexMesh.Cell2DVertices(p2);
          sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
          vector<unsigned int> cell2D2Edges = convexMesh.Cell2DEdges(p2);
          sort(cell2D2Edges.begin(), cell2D2Edges.end());

          Output::Assert(cell2D1Vertices.size() != cell2D2Vertices.size() || !equal(cell2D1Vertices.begin(),
                                                                                    cell2D1Vertices.end(),
                                                                                    cell2D2Vertices.begin()));
          Output::Assert(cell2D1Edges.size() != cell2D2Edges.size() || !equal(cell2D1Edges.begin(),
                                                                              cell2D1Edges.end(),
                                                                              cell2D2Edges.begin()));
        }
      }
    }

    if (configuration.Cell2D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const Eigen::MatrixXd cell2DVertices = convexMesh.Cell2DVerticesCoordinates(p);
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);
        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(convexCell2DUnalignedVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(convexCell2DUnalignedVertices,
                                                                                   convexHull);

        Output::Assert(geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices,
                                                         convexHullVertices));
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        Output::Assert(geometryUtilities.IsValue2DPositive(
                         geometryUtilities.PolygonArea(convexMesh.Cell2DVerticesCoordinates(p))));
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CheckMesh3D(const CheckMesh3DConfiguration& configuration,
                                  const GeometryUtilities& geometryUtilities,
                                  const IMeshDAO& convexMesh) const
  {
    Output::Assert(convexMesh.Dimension() == 3);

    // check Cell0D duplications
    if (configuration.Cell0D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell0DTotalNumber(); p1++)
      {
        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell0DTotalNumber(); p2++)
        {
          Output::Assert(!geometryUtilities.PointsAreCoincident(convexMesh.Cell0DCoordinates(p1),
                                                                convexMesh.Cell0DCoordinates(p2)));
        }
      }
    }

    if (configuration.Cell1D_CheckDuplications)
    {
      for (unsigned int e1 = 0; e1 < convexMesh.Cell1DTotalNumber(); e1++)
      {
        Output::Assert(convexMesh.Cell1DExists(convexMesh.Cell1DOrigin(e1),
                                               convexMesh.Cell1DEnd(e1)));
        Output::Assert(!convexMesh.Cell1DExists(convexMesh.Cell1DEnd(e1),
                                                convexMesh.Cell1DOrigin(e1)));

        for (unsigned int e2 = e1 + 1; e2 < convexMesh.Cell1DTotalNumber(); e2++)
        {
          Output::Assert(!(convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DEnd(e2)));
          Output::Assert(!(convexMesh.Cell1DEnd(e1) == convexMesh.Cell1DOrigin(e2) &&
                           convexMesh.Cell1DOrigin(e1) == convexMesh.Cell1DEnd(e2)));
        }
      }
    }

    if (configuration.Cell1D_CheckMeasure)
    {
      for (unsigned int e = 0; e < convexMesh.Cell1DTotalNumber(); e++)
      {
        Output::Assert(geometryUtilities.IsValue1DPositive(
                         geometryUtilities.SegmentLength(convexMesh.Cell1DOriginCoordinates(e),
                                                         convexMesh.Cell1DEndCoordinates(e))));
      }
    }

    if (configuration.Cell2D_CheckEdges)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        for (unsigned int v = 0; v < convexMesh.Cell2DNumberVertices(p); v++)
        {
          const unsigned int eO = convexMesh.Cell2DVertex(p, v);
          const unsigned int eE = convexMesh.Cell2DVertex(p, (v + 1) % convexMesh.Cell2DNumberVertices(p));
          Output::Assert(convexMesh.Cell1DExists(eO, eE) || convexMesh.Cell1DExists(eE, eO));
          const unsigned int edgeFromVertices = convexMesh.Cell1DExists(eO, eE) ? convexMesh.Cell1DByExtremes(eO, eE) :
                                                                                  convexMesh.Cell1DByExtremes(eE, eO);
          Output::Assert(convexMesh.Cell2DEdge(p, v) == edgeFromVertices);
        }
      }
    }

    if (configuration.Cell2D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell2DTotalNumber(); p1++)
      {
        vector<unsigned int> cell2D1Vertices = convexMesh.Cell2DVertices(p1);
        sort(cell2D1Vertices.begin(), cell2D1Vertices.end());
        vector<unsigned int> cell2D1Edges = convexMesh.Cell2DEdges(p1);
        sort(cell2D1Edges.begin(), cell2D1Edges.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell2DTotalNumber(); p2++)
        {
          vector<unsigned int> cell2D2Vertices = convexMesh.Cell2DVertices(p2);
          sort(cell2D2Vertices.begin(), cell2D2Vertices.end());
          vector<unsigned int> cell2D2Edges = convexMesh.Cell2DEdges(p2);
          sort(cell2D2Edges.begin(), cell2D2Edges.end());

          Output::Assert(cell2D1Vertices.size() != cell2D2Vertices.size() || !equal(cell2D1Vertices.begin(),
                                                                                    cell2D1Vertices.end(),
                                                                                    cell2D2Vertices.begin()));
          Output::Assert(cell2D1Edges.size() != cell2D2Edges.size() || !equal(cell2D1Edges.begin(),
                                                                              cell2D1Edges.end(),
                                                                              cell2D2Edges.begin()));
        }
      }
    }

    if (configuration.Cell2D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const Eigen::MatrixXd cell2DVertices3D = convexMesh.Cell2DVerticesCoordinates(p);
        const Eigen::Vector3d cell2DNormal = geometryUtilities.PolygonNormal(cell2DVertices3D);
        const Eigen::Vector3d cell2DTranslation = geometryUtilities.PolygonTranslation(cell2DVertices3D);
        const Eigen::Matrix3d cell2DRotationMatrix = geometryUtilities.PolygonRotationMatrix(cell2DVertices3D,
                                                                                             cell2DNormal,
                                                                                             cell2DTranslation);
        const Eigen::MatrixXd cell2DVertices2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D,
                                                                                          cell2DRotationMatrix.transpose(),
                                                                                          cell2DTranslation);
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(cell2DVertices2D);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(cell2DVertices2D,
                                                                                              convexCell2DUnalignedVerticesFilter);
        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(convexCell2DUnalignedVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(convexCell2DUnalignedVertices,
                                                                                   convexHull);

        Output::Assert(geometryUtilities.PolygonIsConvex(convexCell2DUnalignedVertices,
                                                         convexHullVertices));
      }
    }

    if (configuration.Cell2D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell2DTotalNumber(); p++)
      {
        const Eigen::MatrixXd cell2DVertices3D = convexMesh.Cell2DVerticesCoordinates(p);
        const Eigen::Vector3d cell2DNormal = geometryUtilities.PolygonNormal(cell2DVertices3D);
        const Eigen::Vector3d cell2DTranslation = geometryUtilities.PolygonTranslation(cell2DVertices3D);
        const Eigen::Matrix3d cell2DRotationMatrix = geometryUtilities.PolygonRotationMatrix(cell2DVertices3D,
                                                                                             cell2DNormal,
                                                                                             cell2DTranslation);
        const Eigen::MatrixXd cell2DVertices2D = geometryUtilities.RotatePointsFrom3DTo2D(cell2DVertices3D,
                                                                                          cell2DRotationMatrix.transpose(),
                                                                                          cell2DTranslation);

        Output::Assert(geometryUtilities.IsValue2DPositive(
                         geometryUtilities.PolygonArea(cell2DVertices2D)));
      }
    }

    if (configuration.Cell3D_CheckDuplications)
    {
      for (unsigned int p1 = 0; p1 < convexMesh.Cell3DTotalNumber(); p1++)
      {
        vector<unsigned int> cell3D1Vertices = convexMesh.Cell3DVertices(p1);
        sort(cell3D1Vertices.begin(), cell3D1Vertices.end());
        vector<unsigned int> cell3D1Edges = convexMesh.Cell3DEdges(p1);
        sort(cell3D1Edges.begin(), cell3D1Edges.end());
        vector<unsigned int> cell3D1Faces = convexMesh.Cell3DFaces(p1);
        sort(cell3D1Faces.begin(), cell3D1Faces.end());

        for (unsigned int p2 = p1 + 1; p2 < convexMesh.Cell3DTotalNumber(); p2++)
        {
          vector<unsigned int> cell3D2Vertices = convexMesh.Cell3DVertices(p2);
          sort(cell3D2Vertices.begin(), cell3D2Vertices.end());
          vector<unsigned int> cell3D2Edges = convexMesh.Cell3DEdges(p2);
          sort(cell3D2Edges.begin(), cell3D2Edges.end());
          vector<unsigned int> cell3D2Faces = convexMesh.Cell3DFaces(p2);
          sort(cell3D2Faces.begin(), cell3D2Faces.end());

          Output::Assert(cell3D1Vertices.size() != cell3D2Vertices.size() || !equal(cell3D1Vertices.begin(),
                                                                                    cell3D1Vertices.end(),
                                                                                    cell3D2Vertices.begin()));
          Output::Assert(cell3D1Edges.size() != cell3D2Edges.size() || !equal(cell3D1Edges.begin(),
                                                                              cell3D1Edges.end(),
                                                                              cell3D2Edges.begin()));
          Output::Assert(cell3D1Faces.size() != cell3D2Faces.size() || !equal(cell3D1Faces.begin(),
                                                                              cell3D1Faces.end(),
                                                                              cell3D2Faces.begin()));
        }
      }
    }


    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);
      }
    }

    if (configuration.Cell3D_CheckConvexity)
    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<Eigen::Vector3d> polyhedronFaceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);

        const GeometryUtilities::PointPolyhedronPositionResult polyhedronBarycenterPosition = geometryUtilities.PointPolyhedronPosition(polyhedronBarycenter,
                                                                                                                                        polyhedron.Faces,
                                                                                                                                        polyhedronFace3DVertices,
                                                                                                                                        polyhedronFace2DVertices,
                                                                                                                                        polyhedronFaceNormals,
                                                                                                                                        polyhedronFaceNormalDirections,
                                                                                                                                        polyhedronFaceTranslations,
                                                                                                                                        polyhedronFaceRotationMatrices);

        Output::Assert(polyhedronBarycenterPosition.Type ==
                       GeometryUtilities::PointPolyhedronPositionResult::Types::Inside);

        Output::Assert(geometryUtilities.PolyhedronIsConvex(polyhedronFace3DVertices,
                                                            polyhedronFace2DVertices,
                                                            polyhedronFaceBarycenters,
                                                            polyhedronFaceNormals,
                                                            polyhedronFaceNormalDirections,
                                                            polyhedronBarycenter));
      }
    }

    if (configuration.Cell3D_CheckMeasure)
    {
      for (unsigned int p = 0; p < convexMesh.Cell3DTotalNumber(); p++)
      {
        GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                          p);

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                          polyhedron.Faces);
        const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                      polyhedronFace3DVertices);

        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
        const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                             polyhedronBarycenter,
                                                                                                             polyhedronFaceNormals);
        const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
        const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                        polyhedronFaceNormals,
                                                                                                                        polyhedronFaceTranslations);

        const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                                 polyhedronFaceTranslations,
                                                                                                                 polyhedronFaceRotationMatrices);

        const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(polyhedronFace2DVertices,
                                                                                                                                                               polyhedronFaceTriangulations);

        Output::Assert(geometryUtilities.IsValue3DPositive(
                         geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                            polyhedronFaceNormals,
                                                            polyhedronFaceNormalDirections,
                                                            polyhedronFaceTranslations,
                                                            polyhedronFaceRotationMatrices)));
      }
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
  void MeshUtilities::Mesh2DFromPolygon(const Eigen::MatrixXd& polygonVertices,
                                        const vector<unsigned int> vertexMarkers,
                                        const vector<unsigned int> edgeMarkers,
                                        IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(2);

    Output::Assert(polygonVertices.rows() == 3 && polygonVertices.cols() > 2);
    const unsigned int numPolygonVertices = polygonVertices.cols();
    Output::Assert(vertexMarkers.size() == numPolygonVertices);
    Output::Assert(edgeMarkers.size() == numPolygonVertices);

    // Create Cell0Ds
    const unsigned int& numCell0Ds = numPolygonVertices;
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polygonVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numPolygonVertices;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                e,
                                (e + 1) % numPolygonVertices);
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, edgeMarkers[e]);
    }

    // Create Cell2Ds
    const unsigned int& numCell2Ds = 1;
    mesh.Cell2DsInitialize(numCell2Ds);

    mesh.Cell2DInitializeVertices(0, numPolygonVertices);
    mesh.Cell2DInitializeEdges(0, numPolygonVertices);

    for (unsigned int v = 0; v < numPolygonVertices; v++)
      mesh.Cell2DInsertVertex(0, v, v);
    for (unsigned int e = 0; e < numPolygonVertices; e++)
      mesh.Cell2DInsertEdge(0, e, e);

    mesh.Cell2DSetState(0, true);
    mesh.Cell2DSetMarker(0, 0);

    // Create Cell1D neighbours
    for (unsigned int e = 0; e < numPolygonVertices; e++)
    {
      mesh.Cell1DInitializeNeighbourCell2Ds(e, 2);
      mesh.Cell1DInsertNeighbourCell2D(e, 1, 0);
    }
  }
  // ***************************************************************************
  void MeshUtilities::Mesh3DFromPolyhedron(const Eigen::MatrixXd& polyhedronVertices,
                                           const Eigen::MatrixXi& polyhedronEdges,
                                           const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                           const vector<unsigned int> vertexMarkers,
                                           const vector<unsigned int> edgeMarkers,
                                           const vector<unsigned int> faceMarkers,
                                           IMeshDAO& mesh) const
  {
    mesh.InitializeDimension(3);

    Output::Assert(polyhedronVertices.rows() == 3 && polyhedronVertices.cols() > 3);
    const unsigned int numVertices = polyhedronVertices.cols();
    const unsigned int numEdges = polyhedronEdges.cols();
    const unsigned int numFaces = polyhedronFaces.size();

    Output::Assert(vertexMarkers.size() == numVertices);
    Output::Assert(edgeMarkers.size() == numEdges);
    Output::Assert(faceMarkers.size() == numFaces);

    // Create Cell0Ds
    const unsigned int& numCell0Ds = numVertices;
    mesh.Cell0DsInitialize(numCell0Ds);
    for (unsigned int v = 0; v < numCell0Ds; v++)
    {
      mesh.Cell0DSetState(v, true);
      mesh.Cell0DInsertCoordinates(v, polyhedronVertices.col(v));
      mesh.Cell0DSetMarker(v, vertexMarkers[v]);
    }

    // Create Cell1Ds
    unsigned int numCell1Ds = numEdges;
    mesh.Cell1DsInitialize(numCell1Ds);
    for (unsigned int e = 0; e < numCell1Ds; e++)
    {
      mesh.Cell1DInsertExtremes(e,
                                polyhedronEdges(0, e),
                                polyhedronEdges(1, e));
      mesh.Cell1DSetState(e, true);
      mesh.Cell1DSetMarker(e, edgeMarkers[e]);
    }


    // Create Cell2Ds
    const unsigned int& numCell2Ds = numFaces;
    mesh.Cell2DsInitialize(numCell2Ds);
    for (unsigned int f = 0; f < numCell2Ds; f++)
    {
      const unsigned int numCell2DVertices = polyhedronFaces.at(f).cols();
      mesh.Cell2DInitializeVertices(f, numCell2DVertices);
      mesh.Cell2DInitializeEdges(f, numCell2DVertices);

      for (unsigned int v = 0; v < numCell2DVertices; v++)
        mesh.Cell2DInsertVertex(f, v, polyhedronFaces.at(f)(0, v));
      for (unsigned int e = 0; e < numCell2DVertices; e++)
        mesh.Cell2DInsertEdge(f, e, polyhedronFaces.at(f)(1, e));

      mesh.Cell2DSetState(f, true);
      mesh.Cell2DSetMarker(f, faceMarkers[f]);
    }

    // Create Cell3Ds
    const unsigned int& numCell3Ds = 1;
    mesh.Cell3DsInitialize(numCell3Ds);

    mesh.Cell3DInitializeVertices(0, numVertices);
    mesh.Cell3DInitializeEdges(0, numEdges);
    mesh.Cell3DInitializeFaces(0, numFaces);

    for (unsigned int v = 0; v < numVertices; v++)
      mesh.Cell3DInsertVertex(0, v, v);
    for (unsigned int e = 0; e < numEdges; e++)
      mesh.Cell3DInsertEdge(0, e, e);
    for (unsigned int f = 0; f < numFaces; f++)
      mesh.Cell3DInsertFace(0, f, f);

    mesh.Cell3DSetState(0, true);
    mesh.Cell3DSetMarker(0, 0);
  }
  // ***************************************************************************
  void MeshUtilities::SetMeshMarkersOnPlane(const GeometryUtilities& geometryUtilities,
                                            const Eigen::Vector3d& planeNormal,
                                            const Eigen::Vector3d& planeOrigin,
                                            const unsigned int& marker,
                                            IMeshDAO& mesh) const
  {
    // set cell0Ds markers
    std::vector<bool> vertices_on_plane(mesh.Cell0DTotalNumber(),
                                        false);
    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      if (geometryUtilities.IsPointOnPlane(mesh.Cell0DCoordinates(v),
                                           planeNormal,
                                           planeOrigin))
      {
        vertices_on_plane[v] = true;
        mesh.Cell0DSetMarker(v,
                             marker);
      }
    }

    // set cell1Ds markers
    for (unsigned int s = 0; s < mesh.Cell1DTotalNumber(); s++)
    {
      const Eigen::VectorXi extremes = mesh.Cell1DExtremes(s);

      if (vertices_on_plane[extremes[0]] &&
          vertices_on_plane[extremes[1]])
      {
        mesh.Cell1DSetMarker(s,
                             marker);
      }
    }

    // set cell2Ds markers
    for (unsigned int p = 0; p < mesh.Cell2DTotalNumber(); p++)
    {
      const vector<unsigned int> extremes = mesh.Cell2DVertices(p);

      bool isOnPlane = true;
      for (unsigned int v = 0; v < extremes.size(); v++)
      {
        if (!vertices_on_plane[extremes[v]])
        {
          isOnPlane = false;
          break;
        }
      }

      if (isOnPlane)
      {
        mesh.Cell2DSetMarker(p,
                             marker);
      }
    }
  }
  // ***************************************************************************
  vector<unsigned int> MeshUtilities::MeshCell2DRoots(const IMeshDAO& mesh) const
  {
    vector<unsigned int> rootCell2Ds(mesh.Cell2DTotalNumber());

    for (unsigned int cc = 0; cc < mesh.Cell2DTotalNumber(); cc++)
    {
      unsigned int rootCell = cc;
      while (mesh.Cell2DHasOriginalCell2D(rootCell))
        rootCell = mesh.Cell2DOriginalCell2D(rootCell);

      rootCell2Ds[cc] = rootCell;
    }

    return rootCell2Ds;
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
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(convexMesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(convexMesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < convexMesh.Cell2DTotalNumber(); c++)
    {
      const unsigned int& domainCell2DIndex = c;

      // Extract original cell2D geometric information
      vector<Eigen::Matrix3d> convexCell2DTriangulationPoints;
      double convexCell2DArea;
      Eigen::Vector3d convexCell2DCentroid;

      const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // compute original cell2D triangulation
      const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
      const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                            convexCell2DUnalignedVerticesFilter);

      const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
      vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
      for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
        convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

      convexCell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                     convexCell2DTriangulation);

      const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();
      unsigned int cell2DTriangulationSize = numConvexCell2DTriangulation;

      // compute original cell2D area and centroids
      Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
      Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
      for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
      {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
      }

      convexCell2DArea = convexCell2DTriangulationAreas.sum();
      convexCell2DCentroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                               convexCell2DTriangulationAreas,
                                                               convexCell2DArea);

      result.Cell2DsVertices[c] = convexCell2DVertices;

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints.size(); cct++)
        result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cct];

      result.Cell2DsAreas[c] = convexCell2DArea;
      result.Cell2DsCentroids[c] = convexCell2DCentroid;
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      result.Cell2DsEdgeDirections[c].resize(convexMesh.Cell2DNumberEdges(domainCell2DIndex));
      for (unsigned int e = 0; e < convexMesh.Cell2DNumberEdges(domainCell2DIndex); e++)
      {
        const unsigned int origin = convexMesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = convexMesh.Cell2DVertex(domainCell2DIndex,
                                                         (e + 1) % convexMesh.Cell2DNumberEdges(domainCell2DIndex));

        result.Cell2DsEdgeDirections[c][e] = convexMesh.Cell1DExists(origin,
                                                                     end);
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData2D MeshUtilities::FillMesh2DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& mesh,
                                                                            const IMeshDAO& convexMesh,
                                                                            const vector<vector<unsigned int>>& meshCell2DToConvexCell2DIndices) const
  {
    MeshGeometricData2D result;

    result.Cell2DsVertices.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsTriangulations.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsAreas.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsCentroids.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsDiameters.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeDirections.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeLengths.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeTangents.resize(mesh.Cell2DTotalNumber());
    result.Cell2DsEdgeNormals.resize(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      const unsigned int& domainCell2DIndex = c;
      const vector<unsigned int>& domainConvexCell2DIndices = meshCell2DToConvexCell2DIndices[domainCell2DIndex];
      const unsigned int& numConvexCells = domainConvexCell2DIndices.size();

      // Get domain cell2D geometry information
      map<unsigned int, unsigned int> cell2DVerticesToPosition;
      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(domainCell2DIndex); v++)
        cell2DVerticesToPosition.insert(pair<unsigned int, unsigned int>(mesh.Cell2DVertex(domainCell2DIndex,
                                                                                           v),
                                                                         v));

      // Extract original cell2D geometric information
      unsigned int cell2DTriangulationSize = 0;
      vector<vector<Eigen::Matrix3d>> convexCell2DTriangulationPoints(numConvexCells);
      Eigen::VectorXd convexCell2DAreas(numConvexCells);
      Eigen::MatrixXd convexCell2DCentroids(3, numConvexCells);

      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        const unsigned int& domainConvexCell2DIndex = domainConvexCell2DIndices[cc];
        const Eigen::MatrixXd convexCell2DVertices = convexMesh.Cell2DVerticesCoordinates(domainConvexCell2DIndex);

        // compute original cell2D triangulation
        const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(convexCell2DVertices);
        const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(convexCell2DVertices,
                                                                                              convexCell2DUnalignedVerticesFilter);

        const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
        vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
        for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
          convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

        convexCell2DTriangulationPoints[cc] = geometryUtilities.ExtractTriangulationPoints(convexCell2DVertices,
                                                                                           convexCell2DTriangulation);

        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints[cc].size();
        cell2DTriangulationSize += numConvexCell2DTriangulation;

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cc][cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cc][cct],
                                                                                          convexCell2DTriangulationAreas[cct]);
        }

        convexCell2DAreas[cc] = convexCell2DTriangulationAreas.sum();
        convexCell2DCentroids.col(cc) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                          convexCell2DTriangulationAreas,
                                                                          convexCell2DAreas[cc]);
      }

      result.Cell2DsVertices[c] = mesh.Cell2DVerticesCoordinates(domainCell2DIndex);

      // Compute cell2D triangulation from original cell2Ds
      result.Cell2DsTriangulations[c].resize(cell2DTriangulationSize);
      unsigned int triangulationCounter = 0;
      for (unsigned int cc = 0; cc < numConvexCells; cc++)
      {
        for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints[cc].size(); cct++)
          result.Cell2DsTriangulations[c][triangulationCounter++] = convexCell2DTriangulationPoints[cc][cct];
      }

      result.Cell2DsAreas[c] = convexCell2DAreas.sum();
      result.Cell2DsCentroids[c] = geometryUtilities.PolygonCentroid(convexCell2DCentroids,
                                                                     convexCell2DAreas,
                                                                     result.Cell2DsAreas[c]);
      result.Cell2DsDiameters[c] = geometryUtilities.PolygonDiameter(result.Cell2DsVertices[c]);

      result.Cell2DsEdgeDirections[c].resize(mesh.Cell2DNumberEdges(domainCell2DIndex));
      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(domainCell2DIndex); e++)
      {
        const unsigned int origin = mesh.Cell2DVertex(domainCell2DIndex, e);
        const unsigned int end = mesh.Cell2DVertex(domainCell2DIndex,
                                                   (e + 1) % mesh.Cell2DNumberEdges(domainCell2DIndex));

        result.Cell2DsEdgeDirections[c][e] = mesh.Cell1DExists(origin,
                                                               end);
      }

      result.Cell2DsEdgeLengths[c] = geometryUtilities.PolygonEdgeLengths(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeTangents[c] = geometryUtilities.PolygonEdgeTangents(result.Cell2DsVertices[c]);
      result.Cell2DsEdgeNormals[c] = geometryUtilities.PolygonEdgeNormals(result.Cell2DsVertices[c]);
    }

    return result;
  }
  // ***************************************************************************
  MeshUtilities::MeshGeometricData3D MeshUtilities::FillMesh3DGeometricData(const GeometryUtilities& geometryUtilities,
                                                                            const IMeshDAO& convexMesh) const
  {
    MeshGeometricData3D result;

    result.Cell3DsVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdges.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsVolumes.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsEdgeDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsTetrahedronPoints.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFacesTranslations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesRotationMatrices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormals.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesNormalDirections.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeDirections.resize(convexMesh.Cell3DTotalNumber());

    result.Cell3DsFaces3DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DVertices.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DTriangulations.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesAreas.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFaces2DCentroids.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesDiameters.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdgeLengths.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DTangents.resize(convexMesh.Cell3DTotalNumber());
    result.Cell3DsFacesEdge2DNormals.resize(convexMesh.Cell3DTotalNumber());

    for(unsigned int c = 0; c <  convexMesh.Cell3DTotalNumber(); c++)
    {
      GeometryUtilities::Polyhedron polyhedron = MeshCell3DToPolyhedron(convexMesh,
                                                                        c);

      result.Cell3DsVertices[c] = polyhedron.Vertices;
      result.Cell3DsEdges[c] = polyhedron.Edges;
      result.Cell3DsFaces[c] = polyhedron.Faces;

      result.Cell3DsEdgeTangents[c] = geometryUtilities.PolyhedronEdgeTangents(result.Cell3DsVertices[c],
                                                                               result.Cell3DsEdges[c]);
      result.Cell3DsEdgeDirections[c].resize(convexMesh.Cell3DNumberEdges(c));
      for (unsigned int e = 0; e < convexMesh.Cell3DNumberEdges(c); e++)
      {
        const unsigned int meshOrigin = convexMesh.Cell3DVertex(c, polyhedron.Edges(0, e));
        const unsigned int meshEnd = convexMesh.Cell3DVertex(c, polyhedron.Edges(1, e));

        result.Cell3DsEdgeDirections[c][e] = convexMesh.Cell1DExists(meshOrigin,
                                                                     meshEnd);
      }

      result.Cell3DsFaces3DVertices[c] = geometryUtilities.PolyhedronFaceVertices(result.Cell3DsVertices[c],
                                                                                  result.Cell3DsFaces[c]);
      result.Cell3DsFacesTranslations[c] = geometryUtilities.PolyhedronFaceTranslations(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesNormals[c] = geometryUtilities.PolyhedronFaceNormals(result.Cell3DsFaces3DVertices[c]);
      result.Cell3DsFacesRotationMatrices[c] = geometryUtilities.PolyhedronFaceRotationMatrices(result.Cell3DsFaces3DVertices[c],
                                                                                                result.Cell3DsFacesNormals[c],
                                                                                                result.Cell3DsFacesTranslations[c]);



      /// TODO: qui secondo me questo non è corretto, occore utilizzare le informazioni della mesh
      const unsigned int numFaces = result.Cell3DsFaces[c].size();
      result.Cell3DsFacesEdgeDirections[c].resize(numFaces);
      for (unsigned int f = 0; f < numFaces; f++)
      {
        const unsigned int numFaceEdges = polyhedron.Faces[f].cols();
        result.Cell3DsFacesEdgeDirections[c][f].resize(numFaceEdges);
        for (unsigned int e = 0; e < numFaceEdges; e++)
        {
          const unsigned int faceEdgeOrigin = polyhedron.Faces[f](0, e);
          const unsigned int faceEdgeEnd = polyhedron.Faces[f](0, (e + 1) % numFaceEdges);

          const unsigned int meshOrigin = convexMesh.Cell3DVertex(c, faceEdgeOrigin);
          const unsigned int meshEnd = convexMesh.Cell3DVertex(c, faceEdgeEnd);

          result.Cell3DsFacesEdgeDirections[c][f][e] = convexMesh.Cell1DExists(meshOrigin,
                                                                               meshEnd);
        }
      }

      result.Cell3DsDiameters[c] = geometryUtilities.PolyhedronDiameter(result.Cell3DsVertices[c]);

      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    result.Cell3DsFaces3DVertices[c]);

      result.Cell3DsFaces2DVertices[c] = geometryUtilities.PolyhedronFaceRotatedVertices(result.Cell3DsFaces3DVertices[c],
                                                                                         result.Cell3DsFacesTranslations[c],
                                                                                         result.Cell3DsFacesRotationMatrices[c]);

      result.Cell3DsFaces2DTriangulations[c] = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(result.Cell3DsFaces2DVertices[c],
                                                                                                                polyhedronFaceTriangulations);


      result.Cell3DsFacesAreas[c].resize(numFaces);
      result.Cell3DsFacesDiameters[c].resize(numFaces);
      result.Cell3DsFaces2DCentroids[c].resize(numFaces);
      result.Cell3DsFacesEdgeLengths[c].resize(numFaces);
      result.Cell3DsFacesEdge2DNormals[c].resize(numFaces);
      result.Cell3DsFacesEdge2DTangents[c].resize(numFaces);

      for(unsigned int f = 0; f < numFaces; f++)
      {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d>& convexCell2DTriangulationPoints = result.Cell3DsFaces2DTriangulations[c][f];
        const unsigned int& numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
        {
          convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
          convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
        }

        result.Cell3DsFacesAreas[c][f] = convexCell2DTriangulationAreas.sum();
        result.Cell3DsFaces2DCentroids[c][f] = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                                 convexCell2DTriangulationAreas,
                                                                                 result.Cell3DsFacesAreas[c][f]);
        result.Cell3DsFacesDiameters[c][f] = geometryUtilities.PolygonDiameter(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdgeLengths[c][f] = geometryUtilities.PolygonEdgeLengths(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DTangents[c][f] = geometryUtilities.PolygonEdgeTangents(result.Cell3DsFaces2DVertices[c][f]);
        result.Cell3DsFacesEdge2DNormals[c][f] = geometryUtilities.PolygonEdgeNormals(result.Cell3DsFaces2DVertices[c][f]);
      }

      result.Cell3DsFacesNormalDirections[c] = geometryUtilities.PolyhedronFaceNormalDirections(result.Cell3DsFaces3DVertices[c],
                                                                                                geometryUtilities.PolyhedronBarycenter(result.Cell3DsVertices[c]),
                                                                                                result.Cell3DsFacesNormals[c]);

      result.Cell3DsVolumes[c] = geometryUtilities.PolyhedronVolume(result.Cell3DsFaces2DTriangulations[c],
                                                                    result.Cell3DsFacesNormals[c],
                                                                    result.Cell3DsFacesNormalDirections[c],
                                                                    result.Cell3DsFacesTranslations[c],
                                                                    result.Cell3DsFacesRotationMatrices[c]);

      result.Cell3DsCentroids[c] = geometryUtilities.PolyhedronCentroid(result.Cell3DsFaces2DTriangulations[c],
                                                                        result.Cell3DsFacesNormals[c],
                                                                        result.Cell3DsFacesNormalDirections[c],
                                                                        result.Cell3DsFacesTranslations[c],
                                                                        result.Cell3DsFacesRotationMatrices[c],
                                                                        result.Cell3DsVolumes[c]);

      vector<unsigned int> polyhedronTetrahedrons = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(result.Cell3DsVertices[c],
                                                                                                                 result.Cell3DsFaces[c],
                                                                                                                 polyhedronFaceTriangulations,
                                                                                                                 result.Cell3DsCentroids[c]);

      result.Cell3DsTetrahedronPoints[c] = geometryUtilities.ExtractTetrahedronPoints(result.Cell3DsVertices[c],
                                                                                      result.Cell3DsCentroids[c],
                                                                                      polyhedronTetrahedrons);
    }

    return result;
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell1DCell2DNeighbours(IMeshDAO& mesh) const
  {
    // Initialize cell1D neighbours
    for (unsigned int c1D = 0; c1D < mesh.Cell1DTotalNumber(); c1D++)
      mesh.Cell1DInitializeNeighbourCell2Ds(c1D, 2);

    // Compute Cell1D neighbours starting from cell2Ds
    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
      const unsigned int numCell2DEdges = mesh.Cell2DNumberEdges(c2D);
      for (unsigned int e = 0; e < numCell2DEdges; e++)
      {
        const unsigned int cell1D = mesh.Cell2DEdge(c2D, e);
        const unsigned int edgeOrigin =  mesh.Cell2DVertex(c2D, e);
        const unsigned int edgeEnd =  mesh.Cell2DVertex(c2D, (e + 1) % numCell2DEdges);

        if (mesh.Cell1DExists(edgeOrigin,
                              edgeEnd)) // left cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           1,
                                           c2D);
        }
        else // right cell
        {
          mesh.Cell1DInsertNeighbourCell2D(cell1D,
                                           0,
                                           c2D);
        }
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::ComputeCell2DCell3DNeighbours(IMeshDAO& mesh) const
  {
    // Compute Cell2D neighbours starting from cell3Ds
    std::vector<std::list<unsigned int>> cell2DsNeighbours(mesh.Cell2DTotalNumber());
    for (unsigned int c3D = 0; c3D < mesh.Cell3DTotalNumber(); c3D++)
    {
      const unsigned int numCell3DFaces = mesh.Cell3DNumberFaces(c3D);
      for (unsigned int f = 0; f < numCell3DFaces; f++)
      {
        const unsigned int cell2D = mesh.Cell3DFace(c3D, f);
        cell2DsNeighbours[cell2D].push_back(c3D);
      }
    }

    for (unsigned int c2D = 0; c2D < mesh.Cell2DTotalNumber(); c2D++)
    {
      mesh.Cell2DInitializeNeighbourCell3Ds(c2D, cell2DsNeighbours[c2D].size());

      unsigned int n = 0;
      for (const auto& cell3DIndex : cell2DsNeighbours[c2D])
        mesh.Cell2DInsertNeighbourCell3D(c2D,
                                         n++,
                                         cell3DIndex);
    }
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
  std::vector<unsigned int> MeshUtilities::SplitCell2D(const unsigned int& cell2DIndex,
                                                       const std::vector<Eigen::MatrixXi> subCell2Ds,
                                                       IMeshDAO& mesh) const
  {
    const unsigned int numSubCells = subCell2Ds.size();
    unsigned int newCell2DsStartingIndex = mesh.Cell2DAppend(numSubCells);

    vector<unsigned int> newCell2DsIndex(numSubCells);

    mesh.Cell2DSetState(cell2DIndex, false);

    for (unsigned int c = 0; c < numSubCells; c++)
    {
      newCell2DsIndex[c] = newCell2DsStartingIndex + c;

      const unsigned int& newCell2DIndex = newCell2DsIndex[c];
      mesh.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                     subCell2Ds[c]);

      mesh.Cell2DSetMarker(newCell2DIndex, mesh.Cell2DMarker(cell2DIndex));
      mesh.Cell2DSetState(newCell2DIndex, true);

      mesh.Cell2DInsertUpdatedCell2D(cell2DIndex, newCell2DIndex);

      for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(newCell2DIndex); e++)
      {
        const unsigned int cell1DIndex = mesh.Cell2DEdge(newCell2DIndex, e);

        for (unsigned int n = 0; n < mesh.Cell1DNumberNeighbourCell2D(cell1DIndex); n++)
        {
          if (!mesh.Cell1DHasNeighbourCell2D(cell1DIndex, n))
            continue;

          if (mesh.Cell1DNeighbourCell2D(cell1DIndex, n) == cell2DIndex)
            mesh.Cell1DInsertNeighbourCell2D(cell1DIndex,
                                             n,
                                             newCell2DIndex);
        }
      }
    }

    return newCell2DsIndex;
  }
  // ***************************************************************************
  void MeshUtilities::ImportOpenVolumeMesh(const std::string& ovmFilePath,
                                           IMeshDAO& mesh,
                                           std::vector<std::vector<bool>>& meshCell3DsFacesOrientation) const
  {
    OpenVolumeMeshInterface openVolumeMeshInterface;
    openVolumeMeshInterface.ImportMeshFromFile(ovmFilePath,
                                               mesh,
                                               meshCell3DsFacesOrientation);
  }
  // ***************************************************************************
  void MeshUtilities::ExportMeshToOpenVolume(const IMeshDAO& mesh,
                                             const std::vector<std::vector<bool>>& meshCell3DsFacesOrientation,
                                             const std::string& ovmFilePath) const
  {
    OpenVolumeMeshInterface openVolumeMeshInterface;
    openVolumeMeshInterface.ExportMeshToFile(mesh,
                                             meshCell3DsFacesOrientation,
                                             ovmFilePath);
  }
  // ***************************************************************************
  void MeshUtilities::ChangePolygonMeshMarkers(const Eigen::MatrixXd& polygonVertices,
                                               const vector<unsigned int>& cell0DMarkers,
                                               const vector<unsigned int>& cell1DMarkers,
                                               IMeshDAO& mesh) const
  {
    Output::Assert(mesh.Dimension() == 2);

    const unsigned int numVertices = polygonVertices.cols();

    for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
    {
      if (mesh.Cell0DMarker(v) == 0)
        continue;

      const unsigned int newMarker = (mesh.Cell0DMarker(v) <= numVertices) ?
                                       cell0DMarkers.at(mesh.Cell0DMarker(v) - 1) :
                                       cell1DMarkers.at(mesh.Cell0DMarker(v) - numVertices - 1);

      mesh.Cell0DSetMarker(v,
                           newMarker);
    }

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
    {
      if (mesh.Cell1DMarker(e) == 0)
        continue;

      const unsigned int newMarker = cell1DMarkers.at(mesh.Cell1DMarker(e) - numVertices - 1);

      mesh.Cell1DSetMarker(e,
                           newMarker);
    }
  }
  // ***************************************************************************
  void MeshUtilities::ExportMeshToVTU(const IMeshDAO& mesh,
                                      const string& exportFolder,
                                      const string& fileName) const
  {
    // Export Cell0Ds
    if (mesh.Cell0DTotalNumber() > 0)
    {
      vector<double> id(mesh.Cell0DTotalNumber());
      vector<double> marker(mesh.Cell0DTotalNumber());
      vector<double> active(mesh.Cell0DTotalNumber());

      for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
      {
        id[g] = g;
        marker[g] = mesh.Cell0DMarker(g);
        active[g] = mesh.Cell0DIsActive(g);
      }

      vector<VTPProperty> properties(3 + mesh.Cell0DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell0DNumberDoubleProperties());

      properties[0] = {
        "Id",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(id.size()),
        id.data()
      };
      properties[1] = {
        "Marker",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(marker.size()),
        marker.data()
      };
      properties[2] = {
        "Active",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(active.size()),
        active.data()
      };

      for (unsigned int p = 0; p < mesh.Cell0DNumberDoubleProperties(); p++)
      {
        propertyValues[p].resize(mesh.Cell0DTotalNumber());
        for (unsigned int g = 0; g < mesh.Cell0DTotalNumber(); g++)
        {
          propertyValues[p][g] = mesh.Cell0DDoublePropertySize(g, p) == 1 ? mesh.Cell0DDoublePropertyValue(g, p, 0) :
                                                                            0.0;
        }

        properties[3 + p] = {
          mesh.Cell0DDoublePropertyId(p),
          Gedim::VTPProperty::Formats::Cells,
          static_cast<unsigned int>(propertyValues[p].size()),
          propertyValues[p].data()
        };
      }

      Gedim::VTKUtilities vtpUtilities;
      vtpUtilities.AddPoints(mesh.Cell0DsCoordinates(),
                             properties);
      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell0Ds.vtu");
    }

    // Export Cell1Ds
    if (mesh.Cell1DTotalNumber() > 0)
    {
      vector<double> id(mesh.Cell1DTotalNumber());
      vector<double> marker(mesh.Cell1DTotalNumber());
      vector<double> active(mesh.Cell1DTotalNumber());

      for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
      {
        id[g] = g;
        marker[g] = mesh.Cell1DMarker(g);
        active[g] = mesh.Cell1DIsActive(g);
      }

      vector<VTPProperty> properties(3 + mesh.Cell1DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell1DNumberDoubleProperties());

      properties[0] = {
        "Id",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(id.size()),
        id.data()
      };
      properties[1] = {
        "Marker",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(marker.size()),
        marker.data()
      };
      properties[2] = {
        "Active",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(active.size()),
        active.data()
      };

      for (unsigned int p = 0; p < mesh.Cell1DNumberDoubleProperties(); p++)
      {
        propertyValues[p].resize(mesh.Cell1DTotalNumber());
        for (unsigned int g = 0; g < mesh.Cell1DTotalNumber(); g++)
        {
          propertyValues[p][g] = mesh.Cell1DDoublePropertySize(g, p) == 1 ? mesh.Cell1DDoublePropertyValue(g, p, 0) :
                                                                            0.0;
        }

        properties[3 + p] = {
          mesh.Cell1DDoublePropertyId(p),
          Gedim::VTPProperty::Formats::Cells,
          static_cast<unsigned int>(propertyValues[p].size()),
          propertyValues[p].data()
        };
      }

      Gedim::VTKUtilities vtpUtilities;
      vtpUtilities.AddSegments(mesh.Cell0DsCoordinates(),
                               mesh.Cell1DsExtremes(),
                               properties);
      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell1Ds.vtu");
    }

    // Export Cell2Ds
    if (mesh.Cell2DTotalNumber() > 0)
    {
      vector<double> id(mesh.Cell2DTotalNumber());
      vector<double> marker(mesh.Cell2DTotalNumber());
      vector<double> active(mesh.Cell2DTotalNumber());

      for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
      {
        id[g] = g;
        marker[g] = mesh.Cell2DMarker(g);
        active[g] = mesh.Cell2DIsActive(g);
      }

      vector<VTPProperty> properties(3 + mesh.Cell2DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell2DNumberDoubleProperties());

      properties[0] = {
        "Id",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(id.size()),
        id.data()
      };
      properties[1] = {
        "Marker",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(marker.size()),
        marker.data()
      };
      properties[2] = {
        "Active",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(active.size()),
        active.data()
      };

      for (unsigned int p = 0; p < mesh.Cell2DNumberDoubleProperties(); p++)
      {
        propertyValues[p].resize(mesh.Cell2DTotalNumber());
        for (unsigned int g = 0; g < mesh.Cell2DTotalNumber(); g++)
        {
          propertyValues[p][g] = mesh.Cell2DDoublePropertySize(g, p) == 1 ? mesh.Cell2DDoublePropertyValue(g, p, 0) :
                                                                            0.0;
        }

        properties[3 + p] = {
          mesh.Cell2DDoublePropertyId(p),
          Gedim::VTPProperty::Formats::Cells,
          static_cast<unsigned int>(propertyValues[p].size()),
          propertyValues[p].data()
        };
      }

      Gedim::VTKUtilities vtpUtilities;
      vtpUtilities.AddPolygons(mesh.Cell0DsCoordinates(),
                               mesh.Cell2DsVertices(),
                               properties);
      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell2Ds.vtu");
    }

    // Export Cell3Ds
    if (mesh.Cell3DTotalNumber() > 0)
    {
      vector<double> id(mesh.Cell3DTotalNumber());
      vector<double> marker(mesh.Cell3DTotalNumber());
      vector<double> active(mesh.Cell3DTotalNumber());

      for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
      {
        id[g] = g;
        marker[g] = mesh.Cell3DMarker(g);
        active[g] = mesh.Cell3DIsActive(g);
      }

      vector<VTPProperty> properties(3 + mesh.Cell3DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell3DNumberDoubleProperties());

      properties[0] = {
        "Id",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(id.size()),
        id.data()
      };
      properties[1] = {
        "Marker",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(marker.size()),
        marker.data()
      };
      properties[2] = {
        "Active",
        Gedim::VTPProperty::Formats::Cells,
        static_cast<unsigned int>(active.size()),
        active.data()
      };

      for (unsigned int p = 0; p < mesh.Cell3DNumberDoubleProperties(); p++)
      {
        propertyValues[p].resize(mesh.Cell3DTotalNumber());
        for (unsigned int g = 0; g < mesh.Cell3DTotalNumber(); g++)
        {
          propertyValues[p][g] = mesh.Cell3DDoublePropertySize(g, p) == 1 ? mesh.Cell3DDoublePropertyValue(g, p, 0) :
                                                                            0.0;
        }

        properties[3 + p] = {
          mesh.Cell3DDoublePropertyId(p),
          Gedim::VTPProperty::Formats::Cells,
          static_cast<unsigned int>(propertyValues[p].size()),
          propertyValues[p].data()
        };
      }

      Gedim::VTKUtilities vtpUtilities;
      vtpUtilities.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                                  mesh.Cell3DsFacesVertices(),
                                  properties);
      vtpUtilities.Export(exportFolder + "/" + fileName + "_Cell3Ds.vtu");
    }
  }
  // ***************************************************************************
  void MeshUtilities::ExportCell2DToVTU(const IMeshDAO&,
                                        const unsigned int& cell2DIndex,
                                        const Eigen::MatrixXd& cell2DVertices,
                                        const vector<Eigen::Matrix3d>& cell2DTriangulations,
                                        const double& cell2DArea,
                                        const Eigen::Vector3d& cell2DCentroid,
                                        const string& exportFolder) const
  {
    {
      Gedim::VTKUtilities vtpUtilities;

      vector<double> id(1, cell2DIndex);
      vector<double> area(1, cell2DArea);

      // Export cell2D
      vtpUtilities.AddPolygon(cell2DVertices,
                              {
                                {
                                  "Id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(id.size()),
                                  id.data()
                                },
                                {
                                  "Area",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(area.size()),
                                  area.data()
                                }
                              });

      vtpUtilities.Export(exportFolder + "/" + "Cell2D_" + to_string(cell2DIndex) + ".vtu");
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      // Export cell2D triangulation
      for (unsigned int t = 0; t < cell2DTriangulations.size(); t++)
      {
        vector<double> id(1, t);
        vtpUtilities.AddPolygon(cell2DTriangulations[t],
                                {
                                  {
                                    "Id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(id.size()),
                                    id.data()
                                  }
                                });
      }

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell2D_" + to_string(cell2DIndex) +
                          "_Triangles" + ".vtu");
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      // Export cell2D centroid
      vtpUtilities.AddPoint(cell2DCentroid);

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell2D_" + to_string(cell2DIndex) +
                          "_Centroid" + ".vtu");
    }
  }
  // ***************************************************************************
  GeometryUtilities::Polyhedron MeshUtilities::MeshCell3DToPolyhedron(const IMeshDAO& mesh,
                                                                      const unsigned int& cell3DIndex) const
  {
    GeometryUtilities::Polyhedron polyhedron;

    unordered_map<unsigned int, unsigned int> cell0DIndexToVertexIndex;
    unordered_map<unsigned int, unsigned int> cell1DIndexToEdgeIndex;

    polyhedron.Vertices = mesh.Cell3DVerticesCoordinates(cell3DIndex);
    polyhedron.Edges.resize(2, mesh.Cell3DNumberEdges(cell3DIndex));
    polyhedron.Faces.resize(mesh.Cell3DNumberFaces(cell3DIndex));

    for (unsigned int v = 0; v < polyhedron.Vertices.cols(); v++)
    {
      cell0DIndexToVertexIndex.insert(make_pair(mesh.Cell3DVertex(cell3DIndex,
                                                                  v),
                                                v));
    }

    for (unsigned int e = 0; e < polyhedron.Edges.cols(); e++)
    {
      const unsigned int cell1DIndex = mesh.Cell3DEdge(cell3DIndex,
                                                       e);
      cell1DIndexToEdgeIndex.insert(make_pair(cell1DIndex,
                                              e));

      polyhedron.Edges(0, e) = cell0DIndexToVertexIndex.at(mesh.Cell1DOrigin(cell1DIndex));
      polyhedron.Edges(1, e) = cell0DIndexToVertexIndex.at(mesh.Cell1DEnd(cell1DIndex));
    }

    for (unsigned int f = 0; f < polyhedron.Faces.size(); f++)
    {
      const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                       f);
      const unsigned int numFaceVertices = mesh.Cell2DNumberVertices(cell2DIndex);

      polyhedron.Faces[f].resize(2, numFaceVertices);
      for (unsigned int v = 0; v < numFaceVertices; v++)
      {
        polyhedron.Faces[f](0, v) = cell0DIndexToVertexIndex.at(mesh.Cell2DVertex(cell2DIndex,
                                                                                  v));
        polyhedron.Faces[f](1, v) = cell1DIndexToEdgeIndex.at(mesh.Cell2DEdge(cell2DIndex ,
                                                                              v));
      }
    }

    return polyhedron;
  }
  // ***************************************************************************
  MeshUtilities::VTPPolyhedron MeshUtilities::MeshCell3DToVTPPolyhedron(const IMeshDAO& mesh,
                                                                        const unsigned int& cell3DIndex) const
  {
    VTPPolyhedron vtpPolyhedron;

    unordered_map<unsigned int, unsigned int> cell0DIndexToVertexIndex;

    vtpPolyhedron.Vertices = mesh.Cell3DVerticesCoordinates(cell3DIndex);
    vtpPolyhedron.PolyhedronFaces.resize(mesh.Cell3DNumberFaces(cell3DIndex));

    for (unsigned int v = 0; v < vtpPolyhedron.Vertices.cols(); v++)
    {
      cell0DIndexToVertexIndex.insert(make_pair(mesh.Cell3DVertex(cell3DIndex,
                                                                  v),
                                                v));
    }

    for (unsigned int f = 0; f < vtpPolyhedron.PolyhedronFaces.size(); f++)
    {
      const unsigned int cell2DIndex = mesh.Cell3DFace(cell3DIndex,
                                                       f);
      const unsigned int numFaceVertices = mesh.Cell2DNumberVertices(cell2DIndex);

      vtpPolyhedron.PolyhedronFaces[f].resize(numFaceVertices);
      for (unsigned int v = 0; v < numFaceVertices; v++)
      {
        vtpPolyhedron.PolyhedronFaces[f][v] = cell0DIndexToVertexIndex.at(mesh.Cell2DVertex(cell2DIndex,
                                                                                            v));
      }
    }

    return vtpPolyhedron;
  }
  // ***************************************************************************
}
