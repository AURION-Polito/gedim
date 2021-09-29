#ifndef __ConformerMeshPolygon_H
#define __ConformerMeshPolygon_H

#include "Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"
#include "MeshMatrices.hpp"
#include "ConformerMeshSegment.hpp"

using namespace std;

namespace Gedim
{
  class ConformerMeshPolygon final
  {
    public:
      struct ConformMesh final {
          vector<list<unsigned int>> Mesh2DCell2DMesh1DCell1DIds; ///< for each cell2D of mesh2D the list of cell1D of mesh1D
      };

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;

      ConformerMeshSegment::ConformMesh::ConformMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                 ConformerMeshSegment::ConformMesh& result,
                                                                                 bool& found);

      void CheckSegmentOriginAndEndIntersections(const Vector3d& segmentOrigin,
                                                 const Vector3d& segmentEnd,
                                                 ConformerMeshSegment::ConformMesh& mesh1D,
                                                 const Gedim::IMeshDAO& mesh2DReader);

      void Cell2DMesh2DListCell1DToListCell0D(const Gedim::IMeshDAO& mesh2DReader,
                                              const unsigned int& cell1DMesh2DId,
                                              const list<unsigned int>& cell1DMesh2DUpdated,
                                              vector<unsigned int>& cell0DMesh2Ds,
                                              vector<unsigned int>& cell1DMesh2Ds);
      void Cell2DMesh2DUpdatedPointsAndEdges(const Gedim::IMeshDAO& mesh2DReader,
                                             const unsigned int& cell2DMesh2DId,
                                             const unsigned int& oldCell1DMesh2DId,
                                             const vector<unsigned int>& cell0DMesh2DUpdated,
                                             const vector<unsigned int>& cell1DMesh2DUpdated,
                                             vector<unsigned int>& cell2DMesh2DNewVertices,
                                             vector<unsigned int>& cell2DMesh2DNewEdges);

      void Cell2DMesh2DToMaps(const Gedim::IMeshDAO& mesh2DReader,
                              const unsigned int& cell2DMesh2DId,
                              map<unsigned int, unsigned int>& verticesMap,
                              map<unsigned int, unsigned int>& edgesMap);

      void Cell2DMesh2DToSplitInput(const list<unsigned int> cell1DMesh1DIds,
                                    const ConformerMeshSegment::ConformMesh& mesh1D,
                                    const Gedim::IMeshDAO& mesh2D,
                                    const unsigned int& cell2DMesh2DId,
                                    const map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                                    const map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                                    Gedim::GeometryUtilities::SplitPolygonInput& splitInput);

      void SplitCell2DMesh2D(const Vector3d& segmentOrigin,
                             const Vector3d& segmentTangent,
                             const list<unsigned int> cell1DMesh1DIds,
                             ConformerMeshSegment::ConformMesh& mesh1D,
                             Gedim::IMeshDAO& mesh2D,
                             const unsigned int& cell2DMesh2DId,
                             map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                             map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                             const Gedim::GeometryUtilities::SplitPolygonResult& splitResult);

      void InsertCell2DMesh2DMiddleEdgesPolygonUpdate(const Vector3d& segmentOrigin,
                                                      const Vector3d& segmentTangent,
                                                      ConformerMeshSegment::ConformMesh& mesh1D,
                                                      Gedim::IMeshDAO& mesh2D,
                                                      const list<unsigned int> cell1DMesh1DIds,
                                                      const unsigned int& cell2DMesh2DId);

      void InsertCell2DMesh2DMiddleEdgesPolygonCreation(const Vector3d& segmentOrigin,
                                                        const Vector3d& segmentTangent,
                                                        ConformerMeshSegment::ConformMesh& mesh1D,
                                                        Gedim::IMeshDAO& mesh2D,
                                                        const list<unsigned int> cell1DMesh1DIds,
                                                        const unsigned int& cell2DMesh2DId);

      void UpdateCell2DNeighbours(ConformerMeshSegment::ConformMesh& mesh1D,
                                  Gedim::IMeshDAO& mesh2D,
                                  const unsigned int& cell2DMesh2DId);

      /// \brief Get the leaf cells id of the mesh2D tree structure
      /// \param fatherCellId the father cellid
      /// \param mesh2DUpdatedCellIds the mesh2D tree structure
      /// \param newCellIds the leaf cell ids
      /// \return true if the father cell is a leaf, false otherwise
      /// \note works for Cell1D and Cell2D
      bool GetLeafCells(const unsigned int& fatherCellId,
                        const map<unsigned int, set<unsigned int> >& mesh2DUpdatedCellIds,
                        list<unsigned int>& newCellIds);

    public:
      ConformerMeshPolygon(const Gedim::GeometryUtilities& geometryUtilities);
      ~ConformerMeshPolygon();

      /// \brief Conformer the input Mesh2D with a linear mesh1D
      /// \param mesh1D the 1D mesh
      /// \param mesh2DConformed the resulting conformed mesh
      void CreateConformMesh(const Vector3d& segmentOrigin,
                             const Vector3d& segmentEnd,
                             ConformerMeshSegment::ConformMesh& mesh1D,
                             Gedim::IMeshDAO& mesh2D,
                             ConformerMeshPolygon::ConformMesh& meshConformedInformation);
  };
}

#endif
