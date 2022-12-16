#ifndef __ConformerMeshPolygon_H
#define __ConformerMeshPolygon_H

#include "Eigen/Eigen"
#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"
#include "MeshMatrices.hpp"
#include "ConformerMeshSegment.hpp"

namespace Gedim
{
  class ConformerMeshPolygon final
  {
    public:
      struct ConformerMeshPolygonConfiguration final
      {
          enum struct Types
          {
            Generalized = 0, // conform checking the intersection types
            OnlyOnEdges = 1 // mesh 2D is already conform on edges, check only vertices
          };

          Types Type = Types::Generalized;
      };

    private:
      const Gedim::GeometryUtilities& _geometryUtilities;
      ConformerMeshPolygonConfiguration configuration;

      ConformerMeshSegment::ConformMesh::ConformMeshPoint& InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                 ConformerMeshSegment::ConformMesh& result,
                                                                                 bool& found);

      void CheckSegmentOriginAndEndIntersections(const Eigen::Vector3d& segmentOrigin,
                                                 const Eigen::Vector3d& segmentEnd,
                                                 ConformerMeshSegment::ConformMesh& mesh1D,
                                                 const Gedim::IMeshDAO& mesh2DReader);

      void Cell2DMesh2DListCell1DToListCell0D(const Gedim::IMeshDAO& mesh2DReader,
                                              const unsigned int& cell1DMesh2DId,
                                              const std::list<unsigned int>& cell1DMesh2DUpdated,
                                              std::vector<unsigned int>& cell0DMesh2Ds,
                                              std::vector<unsigned int>& cell1DMesh2Ds);
      void Cell2DMesh2DUpdatedPointsAndEdges(const Gedim::IMeshDAO& mesh2DReader,
                                             const unsigned int& cell2DMesh2DId,
                                             const unsigned int& oldCell1DMesh2DId,
                                             const std::vector<unsigned int>& cell0DMesh2DUpdated,
                                             const std::vector<unsigned int>& cell1DMesh2DUpdated,
                                             std::vector<unsigned int>& cell2DMesh2DNewVertices,
                                             std::vector<unsigned int>& cell2DMesh2DNewEdges);

      void Cell2DMesh2DToMaps(const Gedim::IMeshDAO& mesh2DReader,
                              const unsigned int& cell2DMesh2DId,
                              std::map<unsigned int, unsigned int>& verticesMap,
                              std::map<unsigned int, unsigned int>& edgesMap);

      void Cell2DMesh2DToSplitInput(const std::list<unsigned int> cell1DMesh1DIds,
                                    const ConformerMeshSegment::ConformMesh& mesh1D,
                                    const Gedim::IMeshDAO& mesh2D,
                                    const unsigned int& cell2DMesh2DId,
                                    const std::map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                                    const std::map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                                    Gedim::GeometryUtilities::SplitPolygonInput& splitInput);

      void SplitCell2DMesh2D(const Eigen::Vector3d& segmentOrigin,
                             const Eigen::Vector3d& segmentTangent,
                             const std::list<unsigned int> cell1DMesh1DIds,
                             ConformerMeshSegment::ConformMesh& mesh1D,
                             Gedim::IMeshDAO& mesh2D,
                             const unsigned int& cell2DMesh2DId,
                             std::map<unsigned int, unsigned int>& cell2DMesh2DVerticesMap,
                             std::map<unsigned int, unsigned int>& cell2DMesh2DEdgesMap,
                             const Gedim::GeometryUtilities::SplitPolygonWithSegmentResult& splitResult);

      void UpdateCell2DMesh2DWithSegmentOnEdges(const Eigen::Vector3d& segmentOrigin,
                                                const Eigen::Vector3d& segmentTangent,
                                                ConformerMeshSegment::ConformMesh& mesh1D,
                                                Gedim::IMeshDAO& mesh2D,
                                                const std::vector<unsigned int>& cell1DMesh1DIds,
                                                const unsigned int& cell2DMesh2DId);

      void InsertCell2DMesh2DMiddleEdgesPolygonUpdate(const Eigen::Vector3d& segmentOrigin,
                                                      const Eigen::Vector3d& segmentTangent,
                                                      ConformerMeshSegment::ConformMesh& mesh1D,
                                                      Gedim::IMeshDAO& mesh2D,
                                                      const std::list<unsigned int>& cell1DMesh1DIds,
                                                      const unsigned int& cell2DMesh2DId);

      void InsertCell2DMesh2DMiddleEdgesPolygonCreation(const Eigen::Vector3d& segmentOrigin,
                                                        const Eigen::Vector3d& segmentTangent,
                                                        ConformerMeshSegment::ConformMesh& mesh1D,
                                                        Gedim::IMeshDAO& mesh2D,
                                                        const std::list<unsigned int> cell1DMesh1DIds,
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
                        const std::map<unsigned int, std::set<unsigned int> >& mesh2DUpdatedCellIds,
                        std::list<unsigned int>& newCellIds);

      void CreateConformMeshGeneralized(const Eigen::Vector3d& segmentOrigin,
                                        const Eigen::Vector3d& segmentEnd,
                                        const Eigen::Vector3d& segmentTangent,
                                        ConformerMeshSegment::ConformMesh& mesh1D,
                                        Gedim::IMeshDAO& mesh2D);

      void CreateConformMeshOnlyOnEdges(const Eigen::Vector3d& segmentOrigin,
                                        const Eigen::Vector3d& segmentEnd,
                                        const Eigen::Vector3d& segmentTangent,
                                        ConformerMeshSegment::ConformMesh& mesh1D,
                                        Gedim::IMeshDAO& mesh2D);

    public:
      ConformerMeshPolygon(const Gedim::GeometryUtilities& geometryUtilities);
      ConformerMeshPolygon(const Gedim::GeometryUtilities& geometryUtilities,
                           const ConformerMeshPolygonConfiguration& configuration);
      ~ConformerMeshPolygon();

      /// \brief Conformer the input Mesh2D with a linear mesh1D
      /// \param mesh1D the 1D mesh
      /// \param mesh2DConformed the resulting conformed mesh
      void CreateConformMesh(const Eigen::Vector3d& segmentOrigin,
                             const Eigen::Vector3d& segmentEnd,
                             const Eigen::Vector3d& segmentTangent,
                             ConformerMeshSegment::ConformMesh& mesh1D,
                             Gedim::IMeshDAO& mesh2D);
  };
}

#endif
