#ifndef __RefinementUtilities_H
#define __RefinementUtilities_H

#include "IMeshDAO.hpp"
#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include <numeric>

namespace Gedim
{
  /// \brief RefinementUtilities
  /// \copyright See top level LICENSE file for details.
  class RefinementUtilities final
  {
    public:
      struct MaxEdgeDirection final
      {
          unsigned int MaxEdgeIndex = 0;
          unsigned int OppositeVertexIndex = 0;
      };

      struct PolygonDirection final
      {
          Eigen::Vector3d LineOrigin;
          Eigen::Vector3d LineTangent;
      };

      struct MeshQuality final
      {
          std::vector<double> Cell2DsQuality;
          std::vector<double> Cell1DsQuality;
      };

      struct SplitCell1D_Result final
      {
          unsigned int NewCell0DIndex = 0;
          std::vector<unsigned int> NewCell1DsIndex = {};
      };

      struct SplitPolygon_Result final
      {
          enum struct Types
          {
            Unknown = 0,
            Split = 1,
            NoSplit = 2
          };

          Types Type = Types::Unknown;
          unsigned int NewCell1DIndex = 0;
          std::vector<unsigned int> NewCell2DsIndex = {};
      };

      struct RefinePolygon_Result final
      {
          struct RefinedCell1D final
          {
              enum struct Types
              {
                Unknown = 0,
                Updated = 1,
                New = 2
              };

              Types Type = Types::Unknown;
              std::vector<unsigned int> NewCell1DsIndex = {};
              unsigned int OriginalCell1DIndex = 0;
              unsigned int NewCell0DIndex = 0;
              unsigned int OriginalCell2DEdgeIndex = 0;
          };

          struct Cell1DToSplit final
          {
              enum struct Types
              {
                Unknown = -1,
                NotInside = 0,
                EdgeLengthNotEnough = 1,
                OnlyLocalQualityNotEnough = 2,
                OnlyNeighQualityNotEnough = 3,
                BothQualityNotEnough = 4,
                OnlyLocalAlignedNotRespect = 5,
                OnlyNeighAlignedNotRespect = 6,
                BothAlignedNotRespect = 7,
                ToSplit = 8
              };

              bool IsIntersectionInside = false;
              bool IsEdgeLengthEnough = false;
              bool IsLocalQualityEnough = false;
              bool IsQualityEnough = false;
              std::vector<bool> IsNeighQualityEnough = {};
              bool IsLocalAlignedRespect = false;
              bool IsAlignedRespect = false;
              std::vector<bool> IsNeighAlignedRespect = {};
              bool IsToSplit = false;
              unsigned int Cell2DEdgeIndex = 0;
              Types Type = Types::Unknown;
          };

          std::vector<Cell1DToSplit> Cell1DsToSplit = {};

          std::vector<unsigned int> NewCell0DsIndex = {};
          std::vector<RefinedCell1D> NewCell1DsIndex = {};
          std::vector<unsigned int> NewCell2DsIndex = {};
      };

      struct RefinePolygon_UpdateNeighbour_Result final
      {
          struct UpdatedCell2D final
          {
              unsigned int OriginalCell2DIndex = 0;
              unsigned int NewCell2DIndex = 0;
          };

          std::vector<UpdatedCell2D> UpdatedCell2Ds = {};
      };

      struct Cell2Ds_GeometricData final
      {
          struct Cell2D_GeometricData final
          {
              std::vector<std::vector<unsigned int>> UnalignedVerticesIndex = {};
              std::vector<Eigen::MatrixXd> Vertices = {};
              std::vector<double> Area = {};
              std::vector<Eigen::Vector3d> Centroid = {};
              std::vector<std::vector<bool>> EdgesDirection = {};
              std::vector<Eigen::MatrixXd> EdgesNormal = {};
              std::vector<Eigen::VectorXd> EdgesLength = {};
              std::vector<std::vector<Eigen::Matrix3d>> Triangulations = {};
              std::vector<Eigen::Matrix3d> Inertia = {};
              std::vector<Eigen::MatrixXd> UnalignedVertices = {};
              std::vector<Eigen::VectorXd> UnalignedEdgesLength = {};
              std::vector<Eigen::VectorXd> CentroidEdgesDistance = {};
              std::vector<double> InRadius = {};
              std::vector<double> Quality = {};
          };

          struct Cell1D_GeometricData final
          {
              std::vector<unsigned int> Aligned = {};
          };

          Cell1D_GeometricData Cell1Ds;
          Cell2D_GeometricData Cell2Ds;
      };

    private:
      const GeometryUtilities& geometryUtilities;
      const MeshUtilities& meshUtilities;

    public:
      RefinementUtilities(const GeometryUtilities& geometryUtilities,
                          const MeshUtilities& meshUtilities);
      ~RefinementUtilities();

      SplitCell1D_Result SplitCell1D_MiddlePoint(const unsigned int& cell1DIndex,
                                                 IMeshDAO& mesh) const;
      bool SplitPolygon_CheckIsToSplit(const RefinePolygon_Result::Cell1DToSplit& cell1DSplitOne,
                                       const RefinePolygon_Result::Cell1DToSplit& cell1DSplitTwo) const;

      bool SplitPolygon_CheckArea(const Eigen::VectorXi& newCell2DVertices,
                                  IMeshDAO& mesh) const;

      SplitPolygon_Result SplitPolygon_NoNewVertices(const unsigned int& cell2DIndex,
                                                     const unsigned int cell2DNumVertices,
                                                     const unsigned int& fromVertex,
                                                     const unsigned int& toVertex,
                                                     IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertexFrom(const unsigned int& cell2DIndex,
                                                     const unsigned int cell2DNumVertices,
                                                     const unsigned int& fromEdge,
                                                     const unsigned int& toVertex,
                                                     const unsigned int& fromNewCell0DIndex,
                                                     const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                     const bool& fromEdgeDirection,
                                                     IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertexTo(const unsigned int& cell2DIndex,
                                                   const unsigned int cell2DNumVertices,
                                                   const unsigned int& fromVertex,
                                                   const unsigned int& toEdge,
                                                   const unsigned int& toNewCell0DIndex,
                                                   const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                   const bool& toEdgeDirection,
                                                   IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertices(const unsigned int& cell2DIndex,
                                                   const unsigned int cell2DNumVertices,
                                                   const unsigned int& fromEdge,
                                                   const unsigned int& toEdge,
                                                   const unsigned int& fromNewCell0DIndex,
                                                   const unsigned int& toNewCell0DIndex,
                                                   const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                   const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                   const bool& fromEdgeDirection,
                                                   const bool& toEdgeDirection,
                                                   IMeshDAO& mesh) const;

      MaxEdgeDirection ComputeTriangleMaxEdgeDirection(const Eigen::VectorXd& edgesLength) const;

      PolygonDirection ComputePolygonMaxDiameterDirection(const Eigen::MatrixXd unalignedVertices,
                                                          const Eigen::Vector3d& centroid) const;
      PolygonDirection ComputePolygonMaxInertiaDirection(const Eigen::MatrixXd& unalignedVertices,
                                                         const Eigen::VectorXd& unalignedEdgesLength,
                                                         const Eigen::Vector3d& centroid,
                                                         const Eigen::Matrix3d& inertia) const;

      /// \brief Refine Triangle Cell2D By Edge
      /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
      /// \param edgeIndex the edge local index to split
      /// \param oppositeVertexIndex the vertex opposite to edge local index
      /// \param mesh the mesh to be updated
      RefinePolygon_Result RefineTriangleCell_ByEdge(const unsigned int& cell2DIndex,
                                                     const unsigned int& edgeIndex,
                                                     const unsigned int& oppositeVertexIndex,
                                                     const std::vector<bool>& cell2DEdgesDirection,
                                                     const double& cell2DArea,
                                                     const Eigen::VectorXd& cell2DEdgesLength,
                                                     IMeshDAO& mesh) const;

      /// \brief Update Cell1D neighbours of refined triangle by edge with refine by edge
      /// \param cell2DIndex the index of Cell2D refined, from 0 to Cell2DTotalNumber()
      /// \param cell1DIndex the index of Cell1D splitted by the refinement, from 0 to Cell1DTotalNumber()
      /// \param newCell0DIndex the index of Cell0D created by the cell1D splitting process, from 0 to Cell0DTotalNumber()
      /// \param splitCell1DsIndex the indices of the new Cell1Ds created by the splitting process, from 0 to Cell1DTotalNumber()
      /// \param cell2DEdgeDirection the direction of the Cell1D splitted in the Cell2D
      /// \param mesh the mesh to be updated
      void RefineTriangleCell_UpdateNeighbours(const unsigned int& cell2DIndex,
                                               const unsigned int& cell1DIndex,
                                               const unsigned int& newCell0DIndex,
                                               const std::vector<unsigned int>& splitCell1DsIndex,
                                               const bool& cell2DEdgeDirection,
                                               IMeshDAO& mesh) const;

      /// \brief Refine Polygon Cell2D By Direction
      RefinePolygon_Result RefinePolygonCell_ByDirection(const unsigned int& cell2DIndex,
                                                         const Eigen::MatrixXd& cell2DVertices,
                                                         const Eigen::Vector3d& lineTangent,
                                                         const Eigen::Vector3d& lineOrigin,
                                                         const std::vector<double>& cell2DsQuality,
                                                         const std::vector<unsigned int>& cell1DsAligned,
                                                         const double& cell1DsQualityWeight,
                                                         const double& cell2DArea,
                                                         const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                         const std::vector<bool>& cell2DEdgesDirection,
                                                         IMeshDAO& mesh) const;

      RefinePolygon_UpdateNeighbour_Result RefinePolygonCell_UpdateNeighbours(const unsigned int& cell2DIndex,
                                                                              const unsigned int& cell1DIndex,
                                                                              const unsigned int& newCell0DIndex,
                                                                              const std::vector<unsigned int>& splitCell1DsIndex,
                                                                              const std::vector<std::vector<bool>>& cell2DsEdgesDirection,
                                                                              IMeshDAO& mesh) const;

      /// Compute the geometric data for all the mesh
      Cell2Ds_GeometricData RefinePolygonCell_InitializeGeometricData(const IMeshDAO& mesh) const;

      /// \brief Update the geometric data for only cell2Ds
      void RefinePolygonCell_UpdateGeometricData(const IMeshDAO& mesh,
                                                 const std::vector<unsigned int>& cell2DsIndex,
                                                 Cell2Ds_GeometricData& geometricData) const;

      RefinePolygon_Result::Cell1DToSplit RefinePolygonCell_IsCell1DToSplit(const unsigned int& cell1DIndex,
                                                                            const unsigned int& cell2DIndex,
                                                                            const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersection,
                                                                            const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                                            const double& cell1DsQualityWeight,
                                                                            const std::vector<double>& cell2DsQuality,
                                                                            const std::vector<unsigned int>& cell1DsAligned,
                                                                            const IMeshDAO& mesh) const;
  };

}

#endif // __RefinementUtilities_H
