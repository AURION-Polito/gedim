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
      struct TriangleMaxEdgeDirection final
      {
          unsigned int MaxEdgeIndex = 0;
          unsigned int OppositeVertexIndex = 0;
      };

      struct TetrahedronMaxEdgeDirection final
      {
          unsigned int MaxEdgeIndex = 0;
          std::array<unsigned int, 2> OppositeVerticesIndex = {};
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

      struct RefinePolygon_CheckResult final
      {
          enum struct ResultTypes
          {
            Unknown = 0,
            Cell2DToBeSplitted = 1,
            Cell2DAlreadySplitted = 2,
            Cell2DSplitUnderTolerance = 3,
            SplitDirectionNotInsideCell2D = 4
          };

          struct Cell1DToSplit final
          {
              enum struct Types
              {
                Unknown = 0,
                NotInside = 1,
                EdgeLengthNotEnough = 2,
                OnlyLocalQualityNotEnough = 3,
                OnlyNeighQualityNotEnough = 4,
                BothQualityNotEnough = 5,
                OnlyLocalAlignedNotRespect = 6,
                OnlyNeighAlignedNotRespect = 7,
                BothAlignedNotRespect = 8,
                ToSplit = 9
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

          std::vector<unsigned int> Cell1DsIndex = {};
          std::vector<GeometryUtilities::LinePolygonPositionResult::EdgeIntersection> Cell1DsIntersection = {};
          std::vector<Cell1DToSplit> Cell1DsToSplit = {};
          ResultTypes ResultType = ResultTypes::Unknown;
      };

      struct CheckSplitType_Result final
      {
          enum struct SplitTypes
          {
            Unknown = 0,
            NoSplit = 1,
            NoNewVertices = 2,
            NewVertexFrom = 3,
            NewVertexTo = 4,
            NewVertices = 5
          };

          std::array<unsigned int, 2> NoNewVerticesIndex = {}; ///< valid only for NoNewVertices type
          SplitTypes Type = SplitTypes::Unknown;
      };

      struct RefinePolyhedron_Result final
      {
          enum struct ResultTypes
          {
            Unknown = 0,
            Successfull = 1,
            Cell3DAlreadySplitted = 2,
            Cell3DSplitUnderTolerance = 3,
            Cell3DSplitNone = 4
          };

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
              unsigned int OriginalCell3DEdgeIndex = 0;
          };

          struct RefinedCell2D final
          {
              enum struct Types
              {
                Unknown = 0,
                Updated = 1,
                New = 2
              };

              Types Type = Types::Unknown;
              std::vector<unsigned int> NewCell2DsIndex = {};
              unsigned int OriginalCell2DIndex = 0;
              unsigned int NewCell1DIndex = 0;
              std::vector<unsigned int> NewCell1DsPosition = {}; ///< Position in NewCell1DsIndex array
              unsigned int OriginalCell3DFaceIndex = 0;
          };

          ResultTypes ResultType = ResultTypes::Unknown;
          std::vector<unsigned int> NewCell0DsIndex = {};
          std::vector<RefinedCell1D> NewCell1DsIndex = {};
          std::vector<RefinedCell2D> NewCell2DsIndex = {};
          std::vector<unsigned int> NewCell3DsIndex = {};
      };

      struct RefinePolygon_Result final
      {
          enum struct ResultTypes
          {
            Unknown = 0,
            Successfull = 1,
            Cell2DAlreadySplitted = 2,
            Cell2DSplitUnderTolerance = 3,
            SplitDirectionNotInsideCell2D = 4,
            SplitQualityCheckCell2DFailed = 5
          };

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

          std::vector<unsigned int> NewCell0DsIndex = {};
          std::vector<RefinedCell1D> NewCell1DsIndex = {};
          std::vector<unsigned int> NewCell2DsIndex = {};

          CheckSplitType_Result::SplitTypes SplitType = CheckSplitType_Result::SplitTypes::Unknown;
          ResultTypes ResultType = ResultTypes::Unknown;
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

      struct RefinePolyhedron_UpdateNeighbour_Result final
      {
          struct UpdatedCell3D final
          {
              unsigned int OriginalCell3DIndex = 0;
              unsigned int NewCell3DIndex = 0;
          };

          std::vector<UpdatedCell3D> UpdatedCell3Ds = {};
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
              std::vector<Eigen::VectorXd> CentroidVerticesDistance = {};
              std::vector<double> InRadius = {};
              std::vector<double> Quality = {};
          };

          struct Cell1D_GeometricData final
          {
              unsigned int MaxAligned = 0;
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

      SplitCell1D_Result SplitCell1D(const unsigned int& cell1DIndex,
                                     const Eigen::Vector3d& newVertexCoordinate,
                                     IMeshDAO& mesh) const;
      /// \brief update cell2DIndex with a new splitted edge cell1DIndex by newCell0DIndex
      unsigned int UpdateCell2D_NewVertex(const unsigned int cell2DIndex,
                                          const bool cell2DEdgeDirection,
                                          const unsigned int cell2DEdgePosition,
                                          const std::vector<unsigned int>& newCell1DsIndex,
                                          const unsigned int newCell0DIndex,
                                          IMeshDAO& mesh) const;

      inline SplitCell1D_Result SplitCell1D_MiddlePoint(const unsigned int& cell1DIndex,
                                                        IMeshDAO& mesh) const
      {
        return SplitCell1D(cell1DIndex,
                           0.5 * (mesh.Cell0DCoordinates(mesh.Cell1DOrigin(cell1DIndex)) +
                                  mesh.Cell0DCoordinates(mesh.Cell1DEnd(cell1DIndex))),
                           mesh);
      }

      bool AreVerticesAligned(const Eigen::MatrixXd& cell2DVertices,
                              const unsigned int fromVertex,
                              const unsigned int toVertex) const;

      CheckSplitType_Result SplitPolygon_CheckSplitType(const Gedim::GeometryUtilities::PolygonTypes& cell2DPolygonType,
                                                        const Gedim::GeometryUtilities::PolygonTypes& cell2DUnalignedPolygonType,
                                                        const Eigen::MatrixXd& cell2DVertices,
                                                        const RefinePolygon_CheckResult& cell2DCheckToRefine) const;

      bool SplitPolygon_CheckIsNotToExtend(const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitOne,
                                           const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitTwo) const;
      bool SplitPolygon_CheckIsToSplit_Relaxed(const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitOne,
                                               const RefinePolygon_CheckResult::Cell1DToSplit& cell1DSplitTwo) const;

      bool SplitPolygon_IsAreaPositive(const Eigen::VectorXi& newCell2D_Indices,
                                       const Eigen::Matrix3d& cell2DRotation,
                                       const Eigen::Vector3d& cell2DTranslation,
                                       IMeshDAO& mesh) const;

      SplitPolygon_Result SplitPolygon_NoNewVertices(const unsigned int cell2DIndex,
                                                     const unsigned int cell2DNumVertices,
                                                     const unsigned int fromVertex,
                                                     const unsigned int toVertex,
                                                     const Eigen::Matrix3d& cell2DRotation,
                                                     const Eigen::Vector3d& cell2DTranslation,
                                                     IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertexFrom(const unsigned int cell2DIndex,
                                                     const unsigned int cell2DNumVertices,
                                                     const unsigned int fromEdge,
                                                     const unsigned int toVertex,
                                                     const Eigen::Matrix3d& cell2DRotation,
                                                     const Eigen::Vector3d& cell2DTranslation,
                                                     const unsigned int fromNewCell0DIndex,
                                                     const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                     const bool& fromEdgeDirection,
                                                     IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertexTo(const unsigned int cell2DIndex,
                                                   const unsigned int cell2DNumVertices,
                                                   const unsigned int fromVertex,
                                                   const unsigned int toEdge,
                                                   const Eigen::Matrix3d& cell2DRotation,
                                                   const Eigen::Vector3d& cell2DTranslation,
                                                   const unsigned int toNewCell0DIndex,
                                                   const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                   const bool& toEdgeDirection,
                                                   IMeshDAO& mesh) const;
      SplitPolygon_Result SplitPolygon_NewVertices(const unsigned int cell2DIndex,
                                                   const unsigned int cell2DNumVertices,
                                                   const unsigned int fromEdge,
                                                   const unsigned int toEdge,
                                                   const Eigen::Matrix3d& cell2DRotation,
                                                   const Eigen::Vector3d& cell2DTranslation,
                                                   const unsigned int fromNewCell0DIndex,
                                                   const unsigned int toNewCell0DIndex,
                                                   const std::vector<unsigned int>& fromSplitCell1DsIndex,
                                                   const std::vector<unsigned int>& toSplitCell1DsIndex,
                                                   const bool& fromEdgeDirection,
                                                   const bool& toEdgeDirection,
                                                   IMeshDAO& mesh) const;

      TriangleMaxEdgeDirection ComputeTriangleMaxEdgeDirection(const Eigen::VectorXd& edgesLength) const;

      PolygonDirection ComputePolygonMaxDiameterDirection(const Eigen::MatrixXd unalignedVertices,
                                                          const Eigen::Vector3d& centroid) const;
      PolygonDirection ComputePolygonMaxInertiaDirection(const Eigen::MatrixXd& unalignedVertices,
                                                         const Eigen::VectorXd& unalignedEdgesLength,
                                                         const Eigen::Vector3d& centroid,
                                                         const Eigen::Matrix3d& inertia) const;

      TetrahedronMaxEdgeDirection ComputeTetrahedronMaxEdgeDirection(const Eigen::MatrixXi& polyhedronEdges,
                                                                     const Eigen::VectorXd& edgesLength) const;

      /// \brief Refine Triangle Cell2D By Edge
      /// \param cell2DIndex the index of Cell2D from 0 to Cell2DTotalNumber()
      /// \param cell2DVertices the cell2D 2D vertices
      /// \param edgeIndex the edge local index to split
      /// \param oppositeVertexIndex the vertex opposite to edge local index
      /// \param mesh the mesh to be updated
      RefinePolygon_Result RefineTriangleCell_ByEdge(const unsigned int& cell2DIndex,
                                                     const unsigned int& edgeIndex,
                                                     const unsigned int& oppositeVertexIndex,
                                                     const std::vector<bool>& cell2DEdgesDirection,
                                                     const double& cell2DArea,
                                                     const Eigen::Matrix3d& cell2DRotation,
                                                     const Eigen::Vector3d& cell2DTranslation,
                                                     const Eigen::VectorXd& cell2DEdgesLength,
                                                     IMeshDAO& mesh) const;

      /// \brief Refine Polyhedral Cell3D By Plane
      RefinePolyhedron_Result RefinePolyhedronCell_ByPlane(const unsigned int& cell3DIndex,
                                                           const Eigen::MatrixXd& cell3DVertices,
                                                           const Eigen::MatrixXi& cell3DEdges,
                                                           const Eigen::VectorXd& cell3DEdgesLength,
                                                           const std::vector<Eigen::MatrixXi>& cell3DFaces,
                                                           const std::vector<Eigen::MatrixXd>& cell3DFaces3DVertices,
                                                           const std::vector<Eigen::MatrixXd>& cell3DFacesEdges3DTangent,
                                                           const std::vector<Eigen::Vector3d>& cell3DFacesTranslation,
                                                           const std::vector<Eigen::Matrix3d>& cell3DFacesRotationMatrix,
                                                           const double& cell3DVolume,
                                                           const Eigen::Vector3d& planeNormal,
                                                           const Eigen::Vector3d& planeOrigin,
                                                           const Eigen::Matrix3d& planeRotationMatrix,
                                                           const Eigen::Vector3d& planeTranslation,
                                                           IMeshDAO& mesh) const;

      RefinePolyhedron_UpdateNeighbour_Result RefinePolyhedronCell_UpdateFaceNeighbours(const unsigned int& cell3DIndex,
                                                                                        const unsigned int& cell2DIndex,
                                                                                        const unsigned int& newCell1DIndex,
                                                                                        const std::vector<unsigned int>& splitCell1DsOriginalIndex,
                                                                                        const std::vector<unsigned int>& splitCell1DsNewCell0DIndex,
                                                                                        const std::vector<std::vector<unsigned int>>& splitCell1DsUpdatedIndices,
                                                                                        const std::vector<unsigned int>& splitCell2DsIndex,
                                                                                        const std::vector<std::vector<std::vector<bool>>>& cell3DsFacesEdgesDirection,
                                                                                        std::map<unsigned int, unsigned int>& updatedCell2Ds,
                                                                                        IMeshDAO& mesh) const;

      RefinePolyhedron_UpdateNeighbour_Result RefinePolyhedronCell_UpdateEdgeNeighbours(const unsigned int& cell3DIndex,
                                                                                        const unsigned int& cell1DIndex,
                                                                                        const std::vector<unsigned int>& newCell1DsIndex,
                                                                                        const unsigned int& newCell0DIndex,
                                                                                        const std::vector<std::vector<std::vector<bool>>>& cell3DsFacesEdgesDirection,
                                                                                        std::map<unsigned int, unsigned int>& updatedCell2Ds,
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
                                               const std::vector<Eigen::Matrix3d>& cell2DsRotation,
                                               const std::vector<Eigen::Vector3d>& cell2DsTranslation,
                                               IMeshDAO& mesh) const;

      /// \brief Refine Polygon Cell2D By Direction
      RefinePolygon_CheckResult RefinePolygonCell_CheckRefinement(const unsigned int& cell2DIndex,
                                                                  const Eigen::MatrixXd& cell2DVertices,
                                                                  const Eigen::Vector3d& lineTangent,
                                                                  const Eigen::Vector3d& lineOrigin,
                                                                  const std::vector<double>& cell2DsQuality,
                                                                  const std::vector<unsigned int>& cell1DsAligned,
                                                                  const double& cell1DsQualityWeight,
                                                                  const double& cell1DsAlignedWeight,
                                                                  const double& cell2DArea,
                                                                  const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                                  const std::vector<bool>& cell2DEdgesDirection,
                                                                  IMeshDAO& mesh) const;

      /// \brief Refine Polygon Cell2D By Direction
      RefinePolygon_Result RefinePolygonCell_ByDirection(const unsigned int& cell2DIndex,
                                                         const Gedim::GeometryUtilities::PolygonTypes& cell2DPolygonType,
                                                         const Gedim::GeometryUtilities::PolygonTypes& cell2DUnalignedPolygonType,
                                                         const Eigen::MatrixXd& cell2DVertices,
                                                         const RefinePolygon_CheckResult& cell2DCheckToRefine,
                                                         const Eigen::Matrix3d& cell2DRotation,
                                                         const Eigen::Vector3d& cell2DTranslation,
                                                         const std::vector<bool>& cell2DEdgesDirection,
                                                         const bool& extendToNeighbours,
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

      RefinePolygon_CheckResult::Cell1DToSplit RefinePolygonCell_IsCell1DToSplit(const unsigned int& cell1DIndex,
                                                                                 const unsigned int& cell2DIndex,
                                                                                 const GeometryUtilities::LinePolygonPositionResult::EdgeIntersection& edgeIntersection,
                                                                                 const std::vector<Eigen::VectorXd>& cell2DsEdgesLength,
                                                                                 const double& cell1DsQualityWeight,
                                                                                 const double& cell1DsAlignedWeight,
                                                                                 const std::vector<double>& cell2DsQuality,
                                                                                 const std::vector<unsigned int>& cell1DsAligned,
                                                                                 const IMeshDAO& mesh) const;
  };

}

#endif // __RefinementUtilities_H
