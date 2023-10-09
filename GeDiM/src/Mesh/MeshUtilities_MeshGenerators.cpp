#include "MeshUtilities.hpp"

#include "TriangleInterface.hpp"
#include "TetgenInterface.hpp"
#include "VTKUtilities.hpp"
#include "MapTetrahedron.hpp"

#include "numeric"

using namespace std;
using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  void MeshUtilities::CreateRectangleMesh(const Eigen::Vector3d& rectangleOrigin,
                                          const Eigen::Vector3d& rectangleBaseTangent,
                                          const Eigen::Vector3d& rectangleHeightTangent,
                                          const vector<double>& baseMeshCurvilinearCoordinates,
                                          const vector<double>& heightMeshCurvilinearCoordinates,
                                          IMeshDAO& mesh) const
  {
    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    const unsigned int numCell0Ds = numBasePoints * numHeightPoints;
    const unsigned int numCell1Ds = numHeightPoints * (numBasePoints - 1) + numBasePoints * (numHeightPoints - 1);
    const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

    mesh.InitializeDimension(2);

    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell2DsInitialize(numCell2Ds);

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const Eigen::Vector3d coordinate = rectangleOrigin +
                                           baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                           heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;
        const unsigned int marker = 1 * (b == 0 && h == 0) +
                                    2 * (b == (numBasePoints - 1) && h == 0) +
                                    4 * (b == 0 && h == (numHeightPoints - 1)) +
                                    3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                    5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                    7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                    8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                    6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

        mesh.Cell0DSetState(cell0DIndex, true);
        mesh.Cell0DInsertCoordinates(cell0DIndex,
                                     coordinate);

        mesh.Cell0DSetMarker(cell0DIndex, marker);
        cell0DIndex++;
      }
    }

    // create cell1Ds
    unsigned int cell1DIndex = 0;

    // create horizontal cell1Ds
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + 1;

        const unsigned int marker = 5 * (h == 0) +
                                    7 * (h == (numHeightPoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create vertical cell1Ds
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

        const unsigned int marker = 8 * (b == 0) +
                                    6 * (b == (numBasePoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create cell2Ds
    unsigned int cell2DIndex = 0;
    mesh.Cell2DsInitializeVertices(4);
    mesh.Cell2DsInitializeEdges(4);
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DHorizontalIndex = b + h * (numBasePoints - 1);
        const unsigned int cell1DVerticalIndex = cell0DIndex + numHeightPoints * (numBasePoints - 1);

        vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                cell0DIndex + 1,
                                                cell0DIndex + numBasePoints + 1,
                                                cell0DIndex + numBasePoints };
        vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                             cell1DVerticalIndex + 1,
                                             cell1DHorizontalIndex + (numBasePoints - 1),
                                             cell1DVerticalIndex
                                           };

        mesh.Cell2DInsertVertices(cell2DIndex, cell2DVertices);
        mesh.Cell2DInsertEdges(cell2DIndex, cell2DEdges);

        mesh.Cell2DSetState(cell2DIndex, true);

        mesh.Cell2DSetMarker(cell2DIndex, 0);

        cell2DIndex++;
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateRectanglePlusHangingNodesMesh(const Eigen::Vector3d& rectangleOrigin,
                                                          const Eigen::Vector3d& rectangleBaseTangent,
                                                          const Eigen::Vector3d& rectangleHeightTangent,
                                                          const vector<double>& baseMeshCurvilinearCoordinates,
                                                          const vector<double>& heightMeshCurvilinearCoordinates,
                                                          const vector<unsigned int>& numberOfAddedVerticesForEachRectangle,
                                                          const GeometryUtilities& geometryUtilities,
                                                          IMeshDAO& mesh) const
  {
    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    unsigned int totalAddedHangNodes = std::accumulate(numberOfAddedVerticesForEachRectangle.begin(), numberOfAddedVerticesForEachRectangle.end(), 0);

    if ( totalAddedHangNodes == 0)
    {
      MeshUtilities::CreateRectangleMesh(rectangleOrigin,
                                         rectangleBaseTangent,
                                         rectangleHeightTangent,
                                         baseMeshCurvilinearCoordinates,
                                         heightMeshCurvilinearCoordinates,
                                         mesh);
    }
    else
    {

      const unsigned int numCell0Ds = numBasePoints * numHeightPoints
                                      + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                      + ceil(numBasePoints * 0.5) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[1]
                                      + (numHeightPoints / 2) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2]
                                      + (numBasePoints / 2) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[3];

      const unsigned int numCell1Ds = ceil((numHeightPoints * 0.5)) *
                                      (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1)
          + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
          * ( numberOfAddedVerticesForEachRectangle[1] + 1)
          + (numHeightPoints / 2) * (numBasePoints - 1) *
          (numberOfAddedVerticesForEachRectangle[2] + 1)
          + (numBasePoints / 2) * (numHeightPoints - 1) *
          (numberOfAddedVerticesForEachRectangle[3] + 1);

      const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

      mesh.InitializeDimension(2);

      mesh.Cell0DsInitialize(numCell0Ds);
      mesh.Cell1DsInitialize(numCell1Ds);
      mesh.Cell2DsInitialize(numCell2Ds);

      // create cell0Ds (# numBasePoints * numHeightPoints)
      unsigned int cell0DIndex = 0;
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b++)
        {
          const Eigen::Vector3d coordinate = rectangleOrigin +
                                             baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                             heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;
          const unsigned int marker = 1 * (b == 0 && h == 0) +
                                      2 * (b == (numBasePoints - 1) && h == 0) +
                                      4 * (b == 0 && h == (numHeightPoints - 1)) +
                                      3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                      5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                      7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                      8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                      6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

          mesh.Cell0DSetState(cell0DIndex, true);
          mesh.Cell0DInsertCoordinates(cell0DIndex,
                                       coordinate);

          mesh.Cell0DSetMarker(cell0DIndex, marker);
          cell0DIndex++;
        }
      }

      // create hanging nodes onto the base of rectangle
      // (# numberOfAddedNodesToFirstNEdges * (numBasePoints - 1) * ceil((numHeightPoints * 0.5)))
      for (unsigned int h = 0; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < (numBasePoints - 1); b++)
        {

          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[0],
              baseMeshCurvilinearCoordinates[b],
              baseMeshCurvilinearCoordinates[b+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               curvilinearPoints[s] * rectangleBaseTangent +
                                               heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;



            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 0; h < (numHeightPoints - 1); h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b = b + 2)
        {

          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[1],
              heightMeshCurvilinearCoordinates[h],
              heightMeshCurvilinearCoordinates[h+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                               curvilinearPoints[s] * rectangleHeightTangent;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 1; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < (numBasePoints - 1); b++)
        {
          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[2],
              baseMeshCurvilinearCoordinates[b],
              baseMeshCurvilinearCoordinates[b+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               curvilinearPoints[s] * rectangleBaseTangent +
                                               heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;



            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      for (unsigned int h = 0; h < (numHeightPoints-1); h++)
      {
        for (unsigned int b = 1; b < numBasePoints; b = b + 2)
        {
          std::vector<double> curvilinearPoints = geometryUtilities.EquispaceCoordinates(numberOfAddedVerticesForEachRectangle[3],
              heightMeshCurvilinearCoordinates[h],
              heightMeshCurvilinearCoordinates[h+1],
              false);
          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
          {
            const Eigen::Vector3d coordinate = rectangleOrigin +
                                               baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                               curvilinearPoints[s] * rectangleHeightTangent;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell0DSetState(cell0DIndex, true);
            mesh.Cell0DInsertCoordinates(cell0DIndex,
                                         coordinate);

            mesh.Cell0DSetMarker(cell0DIndex, marker);
            cell0DIndex++;
          }
        }
      }

      // create cell1Ds
      unsigned int cell1DIndex = 0;

      // create horizontal cell1Ds on b % 2 == 0
      unsigned int cell0DIndexBaseHangsEvenB = numBasePoints * numHeightPoints - 1;
      for (unsigned int h = 0; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {
          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsEvenB + 1;

            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsEvenB++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + 1;

          const unsigned int marker = 5 * (h == 0) +
                                      7 * (h == (numHeightPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create horizontal cell1Ds on b % 2 == 1
      unsigned int cell0DIndexBaseHangsOddB = numBasePoints * numHeightPoints - 1
                                              + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                              + ceil(numBasePoints * 0.5) * (numHeightPoints - 1) * numberOfAddedVerticesForEachRectangle[1];

      for (unsigned int h = 1; h < numHeightPoints; h = h + 2)
      {
        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsOddB + 1;

            const unsigned int marker = 5 * (h == 0) +
                                        7 * (h == (numHeightPoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsOddB++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + 1;

          const unsigned int marker = 5 * (h == 0) +
                                      7 * (h == (numHeightPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create vertical cell1Ds on h % 2 == 0
      unsigned int cell0DIndexBaseHangsEvenH = numBasePoints * numHeightPoints - 1
                                               + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int b = 0; b < numBasePoints; b = b + 2)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsEvenH + 1;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsEvenH++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

          const unsigned int marker = 8 * (b == 0) +
                                      6 * (b == (numBasePoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create vertical cell1Ds on h % 2 == 1
      unsigned int cell0DIndexBaseHangsOddH = numBasePoints * numHeightPoints - 1
                                              + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0]
                                              + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                              * numberOfAddedVerticesForEachRectangle[1]
                                              + (numHeightPoints / 2) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2];

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int b = 1; b < numBasePoints; b = b + 2)
        {

          const unsigned int cell0DIndex = b + h * numBasePoints;
          unsigned int cell1DOrigin = cell0DIndex;

          for (unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
          {
            const unsigned int cell1DEnd = cell0DIndexBaseHangsOddH + 1;

            const unsigned int marker = 8 * (b == 0) +
                                        6 * (b == (numBasePoints - 1));

            mesh.Cell1DInsertExtremes(cell1DIndex,
                                      cell1DOrigin,
                                      cell1DEnd);
            mesh.Cell1DSetState(cell1DIndex, true);
            mesh.Cell1DSetMarker(cell1DIndex, marker);

            cell1DIndex++;
            cell0DIndexBaseHangsOddH++;
            cell1DOrigin = cell1DEnd;
          }

          const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

          const unsigned int marker = 8 * (b == 0) +
                                      6 * (b == (numBasePoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }

      // create cell2Ds
      unsigned int cell2DIndex = 0;

      cell0DIndexBaseHangsEvenB = numBasePoints * numHeightPoints - 1;

      cell0DIndexBaseHangsEvenH = cell0DIndexBaseHangsEvenB
                                  + ceil((numHeightPoints * 0.5)) * (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];

      cell0DIndexBaseHangsOddB = cell0DIndexBaseHangsEvenH
                                 + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                 * (numberOfAddedVerticesForEachRectangle[1]);

      cell0DIndexBaseHangsOddH = cell0DIndexBaseHangsOddB
                                 + (numHeightPoints / 2) * (numBasePoints - 1) *
                                 (numberOfAddedVerticesForEachRectangle[2]);

      unsigned int cell1DHorizontalEvenB = 0;
      unsigned int cell1DHorizontalOddB = ceil((numHeightPoints * 0.5)) *
                                          (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1);
      unsigned int cell1DVerticalEvenH = cell1DHorizontalOddB
                                         + (numHeightPoints / 2) * (numBasePoints - 1) *
                                         (numberOfAddedVerticesForEachRectangle[2] + 1);
      unsigned int cell1DVerticalOddH = cell1DVerticalEvenH
                                        + ceil(numBasePoints * 0.5) * (numHeightPoints - 1)
                                        * ( numberOfAddedVerticesForEachRectangle[1] + 1);

      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        if ( h != 0)
        {
          if (h % 2 == 0)
          {
            cell0DIndexBaseHangsEvenB -= (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[0];
            cell1DHorizontalEvenB -= (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[0] + 1);
          }
          else
          {
            cell0DIndexBaseHangsOddB -= (numBasePoints - 1) * numberOfAddedVerticesForEachRectangle[2];
            cell1DHorizontalOddB -= (numBasePoints - 1) * (numberOfAddedVerticesForEachRectangle[2] + 1);
          }
        }

        for (unsigned int b = 0; b < numBasePoints - 1; b++)
        {

          if ( b != 0)
          {
            if (b % 2 == 0)
            {
              cell0DIndexBaseHangsEvenH -= numberOfAddedVerticesForEachRectangle[1];
              cell1DVerticalEvenH -= numberOfAddedVerticesForEachRectangle[1] + 1;
            }
            else
            {
              cell0DIndexBaseHangsOddH -= numberOfAddedVerticesForEachRectangle[3];
              cell1DVerticalOddH -= numberOfAddedVerticesForEachRectangle[3] + 1;
            }
          }

          const unsigned int cell0DIndex = b + h * numBasePoints;

          vector<unsigned int> cell2DVertices(4 + totalAddedHangNodes);
          vector<unsigned int> cell2DEdges(4 + totalAddedHangNodes);

          unsigned int countV = 0;
          unsigned int countE = 0;
          cell2DVertices[countV] = cell0DIndex;
          countV++;

          //unsigned int numberOfAddedNodesTolastEdges;

          if (h % 2 == 0)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenB + 1;
              cell0DIndexBaseHangsEvenB++;
              cell2DEdges[countE++] = cell1DHorizontalEvenB;
              cell1DHorizontalEvenB++;
            }
            cell2DEdges[countE++] = cell1DHorizontalEvenB;
            cell1DHorizontalEvenB++;
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddB + 1;
              cell0DIndexBaseHangsOddB++;
              cell2DEdges[countE++] = cell1DHorizontalOddB;
              cell1DHorizontalOddB++;
            }

            cell2DEdges[countE++] = cell1DHorizontalOddB;
            cell1DHorizontalOddB++;
          }

          cell2DVertices[countV] = cell0DIndex + 1;
          countV++;

          if(b % 2 == 0)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddH + 1;
              cell0DIndexBaseHangsOddH++;
              cell2DEdges[countE++] = cell1DVerticalOddH;
              cell1DVerticalOddH++;
            }
            cell2DEdges[countE++] = cell1DVerticalOddH;
            cell1DVerticalOddH++;
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenH + 1;
              cell0DIndexBaseHangsEvenH++;
              cell2DEdges[countE++] = cell1DVerticalEvenH;
              cell1DVerticalEvenH++;
            }
            cell2DEdges[countE++] = cell1DVerticalEvenH;
            cell1DVerticalEvenH++;
          }

          cell2DVertices[countV] = cell0DIndex + numBasePoints + 1;
          countV++;

          if (h % 2 == 1)
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[0]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenB + (numberOfAddedVerticesForEachRectangle[0] - s);
              cell2DEdges[countE++] = cell1DHorizontalEvenB + (numberOfAddedVerticesForEachRectangle[0] - s);
            }

            cell2DEdges[countE++] = cell1DHorizontalEvenB;
            cell1DHorizontalEvenB += numberOfAddedVerticesForEachRectangle[0] + 1;
            cell0DIndexBaseHangsEvenB += numberOfAddedVerticesForEachRectangle[0];
          }
          else
          {

            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[2]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddB + (numberOfAddedVerticesForEachRectangle[2] - s);
              cell2DEdges[countE++] = cell1DHorizontalOddB + (numberOfAddedVerticesForEachRectangle[2] - s);
            }

            cell2DEdges[countE++] = cell1DHorizontalOddB;
            cell1DHorizontalOddB += numberOfAddedVerticesForEachRectangle[2] + 1;
            cell0DIndexBaseHangsOddB += numberOfAddedVerticesForEachRectangle[2];
          }

          cell2DVertices[countV] = cell0DIndex + numBasePoints;
          countV++;

          if(b % 2 == 1)
          {

            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[3]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsOddH + (numberOfAddedVerticesForEachRectangle[3] - s);
              cell2DEdges[countE++] = cell1DVerticalOddH + (numberOfAddedVerticesForEachRectangle[3] - s);
            }

            cell2DEdges[countE++] = cell1DVerticalOddH;
            cell1DVerticalOddH += numberOfAddedVerticesForEachRectangle[3] + 1;
            cell0DIndexBaseHangsOddH += numberOfAddedVerticesForEachRectangle[3];
          }
          else
          {
            for(unsigned int s = 0; s < numberOfAddedVerticesForEachRectangle[1]; s++)
            {
              cell2DVertices[countV++] = cell0DIndexBaseHangsEvenH + (numberOfAddedVerticesForEachRectangle[1] - s);
              cell2DEdges[countE++] = cell1DVerticalEvenH + (numberOfAddedVerticesForEachRectangle[1] - s);
            }

            cell2DEdges[countE++] = cell1DVerticalEvenH;
            cell1DVerticalEvenH += numberOfAddedVerticesForEachRectangle[1] + 1;
            cell0DIndexBaseHangsEvenH += numberOfAddedVerticesForEachRectangle[1];
          }

          mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
          mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

          mesh.Cell2DSetState(cell2DIndex, true);

          mesh.Cell2DSetMarker(cell2DIndex, 0);

          cell2DIndex++;
        }
      }
    }
  }

  // ***************************************************************************
  void MeshUtilities::CreateRandomlyDeformedQuadrilaterals(const GeometryUtilities& geometryUtilities,
                                                           const Eigen::Vector3d& rectangleOrigin,
                                                           const Eigen::Vector3d& rectangleBaseTangent,
                                                           const Eigen::Vector3d& rectangleHeightTangent,
                                                           const unsigned int& numQuadrilateralsBaseTangent,
                                                           const unsigned int& numQuadrilateralsHeightTangent,
                                                           const double& maxDeformingPercentageBase,
                                                           const double& maxDeformingPercentageHeight,
                                                           IMeshDAO& mesh) const
  {


    vector<double> baseMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numQuadrilateralsBaseTangent + 1,
                                                                                           0.0, 1.0, 1);


    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numQuadrilateralsHeightTangent + 1,
                                                                                             0.0, 1.0, 1);

    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    const unsigned int numCell0Ds = numBasePoints * numHeightPoints;
    const unsigned int numCell1Ds = numHeightPoints * (numBasePoints - 1) + numBasePoints * (numHeightPoints - 1);
    const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

    mesh.InitializeDimension(2);

    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell2DsInitialize(numCell2Ds);

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        double diffBase = baseMeshCurvilinearCoordinates[1] - baseMeshCurvilinearCoordinates[0];
        double baseRand = diffBase * (static_cast<double>(rand()) / RAND_MAX) - diffBase * 0.5;
        double diffHeight = heightMeshCurvilinearCoordinates[1] - heightMeshCurvilinearCoordinates[0];
        double heightRand = diffHeight * (static_cast<double>(rand()) / RAND_MAX) - diffHeight * 0.5;
        const Eigen::Vector3d coordinate = rectangleOrigin +
                                           baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                           (b!= 0 && b != numBasePoints - 1) * baseRand * maxDeformingPercentageBase * rectangleBaseTangent +
                                           (h!= 0 && h != numHeightPoints - 1) * heightRand * maxDeformingPercentageHeight * rectangleHeightTangent +
                                           heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;
        const unsigned int marker = 1 * (b == 0 && h == 0) +
                                    2 * (b == (numBasePoints - 1) && h == 0) +
                                    4 * (b == 0 && h == (numHeightPoints - 1)) +
                                    3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                    5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                    7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                    8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                    6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

        mesh.Cell0DSetState(cell0DIndex, true);
        mesh.Cell0DInsertCoordinates(cell0DIndex,
                                     coordinate);

        mesh.Cell0DSetMarker(cell0DIndex, marker);
        cell0DIndex++;
      }
    }

    // create cell1Ds
    unsigned int cell1DIndex = 0;

    // create horizontal cell1Ds
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + 1;

        const unsigned int marker = 5 * (h == 0) +
                                    7 * (h == (numHeightPoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create vertical cell1Ds
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

        const unsigned int marker = 8 * (b == 0) +
                                    6 * (b == (numBasePoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create cell2Ds
    unsigned int cell2DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DHorizontalIndex = b + h * (numBasePoints - 1);
        const unsigned int cell1DVerticalIndex = cell0DIndex + numHeightPoints * (numBasePoints - 1);

        vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                cell0DIndex + 1,
                                                cell0DIndex + numBasePoints + 1,
                                                cell0DIndex + numBasePoints };
        vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                             cell1DVerticalIndex + 1,
                                             cell1DHorizontalIndex + (numBasePoints - 1),
                                             cell1DVerticalIndex
                                           };

        mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
        mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

        mesh.Cell2DSetState(cell2DIndex, true);

        mesh.Cell2DSetMarker(cell2DIndex, 0);

        cell2DIndex++;
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateDistortedQuadrilaterals(const GeometryUtilities& geometryUtilities,
                                                    const Eigen::Vector3d& rectangleOrigin,
                                                    const Eigen::Vector3d& rectangleBaseTangent,
                                                    const Eigen::Vector3d& rectangleHeightTangent,
                                                    const unsigned int& numQuadrilateralsBaseTangent,
                                                    const unsigned int& numQuadrilateralsHeightTangent,
                                                    IMeshDAO& mesh) const
  {


    vector<double> baseMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numQuadrilateralsBaseTangent + 1,
                                                                                           0.0, 1.0, 1);


    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numQuadrilateralsHeightTangent + 1,
                                                                                             0.0, 1.0, 1);

    const unsigned int& numBasePoints = baseMeshCurvilinearCoordinates.size();
    const unsigned int& numHeightPoints = heightMeshCurvilinearCoordinates.size();

    const unsigned int numCell0Ds = numBasePoints * numHeightPoints;
    const unsigned int numCell1Ds = numHeightPoints * (numBasePoints - 1) + numBasePoints * (numHeightPoints - 1);
    const unsigned int numCell2Ds = (numBasePoints - 1) * (numHeightPoints - 1);

    mesh.InitializeDimension(2);

    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell2DsInitialize(numCell2Ds);

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        Eigen::Vector3d coordinate = rectangleOrigin +
                                     baseMeshCurvilinearCoordinates[b] * rectangleBaseTangent +
                                     heightMeshCurvilinearCoordinates[h] * rectangleHeightTangent;

        coordinate(0) += 0.1 * sin(2.0 * M_PI * coordinate(0)) * sin(2.0 * M_PI * coordinate(1));
        coordinate(1) += 0.1 * sin(2.0 * M_PI * coordinate(0)) * sin(2.0 * M_PI * coordinate(1));

        const unsigned int marker = 1 * (b == 0 && h == 0) +
                                    2 * (b == (numBasePoints - 1) && h == 0) +
                                    4 * (b == 0 && h == (numHeightPoints - 1)) +
                                    3 * (b == (numBasePoints - 1) && h == (numHeightPoints - 1)) +
                                    5 * (h == 0 && b != 0 && b != (numBasePoints - 1)) +
                                    7 * (h == (numHeightPoints - 1) && b != 0 && b != (numBasePoints - 1)) +
                                    8 * (b == 0 && h != 0 && h != (numHeightPoints - 1)) +
                                    6 * (b == (numBasePoints - 1) && h != 0 && h != (numHeightPoints - 1));

        mesh.Cell0DSetState(cell0DIndex, true);
        mesh.Cell0DInsertCoordinates(cell0DIndex,
                                     coordinate);

        mesh.Cell0DSetMarker(cell0DIndex, marker);
        cell0DIndex++;
      }
    }

    // create cell1Ds
    unsigned int cell1DIndex = 0;

    // create horizontal cell1Ds
    for (unsigned int h = 0; h < numHeightPoints; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + 1;

        const unsigned int marker = 5 * (h == 0) +
                                    7 * (h == (numHeightPoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create vertical cell1Ds
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DOrigin = cell0DIndex;
        const unsigned int cell1DEnd = cell0DIndex + numBasePoints;

        const unsigned int marker = 8 * (b == 0) +
                                    6 * (b == (numBasePoints - 1));

        mesh.Cell1DInsertExtremes(cell1DIndex,
                                  cell1DOrigin,
                                  cell1DEnd);
        mesh.Cell1DSetState(cell1DIndex, true);
        mesh.Cell1DSetMarker(cell1DIndex, marker);

        cell1DIndex++;
      }
    }

    // create cell2Ds
    unsigned int cell2DIndex = 0;
    for (unsigned int h = 0; h < numHeightPoints - 1; h++)
    {
      for (unsigned int b = 0; b < numBasePoints - 1; b++)
      {
        const unsigned int cell0DIndex = b + h * numBasePoints;
        const unsigned int cell1DHorizontalIndex = b + h * (numBasePoints - 1);
        const unsigned int cell1DVerticalIndex = cell0DIndex + numHeightPoints * (numBasePoints - 1);

        vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                cell0DIndex + 1,
                                                cell0DIndex + numBasePoints + 1,
                                                cell0DIndex + numBasePoints };
        vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                             cell1DVerticalIndex + 1,
                                             cell1DHorizontalIndex + (numBasePoints - 1),
                                             cell1DVerticalIndex
                                           };

        mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
        mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

        mesh.Cell2DSetState(cell2DIndex, true);

        mesh.Cell2DSetMarker(cell2DIndex, 0);

        cell2DIndex++;
      }
    }
  }
  // ***************************************************************************
  void MeshUtilities::CreateParallelepipedMesh(const Eigen::Vector3d& parallelepipedOrigin,
                                               const Eigen::Vector3d& parallelepipedLengthTangent,
                                               const Eigen::Vector3d& parallelepipedHeightTangent,
                                               const Eigen::Vector3d& parallelepipedWidthTangent,
                                               const vector<double>& lengthMeshCurvilinearCoordinates,
                                               const vector<double>& heightMeshCurvilinearCoordinates,
                                               const vector<double>& widthMeshCurvilinearCoordinates,
                                               IMeshDAO& mesh) const
  {
    const unsigned int numLengthPoints = lengthMeshCurvilinearCoordinates.size();
    const unsigned int numHeightPoints = heightMeshCurvilinearCoordinates.size();
    const unsigned int numWidthPoints = widthMeshCurvilinearCoordinates.size();

    const unsigned int numCell0Ds = numLengthPoints * numHeightPoints * numWidthPoints;
    const unsigned int numCell1Ds = numHeightPoints * numWidthPoints * (numLengthPoints - 1)
                                    + numLengthPoints * numWidthPoints * (numHeightPoints - 1)
                                    + numLengthPoints * numHeightPoints * (numWidthPoints -1);
    const unsigned int numCell2Ds = numWidthPoints * (numLengthPoints - 1) * (numHeightPoints - 1)
                                    + numLengthPoints * (numWidthPoints - 1) * (numHeightPoints - 1)
                                    + numHeightPoints * (numWidthPoints - 1) * (numLengthPoints - 1);

    const unsigned int numCell3Ds =  (numLengthPoints - 1) * (numHeightPoints - 1) * (numWidthPoints - 1);

    mesh.InitializeDimension(3);

    mesh.Cell0DsInitialize(numCell0Ds);
    mesh.Cell1DsInitialize(numCell1Ds);
    mesh.Cell2DsInitialize(numCell2Ds);
    mesh.Cell3DsInitialize(numCell3Ds);

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints; l++)
        {
          const Eigen::Vector3d coordinate = parallelepipedOrigin +
                                             lengthMeshCurvilinearCoordinates[l] * parallelepipedLengthTangent +
                                             heightMeshCurvilinearCoordinates[h] * parallelepipedHeightTangent +
                                             widthMeshCurvilinearCoordinates[w] * parallelepipedWidthTangent;

          const unsigned int marker = 1 * (l == 0 && h == 0 && w == 0) +
                                      2 * (l == (numLengthPoints - 1) && h == 0 && w == 0) +
                                      3 * (l == (numLengthPoints - 1) && h == (numHeightPoints - 1) && w == 0) +
                                      4 * (l == 0 && h == (numHeightPoints - 1) && w == 0) +
                                      5 * (l == 0 && h == 0 && w == (numWidthPoints - 1)) +
                                      6 * (l == (numLengthPoints - 1) && h == 0 && w == (numWidthPoints - 1)) +
                                      7 * (l == (numLengthPoints - 1) && h == (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                      8 * (l == 0 && h == (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                      9 * (l != 0 && l != (numLengthPoints - 1) && h == 0 && w == 0) +
                                      10 * (l == (numLengthPoints - 1) && h != 0 && h != (numHeightPoints - 1) && w == 0) +
                                      11 * (l != 0 && l != (numLengthPoints - 1) && h == (numHeightPoints - 1) && w == 0) +
                                      12 * (l == 0 && h != 0 && h != (numHeightPoints - 1) && w == 0) +
                                      13 * (l != 0 && l != (numLengthPoints - 1) && h == 0 && w == (numWidthPoints - 1)) +
                                      14 * (l == (numLengthPoints - 1) && h != 0 && h != (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                      15 * (l != 0 && l != (numLengthPoints - 1) && h == (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                      16 * (l == 0 && h != 0 && h != (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                      17 * (l == 0 && h == 0 && w != 0 && w != (numWidthPoints - 1)) +
                                      18 * (l == (numLengthPoints -1) && h == 0 && w != 0 && w != (numWidthPoints - 1)) +
                                      19 * (l == (numLengthPoints -1) && h == (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1)) +
                                      20 * (l == 0 && h == (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1)) +
                                      21 * (l != 0 && l != (numLengthPoints - 1) && h != 0 && h != (numHeightPoints -1) && w == 0) +
                                      22 * (l != 0 && l != (numLengthPoints - 1) && h != 0 && h != (numHeightPoints -1) && w == (numWidthPoints - 1)) +
                                      23 * (l == 0 && h != 0 && h != (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1)) +
                                      24 * (l == (numLengthPoints - 1) && h != 0 && h != (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1)) +
                                      25 * (l != 0 && l  != (numLengthPoints - 1) && h == 0  && w != 0 && w != (numWidthPoints - 1)) +
                                      26 * (l != 0 && l != (numLengthPoints - 1) &&  h == (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1));

          mesh.Cell0DSetState(cell0DIndex, true);
          mesh.Cell0DInsertCoordinates(cell0DIndex,
                                       coordinate);

          mesh.Cell0DSetMarker(cell0DIndex, marker);
          cell0DIndex++;
        }
      }
    }

    // create cell1Ds
    unsigned int cell1DIndex = 0;

    // create length cell1Ds
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints - 1; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;
          const unsigned int cell1DOrigin = cell0DIndex;
          const unsigned int cell1DEnd = cell0DIndex + 1;

          const unsigned int marker =  9 * ( h == 0 && w == 0) +
                                       11 * (h == (numHeightPoints - 1) && w == 0) +
                                       13 * (h == 0 && w == (numWidthPoints - 1)) +
                                       15 * (h == (numHeightPoints - 1) && w == (numWidthPoints - 1)) +
                                       21 * (h != 0 && h != (numHeightPoints -1) && w == 0) +
                                       22 * (h != 0 && h != (numHeightPoints -1) && w == (numWidthPoints - 1)) +
                                       25 * (h == 0  && w != 0 && w != (numWidthPoints - 1)) +
                                       26 * (h == (numHeightPoints -1) && w != 0 && w != (numWidthPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }
    }

    // create height cell1Ds
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;
          const unsigned int cell1DOrigin = cell0DIndex;
          const unsigned int cell1DEnd = cell0DIndex + numLengthPoints;

          const unsigned int marker = 10 * (l == (numLengthPoints - 1) && w == 0) +
                                      12 * (l == 0 && w == 0) +
                                      14 * (l == (numLengthPoints - 1) && w == (numWidthPoints - 1)) +
                                      16 * (l == 0 && w == (numWidthPoints - 1)) +
                                      21 * (l != 0 && l != (numLengthPoints - 1) && w == 0) +
                                      22 * (l != 0 && l != (numLengthPoints - 1) && w == (numWidthPoints - 1)) +
                                      23 * (l == 0 && w != 0 && w != (numWidthPoints - 1)) +
                                      24 * (l == (numLengthPoints - 1) && w != 0 && w != (numWidthPoints - 1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }
    }

    // create width cell1Ds
    for (unsigned int w = 0; w < numWidthPoints - 1; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;
          const unsigned int cell1DOrigin = cell0DIndex;
          const unsigned int cell1DEnd = cell0DIndex + numLengthPoints * numHeightPoints;

          const unsigned int marker = 17 * (l == 0 && h == 0 ) +
                                      18 * (l == (numLengthPoints -1) && h == 0 ) +
                                      19 * (l == (numLengthPoints -1) && h == (numHeightPoints -1) ) +
                                      20 * (l == 0 && h == (numHeightPoints -1) ) +
                                      23 * (l == 0 && h != 0 && h != (numHeightPoints -1)) +
                                      24 * (l == (numLengthPoints - 1) && h != 0 && h != (numHeightPoints -1)) +
                                      25 * (l != 0 && l  != (numLengthPoints - 1) && h == 0) +
                                      26 * (l != 0 && l != (numLengthPoints - 1) &&  h == (numHeightPoints -1));

          mesh.Cell1DInsertExtremes(cell1DIndex,
                                    cell1DOrigin,
                                    cell1DEnd);
          mesh.Cell1DSetState(cell1DIndex, true);
          mesh.Cell1DSetMarker(cell1DIndex, marker);

          cell1DIndex++;
        }
      }
    }


    // create length x height cell2Ds
    unsigned int cell2DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints - 1; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;
          const unsigned int cell1DHorizontalIndex = l + h * (numLengthPoints - 1) + w * (numLengthPoints - 1) * numHeightPoints;
          const unsigned int cell1DVerticalIndex = l + h * numLengthPoints + w * numLengthPoints * (numHeightPoints -1)
                                                   + numHeightPoints * numWidthPoints * (numLengthPoints - 1);

          vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                  cell0DIndex + 1,
                                                  cell0DIndex + numLengthPoints + 1,
                                                  cell0DIndex + numLengthPoints };
          vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                               cell1DVerticalIndex + 1,
                                               cell1DHorizontalIndex + (numLengthPoints - 1),
                                               cell1DVerticalIndex
                                             };

          mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
          mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

          mesh.Cell2DSetState(cell2DIndex, true);

          const unsigned int marker = 21 * (w == 0) +
                                      22 * (w == (numWidthPoints - 1));


          mesh.Cell2DSetMarker(cell2DIndex, marker);

          cell2DIndex++;
        }
      }
    }


    // create length x width cell2Ds
    for (unsigned int w = 0; w < numWidthPoints - 1; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints - 1; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;
          const unsigned int cell1DHorizontalIndex = l + h * (numLengthPoints - 1) + w * (numLengthPoints -1) * numHeightPoints;
          const unsigned int cell1DVerticalIndex = cell0DIndex
                                                   + numHeightPoints * numWidthPoints * (numLengthPoints - 1)
                                                   + numLengthPoints * numWidthPoints * (numHeightPoints - 1);

          vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                  cell0DIndex + 1,
                                                  cell0DIndex + numLengthPoints * numHeightPoints + 1,
                                                  cell0DIndex + numLengthPoints * numHeightPoints };

          vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                               cell1DVerticalIndex + 1,
                                               cell1DHorizontalIndex + (numLengthPoints -1) * numHeightPoints,
                                               cell1DVerticalIndex
                                             };

          mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
          mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

          mesh.Cell2DSetState(cell2DIndex, true);

          const unsigned int marker = 25 * (h == 0) +
                                      26 * (h == (numHeightPoints -1));


          mesh.Cell2DSetMarker(cell2DIndex, marker);

          cell2DIndex++;
        }
      }
    }


    // create Height x width cell2Ds
    for (unsigned int w = 0; w < numWidthPoints - 1; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;

          const unsigned int cell1DHorizontalIndex = l + h * numLengthPoints + w * numLengthPoints * (numHeightPoints - 1)
                                                     + numHeightPoints * numWidthPoints * (numLengthPoints - 1);

          const unsigned int cell1DVerticalIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints
                                                   + numHeightPoints * numWidthPoints * (numLengthPoints - 1)
                                                   + numLengthPoints * numWidthPoints * (numHeightPoints - 1);

          vector<unsigned int> cell2DVertices = { cell0DIndex,
                                                  cell0DIndex + numLengthPoints,
                                                  cell0DIndex + numLengthPoints * numHeightPoints + numLengthPoints ,
                                                  cell0DIndex + numLengthPoints * numHeightPoints
                                                };

          vector<unsigned int> cell2DEdges = { cell1DHorizontalIndex,
                                               cell1DVerticalIndex + numLengthPoints,
                                               cell1DHorizontalIndex + numLengthPoints * (numHeightPoints - 1),
                                               cell1DVerticalIndex
                                             };

          mesh.Cell2DAddVertices(cell2DIndex, cell2DVertices);
          mesh.Cell2DAddEdges(cell2DIndex, cell2DEdges);

          mesh.Cell2DSetState(cell2DIndex, true);

          const unsigned int marker = 23 * (l == 0) +
                                      24 * (l == (numLengthPoints - 1));

          mesh.Cell2DSetMarker(cell2DIndex, marker);

          cell2DIndex++;
        }
      }
    }



    // create cell3Ds
    unsigned int cell3DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints - 1; w++)
    {
      for (unsigned int h = 0; h < numHeightPoints - 1; h++)
      {
        for (unsigned int l = 0; l < numLengthPoints - 1; l++)
        {
          const unsigned int cell0DIndex = l + h * numLengthPoints + w * numLengthPoints * numHeightPoints;

          const unsigned int cell1DHorizontalIndex = l + h * (numLengthPoints - 1) + w * (numLengthPoints -1) * numHeightPoints;

          const unsigned int cell1DVerticalIndex = l + h * numLengthPoints + w * numLengthPoints * (numHeightPoints - 1)
                                                   + numHeightPoints * numWidthPoints * (numLengthPoints - 1);

          const unsigned int cell1DOrthogonalIndex = cell0DIndex
                                                     + numHeightPoints * numWidthPoints * (numLengthPoints - 1)
                                                     + numLengthPoints * numWidthPoints * (numHeightPoints - 1);

          const unsigned int cell2DHorizontalIndex = l + h * (numLengthPoints - 1) + w * (numLengthPoints -1) * (numHeightPoints - 1);

          const unsigned int cell2DVerticalIndex = l + h * (numLengthPoints - 1) + w * (numLengthPoints - 1) * numHeightPoints
                                                   + (numHeightPoints - 1) * numWidthPoints * (numLengthPoints - 1);

          const unsigned int cell2DOrthogonalIndex = l + h * numLengthPoints + w * numLengthPoints * (numHeightPoints - 1)
                                                     + (numHeightPoints - 1) * numWidthPoints * (numLengthPoints - 1)
                                                     + (numLengthPoints - 1) * (numWidthPoints - 1) * numHeightPoints;

          vector<unsigned int> cell3DVertices = { cell0DIndex,
                                                  cell0DIndex + 1,
                                                  cell0DIndex + numLengthPoints + 1,
                                                  cell0DIndex + numLengthPoints,
                                                  cell0DIndex + numLengthPoints * numHeightPoints,
                                                  cell0DIndex + 1 + numLengthPoints * numHeightPoints,
                                                  cell0DIndex + numLengthPoints + 1 + numLengthPoints * numHeightPoints,
                                                  cell0DIndex + numLengthPoints + numLengthPoints * numHeightPoints
                                                };

          vector<unsigned int> cell3DEdges = { cell1DHorizontalIndex,
                                               cell1DHorizontalIndex + (numLengthPoints - 1),
                                               cell1DHorizontalIndex + (numLengthPoints - 1) * numHeightPoints ,
                                               cell1DHorizontalIndex + (numLengthPoints - 1) + (numLengthPoints -1) * numHeightPoints,
                                               cell1DVerticalIndex,
                                               cell1DVerticalIndex + 1,
                                               cell1DVerticalIndex + numLengthPoints * (numHeightPoints - 1),
                                               cell1DVerticalIndex + numLengthPoints * (numHeightPoints - 1) + 1,
                                               cell1DOrthogonalIndex,
                                               cell1DOrthogonalIndex + 1,
                                               cell1DOrthogonalIndex + numLengthPoints,
                                               cell1DOrthogonalIndex + numLengthPoints + 1,
                                             };

          vector<unsigned int> cell3DFaces = { cell2DHorizontalIndex,
                                               cell2DHorizontalIndex + (numLengthPoints - 1) * (numHeightPoints - 1),
                                               cell2DVerticalIndex,
                                               cell2DVerticalIndex + (numLengthPoints - 1),
                                               cell2DOrthogonalIndex,
                                               cell2DOrthogonalIndex + 1
                                             };

          mesh.Cell3DAddVertices(cell3DIndex, cell3DVertices);
          mesh.Cell3DAddEdges(cell3DIndex, cell3DEdges);
          mesh.Cell3DAddFaces(cell3DIndex, cell3DFaces);

          mesh.Cell3DSetState(cell3DIndex, true);

          mesh.Cell3DSetMarker(cell3DIndex, 0);

          cell3DIndex++;
        }
      }
    }




  }
  // ***************************************************************************
  void MeshUtilities::CreateTriangularMesh(const Eigen::MatrixXd& polygonVertices,
                                           const double& maxTriangleArea,
                                           IMeshDAO& mesh,
                                           const string& options) const
  {
    TriangleInterface triangleInterface;

    triangleInterface.CreateMesh(polygonVertices,
                                 maxTriangleArea,
                                 mesh,
                                 options);
  }
  // ***************************************************************************
  void MeshUtilities::CreateTetrahedralMesh(const Eigen::MatrixXd& polyhedronVertices,
                                            const Eigen::MatrixXi& polyhedronEdges,
                                            const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                            const double& maxTetrahedronVolume,
                                            IMeshDAO& mesh,
                                            const string& options) const
  {
    TetgenInterface tetgenInterface;

    tetgenInterface.CreateMesh(polyhedronVertices,
                               polyhedronEdges,
                               polyhedronFaces,
                               maxTetrahedronVolume,
                               mesh,
                               options);
  }
  // ***************************************************************************
}
