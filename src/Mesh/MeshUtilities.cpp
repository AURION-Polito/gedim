#include "MeshUtilities.hpp"
#include <iostream>
#include <fstream>

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
  void MeshUtilities::ExtractActiveMesh(const IMeshDAO& mesh,
                                        IMeshDAO& extractedMesh,
                                        ExtractActiveMeshData& extractionData)
  {
    // extract active Cell0Ds
    unsigned int numNewCell0Ds = 0;
    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
      if (!mesh.Cell0DIsActive(c))
        continue;

      extractionData.NewCell0DToOldCell0D.insert(pair<unsigned int,
                                                 unsigned int>(numNewCell0Ds, c));
      extractionData.OldCell0DToNewCell0D.insert(pair<unsigned int,
                                                 unsigned int>(c, numNewCell0Ds));
      numNewCell0Ds++;
    }
  }
  // ***************************************************************************
}
