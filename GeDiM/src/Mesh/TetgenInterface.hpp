#ifndef __MESHCREATOR3DTETGEN_H
#define __MESHCREATOR3DTETGEN_H

#include "MeshCreator.hpp"
#include "tetgen.h"

using namespace std;
using namespace Eigen;
using namespace MainApplication;

namespace GeDiM
{
  class MeshCreator3DTetgen : public MeshCreator
  {
    protected:
      tetgenio* inputMeshPointer;
      tetgenio* outputMeshPointer;

      string tetgenOptions;

      // http://wias-berlin.de/software/tetgen/files/tetcall.cxx
      void SetInputMesh(tetgenio* _inputMeshPointer) { inputMeshPointer = _inputMeshPointer; }
      void SetOutputMesh(tetgenio* _outputMeshPointer) { outputMeshPointer = _outputMeshPointer; }
      void SetTetgenOptions(const string& _tetgenOptions) { tetgenOptions = _tetgenOptions; }

      Output::ExitCodes CreateTetgenInput(const Polyhedron& domain);
      Output::ExitCodes CreateTetgenOutput(const Polyhedron& domain);

    public:
      MeshCreator3DTetgen();
      virtual ~MeshCreator3DTetgen();

      Output::ExitCodes CreateMesh(const IGeometricObject& domain,
                                   IMesh& mesh);

      Output::ExitCodes ExportTetgenMesh(const string& nameFolder, const string& nameFile) const;
  };

}

#endif // __MESHCREATOR3DTETGEN_H
