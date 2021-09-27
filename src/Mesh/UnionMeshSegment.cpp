#include "UnionMeshSegment.hpp"
#include "IOUtilities.hpp"

using namespace std;

namespace Gedim
{
  // ***************************************************************************
  UnionMeshSegment::UnionMeshSegment(const Gedim::GeometryUtilities& geometryUtilities) :
    _geometryUtilities(geometryUtilities)
  {
  }
  UnionMeshSegment::~UnionMeshSegment()
  {
  }
  // ***************************************************************************
  void UnionMeshSegment::ToCurvilinearCoordinates(const UnionMeshSegment::UnionMesh& unionMesh,
                                                  vector<double>& curvilinearCoordinates)
  {
    curvilinearCoordinates.reserve(unionMesh.Points.size());
    for (std::map<double,
         UnionMeshSegment::UnionMesh::UnionMeshPoint>::const_iterator it = unionMesh.Points.begin();
         it != unionMesh.Points.end(); it++)
    {
      curvilinearCoordinates.push_back(it->first);
    }
  }
  // ***************************************************************************
  void UnionMeshSegment::ToString(const UnionMeshSegment::UnionMesh& unionMesh)
  {
    cerr.precision(16);
    for_each(unionMesh.Points.begin(), unionMesh.Points.end(), [](const std::pair<double,
             UnionMeshSegment::UnionMesh::UnionMeshPoint>& p)
    { cerr<< scientific << "{ Key: " << p.first<< "; Value: T: "<< p.second.Type<< " I: "<< p.second.MeshIndices<< " }\n"; });
  }
  // ***************************************************************************
  UnionMeshSegment::UnionMesh::UnionMeshPoint& UnionMeshSegment::InsertNewIntersection(const double& curvilinearCoordinate,
                                                                                       UnionMeshSegment::UnionMesh& result,
                                                                                       bool& found)
  {
    double foundCoordinate = -1.0;
    for (std::map<double,
         UnionMesh::UnionMeshPoint>::const_iterator it = result.Points.begin();
         it != result.Points.end(); it++)
    {
      GeometryUtilities::CompareTypes result = _geometryUtilities.Compare1DValues(it->first, curvilinearCoordinate);

      if (result == GeometryUtilities::CompareTypes::Coincident)
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

    result.Points.insert(pair<double,
                         UnionMesh::UnionMeshPoint>(curvilinearCoordinate,
                                                    UnionMesh::UnionMeshPoint()));
    found = false;
    return result.Points[curvilinearCoordinate];
  }
  // ***************************************************************************
  void UnionMeshSegment::CreateUnionMesh(const vector<double>& curvilinearCoordinatesMeshOne,
                                         const vector<double>& curvilinearCoordinatesMeshTwo,
                                         UnionMeshSegment::UnionMesh& result)
  {
    // Insert first mesh in union
    for (unsigned int i = 0; i < curvilinearCoordinatesMeshOne.size(); i++)
    {
      const double& curvilinearCoordinate = curvilinearCoordinatesMeshOne[i];
      bool found;
      UnionMesh::UnionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                      result,
                                                                      found);
      intersection.Type = UnionMesh::UnionMeshPoint::First;
      intersection.MeshIndices.resize(2);
      intersection.MeshIndices[0] = i;
    }

    // Insert second mesh in union checking union
    for (unsigned int i = 0; i < curvilinearCoordinatesMeshTwo.size(); i++)
    {
      const double& curvilinearCoordinate = curvilinearCoordinatesMeshTwo[i];
      bool found;
      UnionMesh::UnionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                      result,
                                                                      found);

      if (found)
      {
        intersection.Type = UnionMesh::UnionMeshPoint::Both;
        intersection.MeshIndices[1] = i;
      }
      else
      {
        intersection.Type = UnionMesh::UnionMeshPoint::Second;
        intersection.MeshIndices.resize(2);
        intersection.MeshIndices[1] = i;
      }
    }
  }
  // ***************************************************************************
}
