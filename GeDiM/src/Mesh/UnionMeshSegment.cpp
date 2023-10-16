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
    { cerr<< scientific << "{ Key: " << p.first<< "; Value: T: "<< (unsigned int)p.second.Type<< " I: "<< p.second.MeshIndices<< " }\n"; });
    cerr<< "Segments:"<< endl;
    for_each(unionMesh.Segments.begin(), unionMesh.Segments.end(), [](
             const UnionMeshSegment::UnionMesh::UnionMeshSegment& p)
    { cerr<< scientific << "{ P: "<< p.Points<< " T: "<< (unsigned int)p.Type<< " I: "<< p.MeshIndices<< " }\n"; });
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
      if (!_geometryUtilities.IsValuePositive(abs(it->first - curvilinearCoordinate),
                                              _geometryUtilities.Tolerance1D()))
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
  void UnionMeshSegment::CreateUnionPoints(const vector<double>& curvilinearCoordinatesMeshOne,
                                           const vector<double>& curvilinearCoordinatesMeshTwo,
                                           UnionMesh& result)
  {
    // Insert first mesh in union
    for (unsigned int i = 0; i < curvilinearCoordinatesMeshOne.size(); i++)
    {
      const double& curvilinearCoordinate = curvilinearCoordinatesMeshOne[i];
      bool found;
      UnionMesh::UnionMeshPoint& intersection = InsertNewIntersection(curvilinearCoordinate,
                                                                      result,
                                                                      found);
      intersection.Type = UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::First;
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
        intersection.Type = UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both;
        intersection.MeshIndices[1] = i;
      }
      else
      {
        intersection.Type = UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Second;
        intersection.MeshIndices.resize(2);
        intersection.MeshIndices[1] = i;
      }
    }
  }
  // ***************************************************************************
  void UnionMeshSegment::CreateUnionSegments(const vector<double>& curvilinearCoordinatesMeshOne,
                                             const vector<double>& curvilinearCoordinatesMeshTwo,
                                             UnionMesh& result)
  {
    result.Segments.resize(result.Points.size() - 1);
    vector<vector<int>> meshIndices(result.Points.size() - 1);

    map<double, UnionMesh::UnionMeshPoint>::const_iterator itPoint = result.Points.begin();
    map<double, UnionMesh::UnionMeshPoint>::const_iterator itPointNext = result.Points.begin();
    itPointNext++;
    for (unsigned int p = 0; p < result.Segments.size(); p++)
    {
      const double& curvilinearCoordinatePoint = itPoint->first;
      const double& curvilinearCoordinatePointNext = itPointNext->first;
      const UnionMesh::UnionMeshPoint& intersectionPoint = itPoint->second;

      // fill origin and end of segment
      UnionMesh::UnionMeshSegment& meshSegment = result.Segments[p];
      meshSegment.Points.resize(2);
      meshSegment.Points[0] = curvilinearCoordinatePoint;
      meshSegment.Points[1] = curvilinearCoordinatePointNext;

      meshIndices[p].resize(2);
      meshIndices[p][0] = p == 0 ? -1 : meshIndices[p - 1][0];
      meshIndices[p][1] = p == 0 ? -1 : meshIndices[p - 1][1];

      switch (intersectionPoint.Type)
      {
        case Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Types::First:
          meshIndices[p][0]++;
          break;
        case Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Types::Second:
          meshIndices[p][1]++;
          break;
        case Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Types::Both:
          meshIndices[p][0]++;
          meshIndices[p][1]++;
          break;
        default:
          throw runtime_error("Unmanaged intersectionPoint.Type");
      }

      if ((meshIndices[p][0] + 1) >= static_cast<int>(curvilinearCoordinatesMeshOne.size()))
        meshIndices[p][0] = -1;
      if ((meshIndices[p][1] + 1) >= static_cast<int>(curvilinearCoordinatesMeshTwo.size()))
        meshIndices[p][1] = -1;

      meshSegment.MeshIndices.resize(2);
      meshSegment.MeshIndices[0] = meshIndices[p][0] >= 0 ? meshIndices[p][0] : 0;
      meshSegment.MeshIndices[1] = meshIndices[p][1] >= 0 ? meshIndices[p][1] : 0;

      Output::Assert(meshIndices[p][0] != -1 || meshIndices[p][1] != -1);

      if (meshIndices[p][0] >= 0 && meshIndices[p][1] >= 0)
        meshSegment.Type = Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both;
      else if (meshIndices[p][0] >= 0 && meshIndices[p][1] == -1)
        meshSegment.Type = Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::First;
      else if (meshIndices[p][0] == -1 && meshIndices[p][1] >= 0)
        meshSegment.Type = Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Second;

      itPoint++;
      itPointNext++;
    }
  }
  // ***************************************************************************
  void UnionMeshSegment::CreateUnionMesh(const vector<double>& curvilinearCoordinatesMeshOne,
                                         const vector<double>& curvilinearCoordinatesMeshTwo,
                                         UnionMeshSegment::UnionMesh& result)
  {
    CreateUnionPoints(curvilinearCoordinatesMeshOne,
                      curvilinearCoordinatesMeshTwo,
                      result);

    CreateUnionSegments(curvilinearCoordinatesMeshOne,
                        curvilinearCoordinatesMeshTwo,
                        result);
  }
  // ***************************************************************************
}
