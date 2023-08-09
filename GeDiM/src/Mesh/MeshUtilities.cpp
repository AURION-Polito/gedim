#include "MeshUtilities.hpp"

#include "MeshDAOExporterToCsv.hpp"
#include "VTKUtilities.hpp"
#include "OpenVolumeMeshInterface.hpp"

using namespace std;
using namespace Eigen;

namespace Gedim
{
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
  void MeshUtilities::ExportMeshToCsv(const IMeshDAO& mesh,
                                      const char& separator,
                                      const std::string& exportFolderPath) const
  {
    // export mesh
    Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
    Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
    exportConfiguration.Separator = separator;
    exportConfiguration.Folder = exportFolderPath;
    Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
    exporter.Export(exportConfiguration,
                    mesh);
  }
  // ***************************************************************************
}
