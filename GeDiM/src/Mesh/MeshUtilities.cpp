#include "MeshUtilities.hpp"

#include "MeshDAOExporterToCsv.hpp"
#include "StringsUtilities.hpp"
#include "VTKUtilities.hpp"
#include "UCDUtilities.hpp"
#include "OpenVolumeMeshInterface.hpp"
#include "ObjectFileFormatInterface.hpp"

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
  MeshUtilities::FilterMeshData MeshUtilities::FilterActiveMesh(const IMeshDAO& mesh) const
  {
    std::list<unsigned int> active_cell0Ds;
    std::list<unsigned int> active_cell1Ds;
    std::list<unsigned int> active_cell2Ds;
    std::list<unsigned int> active_cell3Ds;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
    {
      if (mesh.Cell0DIsActive(c))
        active_cell0Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
      if (mesh.Cell1DIsActive(c))
        active_cell1Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (mesh.Cell2DIsActive(c))
        active_cell2Ds.push_back(c);
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
      if (mesh.Cell3DIsActive(c))
        active_cell3Ds.push_back(c);
    }

    return
    {
      std::vector<unsigned int>(active_cell0Ds.begin(),
                                active_cell0Ds.end()),
          std::vector<unsigned int>(active_cell1Ds.begin(),
                                    active_cell1Ds.end()),
          std::vector<unsigned int>(active_cell2Ds.begin(),
                                    active_cell2Ds.end()),
          std::vector<unsigned int>(active_cell3Ds.begin(),
                                    active_cell3Ds.end())
    };
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
  void MeshUtilities::ImportObjectFileFormat(const std::string& offFilePath,
                                             IMeshDAO& mesh) const
  {
    ObjectFileFormatInterface objectFileFormatInterface;
    objectFileFormatInterface.ImportMeshFromFile(offFilePath,
                                                 *this,
                                                 mesh);
  }
  // ***************************************************************************
  void MeshUtilities::ExportMeshToObjectFileFormat(const IMeshDAO& mesh,
                                                   const std::string& offFilePath) const
  {
    ObjectFileFormatInterface objectFileFormatInterface;
    objectFileFormatInterface.ExportMeshToFile(mesh,
                                               offFilePath);
  }
  // ***************************************************************************
  void MeshUtilities::ExportMeshToVTU(const IMeshDAO& mesh,
                                      const string& exportFolder,
                                      const string& fileName,
                                      const bool& separateFile) const
  {
    string cell0DsFolder = exportFolder;
    string cell1DsFolder = exportFolder;
    string cell2DsFolder = exportFolder;
    string cell3DsFolder = exportFolder;

    if (separateFile)
    {
      Gedim::Output::CreateFolder(exportFolder + "/Cell0Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell1Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell2Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell3Ds");

      cell0DsFolder = exportFolder + "/Cell0Ds";
      cell1DsFolder = exportFolder + "/Cell1Ds";
      cell2DsFolder = exportFolder + "/Cell2Ds";
      cell3DsFolder = exportFolder + "/Cell3Ds";
    }

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
      vtpUtilities.Export(cell0DsFolder + "/" + fileName + "_Cell0Ds.vtu");
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
      vtpUtilities.Export(cell1DsFolder + "/" + fileName + "_Cell1Ds.vtu");
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
      vtpUtilities.Export(cell2DsFolder + "/" + fileName + "_Cell2Ds.vtu");
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
      vtpUtilities.Export(cell3DsFolder + "/" + fileName + "_Cell3Ds.vtu");
    }
  }
  // ***************************************************************************
  void MeshUtilities::ExportMeshToUCD(const IMeshDAO& mesh,
                                      const std::string& exportFolder,
                                      const std::string& fileName,
                                      const bool& separateFile) const
  {
    string cell0DsFolder = exportFolder;
    string cell1DsFolder = exportFolder;
    string cell2DsFolder = exportFolder;
    string cell3DsFolder = exportFolder;

    if (separateFile)
    {
      Gedim::Output::CreateFolder(exportFolder + "/Cell0Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell1Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell2Ds");
      Gedim::Output::CreateFolder(exportFolder + "/Cell3Ds");

      cell0DsFolder = exportFolder + "/Cell0Ds";
      cell1DsFolder = exportFolder + "/Cell1Ds";
      cell2DsFolder = exportFolder + "/Cell2Ds";
      cell3DsFolder = exportFolder + "/Cell3Ds";
    }

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

      vector<UCDProperty<double>> properties(3 + mesh.Cell0DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell0DNumberDoubleProperties());

      properties[0] = {
        "Id",
        "",
        static_cast<unsigned int>(id.size()),
        1,
        id.data()
      };
      properties[1] = {
        "Marker",
        "",
        static_cast<unsigned int>(marker.size()),
        1,
        marker.data()
      };
      properties[2] = {
        "Active",
        "",
        static_cast<unsigned int>(active.size()),
        1,
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

        properties[3 + p] =
        {
          mesh.Cell0DDoublePropertyId(p),
          "",
          static_cast<unsigned int>(propertyValues[p].size()),
          1,
          propertyValues[p].data()
        };
      }

      Gedim::UCDUtilities exporter;
      exporter.ExportPoints(cell0DsFolder + "/" + fileName + "_Cell0Ds.inp",
                            mesh.Cell0DsCoordinates(),
                            properties);
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

      vector<UCDProperty<double>> properties(3 + mesh.Cell1DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell1DNumberDoubleProperties());

      properties[0] = {
        "Id",
        "",
        static_cast<unsigned int>(id.size()),
        1,
        id.data()
      };
      properties[1] = {
        "Marker",
        "",
        static_cast<unsigned int>(marker.size()),
        1,
        marker.data()
      };
      properties[2] = {
        "Active",
        "",
        static_cast<unsigned int>(active.size()),
        1,
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

        properties[3 + p] =
        {
          mesh.Cell1DDoublePropertyId(p),
          "",
          static_cast<unsigned int>(propertyValues[p].size()),
          1,
          propertyValues[p].data()
        };
      }


      Gedim::UCDUtilities exporter;
      exporter.ExportSegments(cell1DsFolder + "/" + fileName + "_Cell1Ds.inp",
                              mesh.Cell0DsCoordinates(),
                              mesh.Cell1DsExtremes(),
                              {},
                              properties);
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

      vector<UCDProperty<double>> properties(3 + mesh.Cell2DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell2DNumberDoubleProperties());

      properties[0] = {
        "Id",
        "",
        static_cast<unsigned int>(id.size()),
        1,
        id.data()
      };
      properties[1] = {
        "Marker",
        "",
        static_cast<unsigned int>(marker.size()),
        1,
        marker.data()
      };
      properties[2] = {
        "Active",
        "",
        static_cast<unsigned int>(active.size()),
        1,
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

        properties[3 + p] =
        {
          mesh.Cell2DDoublePropertyId(p),
          "",
          static_cast<unsigned int>(propertyValues[p].size()),
          1,
          propertyValues[p].data()
        };
      }

      Gedim::UCDUtilities exporter;
      exporter.ExportPolygons(cell2DsFolder + "/" + fileName + "_Cell2Ds.inp",
                              mesh.Cell0DsCoordinates(),
                              mesh.Cell2DsVertices(),
                              {},
                              properties);
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

      vector<UCDProperty<double>> properties(3 + mesh.Cell3DNumberDoubleProperties());
      vector<vector<double>> propertyValues(mesh.Cell3DNumberDoubleProperties());

      properties[0] = {
        "Id",
        "",
        static_cast<unsigned int>(id.size()),
        1,
        id.data()
      };
      properties[1] = {
        "Marker",
        "",
        static_cast<unsigned int>(marker.size()),
        1,
        marker.data()
      };
      properties[2] = {
        "Active",
        "",
        static_cast<unsigned int>(active.size()),
        1,
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

        properties[3 + p] =
        {
          mesh.Cell3DDoublePropertyId(p),
          "",
          static_cast<unsigned int>(propertyValues[p].size()),
          1,
          propertyValues[p].data()
        };
      }

      Gedim::UCDUtilities exporter;
      exporter.ExportPolyhedra(cell3DsFolder + "/" + fileName + "_Cell3Ds.inp",
                               mesh.Cell0DsCoordinates(),
                               mesh.Cell3DsVertices(),
                               {},
                               properties);
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
  void MeshUtilities::ExportCell3DToVTU(const GeometryUtilities& geometryUtilities,
                                        const IMeshDAO& mesh,
                                        const unsigned int& cell3DIndex,
                                        const Eigen::MatrixXd& cell3DVertices,
                                        const std::vector<Eigen::MatrixXd>& cell3DTetrahedrons,
                                        const std::vector<std::vector<Eigen::Matrix3d> >& cell3DFaces3DTriangulations,
                                        const double& cell3DVolume,
                                        const Eigen::Vector3d& cell3DCentroid,
                                        const std::vector<Eigen::Vector3d>& cell3DFacesTranslation,
                                        const std::vector<Eigen::Matrix3d>& cell3DFacesRotationMatrix,
                                        const std::vector<double>& cell3DFacesArea,
                                        const std::vector<Eigen::MatrixXd>& cell3DFaces2DVertices,
                                        const std::vector<Eigen::MatrixXd>& cell3DFaces3DVertices,
                                        const std::vector<Eigen::VectorXd>& cell3DFacesEdgeLengths,
                                        const std::vector<std::vector<bool>>& cell3DFacesEdgeDirections,
                                        const std::vector<Eigen::MatrixXd>& cell3DFacesEdges2DTangent,
                                        const std::vector<Eigen::MatrixXd>& cell3DFacesEdges2DNormal,
                                        const std::vector<Eigen::Vector3d>& cell3DFacesNormals,
                                        const std::vector<bool>& cell3DFacesNormalDirections,
                                        const std::vector<Eigen::Vector3d>& cell3DFaces2DCentroids,
                                        const std::string& exportFolder) const
  {
    {
      const GeometryUtilities::Polyhedron cell3D = MeshCell3DToPolyhedron(mesh,
                                                                          cell3DIndex);

      {
        Gedim::VTKUtilities vtpUtilities;
        const vector<double> id(1, cell3DIndex);
        const vector<double> volume(1, cell3DVolume);

        // Export cell3D vertices
        vtpUtilities.AddPoints(cell3DVertices);

        vtpUtilities.Export(exportFolder + "/"
                            + "Cell3D_" + to_string(cell3DIndex) +
                            "_Vertices" +
                            ".vtu");
      }

      {
        Gedim::VTKUtilities vtpUtilities;
        const vector<double> id(1, cell3DIndex);
        const vector<double> volume(1, cell3DVolume);

        // Export cell3D edges
        vtpUtilities.AddSegments(cell3D.Vertices,
                                 cell3D.Edges);

        vtpUtilities.Export(exportFolder + "/"
                            + "Cell3D_" + to_string(cell3DIndex) +
                            "_Edges" +
                            ".vtu");
      }

      {
        Gedim::VTKUtilities vtpUtilities;
        const vector<double> id(1, cell3DIndex);
        const vector<double> volume(1, cell3DVolume);

        // Export cell3D
        vtpUtilities.AddPolyhedron(cell3D.Vertices,
                                   cell3D.Edges,
                                   cell3D.Faces,
                                   {
                                     {
                                       "Id",
                                       Gedim::VTPProperty::Formats::Cells,
                                       static_cast<unsigned int>(id.size()),
                                       id.data()
                                     },
                                     {
                                       "Volume",
                                       Gedim::VTPProperty::Formats::Cells,
                                       static_cast<unsigned int>(volume.size()),
                                       volume.data()
                                     }
                                   });

        vtpUtilities.Export(exportFolder + "/"
                            + "Cell3D_" + to_string(cell3DIndex) +
                            ".vtu");
      }
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      // Export cell3D tetrahedrons
      for (unsigned int t = 0; t < cell3DTetrahedrons.size(); t++)
      {
        const GeometryUtilities::Polyhedron tetra = geometryUtilities.CreateTetrahedronWithVertices(cell3DTetrahedrons[t].col(0),
                                                                                                    cell3DTetrahedrons[t].col(1),
                                                                                                    cell3DTetrahedrons[t].col(2),
                                                                                                    cell3DTetrahedrons[t].col(3));

        const vector<double> id(1, t);
        vtpUtilities.AddPolyhedron(tetra.Vertices,
                                   tetra.Edges,
                                   tetra.Faces,
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
                          "Cell3D_" + to_string(cell3DIndex) +
                          "_Tetra" + ".vtu");
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      unsigned int totT = 0;
      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
      {
        const std::vector<Eigen::Matrix3d>& cell2DTriangulations = cell3DFaces3DTriangulations[f];

        // Export cell2D triangulation
        for (unsigned int t = 0; t < cell2DTriangulations.size(); t++)
        {
          const vector<double> face_local_id(1, f);
          const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
          const vector<double> triangle_local_id(1, t);
          const vector<double> triangle_global_id(1, totT++);
          vtpUtilities.AddPolygon(cell2DTriangulations[t],
                                  {
                                    {
                                      "face_global_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(face_global_id.size()),
                                      face_global_id.data()
                                    },
                                    {
                                      "face_local_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(face_local_id.size()),
                                      face_local_id.data()
                                    },
                                    {
                                      "triangle_global_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(triangle_global_id.size()),
                                      triangle_global_id.data()
                                    },
                                    {
                                      "triangle_local_id",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(triangle_local_id.size()),
                                      triangle_local_id.data()
                                    }
                                  });
        }
      }

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell3D_" + to_string(cell3DIndex) +
                          "_Faces_Triangles" + ".vtu");
    }

    {
      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
      {
        const Eigen::MatrixXd& face2DVertices = cell3DFaces2DVertices[f];
        const Eigen::MatrixXd& face3DVertices = cell3DFaces3DVertices[f];
        const Eigen::VectorXd& faceEdgeLengths = cell3DFacesEdgeLengths[f];
        const std::vector<bool>& faceEdgeDirections = cell3DFacesEdgeDirections[f];
        const Eigen::Matrix3d& faceRotation = cell3DFacesRotationMatrix[f];
        const Eigen::Vector3d& faceTranslation = cell3DFacesTranslation[f];
        const Eigen::MatrixXd& faceEdgesNormal = cell3DFacesEdges2DNormal[f];
        const Eigen::MatrixXd& faceEdgesTangent = cell3DFacesEdges2DTangent[f];

        {
          Gedim::VTKUtilities vtpUtilities;
          vtpUtilities.AddPoints(face3DVertices);
          vtpUtilities.Export(exportFolder + "/" +
                              "Cell3D_" + to_string(cell3DIndex) +
                              "_Face_" + to_string(f) +
                              "_Vertices" + ".vtu");
        }


        {
          Gedim::VTKUtilities vtpUtilities;
          for (unsigned int e = 0; e < faceEdgeDirections.size(); e++)
          {
            // Export cell2D centroids
            const vector<double> edge_length(1, faceEdgeLengths[e]);
            const vector<double> edge_direction(1, faceEdgeDirections[e] ? 1.0 : -1.0);
            const vector<double> face_local_id(1, f);
            const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));

            const Eigen::Vector3d edgeO = face3DVertices.col(e);
            const Eigen::Vector3d edgeE = face3DVertices.col((e + 1) % faceEdgeDirections.size());

            vtpUtilities.AddSegment(edgeO,
                                    edgeE,
                                    {
                                      {
                                        "edge_length",
                                        Gedim::VTPProperty::Formats::Cells,
                                        static_cast<unsigned int>(edge_length.size()),
                                        edge_length.data()
                                      },
                                      {
                                        "edge_direction",
                                        Gedim::VTPProperty::Formats::Cells,
                                        static_cast<unsigned int>(edge_direction.size()),
                                        edge_direction.data()
                                      },
                                      {
                                        "face_global_id",
                                        Gedim::VTPProperty::Formats::Cells,
                                        static_cast<unsigned int>(face_global_id.size()),
                                        face_global_id.data()
                                      },
                                      {
                                        "face_local_id",
                                        Gedim::VTPProperty::Formats::Cells,
                                        static_cast<unsigned int>(face_local_id.size()),
                                        face_local_id.data()
                                      }
                                    });
          }

          vtpUtilities.Export(exportFolder + "/" +
                              "Cell3D_" + to_string(cell3DIndex) +
                              "_Face_" + to_string(f) +
                              "_Edges" + ".vtu");
        }

        {
          Gedim::VTKUtilities vtpUtilitiesNormal, vtpUtilitiesTangent;
          for (unsigned int e = 0; e < faceEdgeDirections.size(); e++)
          {
            // Export cell2D centroids
            const vector<double> edge_direction(1, faceEdgeDirections[e] ? 1.0 : -1.0);
            const vector<double> face_local_id(1, f);
            const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));

            const double sizeNormal = std::max(std::max(std::sqrt(cell3DFacesArea[f]),
                                                        std::cbrt(cell3DVolume)),
                                               faceEdgeLengths[e]);

            const Eigen::Vector3d edgeO2D = face2DVertices.col(e);
            const Eigen::Vector3d edgeE2D = face2DVertices.col((e + 1) % faceEdgeDirections.size());

            const Eigen::Vector3d edgeNormal2DO = 0.5 * (edgeO2D + edgeE2D);
            const Eigen::Vector3d edgeNormal2DE = edgeNormal2DO + sizeNormal * faceEdgesNormal.col(e);

            const Eigen::Vector3d edgeTangent2DO = edgeO2D;
            const Eigen::Vector3d edgeTangent2DE = edgeO2D + faceEdgesTangent.col(e);

            const Eigen::Vector3d edgeNormal3DO = geometryUtilities.RotatePointsFrom2DTo3D(edgeNormal2DO,
                                                                                           faceRotation,
                                                                                           faceTranslation);
            const Eigen::Vector3d edgeNormal3DE = geometryUtilities.RotatePointsFrom2DTo3D(edgeNormal2DE,
                                                                                           faceRotation,
                                                                                           faceTranslation);

            const Eigen::Vector3d edgeTangent3DO = geometryUtilities.RotatePointsFrom2DTo3D(edgeTangent2DO,
                                                                                            faceRotation,
                                                                                            faceTranslation);
            const Eigen::Vector3d edgeTangent3DE = geometryUtilities.RotatePointsFrom2DTo3D(edgeTangent2DE,
                                                                                            faceRotation,
                                                                                            faceTranslation);

            vtpUtilitiesNormal.AddSegment(edgeTangent3DO,
                                          edgeTangent3DE,
                                          {
                                            {
                                              "edge_direction",
                                              Gedim::VTPProperty::Formats::Cells,
                                              static_cast<unsigned int>(edge_direction.size()),
                                              edge_direction.data()
                                            }
                                          });
            vtpUtilitiesTangent.AddSegment(edgeNormal3DO,
                                           edgeNormal3DE);
          }

          vtpUtilitiesNormal.Export(exportFolder + "/" +
                                    "Cell3D_" + to_string(cell3DIndex) +
                                    "_Face_" + to_string(f) +
                                    "_EdgeTangents" + ".vtu");
          vtpUtilitiesTangent.Export(exportFolder + "/" +
                                     "Cell3D_" + to_string(cell3DIndex) +
                                     "_Face_" + to_string(f) +
                                     "_EdgeNormals" + ".vtu");
        }
      }
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
      {
        const Eigen::Vector3d& faceCentroid = geometryUtilities.RotatePointsFrom2DTo3D(cell3DFaces2DCentroids[f],
                                                                                       cell3DFacesRotationMatrix[f],
                                                                                       cell3DFacesTranslation[f]);

        // Export cell2D centroids
        const vector<double> face_local_id(1, f);
        const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
        vtpUtilities.AddPoint(faceCentroid,
                              {
                                {
                                  "face_global_id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(face_global_id.size()),
                                  face_global_id.data()
                                },
                                {
                                  "face_local_id",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(face_local_id.size()),
                                  face_local_id.data()
                                }
                              });
      }

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell3D_" + to_string(cell3DIndex) +
                          "_Faces_Centroid" + ".vtu");
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
      {
        const Eigen::Vector3d& faceCentroid = geometryUtilities.RotatePointsFrom2DTo3D(cell3DFaces2DCentroids[f],
                                                                                       cell3DFacesRotationMatrix[f],
                                                                                       cell3DFacesTranslation[f]);
        const double sizeNormal = std::max(std::sqrt(cell3DFacesArea[f]),
                                           std::cbrt(cell3DVolume));
        const Eigen::Vector3d& faceNormal = cell3DFacesNormals[f];
        const double faceDirection = cell3DFacesNormalDirections[f] ? 1.0 : -1.0;

        // Export cell2D centroids
        const vector<double> face_local_id(1, f);
        const vector<double> face_global_id(1, mesh.Cell3DFace(cell3DIndex, f));
        vtpUtilities.AddSegment(faceCentroid,
                                faceCentroid + faceDirection * sizeNormal * faceNormal,
                                {
                                  {
                                    "face_global_id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(face_global_id.size()),
                                    face_global_id.data()
                                  },
                                  {
                                    "face_local_id",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(face_local_id.size()),
                                    face_local_id.data()
                                  }
                                });
      }

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell3D_" + to_string(cell3DIndex) +
                          "_Faces_Normal" + ".vtu");
    }

    {
      Gedim::VTKUtilities vtpUtilities;

      // Export cell3D centroid
      vtpUtilities.AddPoint(cell3DCentroid);

      vtpUtilities.Export(exportFolder + "/" +
                          "Cell3D_" + to_string(cell3DIndex) +
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
  std::vector<unsigned int> MeshUtilities::MarkCells(const std::function<Eigen::VectorXi(const Eigen::MatrixXd&)>& marking_function,
                                                     const std::vector<Eigen::MatrixXd>& cells_points,
                                                     const unsigned int default_mark) const
  {
    std::vector<unsigned int> cells_marks(cells_points.size());
    for (unsigned int c = 0; c < cells_points.size(); c++)
    {
      if (cells_points.at(c).size() == 0)
      {
        cells_marks[c] = default_mark;
        continue;
      }

      const auto cell_marks = marking_function(cells_points.at(c));

      cells_marks[c] = cell_marks.maxCoeff();
    }

    return cells_marks;
  }
  // ***************************************************************************
}
