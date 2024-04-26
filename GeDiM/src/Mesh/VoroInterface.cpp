#include "VoroInterface.hpp"

#include "Macro.hpp"
#include "Eigen/Eigen"
#include <unistd.h>

#if ENABLE_VORO == 0
#include "CommonUtilities.hpp"
#endif

namespace Gedim
{
// ************************************************************************* //
VoroInterface::VoroInterface(const Gedim::GeometryUtilities& geometryUtilities):
    geometryUtilities(geometryUtilities)
{

}
// ************************************************************************* //
#if ENABLE_VORO == 0
void VoroInterface::GenerateVoronoiTassellations3D(const Eigen::MatrixXd& polyhedronVertices,
                                                   const Eigen::MatrixXi& polyhedronEdges,
                                                   const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                   const unsigned int &numPoints,
                                                   const unsigned int& numIterations,
                                                   Gedim::IMeshDAO& mesh)
{
    Gedim::Utilities::Unused(polyhedronVertices);
    Gedim::Utilities::Unused(polyhedronEdges);
    Gedim::Utilities::Unused(polyhedronFaces);
    Gedim::Utilities::Unused(numPoints);
    Gedim::Utilities::Unused(numIterations);
    Gedim::Utilities::Unused(mesh);
    throw runtime_error("Not active module VORO");
}

void VoroInterface::GenerateVoronoiTassellations2D(const Eigen::MatrixXd& polygonVertices,
                                                   const unsigned int& numPoints,
                                                   const unsigned int& numIterations,
                                                   Gedim::IMeshDAO& mesh)
{
    Gedim::Utilities::Unused(polygonVertices);
    Gedim::Utilities::Unused(numPoints);
    Gedim::Utilities::Unused(numIterations);
    Gedim::Utilities::Unused(mesh);
    throw runtime_error("Not active module VORO");
}
#else
// ************************************************************************* //
bool VoroInterface::InsertNewPoints(Cell0D& cell0D,
                                    list<Cell0D>& cell0Ds)
{
    if(cell0Ds.size() == 0)
    {
        cell0D.id = 0;
        cell0Ds.push_back(cell0D);
        return true;
    }
    else{
        for(auto it = cell0Ds.begin(); it != cell0Ds.end(); it++)
        {
            if( geometryUtilities.IsValueGreater((*it).x, cell0D.x, geometryUtilities.Tolerance1D()))
            {
                cell0D.id = cell0Ds.size();
                cell0Ds.insert(it, cell0D);
                return true;
            }
            else if(geometryUtilities.IsValueZero((*it).x - cell0D.x, geometryUtilities.Tolerance1D()))
            {
                if( geometryUtilities.IsValueGreater((*it).y, cell0D.y, geometryUtilities.Tolerance1D()))
                {
                    cell0D.id = cell0Ds.size();
                    cell0Ds.insert(it, cell0D);
                    return true;
                }
                else if(geometryUtilities.IsValueZero((*it).y - cell0D.y, geometryUtilities.Tolerance1D()))
                {
                    if( geometryUtilities.IsValueGreater((*it).z, cell0D.z, geometryUtilities.Tolerance1D()))
                    {
                        cell0D.id = cell0Ds.size();
                        cell0Ds.insert(it, cell0D);
                        return true;
                    }
                    else if(geometryUtilities.IsValueZero((*it).z - cell0D.z, geometryUtilities.Tolerance1D()))
                    {
                        cell0D.id = (*it).id;
                        return false;
                    }
                }
            }
        }

        cell0D.id = cell0Ds.size();
        cell0Ds.push_back(cell0D);
        return true;
    }
}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations2D(const Eigen::MatrixXd& polygonVertices,
                                                   const unsigned int& numPoints,
                                                   const unsigned int& numIterations,
                                                   Gedim::IMeshDAO& mesh)
{
    const unsigned int numVerticesDomain = polygonVertices.cols();
    const unsigned int numEdgesDomain = numVerticesDomain;

    vector<Eigen::Vector3d> edgeDomainTangent(numEdgesDomain);
    vector<double> edgeDomainTangentSquaredNorm(numEdgesDomain);
    vector<Eigen::Vector3d> edgeDomainLineOrigin(numEdgesDomain);
    for(unsigned int e = 0; e < numEdgesDomain; e++)
    {
        edgeDomainLineOrigin[e] = polygonVertices.col(e);
        const Eigen::Vector3d end = polygonVertices.col((e + 1) % numEdgesDomain);

        Eigen::Vector3d tangent = end - edgeDomainLineOrigin[e];

        edgeDomainTangent[e] = tangent;
        edgeDomainTangentSquaredNorm[e] = tangent.squaredNorm();
    }


    //    vector<unsigned int> map_marker(numFacesDomain, 0);

    const double x_min = polygonVertices.row(0).minCoeff(), x_max = polygonVertices.row(0).maxCoeff();
    const double y_min = polygonVertices.row(1).minCoeff(), y_max = polygonVertices.row(1).maxCoeff();
    const double cvol = (x_max - x_min) * (y_max - y_min);

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = 6, n_z = 6;

    Eigen::MatrixXd VoronoiPoints;
    GenerateRandomPoints(polygonVertices,
                         numPoints,
                         VoronoiPoints);


    // Loop on Voronoi cells
    for(unsigned int it = 0; it < numIterations; it++)
    {
        double vvol = 0.0;
        // Initialize Voronoi mesh
        voro::container con
            (
                x_min, x_max, y_min, y_max, -0.5, 0.5,
                n_x, n_y, n_z,
                false, false, false, 8
                );

        for(unsigned int i=0; i < VoronoiPoints.cols(); i++) {
            con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));
        }

        unsigned int countPoints = 0;

        voro::c_loop_all vl(con);
        int ijk, q; double *pp;
        voro::voronoicell_neighbor c;
        if (vl.start())
        {
            do{
                if (con.compute_cell(c, vl))
                {
                    ijk = vl.ijk;
                    q = vl.q;
                    pp = con.p[ijk] + con.ps*q;

                    // Cell barycenter
                    const double cx = *pp;
                    const double cy = pp[1];
                    const double cz = pp[2];

                    // Get cell vertices
                    std::vector<double> vertices;
                    c.vertices(cx, cy, cz, vertices);

                    // Get cell face cardinality
                    std::vector<int> faceOrders;
                    c.face_orders(faceOrders);

                    // Get face-vertices connectivity
                    std::vector<int> localFaceVertices;
                    c.face_vertices(localFaceVertices);

                    // Get cell neighbors IDs
                    std::vector<int> faceNeighs;
                    c.neighbors(faceNeighs);


                    // Fill mesh connectivity structures
                    // cell vertices
                    const unsigned int nCoords = vertices.size();
                    const unsigned int numVertices = nCoords/3;
                    std::vector<int> cell3DVertices(numVertices);

                    Eigen::MatrixXd cellVertCoordinates = Eigen::MatrixXd::Zero(3, numVertices);

                    unsigned int v = 0;
                    unsigned int i = 0;
                    while (i < nCoords)
                    {
                        double x = vertices[i++];
                        double y = vertices[i++];
                        double z = vertices[i++];

                        if(!geometryUtilities.IsValueZero(z + 0.5,geometryUtilities.Tolerance1D()))
                        {
                            cell3DVertices[v] = -1;
                            v++;
                        }
                        else{
                            cellVertCoordinates.col(v) << x, y, 0.0;
                            v++;
                        }
                    }

                    // cell faces
                    unsigned int nFaces = faceOrders.size();

                    int count = 0;
                    for (unsigned int i = 0; i < nFaces; i++)
                    {
                        const unsigned numFaceVertices = faceOrders[i];
                        vector<unsigned int> faceVertices(numFaceVertices);
                        Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, numFaceVertices);
                        const int faceNeigh = faceNeighs[i];

                        count++;
                        unsigned int numMinus1Vertices = 0;
                        for (unsigned int j = 0; j < numFaceVertices; j++)
                        {
                            if(cell3DVertices[localFaceVertices[count]] == -1)
                                numMinus1Vertices++;
                            else
                            {
                                faceVertices[j] = cell3DVertices[localFaceVertices[count]];
                                faceVertCoordinates.col(j) = cellVertCoordinates.col(localFaceVertices[count]);
                            }
                            count++;
                        }

                        if(numMinus1Vertices != 0)
                            continue;

                        if(faceNeigh < 0) // if it is a boundary face
                        {
                            // compute original cell2D triangulation
                            const vector<unsigned int> convexCell2DUnalignedVerticesFilter = geometryUtilities.UnalignedPoints(faceVertCoordinates);
                            const Eigen::MatrixXd convexCell2DUnalignedVertices = geometryUtilities.ExtractPoints(faceVertCoordinates,
                                                                                                                  convexCell2DUnalignedVerticesFilter);

                            const vector<unsigned int> convexCell2DTriangulationFiltered = geometryUtilities.PolygonTriangulationByFirstVertex(convexCell2DUnalignedVertices);
                            vector<unsigned int> convexCell2DTriangulation(convexCell2DTriangulationFiltered.size());
                            for (unsigned int ocf = 0; ocf < convexCell2DTriangulationFiltered.size(); ocf++)
                                convexCell2DTriangulation[ocf] = convexCell2DUnalignedVerticesFilter[convexCell2DTriangulationFiltered[ocf]];

                            const vector<Eigen::Matrix3d> convexCell2DTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(faceVertCoordinates,
                                                                                                                                         convexCell2DTriangulation);

                            const unsigned int numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();

                            // compute original cell2D area and centroids
                            Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
                            Eigen::MatrixXd convexCell2DTriangulationCentroids(3, numConvexCell2DTriangulation);
                            for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
                            {
                                convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
                                if(convexCell2DTriangulationAreas[cct] < 0)
                                    throw runtime_error("uffa");
                                convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonBarycenter(convexCell2DTriangulationPoints[cct]);
                            }

                            const double convexCell2DArea = convexCell2DTriangulationAreas.sum();
                            VoronoiPoints.col(countPoints++) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                                                 convexCell2DTriangulationAreas,
                                                                                                 convexCell2DArea);
                        }
                        else
                            throw runtime_error("wrong face");
                    }
                    vvol += c.volume();

                }
            }
            while(vl.inc());
        }

        if(abs(cvol - vvol) > 1.0e-12)
            throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");
    }

    double vvol = 0.0;
    // Get cell number
    vector<Cell2D> cell2Ds;
    // chiave: origine-fine, valore: id-marker
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> cell1Ds;
    list<Cell0D> cell0Ds;

    // Initialize Voronoi mesh
    voro::container con
        (
            x_min, x_max, y_min, y_max, -0.5, 0.5,
            n_x, n_y, n_z,
            false, false, false, 8
            );

    for(unsigned int i=0; i < VoronoiPoints.cols(); i++) {
        con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));
    }

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    int ijk, q; double *pp;
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do{
            if (con.compute_cell(c, vl))
            {
                ijk = vl.ijk;
                q = vl.q;
                pp = con.p[ijk] + con.ps*q;

                // Cell barycenter
                const double cx = *pp;
                const double cy = pp[1];
                const double cz = pp[2];

                // Get cell vertices
                std::vector<double> vertices;
                c.vertices(cx, cy, cz, vertices);

                // Get cell face cardinality
                std::vector<int> faceOrders;
                c.face_orders(faceOrders);

                // Get face-vertices connectivity
                std::vector<int> localFaceVertices;
                c.face_vertices(localFaceVertices);

                // Get cell neighbors IDs
                std::vector<int> faceNeighs;
                c.neighbors(faceNeighs);


                // Fill mesh connectivity structures
                // cell vertices
                const unsigned int nCoords = vertices.size();
                const unsigned int numVertices = nCoords/3;
                std::vector<int> cell3DVertices(numVertices);

                Eigen::MatrixXd cellVertCoordinates = Eigen::MatrixXd::Zero(3, numVertices);

                unsigned int v = 0;
                unsigned int i = 0;
                while (i < nCoords)
                {
                    double x = vertices[i++];
                    double y = vertices[i++];
                    double z = vertices[i++];

                    if(!geometryUtilities.IsValueZero(z + 0.5,geometryUtilities.Tolerance1D()))
                    {
                        cell3DVertices[v] = -1;
                        v++;
                    }
                    else{

                        Cell0D cell0D(x, y, 0.0);
                        InsertNewPoints(cell0D, cell0Ds);

                        cell3DVertices[v] = cell0D.id;
                        cellVertCoordinates.col(v) << x, y, 0.0;
                        v++;
                    }
                }

                // cell faces
                unsigned int nFaces = faceOrders.size();

                int count = 0;
                for (unsigned int i = 0; i < nFaces; i++)
                {
                    const unsigned numFaceVertices = faceOrders[i];
                    vector<unsigned int> faceVertices(numFaceVertices);
                    Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, numFaceVertices);
                    const int faceNeigh = faceNeighs[i];

                    count++;
                    unsigned int numMinus1Vertices = 0;
                    for (unsigned int j = 0; j < numFaceVertices; j++)
                    {
                        if(cell3DVertices[localFaceVertices[count]] == -1)
                            numMinus1Vertices++;
                        else
                        {
                            faceVertices[j] = cell3DVertices[localFaceVertices[count]];
                            faceVertCoordinates.col(j) = cellVertCoordinates.col(localFaceVertices[count]);
                        }
                        count++;
                    }

                    if(numMinus1Vertices != 0)
                        continue;

                    if(faceNeigh < 0) // if it is a boundary face
                    {
                        Cell2D cell2D;
                        cell2D.id = cell2Ds.size();
                        cell2D.vertices = faceVertices;
                        cell2D.marker = 0;

                        cell2D.edges.resize(numFaceVertices);

                        for(unsigned int v = 0; v < numFaceVertices; v++)
                        {
                            const unsigned int id_min = min(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);
                            const unsigned int id_max = max(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);

                            auto it = cell1Ds.find({id_min,id_max});

                            if( it != cell1Ds.end())
                                cell2D.edges[v] = it -> second[0];
                            else
                            {
                                unsigned int marker = 0;

                                for(unsigned int e = 0; e < numEdgesDomain; e++)
                                {

                                    if(geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]) &&
                                        geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v+1)%numFaceVertices),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]))
                                    {
                                        marker = numVerticesDomain + e + 1;
                                        break;
                                    }
                                }

                                vector<unsigned int> cell1D(2);
                                cell1D[0] = cell1Ds.size();
                                cell1D[1] = marker;
                                cell1Ds.insert({{id_min, id_max}, cell1D});

                                cell2D.edges[v] = cell1Ds.size() - 1;
                            }
                        }

                        cell2Ds.push_back(cell2D);
                    }
                    else
                        throw runtime_error("wrong face");
                }
                vvol += c.volume();

            }
        }
        while(vl.inc());
    }

    //    if(!geometryUtilities.IsValue3DZero(cvol - vvol))
    //        throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    if(abs(cvol - vvol) > 1.0e-12)
        throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    /// <li> Set Cell0Ds
    const unsigned int numberOfPointsMesh = cell0Ds.size();
    const unsigned int numberOfEdgesMesh = cell1Ds.size();
    const unsigned int numberOfCellsMesh = cell2Ds.size();


    mesh.InitializeDimension(2);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    mesh.Cell2DsInitialize(numberOfCellsMesh);
    mesh.Cell3DsInitialize(0);

    for (const Cell0D& cell0D : cell0Ds)
    {
        const unsigned int id = cell0D.id;

        mesh.Cell0DSetState(id, true);
        mesh.Cell0DInsertCoordinates(id,
                                     Eigen::Vector3d(cell0D.x,
                                                     cell0D.y,
                                                     cell0D.z));


        const bool is_xmin = (geometryUtilities.IsValueZero(cell0D.x - x_min, geometryUtilities.Tolerance1D())); // l == 0
        const bool is_xmax = (geometryUtilities.IsValueZero(cell0D.x - x_max, geometryUtilities.Tolerance1D())); // l == (numLengthPoints - 1)
        const bool is_ymin = (geometryUtilities.IsValueZero(cell0D.y - y_min, geometryUtilities.Tolerance1D())); // h == 0
        const bool is_ymax = (geometryUtilities.IsValueZero(cell0D.y - y_max, geometryUtilities.Tolerance1D())); // h == (numHeightPoints - 1)


        const unsigned int marker = 1 * (is_xmin && is_ymin) +
                                    2 * (is_xmax && is_ymin) +
                                    4 * (is_xmin && is_ymax) +
                                    3 * (is_xmax && is_ymax) +
                                    5 * (is_ymin && !is_xmin && !is_xmax) +
                                    7 * (is_ymax && !is_xmin && !is_xmax) +
                                    8 * (is_xmin && !is_ymin && !is_ymax) +
                                    6 * (is_xmax && !is_ymin && !is_ymax);

        mesh.Cell0DSetMarker(id, marker);
    }

    /// <li> Set Edges
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    for (const auto& cell1D : cell1Ds)
    {
        const unsigned int e = cell1D.second[0];

        mesh.Cell1DSetState(e, true);
        mesh.Cell1DSetMarker(e, cell1D.second[1]);
        mesh.Cell1DInsertExtremes(e,
                                  cell1D.first.first,
                                  cell1D.first.second);
    }

    /// <li> Set cells
    for (const Cell2D& cell2D : cell2Ds)
    {
        const unsigned int f = cell2D.id;

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, cell2D.marker);

        const unsigned int numCellVertices = cell2D.vertices.size();
        const unsigned int numCellEdges = cell2D.edges.size();

        Gedim::Output::Assert(numCellVertices == numCellEdges);

        bool oriented = false;
        for(unsigned int v = 0; v < numCellVertices; v++)
        {
            Eigen::MatrixXd vert = Eigen::MatrixXd::Zero(3, 3);
            vert.col(0) = mesh.Cell0DCoordinates(cell2D.vertices[v]);
            vert.col(1) = mesh.Cell0DCoordinates(cell2D.vertices[(v + 1) % numCellVertices]);
            vert.col(2) = mesh.Cell0DCoordinates(cell2D.vertices[(v + 2) % numCellVertices]);

            oriented = geometryUtilities.IsValuePositive(
                geometryUtilities.PolygonArea(vert), geometryUtilities.Tolerance2D());

            if(oriented)
                break;
        }

        mesh.Cell2DInitializeVertices(f, numCellVertices);
        mesh.Cell2DInitializeEdges(f, numCellEdges);

        if(oriented)
        {
            for (unsigned int v = 0; v < numCellVertices; v++)
                mesh.Cell2DInsertVertex(f, v, cell2D.vertices[v]);

            for (unsigned int e = 0; e < numCellEdges; e++)
                mesh.Cell2DInsertEdge(f, e, cell2D.edges[e]);
        }
        else
        {
            for (unsigned int v = 0; v < numCellVertices; v++)
                mesh.Cell2DInsertVertex(f, v, cell2D.vertices[(numCellVertices - 1) - v]);

            for (unsigned int e = 0; e < numCellEdges - 1; e++)
                mesh.Cell2DInsertEdge(f, e, cell2D.edges[(numCellEdges - 2) - e]);

            mesh.Cell2DInsertEdge(f, numCellEdges - 1, cell2D.edges[numCellEdges - 1]);
        }
    }

}
// ************************************************************************* //
void VoroInterface::GenerateVoronoiTassellations3D(const Eigen::MatrixXd& polyhedronVertices,
                                                   const Eigen::MatrixXi& polyhedronEdges,
                                                   const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                   const unsigned int& numPoints,
                                                   const unsigned int& numIterations,
                                                   Gedim::IMeshDAO& mesh)
{
    const unsigned int numVerticesDomain = polyhedronVertices.cols();
    const unsigned int numEdgesDomain = polyhedronEdges.cols();
    const unsigned int numFacesDomain = polyhedronFaces.size();


    vector<Eigen::Vector3d> edgeDomainTangent(numEdgesDomain);
    vector<double> edgeDomainTangentSquaredNorm(numEdgesDomain);
    vector<Eigen::Vector3d> edgeDomainLineOrigin(numEdgesDomain);
    for(unsigned int e = 0; e < numEdgesDomain; e++)
    {
        edgeDomainLineOrigin[e] = polyhedronVertices.col(polyhedronEdges(0, e));
        const Eigen::Vector3d end = polyhedronVertices.col(polyhedronEdges(1, e));

        Eigen::Vector3d tangent = end - edgeDomainLineOrigin[e];

        edgeDomainTangent[e] = tangent;
        edgeDomainTangentSquaredNorm[e] = tangent.squaredNorm();
    }


    vector<Eigen::Vector3d> faceDomainNormal(numFacesDomain);
    for(unsigned int f = 0; f < numFacesDomain; f++)
    {
        const unsigned int numVerticesFaceDomain = polyhedronFaces[f].cols();
        Eigen::MatrixXd vert = Eigen::MatrixXd::Zero(3, numVerticesFaceDomain);

        for(unsigned int v = 0; v < numVerticesFaceDomain; v++)
            vert.col(v) = polyhedronVertices.col(polyhedronFaces[f](0,v));

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(vert);

        faceDomainNormal[f] = normal;
    }

    vector<unsigned int> map_marker(numFacesDomain, 0);

    const double x_min = polyhedronVertices.row(0).minCoeff(), x_max = polyhedronVertices.row(0).maxCoeff();
    const double y_min = polyhedronVertices.row(1).minCoeff(), y_max = polyhedronVertices.row(1).maxCoeff();
    const double z_min = polyhedronVertices.row(2).minCoeff(), z_max = polyhedronVertices.row(2).maxCoeff();
    const double cvol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

    // Set up the number of blocks that the container is divided into
    const int n_x = 6, n_y = 6, n_z = 6;

    Eigen::MatrixXd VoronoiPoints;
    GenerateRandomPoints(polyhedronVertices,
                         numPoints,
                         VoronoiPoints);

    // Loop on Voronoi cells
    for(unsigned int it = 0; it < numIterations; it++)
    {
        double vvol = 0.0;
        // Initialize Voronoi mesh
        voro::container con
            (
                x_min, x_max, y_min, y_max, z_min, z_max,
                n_x, n_y, n_z,
                false, false, false, 8
                );

        for(unsigned int i=0; i < VoronoiPoints.cols(); i++) {
            con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));
        }

        unsigned int countPoints = 0;

        // Loop on Voronoi cells
        voro::c_loop_all vl(con);
        int ijk, q; double *pp;
        voro::voronoicell_neighbor c;
        if (vl.start())
        {
            do{
                if (con.compute_cell(c, vl))
                {
                    ijk = vl.ijk;
                    q = vl.q;
                    pp = con.p[ijk] + con.ps*q;

                    // Cell barycenter
                    const double cx = *pp;
                    const double cy = pp[1];
                    const double cz = pp[2];

                    // Get cell vertices
                    std::vector<double> vertices;
                    c.vertices(cx, cy, cz, vertices);

                    // Get cell face cardinality
                    std::vector<int> faceOrders;
                    c.face_orders(faceOrders);

                    // Get face-vertices connectivity
                    std::vector<int> localFaceVertices;
                    c.face_vertices(localFaceVertices);

                    // Get cell neighbors IDs
                    std::vector<int> faceNeighs;
                    c.neighbors(faceNeighs);


                    // Fill mesh connectivity structures
                    // cell vertices
                    const unsigned int nCoords = vertices.size();
                    const unsigned int numVertices = nCoords/3;
                    Eigen::MatrixXd Cell3DsVertices = Eigen::MatrixXd::Zero(3, numVertices);

                    unsigned int v = 0;
                    unsigned int i = 0;
                    while (i < nCoords)
                    {
                        double x = vertices[i++];
                        double y = vertices[i++];
                        double z = vertices[i++];

                        Cell3DsVertices.col(v) << x, y, z;
                        v++;
                    }

                    // cell faces
                    unsigned int nFaces = faceOrders.size();
                    vector<Eigen::MatrixXi> Cell3DsFaces(nFaces);

                    int count = 0;
                    for (unsigned int i = 0; i < nFaces; i++)
                    {
                        const unsigned numFaceVertices = faceOrders[i];
                        Cell3DsFaces[i] = Eigen::MatrixXi::Zero(2, numFaceVertices);

                        count++;
                        for (unsigned int j = 0; j < numFaceVertices; j++)
                        {
                            Cell3DsFaces[i](0, j) = localFaceVertices[count];
                            count++;
                        }

                    }

                    const std::vector<Eigen::MatrixXd> Cell3DsFaces3DVertices = geometryUtilities.PolyhedronFaceVertices(Cell3DsVertices,
                                                                                                                         Cell3DsFaces);
                    const std::vector<Eigen::Vector3d> Cell3DsFacesTranslations = geometryUtilities.PolyhedronFaceTranslations(Cell3DsFaces3DVertices);
                    const std::vector<Eigen::Vector3d> Cell3DsFacesNormals = geometryUtilities.PolyhedronFaceNormals(Cell3DsFaces3DVertices);
                    const std::vector<Eigen::Matrix3d> Cell3DsFacesRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(Cell3DsFaces3DVertices,
                                                                                                                                       Cell3DsFacesNormals,
                                                                                                                                       Cell3DsFacesTranslations);

                    const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(Cell3DsFaces,
                                                                                                                                                   Cell3DsFaces3DVertices);

                    const std::vector<Eigen::MatrixXd> Cell3DsFaces2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(Cell3DsFaces3DVertices,
                                                                                                                                Cell3DsFacesTranslations,
                                                                                                                                Cell3DsFacesRotationMatrices);


                    const std::vector<std::vector<Eigen::Matrix3d>> Cell3DsFaces2DTriangulations = geometryUtilities.PolyhedronFaceExtractTriangulationPoints(Cell3DsFaces2DVertices,
                                                                                                                                                               polyhedronFaceTriangulations);



                    const std::vector<bool> Cell3DsFacesNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(Cell3DsFaces3DVertices,
                                                                                                                            geometryUtilities.PolyhedronBarycenter(Cell3DsVertices),
                                                                                                                            Cell3DsFacesNormals);


                    const double Cell3DsVolumes = geometryUtilities.PolyhedronVolumeByBoundaryIntegral(Cell3DsFaces2DTriangulations,
                                                                                                       Cell3DsFacesNormals,
                                                                                                       Cell3DsFacesNormalDirections,
                                                                                                       Cell3DsFacesTranslations,
                                                                                                       Cell3DsFacesRotationMatrices);

                    VoronoiPoints.col(countPoints++) = geometryUtilities.PolyhedronCentroid(Cell3DsFaces2DTriangulations,
                                                                                            Cell3DsFacesNormals,
                                                                                            Cell3DsFacesNormalDirections,
                                                                                            Cell3DsFacesTranslations,
                                                                                            Cell3DsFacesRotationMatrices,
                                                                                            Cell3DsVolumes);

                    vvol += c.volume();

                }
            }
            while(vl.inc());
        }

        if(abs(cvol - vvol) > 1.0e-12)
            throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    }



    double vvol = 0.0;

    // Initialize Voronoi mesh
    voro::container con
        (
            x_min, x_max, y_min, y_max, z_min, z_max,
            n_x, n_y, n_z,
            false, false, false, 8
            );

    for(unsigned int i=0; i < VoronoiPoints.cols(); i++) {
        con.put(i, VoronoiPoints(0, i), VoronoiPoints(1, i), VoronoiPoints(2, i));
    }

    // Get cell number
    unsigned int numberOfCells3D = con.total_particles();
    unsigned int numberOfCellsMesh = 0;

    vector<vector<unsigned int>> cell3DVertices(numberOfCells3D);


    vector<Cell2D> cell2Ds;
    vector<vector<unsigned int>> cellEdges(numberOfCells3D);
    vector<vector<unsigned int>> cellFaces(numberOfCells3D);
    vector<vector<int>> cellFaceNeighbour(numberOfCells3D);
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> cell1Ds;


    list<Cell0D> cell0Ds;

    // Loop on Voronoi cells
    voro::c_loop_all vl(con);
    int ijk, q; double *pp;
    voro::voronoicell_neighbor c;
    if (vl.start())
    {
        do{
            if (con.compute_cell(c, vl))
            {
                ijk = vl.ijk;
                q = vl.q;
                pp = con.p[ijk] + con.ps*q;

                // Cell barycenter
                const double cx = *pp;
                const double cy = pp[1];
                const double cz = pp[2];

                // Get cell3D ID
                const int ID = con.id[ijk][q];
                numberOfCellsMesh++;

                // Get cell vertices
                std::vector<double> vertices;
                c.vertices(cx, cy, cz, vertices);

                // Get cell face cardinality
                std::vector<int> faceOrders;
                c.face_orders(faceOrders);

                // Get face-vertices connectivity
                std::vector<int> localFaceVertices;
                c.face_vertices(localFaceVertices);

                // Get cell neighbors IDs
                std::vector<int> faceNeighs;
                c.neighbors(faceNeighs);


                // conservative resize
                unsigned int ID_max = numberOfCells3D - 1;
                for(unsigned int id_max = 0; id_max < faceNeighs.size(); id_max++)
                    if(faceNeighs[id_max] >= 0)
                        if(static_cast<unsigned int>(faceNeighs[id_max]) > ID_max)
                            ID_max = faceNeighs[id_max];
                ID_max = max(ID_max, static_cast<unsigned int>(ID));
                if(ID_max > numberOfCells3D - 1)
                {
                    numberOfCells3D = ID_max + 1;
                    cell3DVertices.resize(numberOfCells3D);
                    cellEdges.resize(numberOfCells3D);
                    cellFaces.resize(numberOfCells3D);
                    cellFaceNeighbour.resize(numberOfCells3D);
                }

                // Fill mesh connectivity structures
                // cell vertices
                const unsigned int nCoords = vertices.size();
                const unsigned int numVertices = nCoords/3;
                cell3DVertices[ID].resize(numVertices);

                Eigen::MatrixXd cellVertCoordinates = Eigen::MatrixXd::Zero(3, numVertices);

                unsigned int v = 0;
                unsigned int i = 0;
                while (i < nCoords)
                {
                    double x = vertices[i++];
                    double y = vertices[i++];
                    double z = vertices[i++];

                    Cell0D cell0D(x,y,z);
                    InsertNewPoints(cell0D, cell0Ds);

                    cell3DVertices[ID][v] = cell0D.id;
                    cellVertCoordinates.col(v) << x, y, z;
                    v++;
                }

                // cell faces
                unsigned int nFaces = faceOrders.size();
                cellFaceNeighbour[ID].resize(nFaces);
                cellFaces[ID].resize(nFaces);

                int count = 0;
                for (unsigned int i = 0; i < nFaces; i++)
                {
                    const unsigned numFaceVertices = faceOrders[i];
                    vector<unsigned int> faceVertices(numFaceVertices);
                    Eigen::MatrixXd faceVertCoordinates = Eigen::MatrixXd::Zero(3, numFaceVertices);
                    cellFaceNeighbour[ID][i] = faceNeighs[i];

                    count++;
                    for (unsigned int j = 0; j < numFaceVertices; j++)
                    {
                        faceVertices[j] = cell3DVertices[ID][localFaceVertices[count]];
                        faceVertCoordinates.col(j) = cellVertCoordinates.col(localFaceVertices[count]);
                        count++;
                    }

                    if(cellFaceNeighbour[ID][i] < 0) // if it is a boundary face
                    {
                        // prima volta che incontro la faccia  -  nuovi lati possibili - faccia di bordo

                        Cell2D cell2D;
                        cell2D.id = cell2Ds.size();;
                        cell2D.vertices = faceVertices;

                        if( map_marker[abs(faceNeighs[i]) - 1] == 0)
                        {
                            const Eigen::Vector3d barycenter = faceVertCoordinates.rowwise().mean();

                            for(unsigned int f = 0; f < numFacesDomain; f++)
                            {

                                if(geometryUtilities.IsPointOnPlane(barycenter,
                                                                     faceDomainNormal[f],
                                                                     polyhedronVertices.col(polyhedronFaces[f](0,0))))
                                {
                                    map_marker[abs(faceNeighs[i]) - 1] = numVerticesDomain + numEdgesDomain + f + 1;
                                    cell2D.marker = map_marker[abs(faceNeighs[i]) - 1];
                                    break;
                                }

                            }

                            if( map_marker[abs(faceNeighs[i]) - 1] == 0)
                                throw runtime_error("Marker face wrong");

                        }
                        else
                            cell2D.marker = map_marker[abs(faceNeighs[i]) - 1];

                        cell2D.edges.resize(numFaceVertices);
                        for(unsigned int v = 0; v < numFaceVertices; v++)
                        {
                            const unsigned int id_min = min(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);
                            const unsigned int id_max = max(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);

                            auto it = cell1Ds.find({id_min,id_max});

                            if( it != cell1Ds.end())
                            {
                                cell2D.edges[v] = it -> second[0];

                                if((it -> second[1]) == 0)
                                {
                                    for(unsigned int e = 0; e < numEdgesDomain; e++)
                                    {

                                        if(geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                            edgeDomainLineOrigin[e],
                                                                            edgeDomainTangent[e],
                                                                            edgeDomainTangentSquaredNorm[e]) &&
                                            geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v+1)%numFaceVertices),
                                                                            edgeDomainLineOrigin[e],
                                                                            edgeDomainTangent[e],
                                                                            edgeDomainTangentSquaredNorm[e]))
                                        {
                                            (it -> second[1]) = numVerticesDomain + e + 1;
                                            break;
                                        }

                                        if((it -> second[1]) == 0)
                                            (it -> second[1]) = map_marker[abs(faceNeighs[i]) - 1];

                                    }
                                }

                            }
                            else
                            {
                                unsigned int marker = 0;

                                for(unsigned int e = 0; e < numEdgesDomain; e++)
                                {

                                    if(geometryUtilities.IsPointOnLine(faceVertCoordinates.col(v),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]) &&
                                        geometryUtilities.IsPointOnLine(faceVertCoordinates.col((v+1)%numFaceVertices),
                                                                        edgeDomainLineOrigin[e],
                                                                        edgeDomainTangent[e],
                                                                        edgeDomainTangentSquaredNorm[e]))
                                    {
                                        marker = numVerticesDomain + e + 1;
                                        break;
                                    }

                                    if(marker == 0)
                                        marker = map_marker[abs(faceNeighs[i]) - 1];

                                }
                                vector<unsigned int> cell1D(2);
                                cell1D[0] = cell1Ds.size();
                                cell1D[1] = marker;
                                cell1Ds.insert({{id_min, id_max}, cell1D});

                                cell2D.edges[v] = cell1Ds.size() - 1;
                            }
                        }

                        cell2Ds.push_back(cell2D);

                        cellFaces[ID][i] = cell2D.id;
                        for(unsigned int e = 0; e < cell2D.edges.size(); e++)
                            if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2D.edges[e]) == cellEdges[ID].end())
                                cellEdges[ID].push_back(cell2D.edges[e]);

                    }
                    else
                    {
                        if(cellFaceNeighbour[faceNeighs[i]].size()==0)
                        {
                            // prima volta che incontro la faccia -  nuovi lati possibili

                            Cell2D cell2D;
                            cell2D.id = cell2Ds.size();
                            cell2D.marker = 0;
                            cell2D.vertices = faceVertices;

                            cell2D.edges.resize(numFaceVertices);
                            for(unsigned int v = 0; v < numFaceVertices; v++)
                            {
                                const unsigned int id_min = min(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);
                                const unsigned int id_max = max(faceVertices[v],faceVertices[(v+1)%numFaceVertices]);

                                auto it = cell1Ds.find({id_min,id_max});

                                if( it != cell1Ds.end())
                                    cell2D.edges[v] = it -> second[0];
                                else{
                                    vector<unsigned int> cell1D(2);
                                    cell1D[0] = cell1Ds.size();
                                    cell1D[1] = 0;
                                    cell1Ds.insert({{id_min, id_max}, cell1D});
                                    cell2D.edges[v] = cell1Ds.size() - 1;
                                }
                            }

                            cell2Ds.push_back(cell2D);

                            cellFaces[ID][i] = cell2D.id;
                            for(unsigned int e = 0; e < cell2D.edges.size(); e++)
                                if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2D.edges[e]) == cellEdges[ID].end())
                                    cellEdges[ID].push_back(cell2D.edges[e]);

                        }
                        else{
                            for(unsigned int f = 0; f < cellFaceNeighbour[faceNeighs[i]].size(); f++)
                            {

                                // faccia già esistente -  nessun nuovo lato

                                if(cellFaceNeighbour[faceNeighs[i]][f] == ID)
                                {
                                    cellFaces[ID][i] = cellFaces[faceNeighs[i]][f];

                                    for(unsigned int e = 0; e < cell2Ds[cellFaces[ID][i]].edges.size(); e++)
                                        if(find(cellEdges[ID].begin(), cellEdges[ID].end(), cell2Ds[cellFaces[ID][i]].edges[e]) == cellEdges[ID].end())
                                            cellEdges[ID].push_back(cell2Ds[cellFaces[ID][i]].edges[e]);

                                    break;
                                }

                                if(f == cellFaceNeighbour[faceNeighs[i]].size() - 1)
                                    throw runtime_error("ID face not found");
                            }

                        }
                    }
                }
                vvol += c.volume();

            }
        }
        while(vl.inc());
    }

    if(abs(cvol - vvol) > 1.0e-12)
        throw runtime_error("Error generating Voronoi cells: volumes do not mathc each other");

    /// <li> Set Cell0Ds
    const unsigned int numberOfPointsMesh = cell0Ds.size();
    const unsigned int numberOfEdgesMesh = cell1Ds.size();
    const unsigned int numberOfFacesMesh = cell2Ds.size();


    mesh.InitializeDimension(3);
    mesh.Cell0DsInitialize(numberOfPointsMesh);
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    mesh.Cell2DsInitialize(numberOfFacesMesh);
    mesh.Cell3DsInitialize(numberOfCellsMesh);

    for (const Cell0D& cell0D : cell0Ds)
    {
        const unsigned int id = cell0D.id;

        mesh.Cell0DSetState(id, true);
        mesh.Cell0DInsertCoordinates(id,
                                     Eigen::Vector3d(cell0D.x,
                                                     cell0D.y,
                                                     cell0D.z));


        const bool is_xmin = (geometryUtilities.IsValueZero(cell0D.x - x_min, geometryUtilities.Tolerance1D())); // l == 0
        const bool is_xmax = (geometryUtilities.IsValueZero(cell0D.x - x_max, geometryUtilities.Tolerance1D())); // l == (numLengthPoints - 1)
        const bool is_ymin = (geometryUtilities.IsValueZero(cell0D.y - y_min, geometryUtilities.Tolerance1D())); // h == 0
        const bool is_ymax = (geometryUtilities.IsValueZero(cell0D.y - y_max, geometryUtilities.Tolerance1D())); // h == (numHeightPoints - 1)
        const bool is_zmin = (geometryUtilities.IsValueZero(cell0D.z - z_min, geometryUtilities.Tolerance1D())); // w == 0
        const bool is_zmax = (geometryUtilities.IsValueZero(cell0D.z - z_max, geometryUtilities.Tolerance1D())); // w == (numWidthPoints - 1)

        const unsigned int marker = 1 * (is_xmin && is_ymin && is_zmin) +
                                    2 * (is_xmax && is_ymin && is_zmin) +
                                    3 * (is_xmax && is_ymax && is_zmin) +
                                    4 * (is_xmin && is_ymax && is_zmin) +
                                    5 * (is_xmin && is_ymin && is_zmax) +
                                    6 * (is_xmax && is_ymin && is_zmax) +
                                    7 * (is_xmax && is_ymax && is_zmax) +
                                    8 * (is_xmin && is_ymax && is_zmax) +
                                    9 * (!is_xmin && !is_xmax && is_ymin && is_zmin) +
                                    10 * (is_xmax && !is_ymin && !is_ymax && is_zmin) +
                                    11 * (!is_xmin && !is_xmax && is_ymax && is_zmin) +
                                    12 * (is_xmin && !is_ymin && !is_ymax && is_zmin) +
                                    13 * (!is_xmin && !is_xmax && is_ymin && is_zmax) +
                                    14 * (is_xmax && !is_ymin && !is_ymax && is_zmax) +
                                    15 * (!is_xmin && !is_xmax && is_ymax && is_zmax) +
                                    16 * (is_xmin && !is_ymin && !is_ymax && is_zmax) +
                                    17 * (is_xmin && is_ymin && !is_zmin && !is_zmax) +
                                    18 * (is_xmax && is_ymin && !is_zmin && !is_zmax) +
                                    19 * (is_xmax && is_ymax && !is_zmin && !is_zmax) +
                                    20 * (is_xmin && is_ymax && !is_zmin && !is_zmax) +
                                    21 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmin) +
                                    22 * (!is_xmin && !is_xmax && !is_ymin && !is_ymax && is_zmax) +
                                    23 * (is_xmin && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
                                    24 * (is_xmax && !is_ymin && !is_ymax && !is_zmin && !is_zmax) +
                                    25 * (!is_xmin && !is_xmax && is_ymin  && !is_zmin && !is_zmax) +
                                    26 * (!is_xmin && !is_xmax &&  is_ymax && !is_zmin && !is_zmax);

        mesh.Cell0DSetMarker(id, marker);
    }

    /// <li> Set Edges
    mesh.Cell1DsInitialize(numberOfEdgesMesh);
    for (const auto& cell1D : cell1Ds)
    {
        const unsigned int e = cell1D.second[0];

        mesh.Cell1DSetState(e, true);
        mesh.Cell1DSetMarker(e, cell1D.second[1]);
        mesh.Cell1DInsertExtremes(e,
                                  cell1D.first.first,
                                  cell1D.first.second);
    }

    /// <li> Set Faces
    for (const Cell2D& cell2D : cell2Ds)
    {
        const unsigned int f = cell2D.id;

        mesh.Cell2DSetState(f, true);
        mesh.Cell2DSetMarker(f, cell2D.marker);

        const unsigned int numCellVertices = cell2D.vertices.size();
        const unsigned int numCellEdges = cell2D.edges.size();

        Gedim::Output::Assert(numCellVertices == numCellEdges);

        mesh.Cell2DInitializeVertices(f, numCellVertices);
        mesh.Cell2DInitializeEdges(f, numCellEdges);

        for (unsigned int v = 0; v < numCellVertices; v++)
            mesh.Cell2DInsertVertex(f, v, cell2D.vertices[v]);

        for (unsigned int e = 0; e < numCellEdges; e++)
            mesh.Cell2DInsertEdge(f, e, cell2D.edges[e]);
    }

    /// <li> Set Cells
    unsigned int cells3DIndex = 0;
    for (unsigned int c = 0; c < numberOfCells3D; c++)
    {
        const unsigned int numVertices = cell3DVertices[c].size();
        const unsigned int numEdges = cellEdges[c].size();
        const unsigned int numFaces = cellFaces[c].size();

        if(numVertices == 0)
            continue;

        mesh.Cell3DSetState(cells3DIndex, true);
        mesh.Cell3DSetMarker(cells3DIndex, 0);

        mesh.Cell3DInitializeVertices(cells3DIndex, numVertices);
        mesh.Cell3DInitializeEdges(cells3DIndex, numEdges);
        mesh.Cell3DInitializeFaces(cells3DIndex, numFaces);

        for (unsigned int v = 0; v < numVertices; v++)
            mesh.Cell3DInsertVertex(cells3DIndex, v, cell3DVertices[c][v]);

        for (unsigned int e = 0; e < numEdges; e++)
            mesh.Cell3DInsertEdge(cells3DIndex, e, cellEdges[c][e]);

        for (unsigned int f = 0; f < numFaces; f++)
            mesh.Cell3DInsertFace(cells3DIndex, f, cellFaces[c][f]);

        cells3DIndex++;

    }
}
// ************************************************************************* //
void VoroInterface::GenerateRandomPoints(const Eigen::MatrixXd& domainVertices,
                                         const unsigned int& numPoints,
                                         Eigen::MatrixXd& VoronoiPoints)
{
    const double x_min = domainVertices.row(0).minCoeff(), x_max = domainVertices.row(0).maxCoeff();
    const double y_min = domainVertices.row(1).minCoeff(), y_max = domainVertices.row(1).maxCoeff();
    const double z_min = domainVertices.row(2).minCoeff(), z_max = domainVertices.row(2).maxCoeff();


    sleep(1);
    time_t t = time(nullptr);
    srand((unsigned int) t);
    cout << "time in seconds: " << t << endl;
    VoronoiPoints = Eigen::MatrixXd::Zero(3, numPoints);
    for(unsigned int i=0; i < numPoints; i++)
    {
        VoronoiPoints(0, i) = x_min + rnd() * ( x_max - x_min );
        VoronoiPoints(1, i) = y_min + rnd() * ( y_max - y_min );
        VoronoiPoints(2, i) = z_min + rnd() * ( z_max - z_min );
    }
}
// ************************************************************************* //
void VoroInterface::GenerateRandomPoints(const Eigen::MatrixXd& domainVertices,
                                         const unsigned int& numPoints,
                                         voro::container& con)
{
    const double x_min = domainVertices.row(0).minCoeff(), x_max = domainVertices.row(0).maxCoeff();
    const double y_min = domainVertices.row(1).minCoeff(), y_max = domainVertices.row(1).maxCoeff();
    const double z_min = domainVertices.row(2).minCoeff(), z_max = domainVertices.row(2).maxCoeff();

    srand((unsigned int) time(0));
    for(unsigned int i=0; i < numPoints; i++)
    {
        double x = x_min + rnd() * ( x_max - x_min );
        double y = y_min + rnd() * ( y_max - y_min );
        double z = z_min + rnd() * ( z_max - z_min );
        con.put(i, x, y, z);
    }
}
// ************************************************************************* //
void VoroInterface::GenerateCartesianPoints3D(const Eigen::MatrixXd& polyhedronVertices,
                                              const unsigned int& numPoints,
                                              voro::container& con)
{
    Eigen::Vector3d parallelepipedLengthTangent = polyhedronVertices.col(1) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedHeightTangent = polyhedronVertices.col(3) - polyhedronVertices.col(0);
    Eigen::Vector3d parallelepipedWidthTangent = polyhedronVertices.col(4) - polyhedronVertices.col(0);

    vector<double> lengthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                             0.0, 1.0, 1);
    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                             0.0, 1.0, 1);
    vector<double> widthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numPoints + 1,
                                                                                            0.0, 1.0, 1);

    const unsigned int numLengthPoints = lengthMeshCurvilinearCoordinates.size();
    const unsigned int numHeightPoints = heightMeshCurvilinearCoordinates.size();
    const unsigned int numWidthPoints = widthMeshCurvilinearCoordinates.size();

    // create cell0Ds
    unsigned int cell0DIndex = 0;
    for (unsigned int w = 0; w < numWidthPoints; w++)
    {
        for (unsigned int h = 0; h < numHeightPoints; h++)
        {
            for (unsigned int l = 0; l < numLengthPoints; l++)
            {
                const Eigen::Vector3d coordinate = polyhedronVertices.col(0) +
                                                   lengthMeshCurvilinearCoordinates[l] * parallelepipedLengthTangent +
                                                   heightMeshCurvilinearCoordinates[h] * parallelepipedHeightTangent +
                                                   widthMeshCurvilinearCoordinates[w] * parallelepipedWidthTangent;
                con.put(cell0DIndex, coordinate(0), coordinate(1), coordinate(2));
                cell0DIndex++;
            }
        }
    }
}
#endif
// ************************************************************************* //
}

