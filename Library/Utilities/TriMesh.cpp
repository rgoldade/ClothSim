#include "TriMesh.h"

namespace ClothSim
{

TriMesh::TriMesh(const VecVec3d& vertices, const VecVec3i& triangles)
: myVertices(vertices)
, myTriangles(triangles)
{
    // Build vertex to triangle list
    myAdjacentTriangles.resize(vertices.size());

    for (int triIndex = 0; triIndex < myTriangles.size(); ++triIndex)
    {
        const Vec3i& tri = myTriangles[triIndex];

        for (int localIndex : {0, 1, 2})
            myAdjacentTriangles[tri[localIndex]].push_back(triIndex);
    }
}

const VecVec3i& TriMesh::triangles() const
{
    return myTriangles;
}

const VecVec3d& TriMesh::vertices() const
{
    return myVertices;
}

int TriMesh::triangleCount() const
{
    return int(myTriangles.size());
}

int TriMesh::vertexCount() const
{
    return int(myVertices.size());
}

const std::vector<std::vector<int>>& TriMesh::adjacentTriangles() const
{
    return myAdjacentTriangles;
}

void TriMesh::drawMesh(Renderer& renderer,
                    const bool renderFaces,
                    const Vec3d& faceColor,
                    const bool renderEdges,
                    const Vec3d& edgeColor,
                    const double edgeWidth,
                    const bool renderVertices,
                    const Vec3d& vertex_color) const
{
    if (renderFaces)
    {
        // Build normals per triangle
        VecVec3d triNormals(myTriangles.size());
        for (int triIndex = 0; triIndex < myTriangles.size(); ++triIndex)
            triNormals[triIndex] = normal(triIndex);

        VecVec3d normals(myVertices.size(), Vec3d::Zero());

        for (int vertIndex = 0; vertIndex < myVertices.size(); ++vertIndex)
        {
            for (int triIndex : myAdjacentTriangles[vertIndex])
                normals[vertIndex] += triNormals[triIndex];

            assert(myAdjacentTriangles[vertIndex].size() > 0);
            normals[vertIndex].array() /= double(myAdjacentTriangles[vertIndex].size());
        }

        renderer.addTriangles(myVertices, normals, myTriangles, faceColor);
    }

    if (renderEdges)
    {
        VecVec3d startPoints;
        VecVec3d endPoints;

        for (const Vec3i& triangle : myTriangles)
        {
            for (int edge_index : {0,1,2})
            {
                startPoints.push_back(myVertices[triangle[edge_index]]);
                endPoints.push_back(myVertices[triangle[(edge_index + 1) % 3]]);
            }
        }

        renderer.addLines(startPoints, endPoints, edgeColor, edgeWidth);
    }

    if (renderVertices)
    {
        renderer.addPoints(myVertices, vertex_color, 5);
    }
}

}