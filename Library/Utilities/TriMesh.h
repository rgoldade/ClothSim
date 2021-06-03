#ifndef CLOTH_SIM_TRI_MESH_H
#define CLOTH_SIM_TRI_MESH_H

#include "Utilities.h"
#include "Renderer.h"

namespace ClothSim
{

class TriMesh
{
public:

    TriMesh() = default;

    TriMesh(const VecVec3d& vertices, const VecVec3i& triangles);

    const VecVec3i& triangles() const;
    const VecVec3d& vertices() const;

    FORCE_INLINE const Vec3i& triangle(int triIndex) const
    {
        return myTriangles[triIndex];
    }

    FORCE_INLINE const Vec3d& vertex(int vertIndex) const
    {
        return myVertices[vertIndex];
    }

    FORCE_INLINE virtual void setVertex(int vertIndex, const Vec3d& vertex)
    {
        myVertices[vertIndex] = vertex;
    }

    FORCE_INLINE Vec3d normal(int triIndex) const
    {
        const Vec3i& triangle = myTriangles[triIndex];
        Vec3d vec0 = myVertices[triangle[1]] - myVertices[triangle[0]];
        Vec3d vec1 = myVertices[triangle[2]] - myVertices[triangle[1]];

        Vec3d normal = vec0.cross(vec1);

        double norm = normal.norm();

        if (norm == 0)
            return Vec3d::Zero();

        return normal / norm;
    }

    int triangleCount() const;
    int vertexCount() const;

    const std::vector<std::vector<int>>& adjacentTriangles() const;

    void drawMesh(Renderer& renderer,
                    const bool renderFaces,
                    const Vec3d& faceColor,
                    const bool renderEdges,
                    const Vec3d& edgeColor,
                    const double edgeWidth,
                    const bool renderVertices,
                    const Vec3d& vertexColor) const;

private:

    VecVec3d myVertices;
    VecVec3i myTriangles;

    // List of adjacent triangles to mesh vertices
    std::vector<std::vector<int>> myAdjacentTriangles;
};

}

#endif