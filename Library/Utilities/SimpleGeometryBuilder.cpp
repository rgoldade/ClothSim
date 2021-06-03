#include "SimpleGeometryBuilder.h"

namespace ClothSim
{

TriMesh makeSquareCloth(const Vec3d& center, const Vec3d& scale, const int divisions)
{
    Vec3d bottomLeftCorner = center;

    bottomLeftCorner[0] -= scale[0];
    bottomLeftCorner[1] -= scale[1];

    Vec3d dx = 2. * scale / double(divisions);

    VecVec3d vertices;
    // Build square sheet of cloth
    for (int x = 0; x < divisions + 1; ++x)
        for (int y = 0; y < divisions + 1; ++y)
        {
            Vec3d point = bottomLeftCorner;
            point[0] += double(x) * dx[0];
            point[1] += double(y) * dx[1];

            vertices.push_back(point);
        }

    // Add mid point
    for (int x = 0; x < divisions; ++x)
        for (int y = 0; y < divisions; ++y)
        {
            Vec3d point = bottomLeftCorner;
            point[0] += (double(x) + .5) * dx[0];
            point[1] += (double(y) + .5) * dx[1];

            vertices.push_back(point);
        }

    VecVec3i triangles;
    // Build square sheet of cloth

    auto flattenCorners = [&divisions](int x, int y) { return x + (divisions + 1) * y; };
    auto flattenMidPoint = [&divisions](int x, int y) { return (divisions + 1) * (divisions + 1) + x + divisions * y; };

    for (int x = 0; x < divisions; ++x)
        for (int y = 0; y < divisions; ++y)
        {
            triangles.emplace_back(flattenMidPoint(x, y), flattenCorners(x, y), flattenCorners(x + 1, y));
            triangles.emplace_back(flattenMidPoint(x, y), flattenCorners(x + 1, y), flattenCorners(x + 1, y + 1));
            triangles.emplace_back(flattenMidPoint(x, y), flattenCorners(x + 1, y + 1), flattenCorners(x, y + 1));
            triangles.emplace_back(flattenMidPoint(x, y), flattenCorners(x, y + 1), flattenCorners(x, y));
        }

    return TriMesh(vertices, triangles);
}

}