#include "gtest/gtest.h"

#include "SimpleGeometryBuilder.h"
#include "TriMesh.h"
#include "Utilities.h"

using namespace ClothSim;

static void testTriMesh(const TriMesh& mesh)
{
	// Build local adjacent triangle list
	std::vector<std::vector<int>> adjacentTriangles(mesh.vertexCount());

	for (int triIndex = 0; triIndex < mesh.triangleCount(); ++triIndex)
	{
		const Vec3i& tri = mesh.triangle(triIndex);
		for (int vertIndex : {0, 1, 2})
			adjacentTriangles[tri[vertIndex]].push_back(triIndex);
	}

	// Verify match against pre-build adjacency map
	ASSERT_EQ(adjacentTriangles.size(), mesh.adjacentTriangles().size());

	for (int vertIndex = 0; vertIndex < mesh.vertexCount(); ++vertIndex)
	{
		ASSERT_EQ(adjacentTriangles[vertIndex].size(), mesh.adjacentTriangles()[vertIndex].size());
		for (int triIndex = 0; triIndex < adjacentTriangles[vertIndex].size(); ++triIndex)
			EXPECT_EQ(adjacentTriangles[vertIndex][triIndex], mesh.adjacentTriangles()[vertIndex][triIndex]);
	}

	// Verify triangle's adjacent vertex reciprocates
	for (int vertIndex = 0; vertIndex < adjacentTriangles.size(); ++vertIndex)
	{
		for (int triIndex : adjacentTriangles[vertIndex])
		{
			const Vec3i& tri = mesh.triangle(triIndex);
			EXPECT_TRUE(tri[0] == vertIndex || tri[1] == vertIndex || tri[2] == vertIndex);
		}
	}

	// Verify face's adjacent vertex reciprocates
	for (int triIndex = 0; triIndex < mesh.triangleCount(); ++triIndex)
	{
		const Vec3i& tri = mesh.triangle(triIndex);
		for (int localVertIndex : {0, 1, 2})
		{
			int vertIndex = tri[localVertIndex];
			EXPECT_TRUE(std::find(adjacentTriangles[vertIndex].begin(), adjacentTriangles[vertIndex].end(), triIndex) != adjacentTriangles[vertIndex].end());
		}
	}

	// Verify no dangling vertices
	for (int vertIndex = 0; vertIndex < adjacentTriangles.size(); ++vertIndex)
	{
		EXPECT_TRUE(adjacentTriangles[vertIndex].size() >= 2);
	}
}

TEST(TRI_MESH_TESTS, SQUARE_SHEET_TEST)
{
	Vec3d center = Vec3d::Random();
	Vec3d scale = 1.5 * Vec3d::Ones();

	TriMesh mesh = makeSquareCloth(center, scale, 20);
}