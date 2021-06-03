#include "ExplicitFEMSolver.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

namespace ClothSim
{

ExplicitFEMSolver::ExplicitFEMSolver(FEMClothMesh& mesh)
	: myMesh(mesh)
	, myCachedAcceleration(mesh.vertexCount(), Vec3d::Zero())
	, myAcceleration(mesh.vertexCount(), Vec3d::Zero())
	, myCachedTriangleForces(mesh.triangleCount())
{}

void ExplicitFEMSolver::computeAcceleration()
{
	// Apply gravity
	Vec3d gravity(0., -9.8, 0.);
	tbb::parallel_for(tbb::blocked_range<int>(0, myMesh.vertexCount()), [&](const tbb::blocked_range<int>& range)
	{
		for (int vertIndex = range.begin(); vertIndex != range.end(); ++vertIndex)
		{
			if (myMesh.isVertexFixed(vertIndex))
				myAcceleration[vertIndex] = Vec3d::Zero();
			else
			{
				myAcceleration[vertIndex] = gravity;
			}
		}
	});

	// Apply elastic forces
	tbb::parallel_for(tbb::blocked_range<int>(0, myMesh.triangleCount()), [&](const tbb::blocked_range<int>& range)
	{
		Mat3x2d pk1;
		Vec6d vectorizedPk1;
		Vec9d stretchForce;
		Vec9d shearForce;
		for (int triIndex = range.begin(); triIndex != range.end(); ++triIndex)
		{
			{
				pk1 = myMesh.computeStretchPk1(triIndex, false);

				vectorizedPk1.block(0, 0, 3, 1) = pk1.col(0);
				vectorizedPk1.block(3, 0, 3, 1) = pk1.col(1);

				stretchForce = -myMesh.stretchStiffness() * myMesh.triangleRestArea(triIndex) * myMesh.dFdX(triIndex).transpose() * vectorizedPk1;
			}
			
			{
				pk1 = myMesh.computeShearPk1(triIndex, false);

				vectorizedPk1.block(0, 0, 3, 1) = pk1.col(0);
				vectorizedPk1.block(3, 0, 3, 1) = pk1.col(1);

				shearForce = -myMesh.shearStiffness() * myMesh.triangleRestArea(triIndex) * myMesh.dFdX(triIndex).transpose() * vectorizedPk1;
			}
		
			myCachedTriangleForces[triIndex] = stretchForce + shearForce;

#if !defined(NDEBUG)
			assert(std::fabs(stretchForce[0] + stretchForce[3] + stretchForce[6]) < 1e-6);
			assert(std::fabs(stretchForce[1] + stretchForce[4] + stretchForce[7]) < 1e-6);
			assert(std::fabs(stretchForce[2] + stretchForce[5] + stretchForce[8]) < 1e-6);
#endif

#if !defined(NDEBUG)
			assert(std::fabs(shearForce[0] + shearForce[3] + shearForce[6]) < 1e-6);
			assert(std::fabs(shearForce[1] + shearForce[4] + shearForce[7]) < 1e-6);
			assert(std::fabs(shearForce[2] + shearForce[5] + shearForce[8]) < 1e-6);
#endif

		}
	});

	tbb::parallel_for(tbb::blocked_range<int>(0, myMesh.vertexCount()), [&](const tbb::blocked_range<int>& range)
	{
		for (int vertIndex = range.begin(); vertIndex != range.end(); ++vertIndex)
		{
			if (myMesh.isVertexFixed(vertIndex))
				continue;
			else
			{
				for (int triIndex : myMesh.adjacentTriangles()[vertIndex])
				{
					const Vec3i& tri = myMesh.triangle(triIndex);

					for (int localVertIndex : {0, 1, 2})
					{
						if (tri[localVertIndex] == vertIndex)
						{
							myAcceleration[vertIndex][0] += myCachedTriangleForces[triIndex][localVertIndex * 3] / myMesh.vertexMass(vertIndex);
							myAcceleration[vertIndex][1] += myCachedTriangleForces[triIndex][localVertIndex * 3 + 1] / myMesh.vertexMass(vertIndex);
							myAcceleration[vertIndex][2] += myCachedTriangleForces[triIndex][localVertIndex * 3 + 2] / myMesh.vertexMass(vertIndex);
						}
					}

				}
			}
		}
	});
}

void ExplicitFEMSolver::solveTimestep(double dt)
{
	// Get half-step
	tbb::parallel_for(tbb::blocked_range<int>(0, myMesh.vertexCount()), [&](const tbb::blocked_range<int>& range)
	{
		for (int vertIndex = range.begin(); vertIndex != range.end(); ++vertIndex)
		{
			if (myMesh.isVertexFixed(vertIndex))
				continue;
			else
			{
				Vec3d velocity = myMesh.vertexVelocity(vertIndex) + dt / 2. * myCachedAcceleration[vertIndex];
				myMesh.setVertex(vertIndex, myMesh.vertex(vertIndex) + dt * velocity);
			}
		}
	});

	myMesh.computeState();

	// Update acceleration
	computeAcceleration();

	// Update velocity
	tbb::parallel_for(tbb::blocked_range<int>(0, myMesh.vertexCount()), [&](const tbb::blocked_range<int>& range)
	{
		for (int vertIndex = range.begin(); vertIndex != range.end(); ++vertIndex)
		{
			if (myMesh.isVertexFixed(vertIndex))
				continue;
			else
			{
				myMesh.setVertexVelocity(vertIndex, myMesh.vertexVelocity(vertIndex) + dt / 2. * (myAcceleration[vertIndex] + myCachedAcceleration[vertIndex]));
			}
		}
	});

	std::swap(myAcceleration, myCachedAcceleration);
}

}