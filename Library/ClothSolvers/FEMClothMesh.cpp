#include "FEMClothMesh.h"

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

namespace ClothSim
{

static std::vector<double> buildTriangleAreas(const TriMesh& mesh, const VecVec2d& restUVs)
{
	std::vector<double> triangleAreas(mesh.triangleCount(), 0);
	tbb::parallel_for(tbb::blocked_range<int>(0, mesh.triangleCount()), [&](tbb::blocked_range<int>& range)
		{
			for (int triIndex = range.begin(); triIndex != range.end(); ++triIndex)
			{
				const Vec3i& tri = mesh.triangle(triIndex);
				const Vec2d& v0 = restUVs[tri[0]];
				const Vec2d& v1 = restUVs[tri[1]];
				const Vec2d& v2 = restUVs[tri[2]];

				triangleAreas[triIndex] = .5 * std::fabs(v0[0] * (v1[1] - v2[1]) + v1[0] * (v2[1] - v0[1]) + v2[0] * (v0[1] - v1[1]));
			}
		});

	return triangleAreas;
}

static std::vector<double> buildVertexMasses(const TriMesh& mesh, const std::vector<double>& triangleAreas, double density)
{
	std::vector<double> vertexMasses(mesh.vertexCount(), 0);
	assert(mesh.vertexCount() == mesh.adjacentTriangles().size());
	assert(triangleAreas.size() == mesh.triangleCount());
	tbb::parallel_for(tbb::blocked_range<int>(0, mesh.vertexCount()), [&](tbb::blocked_range<int>& range)
	{
		for (int vertexIndex = range.begin(); vertexIndex != range.end(); ++vertexIndex)
		{
			for (int triIndex : mesh.adjacentTriangles()[vertexIndex])
				vertexMasses[vertexIndex] += 1. / 3. * density * triangleAreas[triIndex];
		}
	});

	return vertexMasses;
}

static VecMat2x2d buildDmInv(const TriMesh& mesh, const VecVec2d& restUVs)
{
	assert(mesh.vertexCount() == restUVs.size());
	VecMat2x2d dmInv(mesh.triangleCount());
	tbb::parallel_for(tbb::blocked_range<int>(0, mesh.triangleCount()), [&](tbb::blocked_range<int>& range)
	{
		Mat2x2d localD;
		for (int triIndex = range.begin(); triIndex != range.end(); ++triIndex)
		{
			const Vec3i& tri = mesh.triangle(triIndex);

			localD.col(0) = restUVs[tri[1]] - restUVs[tri[0]];
			localD.col(1) = restUVs[tri[2]] - restUVs[tri[0]];

			dmInv[triIndex] = localD.inverse();
		}
	});

	return dmInv;
}

template<typename Scalar>
Vec6t<Scalar> computeF(const Vec6t<Scalar>& D, const Vec4t<Scalar>& Dinv)
{
	Mat3x2t<Scalar> Dmat;
	Dmat.block(0, 0, 3, 1) = D.block(0, 0, 3, 1);
	Dmat.block(0, 1, 3, 1) = D.block(3, 0, 3, 1);

	Mat2x2t<Scalar> DinvMat;
	DinvMat.block(0, 0, 2, 1) = Dinv.block(0, 0, 2, 1);
	DinvMat.block(0, 1, 2, 1) = Dinv.block(2, 0, 2, 1);

	Mat3x2t<Scalar> Fmat = Dmat * DinvMat;

	Vec6t<Scalar> Fvec;
	Fvec.block(0, 0, 3, 1) = Fmat.block(0, 0, 3, 1);
	Fvec.block(3, 0, 3, 1) = Fmat.block(0, 1, 3, 1);

	return Fvec;
}

static Mat3x2d computeF(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Mat2x2d& DmInv)
{
	Mat3x2d Ds;
	Ds.col(0) = v1 - v0;
	Ds.col(1) = v2 - v0;

	return Ds * DmInv;
}

template<typename Scalar>
static Vec6t<Scalar> computeFTemplated(const Vec3t<Scalar>& v0, const Vec3t<Scalar>& v1, const Vec3t<Scalar>& v2, const Mat2x2t<Scalar>& DmInv)
{
	Mat3x2t<Scalar> Ds;
	Ds.col(0) = v1 - v0;
	Ds.col(1) = v2 - v0;

	Mat3x2t<Scalar> Fmat = Ds * DmInv;

	Vec6t<Scalar> Fvec;
	Fvec.block(0, 0, 3, 1) = Fmat.col(0);
	Fvec.block(3, 0, 3, 1) = Fmat.col(1);

	return Fvec;
}

static Mat6x9d builddFdxAutodiff(const TriMesh& mesh, const VecMat2x2d& DmInv, int triIndex)
{
	using autodiff::dual;
	using autodiff::forward::jacobian;
	using autodiff::forward::at;
	using autodiff::forward::wrtpack;

	const Vec3i& tri = mesh.triangle(triIndex);
	const Vec3d& v0 = mesh.vertex(tri[0]);
	const Vec3d& v1 = mesh.vertex(tri[1]);
	const Vec3d& v2 = mesh.vertex(tri[2]);

	Vec3t<dual> dualV0 = v0.cast<dual>();
	Vec3t<dual> dualV1 = v1.cast<dual>();
	Vec3t<dual> dualV2 = v2.cast<dual>();

	Mat2x2t<dual> dualDmInv = DmInv[triIndex].cast<dual>();

	Vec6t<dual> dualF;
	Mat6x9d dFdx = jacobian(computeFTemplated<dual>, wrtpack(dualV0, dualV1, dualV2), at(dualV0, dualV1, dualV2, dualDmInv), dualF);

#if !defined(NDEBUG)
	Mat3x2d F = computeF(v0, v1, v2, DmInv[triIndex]);

	Mat3x2d Fautodiff;
	Fautodiff.col(0) = dualF.block(0, 0, 3, 1).cast<double>();
	Fautodiff.col(1) = dualF.block(3, 0, 3, 1).cast<double>();

	assert((Fautodiff - F).lpNorm<Eigen::Infinity>() < 1e-10);
#endif

	return dFdx;
}

static VecMat6x9d builddFdx(const TriMesh& mesh, const VecMat2x2d& DmInv)
{
	VecMat6x9d dFdx(mesh.triangleCount(), Mat6x9d::Zero());
	tbb::parallel_for(tbb::blocked_range<int>(0, mesh.triangleCount()), [&](tbb::blocked_range<int>& range)
	{
		for (int triIndex = range.begin(); triIndex != range.end(); ++triIndex)
		{
			double s0 = DmInv[triIndex].col(0).sum();
			double s1 = DmInv[triIndex].col(1).sum();

			double d0 = DmInv[triIndex](0, 0);
			double d1 = DmInv[triIndex](1, 0);
			double d2 = DmInv[triIndex](0, 1);
			double d3 = DmInv[triIndex](1, 1);

			// dF / dx0
			dFdx[triIndex](0, 0) = -s0;
			dFdx[triIndex](3, 0) = -s1;

			// dF / dy0
			dFdx[triIndex](1, 1) = -s0;
			dFdx[triIndex](4, 1) = -s1;

			// dF / dz0
			dFdx[triIndex](2, 2) = -s0;
			dFdx[triIndex](5, 2) = -s1;

			// dF / dx1
			dFdx[triIndex](0, 3) = d0;
			dFdx[triIndex](3, 3) = d2;

			// dF / dy1
			dFdx[triIndex](1, 4) = d0;
			dFdx[triIndex](4, 4) = d2;

			// dF / dz1
			dFdx[triIndex](2, 5) = d0;
			dFdx[triIndex](5, 5) = d2;

			// dF / dx2
			dFdx[triIndex](0, 6) = d1;
			dFdx[triIndex](3, 6) = d3;

			// dF / dy2
			dFdx[triIndex](1, 7) = d1;
			dFdx[triIndex](4, 7) = d3;

			// dF / dz2
			dFdx[triIndex](2, 8) = d1;
			dFdx[triIndex](5, 8) = d3;

#if !defined(NDEBUG)
			// Autodiff verification
			Mat6x9d dFdxAutodiff = builddFdxAutodiff(mesh, DmInv, triIndex);
			assert((dFdx[triIndex] - dFdxAutodiff).lpNorm<Eigen::Infinity>() < 1e-10);
#endif
		}
	});

	return dFdx;
}

FEMClothMesh::FEMClothMesh(const std::vector<int>& fixedVertices,
	const double density,
	const double stretchStiffness,
	const double shearStiffness,
	const VecVec2d& restUVs,
	const TriMesh& inputMesh)
	: TriMesh(inputMesh)
	, myVertexVelocities(this->vertexCount(), Vec3d::Zero())
	, myFixedVertices(this->vertexCount(), false)
	, myDensity(density)
	, myStretchStiffness(stretchStiffness)
	, myShearStiffness(shearStiffness)
	, myRestTriangleAreas(buildTriangleAreas(inputMesh, restUVs))
	, myVertexMasses(buildVertexMasses(inputMesh, myRestTriangleAreas, myDensity))
	, myDmInv(buildDmInv(inputMesh, restUVs))
	, mydFdX(builddFdx(inputMesh, myDmInv))
	, myF(this->triangleCount())
	, myStateDirty(true)
{
	for (int vertIndex : fixedVertices)
		myFixedVertices[vertIndex] = true;

	computeState();
}

void FEMClothMesh::buildF()
{
	tbb::parallel_for(tbb::blocked_range<int>(0, this->triangleCount()), [&](const tbb::blocked_range<int>& range)
	{
		for (int triIndex = range.begin(); triIndex != range.end(); ++triIndex)
		{
			const Vec3i& tri = this->triangle(triIndex);
			const Vec3d& v0 = this->vertex(tri[0]);
			const Vec3d& v1 = this->vertex(tri[1]);
			const Vec3d& v2 = this->vertex(tri[2]);

			myF[triIndex] = computeF(v0, v1, v2, myDmInv[triIndex]);
		}
	});
}

void FEMClothMesh::computeState()
{
	if (myStateDirty)
	{
		buildF();
		myStateDirty = false;
	}
}

}