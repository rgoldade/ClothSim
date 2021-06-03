#ifndef CLOTH_SIM_FEM_CLOTH_MESH_H
#define CLOTH_SIM_FEM_CLOTH_MESH_H

#include "TK20Energy.h"
#include "TriMesh.h"
#include "Utilities.h"

namespace ClothSim
{

class FEMClothMesh : public TriMesh
{

public:
	FEMClothMesh(const std::vector<int>& fixedVertices,
					const double density,
					const double stretchStiffness,
					const double shearStiffness,
					const VecVec2d& restUVs,
					const TriMesh& inputMesh);

	FORCE_INLINE void setVertex(int vertIndex, const Vec3d& vertex) override
	{
		myStateDirty = true;
		TriMesh::setVertex(vertIndex, vertex);
	}

	FORCE_INLINE bool isVertexFixed(int vertIndex) const
	{
		return myFixedVertices[vertIndex];
	}

	FORCE_INLINE void setVertexVelocity(int vertIndex, const Vec3d& velocity)
	{
		myVertexVelocities[vertIndex] = velocity;
	}

	FORCE_INLINE const Vec3d& vertexVelocity(int vertIndex) const
	{
		return myVertexVelocities[vertIndex];
	}

	FORCE_INLINE double vertexMass(int vertIndex) const
	{
		return myVertexMasses[vertIndex];
	}

	FORCE_INLINE double stretchStiffness() const
	{
		return myStretchStiffness;
	}

	FORCE_INLINE double shearStiffness() const
	{
		return myShearStiffness;
	}

	FORCE_INLINE double triangleRestArea(int triIndex) const
	{
		return myRestTriangleAreas[triIndex];
	}

	FORCE_INLINE const Mat6x9d& dFdX(int triIndex) const
	{
		return mydFdX[triIndex];
	}

	FORCE_INLINE Mat3x2d computeStretchPk1(int triIndex, bool useAutodiff) const
	{
		assert(!myStateDirty);
		return TK20Energy::computeStretchPK1(myF[triIndex], useAutodiff);
	}

	FORCE_INLINE Mat3x2d computeShearPk1(int triIndex, bool useAutodiff) const
	{
		assert(!myStateDirty);
		return TK20Energy::computeShearPK1(myF[triIndex], useAutodiff);
	}

	void buildF();
	void computeState();

private:


	VecVec3d myVertexVelocities;
	std::vector<bool> myFixedVertices;

	const double myDensity;
	const double myStretchStiffness;
	const double myShearStiffness;

	std::vector<double> myRestTriangleAreas;
	std::vector<double> myVertexMasses;

	VecMat2x2d myDmInv;
	VecMat6x9d mydFdX;

	VecMat3x2d myF;

	bool myStateDirty;
};

}

#endif
