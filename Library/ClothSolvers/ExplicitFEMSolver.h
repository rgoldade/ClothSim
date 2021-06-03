#ifndef CLOTH_SIM_EXPLICIT_FEM_SOLVER_H
#define CLOTH_SIM_EXPLICIT_FEM_SOLVER_H

#include "FEMClothMesh.h"
#include "Utilities.h"

namespace ClothSim
{

class ExplicitFEMSolver
{
public:

	ExplicitFEMSolver(FEMClothMesh& mesh);

	void solveTimestep(double dt);

private:

	void computeAcceleration();

	FEMClothMesh& myMesh;

	VecVec3d myCachedAcceleration;
	VecVec3d myAcceleration;
	VecVec9d myCachedTriangleForces;
};

}

#endif
