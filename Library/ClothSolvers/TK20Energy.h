#ifndef CLOTH_SIM_TK20_ENERGY_H
#define CLOTH_SIM_TK20_ENERGY_H

#include "Utilities.h"

// Stretching and shearing energy based on the FEM
// formulation from "A finite element formulation of Baraff-Witkin cloth"
// by Ted Kim, 2020

namespace ClothSim
{

namespace TK20Energy
{
	double computeStretchEnergy(const Mat3x2d& F);

	Mat3x2d computeStretchPK1Autodiff(const Mat3x2d& F);
	Mat3x2d computeStretchPK1Analytical(const Mat3x2d& F);
	Mat3x2d computeStretchPK1(const Mat3x2d& F, bool useAutodiff);

	double computeShearEnergy(const Mat3x2d& F);

	Mat3x2d computeShearPK1Autodiff(const Mat3x2d& F);
	Mat3x2d computeShearPK1Analytical(const Mat3x2d& F);
	Mat3x2d computeShearPK1(const Mat3x2d& F, bool useAutodiff);
}

}
#endif