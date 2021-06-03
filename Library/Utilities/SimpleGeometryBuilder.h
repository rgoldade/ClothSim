#ifndef CLOTH_SIM_SIMPLE_GEOMETRY_BUILDER_H
#define CLOTH_SIM_SIMPLE_GEOMETRY_BUILDER_H

#include "Utilities.h"
#include "TriMesh.h"

namespace ClothSim
{

TriMesh makeSquareCloth(const Vec3d& center, const Vec3d& scale, const int divisions);

}

#endif