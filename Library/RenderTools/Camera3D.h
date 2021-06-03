#ifndef CLOTH_SIM_CAMERA_3D_H
#define CLOTH_SIM_CAMERA_3D_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Utilities.h"

namespace ClothSim
{

class Camera3D
{
public:
    Camera3D(const Vec3d& target = Vec3d::Zero(), double targetDistance = 1, double heading = 0, double pitch = 0,
             double fieldOfView = 45, double nearClippingDistance = 0.01, double farClippingDistance = 100.);

    void mouse(int button, int state, int x, int y);
    void drag(int x, int y);

    void reset();
    void transform(const Vec2i& windowSize);

private:
    enum class MouseAction
    {
        INACTIVE,
        ROTATE,
        TRUCK,
        DOLLY
    };

    Vec3d myTarget, myDefaultTarget;
    double myDistance, myDefaultDistance;
    double myHeading, myDefaultHeading;
    double myPitch, myDefaultPitch;
    double myFieldOfView;

    double myNearClippingDistance, myFaceClippingDistance;

    MouseAction myMouseAction;
    Vec2i myOldMouseCoord;
};

}
#endif