#include <memory>

#include "Camera3D.h"
#include "ExplicitFEMSolver.h"
#include "FEMClothMesh.h"
#include "Renderer.h"
#include "SimpleGeometryBuilder.h"
#include "TriMesh.h"
#include "Utilities.h"

using namespace ClothSim;

static std::unique_ptr<Renderer> gRenderer;
static std::unique_ptr<Camera3D> gCamera;

static std::unique_ptr<FEMClothMesh> gMesh;
static std::unique_ptr<ExplicitFEMSolver> gSolver;

static constexpr double gDt = .001;
static constexpr double gDrawDt = 1. / 30.;

static constexpr double gDensity = .01;
static constexpr double gStretchStiffness = 1.;
static constexpr double gShearStiffness = 1.;

static bool isDisplayDirty = true;
static bool runSimulation = false;
static bool runSingleFrame = false;

static double gAccumDt = 0;

void keyboard(unsigned char key, int x, int y)
{
    if (key == ' ')
        runSimulation = !runSimulation;
    else if (key == 'n')
        runSingleFrame = true;

    isDisplayDirty = true;
}

void display()
{
    if (runSimulation || runSingleFrame)
    {
        gSolver->solveTimestep(gDt);
        gAccumDt += gDt;

        if (gAccumDt >= gDrawDt)
        {
            gAccumDt = 0;
            isDisplayDirty = true;
        }
        runSingleFrame = false;
    }

    if (isDisplayDirty)
    {
        gRenderer->clear();

        gMesh->drawMesh(*gRenderer,
                        true /* render faces */, .5 * Vec3d::Ones(),
                        true /* render edges */, Vec3d::Zero(), 5.,
                        true /* render vertices */, Vec3d(1, 0, 0));
        glutPostRedisplay();

        isDisplayDirty = false;
    }
}

int main(int argc, char** argv)
{
    Vec3d center = Vec3d::Zero();
    Vec3d scale = Vec3d::Ones();
    Vec3d topRightCorner = center + scale;
    Vec3d bottomLeftCorner = center - scale;

    gRenderer = std::make_unique<Renderer>("FEM cloth sim", Vec2i(1000, 1000), Vec2d(bottomLeftCorner[0], bottomLeftCorner[1]),
                                   topRightCorner[1] - bottomLeftCorner[1], &argc, argv);

    gCamera = std::make_unique<Camera3D>(.5 * (topRightCorner + bottomLeftCorner), 2., 0., 0.);
    gRenderer->setCamera(gCamera.get());

    TriMesh mesh = makeSquareCloth(center, scale, 20);

    // Build rest UVs
    VecVec2d restUVs(mesh.vertexCount());

    for (int vertIndex = 0; vertIndex < mesh.vertexCount(); ++vertIndex)
    {
        restUVs[vertIndex] = Vec2d(mesh.vertex(vertIndex)[0], mesh.vertex(vertIndex)[1]);
        assert(mesh.vertex(vertIndex)[2] == 0);
    }
    
    // Build fixed vertices
    std::vector<int> fixed_vertices;
    for (int vertIndex = 0; vertIndex < mesh.vertexCount(); ++vertIndex)
    {
        if (mesh.vertex(vertIndex)[1] > .99 && (mesh.vertex(vertIndex)[0] > .99 || mesh.vertex(vertIndex)[0] < -.99))
            fixed_vertices.push_back(vertIndex);
    }

    gMesh = std::make_unique<FEMClothMesh>(fixed_vertices, gDensity, gStretchStiffness, gShearStiffness, restUVs, mesh);

    // Move mesh vertices to create "swing"
    Vec3d offset = Vec3d::Zero();
    offset[2] = .5;
    for (int vertIndex = 0; vertIndex < mesh.vertexCount(); ++vertIndex)
    {
        if (!gMesh->isVertexFixed(vertIndex))
        {
            double dx = (1 - gMesh->vertex(vertIndex)[1]) / 2.;
            gMesh->setVertex(vertIndex, gMesh->vertex(vertIndex) + dx * offset);
        }
    }

    gMesh->computeState();

    gSolver = std::make_unique<ExplicitFEMSolver>(*gMesh);

    std::function<void()> displayFunc = display;
    gRenderer->setUserDisplay(displayFunc);

    std::function<void(unsigned char, int, int)> keyboardFunc = keyboard;
    gRenderer->setUserKeyboard(keyboardFunc);

    gRenderer->run();
}