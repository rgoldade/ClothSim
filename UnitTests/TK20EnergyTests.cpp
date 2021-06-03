#include "gtest/gtest.h"

#include "SimpleGeometryBuilder.h"
#include "TK20Energy.h"
#include "TriMesh.h"
#include "Utilities.h"

using namespace ClothSim;

TEST(TK20_ENERGY_TESTS, ZERO_XY_DEFORM_TEST)
{
	Mat3x2d F = Mat3x2d::Zero();
	F(0, 0) = 1.;
	F(1, 1) = 1.;

	EXPECT_EQ(TK20Energy::computeStretchEnergy(F), 0.);
	EXPECT_EQ(TK20Energy::computeShearEnergy(F), 0.);
	
	EXPECT_EQ(TK20Energy::computeStretchPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeStretchPK1Autodiff(F).norm(), 0.);

	EXPECT_EQ(TK20Energy::computeShearPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeShearPK1Autodiff(F).norm(), 0.);
}

TEST(TK20_ENERGY_TESTS, ZERO_YZ_DEFORM_TEST)
{
	Mat3x2d F = Mat3x2d::Zero();
	F(1, 0) = 1.;
	F(2, 1) = 1.;

	EXPECT_EQ(TK20Energy::computeStretchEnergy(F), 0.);
	EXPECT_EQ(TK20Energy::computeShearEnergy(F), 0.);

	EXPECT_EQ(TK20Energy::computeStretchPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeStretchPK1Autodiff(F).norm(), 0.);

	EXPECT_EQ(TK20Energy::computeShearPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeShearPK1Autodiff(F).norm(), 0.);
}

TEST(TK20_ENERGY_TESTS, ZERO_XZ_DEFORM_TEST)
{
	Mat3x2d F = Mat3x2d::Zero();
	F(0, 0) = 1.;
	F(2, 1) = 1.;

	EXPECT_EQ(TK20Energy::computeStretchEnergy(F), 0.);
	EXPECT_EQ(TK20Energy::computeShearEnergy(F), 0.);

	EXPECT_EQ(TK20Energy::computeStretchPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeStretchPK1Autodiff(F).norm(), 0.);

	EXPECT_EQ(TK20Energy::computeShearPK1Analytical(F).norm(), 0.);
	EXPECT_EQ(TK20Energy::computeShearPK1Autodiff(F).norm(), 0.);
}

TEST(TK20_ENERGY_TESTS, RANDOM_DEFORM_TEST)
{
	for (int i = 0; i < 1000; ++i)
	{
		Mat3x2d F = .5 * Mat3x2d::Random() + Mat3x2d::Ones();

		Mat3x2d pk1StretchAnalytical = TK20Energy::computeStretchPK1Analytical(F);
		Mat3x2d pk1StretchAutodiff = TK20Energy::computeStretchPK1Autodiff(F);

		EXPECT_TRUE(isNearlyEqual((pk1StretchAnalytical - pk1StretchAutodiff).lpNorm<Eigen::Infinity>(), 0., 1e-10, false));

		Mat3x2d pk1ShearAnalytical = TK20Energy::computeShearPK1Analytical(F);
		Mat3x2d pk1ShearAutodiff = TK20Energy::computeShearPK1Autodiff(F);

		EXPECT_TRUE(isNearlyEqual((pk1StretchAnalytical - pk1StretchAutodiff).lpNorm<Eigen::Infinity>(), 0., 1e-10, false));
	}
}