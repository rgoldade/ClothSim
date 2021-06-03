#include "TK20Energy.h"

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

namespace ClothSim
{

template<typename Scalar>
static Scalar stretchEnergy(const Vec6t<Scalar>& F)
{
	Vec2t<Scalar> au(1., 0.);
	Vec2t<Scalar> av(0., 1.);

	Mat3x2t<Scalar> Fmat;
	Fmat.col(0) = F.block(0, 0, 3, 1);
	Fmat.col(1) = F.block(3, 0, 3, 1);

	Scalar i5u = au.transpose() * Fmat.transpose() * Fmat * au;
	Scalar i5v = av.transpose() * Fmat.transpose() * Fmat * av;

	assert(i5u >= 0);
	assert(i5v >= 0);

	return pow(sqrt(i5u) - Scalar(1), 2) + pow(sqrt(i5v) - Scalar(1), 2);
}

double TK20Energy::computeStretchEnergy(const Mat3x2d& F)
{
	Vec6d Fvec;
	Fvec.block(0, 0, 3, 1) = F.col(0);
	Fvec.block(3, 0, 3, 1) = F.col(1);
	return stretchEnergy<double>(Fvec);
}

Mat3x2d TK20Energy::computeStretchPK1Autodiff(const Mat3x2d& F)
{
	using autodiff::dual;
	using autodiff::forward::gradient;
	using autodiff::forward::at;
	using autodiff::forward::wrtpack;

	Vec6t<dual> Fvec;
	Fvec.block(0, 0, 3, 1) = F.col(0).cast<dual>();
	Fvec.block(3, 0, 3, 1) = F.col(1).cast<dual>();

	Vec6t<dual> pk1Vec = gradient(stretchEnergy<dual>, wrtpack(Fvec), at(Fvec));

	Mat3x2d pk1;
	pk1.col(0) = pk1Vec.block(0, 0, 3, 1).cast<double>();
	pk1.col(1) = pk1Vec.block(3, 0, 3, 1).cast<double>();

	return pk1;
}

Mat3x2d TK20Energy::computeStretchPK1Analytical(const Mat3x2d& F)
{
	Mat3x2d pk1 = 2. * F;

	double i5u = F.col(0).array().pow(2).sum();
	double i5v = F.col(1).array().pow(2).sum();

	double uCoeff = 1. - 1. / std::sqrt(i5u);
	double vCoeff = 1. - 1. / std::sqrt(i5v);

	pk1.col(0).array() *= uCoeff;
	pk1.col(1).array() *= vCoeff;

#if !defined(NDEBUG)
	Mat3x2d pk1Autodiff = computeStretchPK1Autodiff(F);
	double diff = (pk1 - pk1Autodiff).lpNorm<Eigen::Infinity>();
	assert(diff < 1e-10);
#endif

	return pk1;
}

Mat3x2d TK20Energy::computeStretchPK1(const Mat3x2d& F, bool useAutodiff)
{
	if (useAutodiff)
		return computeStretchPK1Autodiff(F);
	else
		return computeStretchPK1Analytical(F);
}

template<typename Scalar>
static Scalar shearEnergy(const Vec6t<Scalar>& F)
{
	Vec2t<Scalar> au(1., 0.);
	Vec2t<Scalar> av(0., 1.);

	Mat3x2t<Scalar> Fmat;
	Fmat.col(0) = F.block(0, 0, 3, 1);
	Fmat.col(1) = F.block(3, 0, 3, 1);

	Scalar i6 = au.transpose() * Fmat.transpose() * Fmat * av;

	return pow(i6, 2);
}

double TK20Energy::computeShearEnergy(const Mat3x2d& F)
{
	Vec6d Fvec;
	Fvec.block(0, 0, 3, 1) = F.col(0);
	Fvec.block(3, 0, 3, 1) = F.col(1);
	return shearEnergy<double>(Fvec);
}

Mat3x2d TK20Energy::computeShearPK1Autodiff(const Mat3x2d& F)
{
	using autodiff::dual;
	using autodiff::forward::gradient;
	using autodiff::forward::at;
	using autodiff::forward::wrtpack;

	Vec6t<dual> Fvec;
	Fvec.block(0, 0, 3, 1) = F.col(0).cast<dual>();
	Fvec.block(3, 0, 3, 1) = F.col(1).cast<dual>();

	Vec6t<dual> pk1Vec = gradient(shearEnergy<dual>, wrtpack(Fvec), at(Fvec));

	Mat3x2d pk1;
	pk1.col(0) = pk1Vec.block(0, 0, 3, 1).cast<double>();
	pk1.col(1) = pk1Vec.block(3, 0, 3, 1).cast<double>();

	return pk1;
}

Mat3x2d TK20Energy::computeShearPK1Analytical(const Mat3x2d& F)
{
	Mat3x2d pk1;
	pk1.col(0) = 2. * F.col(1);
	pk1.col(1) = 2. * F.col(0);

	pk1.array() *= F.col(0).dot(F.col(1));

#if !defined(NDEBUG)
	Mat3x2d pk1Autodiff = computeShearPK1Autodiff(F);
	double diff = (pk1 - pk1Autodiff).lpNorm<Eigen::Infinity>();
	assert(diff < 1e-10);
#endif

	return pk1;
}

Mat3x2d TK20Energy::computeShearPK1(const Mat3x2d& F, bool useAutodiff)
{
	if (useAutodiff)
		return computeShearPK1Autodiff(F);
	else
		return computeShearPK1Analytical(F);
}

}