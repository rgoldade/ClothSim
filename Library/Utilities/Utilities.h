#ifndef CLOTH_SIM_UTILITIES_H
#define CLOTH_SIM_UTILITIES_H

#include <vector>

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

namespace ClothSim
{

constexpr double PI = 3.1415926535897932384626433;
enum class Axis { XAXIS = 0, YAXIS = 1, ZAXIS = 2 };

// Borrowed from https://github.com/sideeffects/WindingNumber/blob/master/SYS_Types.h
#if defined(__GNUC__) || defined(__clang__)
#define FORCE_INLINE	__attribute__ ((always_inline)) inline
#elif defined(_MSC_VER)
#define FORCE_INLINE	__forceinline
#else
#define FORCE_INLINE	inline
#endif

#ifdef _MSC_VER
#undef min
#undef max
#endif

using Vec2d = Eigen::Matrix<double, 2, 1>;
using Vec2i = Eigen::Matrix<int, 2, 1>;

template<typename T>
using Vec2t = Eigen::Matrix<T, 2, 1>;

using VecVec2d = std::vector<Vec2d, Eigen::aligned_allocator<Vec2d>>;
using VecVec2i = std::vector<Vec2i, Eigen::aligned_allocator<Vec2i>>;

using Vec3f = Eigen::Matrix<float, 3, 1>;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec3i = Eigen::Matrix<int, 3, 1>;

template<typename T>
using Vec3t = Eigen::Matrix<T, 3, 1>;

using VecVec3d = std::vector<Vec3d, Eigen::aligned_allocator<Vec3d>>;
using VecVec3i = std::vector<Vec3i, Eigen::aligned_allocator<Vec3i>>;

template<typename T>
using Vec4t = Eigen::Matrix<T, 4, 1>;

using Vec6d = Eigen::Matrix<double, 6, 1>;

template<typename T>
using Vec6t = Eigen::Matrix<T, 6, 1>;

using Vec9d = Eigen::Matrix<double, 9, 1>;

using VecVec9d = std::vector<Vec9d, Eigen::aligned_allocator<Vec9d>>;

using Mat2x2d = Eigen::Matrix<double, 2, 2>;

template<typename T>
using Mat2x2t = Eigen::Matrix<T, 2, 2>;

using VecMat2x2d = std::vector<Mat2x2d, Eigen::aligned_allocator<Mat2x2d>>;

using Mat3x2d = Eigen::Matrix<double, 3, 2>;

template<typename T>
using Mat3x2t = Eigen::Matrix<T, 3, 2>;

using VecMat3x2d = std::vector<Mat3x2d, Eigen::aligned_allocator<Mat3x2d>>;

using Mat6x9d = Eigen::Matrix<double, 6, 9>;
using VecMat6x9d = std::vector<Mat6x9d, Eigen::aligned_allocator<Mat6x9d>>;

template<typename T, int N>
using VecXt = Eigen::Matrix<T, N, 1>;

template<typename VecType>
using VecVecT = std::vector<VecType, Eigen::aligned_allocator<VecType>>;

using VectorXd = Eigen::VectorXd;

template<typename T>
using VectorXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename RealType>
bool isNearlyEqual(const RealType a, const RealType b, const RealType tolerance = 1e-5, const bool useRelative = true)
{
	if (a == b)
		return true;

	RealType absDiff = std::fabs(a - b);

	RealType avgMag(1);

	if (useRelative)
	{
		avgMag = (std::fabs(a) + std::fabs(b)) / RealType(2);
	}

	return absDiff < tolerance* avgMag;
}

}

#endif