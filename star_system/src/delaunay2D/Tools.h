#ifndef _TOOLS_
#define _TOOLS_

#include <vector>
#include <array>
#include <map>
#include <algorithm>

struct float2 { float x; float y; };
struct float3 { float x; float y; float z; };
struct uint2 { unsigned int i; unsigned int j; };

namespace Delaunay2D {

	/*! \brief Definition of Face and edge Ids.
	*/
	using FaceId = unsigned int;
	using EdgeId = unsigned int;

	/*! \brief Floating point geometric epsilon.
	*/
	static const float GEOM_EPSILON = float(1e-5);

	struct Tools {

		/*! \brief Adding two 3D points, in place function, i.e. p <- p + alpha q
		*/
		static void addIn(float3& p, float alpha, const float3& q)
		{
			p.x += alpha * q.x;
			p.y += alpha * q.y;
			p.z += alpha * q.z;
		}

		/*! \brief Creating PQ vector.
		*/
		static float3 makeVec(const float3& p, const float3& q)
		{
			return { q.x - p.x,  q.y - p.y,  q.z - p.z };
		}

		/*! \brief Computing cross product from two input vectors.
		*/
		static float3 cross(const float3& u, const float3& v)
		{
			return { u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x };
		}

		/*! \brief Computing square distance between two points.
		*/
		static float squareDistance(const float3& p, const float3& q)
		{
			const float dx = q.x - p.x;
			const float dy = q.y - p.y;
			const float dz = q.z - p.z;
			return dx * dx + dy * dy + dz * dz;
		}

	};
}

#endif
