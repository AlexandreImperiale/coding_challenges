#ifndef _CONVEX_HULL_3D_
#define _CONVEX_HULL_3D_

#include <vector>
#include <array>
#include <algorithm>

struct float3 { float x; float y; float z; };

using uint = size_t;
struct uint2 { uint i; uint j; };

namespace ConvexHull3DTools {

	/*! \brief Floating point geometric epsilon.
	*/
	static const float GEOM_EPSILON = float(1e-5);

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
	static auto makeVec(const float3& p, const float3& q) -> float3
	{
		return { q.x - p.x,  q.y - p.y,  q.z - p.z };
	}

	/*! \brief Computing cross product from two input vectors.
	*/
	static auto cross(const float3& u, const float3& v) -> float3
	{
		return { u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x };
	}
}

namespace ConvexHull3D {

	/*! \brief Definition of a face in a 3D convex hull.
	*/
	struct Face3D {
		/*! \brief Associated point indexes.
		*/
		uint p0, p1, p2;

		/*! \brief Associated index of neighbouring faces.
		*/
		std::vector<uint> neighbours;

		/*! \brief Associated plane normal.
		*/
		float3 normal;

		/*! \brief Associated center.
		*/
		float3 center;
	};

	/*! 
		\brief Creating a face from three points.
		\params pnts are the input 3D points.
		\params p0, p1, p2 are the input coordinate points.
	*/
	static auto makeFace(const std::vector<float3>& pnts, uint p0, uint p1, uint p2) -> Face3D
	{
		Face3D face;
		face.p0 = p0; face.p1 = p1; face.p2 = p2;

		// Computing center.
		const float coef = float(1.0 / 3.0);
		face.center.x = face.center.y = face.center.z = 0.;
		ConvexHull3DTools::addIn(face.center, coef, pnts[p0]);
		ConvexHull3DTools::addIn(face.center, coef, pnts[p1]);
		ConvexHull3DTools::addIn(face.center, coef, pnts[p2]);

		// Computing normal.
		const auto du = ConvexHull3DTools::makeVec(pnts[p0], pnts[p1]);
		const auto dv = ConvexHull3DTools::makeVec(pnts[p0], pnts[p2]);
		face.normal = ConvexHull3DTools::cross(du, dv);

		return face;
	}

	/*!
		\brief Extracting every edges associated to a face.
		\params face is the input face.
		\returns the array of three edges forming the face.
	*/
	static auto getEdges(const Face3D& face) -> std::array<uint2, 3>
	{
		return {{ { face.p0, face.p1}, { face.p1, face.p2}, { face.p2, face.p0} }};
	}

	/*!
		\brief Computing center of a convex hull.
		\params pnts are the points defining a set of points.
		\params hull is the convex hull associated to a set of points.
		\returns returns the center of the convex hull.
	*/
	static auto getCenter(const std::vector<float3>& pnts, const std::vector<Face3D>& hull)
	{
		const float coef = float(1.0) / (3 * hull.size());

		float3 center;
		center.x = center.y = center.z = float(0.0);

		for(const auto& face : hull)
		{
			ConvexHull3DTools::addIn(center, coef, pnts[face.p0]);
			ConvexHull3DTools::addIn(center, coef, pnts[face.p1]);
			ConvexHull3DTools::addIn(center, coef, pnts[face.p2]);
		}

		return center;
	}

	/*!
		\brief Testing if a face in a convex hull is a downward facing face.
		\params center is the center of the convex hull.
		\params pnts are the set of points used to build the convex hull.
		\params face is the input face in the convex hull.
		\returns true if the face is facing downward.
	*/
	static auto isDownward(const float3& center, const Face3D& face) -> bool
	{
		const float dot =
			(center.x - face.center.x) * face.normal.x +
			(center.y - face.center.y) * face.normal.y +
			(center.z - face.center.z) * face.normal.z;

		if( dot > float(0.0) ) return face.normal.z > -ConvexHull3DTools::GEOM_EPSILON;
		else return face.normal.z < ConvexHull3DTools::GEOM_EPSILON;
	}

	/*!
		\brief Computing convex hull of a 3D point set.
		\param pnts are the input points.
		\returns a set of 3D faces corresponding to the convex hull of the point set.
	*/
	static auto makeConvexHull(const std::vector<float3>& /*pnts*/) -> std::vector<Face3D>
	{
		return {};
	}
}

#endif
