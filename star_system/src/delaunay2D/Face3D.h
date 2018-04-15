#ifndef _FACE_3D_
#define _FACE_3D_

#include <vector>
#include "Tools.h"

namespace Delaunay2D {

	/*! \brief Definition of a face in a 3D convex hull.
	*/
	struct Face3D {

		/*! \brief Associated point indexes.
		*/
		uint i, j, k;

		/*! \brief Associated plane normal.
		*/
		float3 normal;

		/*! \brief Associated center.
		*/
		float3 center;

		/*!
			\brief Creating a face from three points.
			\params pnts are the input 3D points.
			\params p0, p1, p2 are the input coordinate points.
		*/
		static Face3D makeFace(const std::vector<float3>& pnts, uint p0, uint p1, uint p2)
		{
			Face3D face;
			face.i = p0; face.j = p1; face.k = p2;

			// Computing center.
			const float coef = float(1.0 / 3.0);
			face.center.x = face.center.y = face.center.z = 0.;
			Tools::addIn(face.center, coef, pnts[p0]);
			Tools::addIn(face.center, coef, pnts[p1]);
			Tools::addIn(face.center, coef, pnts[p2]);

			// Computing normal.
			const auto du = Tools::makeVec(pnts[p0], pnts[p1]);
			const auto dv = Tools::makeVec(pnts[p0], pnts[p2]);
			face.normal = Tools::cross(du, dv);

			return face;
		}

		/*!
			\brief Extracting every edges associated to a face.
			\params face is the input face.
			\returns the array of three edges forming the face.
		*/
		static std::array<uint2, 3> getEdges(const Face3D& face)
		{
			return {{ { face.i, face.j}, { face.j, face.k}, { face.k, face.i} }};
		}

		/*!
			\brief Testing if a face in a convex hull is a downward facing face.
			\params center is the center of the convex hull.
			\params pnts are the set of points used to build the convex hull.
			\params face is the input face in the convex hull.
			\returns true if the face is facing downward.
		*/
		static bool isDownward(const float3& center, const Face3D& face)
		{
			const float dot =
				(center.x - face.center.x) * face.normal.x +
				(center.y - face.center.y) * face.normal.y +
				(center.z - face.center.z) * face.normal.z;

			if( dot > float(0.0) ) return face.normal.z > -GEOM_EPSILON;
			else return face.normal.z < GEOM_EPSILON;
		}

		/*!
			\brief Testing if a point is above a face
			\params pnt is the input point to be tested.
			\params face is the input face
			\returns true if the point is above the face.
		*/
		static bool isAbove(const float3& pnt, const Face3D& face)
		{
			const auto v = Tools::makeVec(pnt, face.center);
			return v.x * face.normal.x + v.y * face.normal.y + v.z * face.normal.z > -GEOM_EPSILON;
		}
	};
}

#endif
