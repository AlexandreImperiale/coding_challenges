#ifndef _FACE_3D_
#define _FACE_3D_

#include <vector>
#include "Tools.h"

namespace Delaunay2D {

	/*! \brief Definition of a face in a 3D convex hull.
	*/
	struct Face3D {

		/*! \brief Associated edge ids.
		*/
		EdgeId e0, e1, e2;

		/*! \brief Associated plane normal.
		*/
		float3 normal;

		/*! \brief Associated center.
		*/
		float3 center;
	};

	struct Face3DTools {

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
			const auto v = Tools::makeVec(face.center, pnt);
			return v.x * face.normal.x + v.y * face.normal.y + v.z * face.normal.z > -GEOM_EPSILON;
		}
	};
}

#endif
