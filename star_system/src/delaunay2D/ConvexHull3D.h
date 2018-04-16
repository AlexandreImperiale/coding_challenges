#ifndef _CONVEX_HULL_3D_
#define _CONVEX_HULL_3D_

#include <unordered_map>
#include <unordered_set>


#include "Face3D.h"
#include "Tools.h"

namespace Delaunay2D {

	struct ConvexHull3D {

		/*! \brief Associated faces.
		*/
		std::unordered_map<FaceId, Face3D> faces;

		/*! \brief Associated edges.
		*/
		std::unordered_map<EdgeId, uint2> edges;

		/*! \brief Mapping from face id to face neighbours.
		*/
		std::unordered_map<FaceId, std::unordered_set<FaceId>> faceNeighbours;

		/*! \brief Mapping from edge id to incident faces.
		*/
		std::unordered_map<EdgeId, std::unordered_set<FaceId>> edgeIncidentFaces;
	};

	struct ConvexHull3DTools {

		/*!
			\brief Computing center of a convex hull.
			\params pnts are the points defining a set of points.
			\params hull is the convex hull associated to a set of points.
			\returns returns the center of the convex hull.
		*/
		static float3 getCenter(const std::vector<float3>& pnts, const ConvexHull3D& hull)
		{
			const float coef = float(1.0) / (3 * hull.faces.size());

			float3 center;
			center.x = center.y = center.z = float(0.0);

			for(const auto& face : hull.faces)
			{
				Tools::addIn(center, coef, pnts[hull.edges.at(face.second.e0).i]);
				Tools::addIn(center, coef, pnts[hull.edges.at(face.second.e1).i]);
				Tools::addIn(center, coef, pnts[hull.edges.at(face.second.e2).i]);
			}

			return center;
		}
	};
}

#endif
