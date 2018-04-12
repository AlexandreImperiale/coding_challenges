#ifndef _DELAUNAY_2D_
#define _DELAUNAY_2D_

#include <vector>
#include "ConvexHull3D.h"

struct float2 { float x; float y; };

namespace Delaunay2D {

	/*!
		\brief Creating parabolic lifted 3D point set from input 2D point set.
		\params pnts2D are input 2D points.
		\returns lifted 3D points.
	*/
	static auto lift(const std::vector<float2>& pnts2D) -> std::vector<float3>
	{
		std::vector<float3> pnts3D;

		pnts3D.reserve(pnts2D.size());
		for(const auto& p : pnts2D)
			pnts3D.push_back({ p.x, p.y, p.x * p.x + p.y * p.y });

		return pnts3D;
	}

	/*!
		\brief Projection of underside of a convex hull onto the xy-plane.
		\params pnts are the set of input points included within the convex hull.
		\params hull is the convex hull associated to the input points.
		\returns the set of edges corresponding to the result of the projection onto the xy-plane. WARNING: The set of edges is not unique!
	*/
	static auto projectUnderside(const std::vector<float3>& pnts, const std::vector<ConvexHull3D::Face3D>& hull) -> std::vector<uint2>
	{
		std::vector<uint2> edges;

		// Extracting center of hull.
		const auto center = ConvexHull3D::getCenter(pnts, hull);

		// Applying projection.
		edges.reserve(3 * hull.size());
		for(const auto& face : hull)
			if(ConvexHull3D::isDownward(center, face))
				for(const auto& e : ConvexHull3D::getEdges(face))
					edges.push_back(e);
		edges.shrink_to_fit();

		return edges;
	}

	/*!
		\brief Computing Delaunay triangulation associated to an input set of 2D points.
		\params pnts is the set of 2D points.
		\returns the set of edges forming the Delaunay triangulation. WARNING: The set of edges is not unique!
	*/
	static auto makeTriangulation(const std::vector<float2>& pnts) -> std::vector<uint2>
	{
		const auto pnts3D = lift(pnts);
		const auto hull = ConvexHull3D::makeConvexHull(pnts3D);
		return projectUnderside(pnts3D, hull);
	}
}

#endif
