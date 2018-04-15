#ifndef _DELAUNAY_2D_
#define _DELAUNAY_2D_

#include <vector>

#include "ConvexHull3D.h"
#include "Face3D.h"
#include "Tools.h"

namespace Delaunay2D {

	struct DelaunayBuilder2D {

		/*!
			\brief Creating parabolic lifted 3D point set from input 2D point set.
			\params pnts2D are input 2D points.
			\returns lifted 3D points.
		*/
		static std::vector<float3> lift(const std::vector<float2>& pnts2D)
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
			\returns the set of edges corresponding to the result of the projection onto the xy-plane.
				WARNING: The set of edges is not unique!
		*/
		static std::vector<uint2> projectUnderside(const std::vector<float3>& pnts, const std::vector<Face3D>& hull)
		{
			std::vector<uint2> edges;

			// Extracting center of hull.
			const auto center = ConvexHull3D::getCenter(pnts, hull);

			// Applying projection.
			edges.reserve(3 * hull.size());
			for(const auto& face : hull)
				if(Face3D::isDownward(center, face))
					for(const auto& e : Face3D::getEdges(face))
						edges.push_back(e);
			edges.shrink_to_fit();

			return edges;
		}

		/*!
			\brief Computing Delaunay triangulation associated to an input set of 2D points.
			\params pnts is the set of 2D points.
			\returns the set of edges forming the Delaunay triangulation.
				WARNING: The set of edges is not unique!
		*/
		static std::vector<uint2> makeTriangulation(const std::vector<float2>& pnts)
		{
			const auto pnts3D = lift(pnts);
			const auto hull = ConvexHull3D::makeConvexHull(pnts3D);
			return projectUnderside(pnts3D, hull);
		}

		/*!
			\brief Removing doublons in set of edges. This operation is performed in O(nlog(n)), where
			n is the number of edges. After this operation, the edges are all oriented toward the point
			with bigger index, i.e. edge.i < edge.j.
			\params edges is the set of edges with potential doublons.
			\params maxPntIdx is the maximal index of points linked to an edge in the set of edges.
		*/
		static void removeDoublons(std::vector<uint2>& edges, uint maxPntIdx)
		{
			// Imposing same directions to every edges.
			for(auto& e : edges)
				if(e.i > e.j)
					std::swap(e.i, e.j);

			// Sorting edges.
			std::sort(edges.begin(), edges.end(), [maxPntIdx](const uint2& e0, const uint2& e1)
			{
					return (e0.j * maxPntIdx + e0.i) < (e1.j * maxPntIdx + e1.i);
			});

			// Applying unique.
			edges.erase(std::unique(edges.begin(), edges.end(), [](const uint2& e0, const uint2& e1)
			{
				return e0.i == e1.i && e0.j == e1.j;
			}), edges.end());
		}

	};
}

#endif
