#ifndef _CONVEX_HULL_3D_
#define _CONVEX_HULL_3D_

#include <algorithm>
#include <array>
#include <stack>
#include <set>
#include <vector>

#include "ConvexHullBuilder3D.h"
#include "OutsideSet.h"
#include "Tools.h"

namespace Delaunay2D {

	struct ConvexHull3D {

		/*!
			\brief Computing center of a convex hull.
			\params pnts are the points defining a set of points.
			\params hull is the convex hull associated to a set of points.
			\returns returns the center of the convex hull.
		*/
		static float3 getCenter(const std::vector<float3>& pnts, const std::vector<Face3D>& hull)
		{
			const float coef = float(1.0) / (3 * hull.size());

			float3 center;
			center.x = center.y = center.z = float(0.0);

			for(const auto& face : hull)
			{
				Tools::addIn(center, coef, pnts[face.i]);
				Tools::addIn(center, coef, pnts[face.j]);
				Tools::addIn(center, coef, pnts[face.k]);
			}

			return center;
		}

		/*!
			\brief Computing convex hull of a 3D point set.
			\param pnts are the input points.
			\returns a set of 3D faces corresponding to the convex hull of the point set.
		*/
		static std::vector<Face3D> makeConvexHull(const std::vector<float3>& pnts)
		{
			ConvexHullBuilder3D builder;

			// Creating first simplex.
			const auto initFaces = builder.init(pnts);

			// Creating outside sets associated to created faces.
			std::array<OutsideSet, 4> initOutsideSets;
			for(size_t i = 0, ni = pnts.size(); i < ni; ++i)
				for(size_t j = 0; j < 4; ++j)
				{
					const auto& face = builder.getFace(initFaces[j]);
					if(Face3D::isAbove(pnts[i], face))
						OutsideSet::addPointIndex(initOutsideSets[j], face, pnts, i);
				}

			// Initializing stack of faces to be processed
			std::stack<FaceId> facesInProcess;
			for(size_t i = 0; i < 4; ++i)
			{
				if(!OutsideSet::isEmpty(initOutsideSets[i])) facesInProcess.push(initFaces[i]);
				builder.setFaceOutsideSet(initFaces[i], std::move(initOutsideSets[i]));
			}

			// Core loop of quickhull algo.
			while(!facesInProcess.empty())
			{
				// Poping top face id.
				const auto currentFace = facesInProcess.top();
				facesInProcess.pop();

				// Accessing outside set.
				const auto& outsideSet = builder.getFaceOutsideSet(currentFace);

				// Extracting furthest points.
				const auto& furthestPnt = pnts[outsideSet.furthestPntIdx];

				// Building set of visible faces from the furthest point.
				std::vector<FaceId> visibleFaces;
				visibleFaces.push_back(currentFace);

				std::set<FaceId> visitedCandidates;
				visitedCandidates.insert(currentFace);

				std::stack<FaceId> candidates;
				for(const auto& neighbour : builder.getNeighbours(currentFace))
					candidates.push(neighbour);

				while(!candidates.empty())
				{
					const auto candidateFace = candidates.top();
					candidates.pop();

					visitedCandidates.insert(candidateFace);
					if(Face3D::isAbove(furthestPnt, builder.getFace(candidateFace)))
						visibleFaces.push_back(candidateFace);

					for(const auto& neighbour : builder.getNeighbours(candidateFace))
						if(visitedCandidates.find(neighbour) == visitedCandidates.end())
							candidates.push(neighbour);
				}

				// Creating new faces in the hull builder from the boundary of the visible faces.
				const auto newFaces = builder.addFacesFromBoundary(pnts, visibleFaces, outsideSet.furthestPntIdx);

				// Creating outside sets associated to new faces.
				const size_t nNewFaces = newFaces.size();
				std::vector<OutsideSet> newOutsideSets(nNewFaces);
				for(const auto& visibleFace : visibleFaces)
					for(const auto& outsidePntIdx : builder.getFaceOutsideSet(visibleFace).pntIdxs)
						for(size_t i = 0; i < nNewFaces; ++i)
						{
							const auto& face = builder.getFace(newFaces[i]);
							if(Face3D::isAbove(pnts[outsidePntIdx], face))
								OutsideSet::addPointIndex(newOutsideSets[i], face, pnts, outsidePntIdx);
						}

				// Setting outside set of new faces and updating faces to be processed.
				for(size_t i = 0, ni = newFaces.size(); i < ni; ++i)
				{
					if(!OutsideSet::isEmpty(newOutsideSets[i])) facesInProcess.push(newFaces[i]);
					builder.setFaceOutsideSet(newFaces[i], std::move(newOutsideSets[i]));
				}

				// Deleting visible faces.
				for(const auto& face : visibleFaces) builder.deleteFace(face);
			};

			return builder.get();
		}
	};
}

#endif
