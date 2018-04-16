#ifndef _CONVEX_HULL_BUILDER_3D_
#define _CONVEX_HULL_BUILDER_3D_

#include <algorithm>
#include <array>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "ConvexHull3D.h"
#include "OutsideSet.h"

namespace Delaunay2D {

	class ConvexHullBuilder3D {

	public:

		/*! \brief Builder.
		*/
		ConvexHullBuilder3D() : maxEdgeId_(0), maxFaceId_(0) {}

		/*!
			\brief Computing convex hull of a 3D point set.
			\param pnts are the input points.
			\returns a set of 3D faces corresponding to the convex hull of the point set.
		*/
		ConvexHull3D make(const std::vector<float3>& pnts)
		{
			// Declaration of convex hull to be built.
			ConvexHull3D hull;

			// Creating first simplex.
			init(pnts, hull);

			// Creating outside sets associated to created faces.
			std::unordered_map<FaceId, OutsideSet> faceOutsideSets;
			for(const auto& face : hull.faces)
				faceOutsideSets[face.first] = OutsideSet{};

			// Creating outside sets associated to created faces.
			for(unsigned int i = 0, ni = pnts.size(); i < ni; ++i)
				for(const auto& face : hull.faces)
					if(Face3DTools::isAbove(pnts[i], face.second))
						OutsideSetTools::addPointIndex(faceOutsideSets[face.first], face.second, pnts, i);

			// Initializing stack of faces to be processed
			std::stack<FaceId> facesInProcess;
			for(const auto& face : hull.faces)
				if(!OutsideSetTools::isEmpty(faceOutsideSets[face.first]))
				 	facesInProcess.push(face.first);

			// Core loop of quickhull algo.
			while(!facesInProcess.empty())
			{
				// Poping top face id.
				const auto currentFaceId = facesInProcess.top();
				facesInProcess.pop();

				// Accessing outside set.
				const auto& outsideSet = faceOutsideSets[currentFaceId];

				// Building set of visible faces from the furthest point.
				const auto visibleFacesId = makeConeOfVisibleFaces(hull, currentFaceId, pnts[outsideSet.furthestPntIdx]);

				// Creating new faces in the hull builder from the boundary of the visible faces.
				const auto newFacesId = addFacesFromVisibleFacesBoundary(hull, visibleFacesId, outsideSet.furthestPntIdx);

				// Creating outside sets associated to new faces.
				for(const auto& newFaceId : newFacesId)
					faceOutsideSets[newFaceId] = OutsideSet{};

				// Creating outside sets associated to new faces.
				for(const auto& visibleFaceId : visibleFacesId)
					for(const auto& outsidePntIdx : faceOutsideSets[visibleFaceId].pntIdxs)
						for(const auto& newFaceId : newFacesId)
						{
							const auto& newFace = hull.faces[newFaceId];
							if(Face3DTools::isAbove(pnts[outsidePntIdx], newFace))
									OutsideSetTools::addPointIndex(faceOutsideSets[newFaceId], newFace, pnts, outsidePntIdx);
						}

				// Setting outside set of new faces and updating faces to be processed.
				for(const auto& newFaceId : newFacesId)
					if(!OutsideSetTools::isEmpty(faceOutsideSets[newFaceId]))
						facesInProcess.push(newFaceId);

				// Deleting visible faces & associated outside sets.
				for(const auto& visibleFaceId : visibleFacesId)
				{
					deleteFace(hull, visibleFaceId);
					faceOutsideSets.erase(faceOutsideSets.find(visibleFaceId));
				}
			};

			return hull;
		}

	private:

		/*! \brief Associated max face and edge id;
		*/
		unsigned int maxFaceId_, maxEdgeId_;

		/*! \brief Initializing first simplex.
		*/
		void init(const std::vector<float3>& pnts, ConvexHull3D& hull);

		/*! \brief Extracting from an outside point coordinate and an initial face id the cone of visible faces.
		*/
		std::unordered_set<FaceId>  makeConeOfVisibleFaces(const ConvexHull3D& hull, FaceId currentFace, const float3& p)
		{
			// Initializing set.
			std::unordered_set<FaceId> visibleFaces;
			visibleFaces.insert(currentFace);

			// Initializing set of visited candidates.
			std::unordered_set<FaceId> visitedCandidates;
			visitedCandidates.insert(currentFace);

			std::stack<FaceId> candidates;
			for(const auto& neighbour : hull.faceNeighbours.at(currentFace))
				candidates.push(neighbour);

			while(!candidates.empty())
			{
				const auto candidateFace = candidates.top();
				candidates.pop();

				visitedCandidates.insert(candidateFace);
				if(Face3DTools::isAbove(p, hull.faces.at(candidateFace)))
					visibleFaces.insert(candidateFace);

				for(const auto& neighbour : hull.faceNeighbours.at(candidateFace))
					if(visitedCandidates.find(neighbour) == visitedCandidates.end())
						candidates.push(neighbour);
			}

			return visibleFaces;
		}

		/*! \brief Building new faces from every edges in the input set of faces and a point.
		*/
		std::unordered_set<FaceId> addFacesFromVisibleFacesBoundary(ConvexHull3D& hull, const std::unordered_set<FaceId>& faces, unsigned int /*pntIdx*/)
		{
			std::unordered_set<FaceId> newFaces;

			// Extracting edge boundaries from faces.
			std::unordered_set<EdgeId> boundary;
			for(const auto& faceId : faces)
			{
				const auto& face = hull.faces[faceId];
				boundary.insert(face.e0);
				boundary.insert(face.e1);
				boundary.insert(face.e2);
			}

			// Creating new edges.
			// To Do !

			// Creating new faces.
			// To Do !

			return newFaces;
		}

		/*!
			\brief Creating a face from three points.
			\params pnts are the input 3D points.
		*/
		/*static Face3D makeFace(const std::vector<float3>& pnts, const std::EdgeId e0, EdgeId e1, EdgeId e2)
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
		}*/

		/*! \brief Deleting a face from its face id while maintaining the adjacency table.
		*/
		void deleteFace(ConvexHull3D& hull, FaceId faceId)
		{
			// Removing face from neighbours.
			for(const auto& neighbourFaceId : hull.faceNeighbours[faceId])
			{
				auto& neighbours = hull.faceNeighbours[neighbourFaceId];
				neighbours.erase(neighbours.find(faceId));
			}

			// Removing associated neighbours.
			hull.faceNeighbours.erase(hull.faceNeighbours.find(faceId));

			// Removing edges only if they did not have any incident faces.
			const auto& face = hull.faces.at(faceId);
			if(hull.edgeIncidentFaces[face.e0].empty()) deleteEdge(hull, face.e0);
			if(hull.edgeIncidentFaces[face.e1].empty()) deleteEdge(hull, face.e1);
			if(hull.edgeIncidentFaces[face.e2].empty()) deleteEdge(hull, face.e2);
			
			// Removing face.
			hull.faces.erase(hull.faces.find(faceId));
		}

		/*! \brief Deleting a edge from its id.
		*/
		void deleteEdge(ConvexHull3D& hull, EdgeId edgeId)
		{
			hull.edgeIncidentFaces.erase(hull.edgeIncidentFaces.find(edgeId));
			hull.edges.erase(hull.edges.find(edgeId));
		}
	};
}

#endif
