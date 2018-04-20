#ifndef _CONVEX_HULL_BUILDER_3D_
#define _CONVEX_HULL_BUILDER_3D_

#include <algorithm>
#include <array>
#include <stack>
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
					if (Face3DTools::isAbove(pnts[i], face.second))
					{
						OutsideSetTools::addPointIndex(faceOutsideSets[face.first], face.second, pnts, i);
						break;
					}

			// Initializing stack of faces to be processed
			std::vector<FaceId> facesInProcess;
			for(const auto& face : hull.faces)
				if(!OutsideSetTools::isEmpty(faceOutsideSets[face.first]))
				 	facesInProcess.push_back(face.first);

			// FOR DEBUG.
			size_t idbg = 0;

			// Core loop of quickhull algo.
			while(!facesInProcess.empty())
			{
				// FOR DEBUG.
				ConvexHull3DTools::write(pnts, hull, "Blank/dbg_" + std::to_string(idbg) + ".vtk");
				++idbg;

				// Poping first face id to process.
				const auto currentFaceId = facesInProcess.back();
				facesInProcess.pop_back();

				// Accessing outside set.
				const auto& outsideSet = faceOutsideSets[currentFaceId];

				// Building set of visible faces from the furthest point.
				const auto visibleFacesId = makeConeOfVisibleFaces(hull, currentFaceId, pnts[outsideSet.furthestPntIdx]);

				// Creating new faces in the hull builder from the boundary of the visible faces.
				const auto newFacesId = addFacesFromVisibleFacesBoundary(hull, pnts, visibleFacesId, outsideSet.furthestPntIdx);

				// Creating outside sets associated to new faces.
				for (const auto& newFaceId : newFacesId)
					faceOutsideSets[newFaceId] = OutsideSet{};

				// Creating outside sets associated to new faces.
				{
					for (const auto& visibleFaceId : visibleFacesId)
						for (const auto& outsidePntIdx : faceOutsideSets[visibleFaceId].pntIdxs)
							for (const auto& newFaceId : newFacesId)
							{
								const auto& newFace = hull.faces[newFaceId];
								if (Face3DTools::isAbove(pnts[outsidePntIdx], newFace))
								{
									OutsideSetTools::addPointIndex(faceOutsideSets[newFaceId], newFace, pnts, outsidePntIdx);
									break;
								}
							}
				}

				// Setting outside set of new faces and updating faces to be processed.
				for (const auto& newFaceId : newFacesId)
					if (!OutsideSetTools::isEmpty(faceOutsideSets[newFaceId]))
						facesInProcess.push_back(newFaceId);

				// Deleting visible faces & associated outside sets.
				for (const auto& visibleFaceId : visibleFacesId)
				{
					// Deleting faces.
					deleteFace(hull, visibleFaceId);

					// Deleting associated outside face.
					faceOutsideSets.erase(faceOutsideSets.find(visibleFaceId));

					// Removing face from faces to be processed.
					const auto it = std::find(facesInProcess.begin(), facesInProcess.end(), visibleFaceId);
					if(it != facesInProcess.end()) facesInProcess.erase(it);
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
		void init(const std::vector<float3>& pnts, ConvexHull3D& hull)
		{
			const auto begin = pnts.begin(); const auto end = pnts.end();
			const auto p0 = 0; // std::distance(begin, std::min_element(begin, end, [](const float3& p, const float3& q) { return p.x < q.x; }));
			const auto p1 = 1; // std::distance(begin, std::max_element(begin, end, [](const float3& p, const float3& q) { return p.x < q.x; }));
			const auto p2 = 2; // std::distance(begin, std::max_element(begin, end, [](const float3& p, const float3& q) { return p.y < q.y; }));
			const auto p3 = 3; // std::distance(begin, std::max_element(begin, end, [](const float3& p, const float3& q) { return p.z < q.z; }));

			// Creating edges.
			const auto e01 = addEdge(hull, p0, p1);
			const auto e12 = addEdge(hull, p1, p2);
			const auto e02 = addEdge(hull, p0, p2);
			const auto e03 = addEdge(hull, p0, p3);
			const auto e13 = addEdge(hull, p1, p3);
			const auto e23 = addEdge(hull, p2, p3);

			// Creating faces.
			const auto f0 = addFace(pnts, hull, e02, e01, e12);
			const auto f1 = addFace(pnts, hull, e12, e13, e23);
			const auto f2 = addFace(pnts, hull, e03, e02, e23);
			const auto f3 = addFace(pnts, hull, e01, e03, e13);

			// Creating face neighbours.
			hull.faceNeighbours[f0] = { f1, f2, f3 };
			hull.faceNeighbours[f1] = { f0, f2, f3 };
			hull.faceNeighbours[f2] = { f0, f1, f3 };
			hull.faceNeighbours[f3] = { f0, f1, f2 };

			// Creating edge incident faces.
			hull.edgeIncidentFaces[e01] = { f0, f3 };
			hull.edgeIncidentFaces[e12] = { f0, f1 };
			hull.edgeIncidentFaces[e02] = { f0, f2 };
			hull.edgeIncidentFaces[e03] = { f2, f3 };
			hull.edgeIncidentFaces[e13] = { f1, f3 };
			hull.edgeIncidentFaces[e23] = { f1, f2 };
		}

		/*! \brief Extracting from an outside point coordinate and an initial face id the cone of visible faces.
		*/
		std::unordered_set<FaceId> makeConeOfVisibleFaces(const ConvexHull3D& hull, FaceId currentFace, const float3& p)
		{
			// Initializing set.
			std::unordered_set<FaceId> visibleFaces;
			visibleFaces.insert(currentFace);

			// Initializing set of visited candidates.
			std::unordered_set<FaceId> visitedCandidates;
			visitedCandidates.insert(currentFace);

			std::stack<FaceId> candidates;
			for (const auto& neighbour : hull.faceNeighbours.at(currentFace))
				candidates.push(neighbour);

			while (!candidates.empty())
			{
				const auto candidateFace = candidates.top();
				candidates.pop();

				visitedCandidates.insert(candidateFace);
				if (Face3DTools::isAbove(p, hull.faces.at(candidateFace)))
					visibleFaces.insert(candidateFace);

				for (const auto& neighbour : hull.faceNeighbours.at(candidateFace))
					if (visitedCandidates.find(neighbour) == visitedCandidates.end())
						candidates.push(neighbour);
			}

			return visibleFaces;
		}

		/*! \brief Building new faces from every edges in the input set of faces and a point.
		*/
		std::vector<FaceId> addFacesFromVisibleFacesBoundary(ConvexHull3D& hull, const std::vector<float3>& pnts, const std::unordered_set<FaceId>& faces, unsigned int pntIdx)
		{
			// Extracting edge boundaries from faces.
			const auto boundary = getBoundary(hull, faces);
			const size_t boundarySz = boundary.size();

			// Creating loop of connected edges.
			std::vector<EdgeId> connectedEdges;
			connectedEdges.reserve(boundary.size());
			connectedEdges.push_back(*boundary.begin());
			for (size_t i = 1; i < boundarySz; ++i)
			{
				// Extracting previous edge id.
				const auto previousEdgeId = connectedEdges.back();

				// Extracting last point in previous edge.
				const auto p1 = hull.edges.at(previousEdgeId).j;

				// Find edge in boundary with first point as p1.
				auto it = std::find_if(boundary.begin(), boundary.end(), [&hull, previousEdgeId, p1](const EdgeId& id)
				{
					if (id != previousEdgeId)
					{
						const auto& edge = hull.edges.at(id);
						return (edge.i == p1 || edge.j == p1);
					}
					else return false;
				});

				// If edge is opposite then reverse edge.
				auto& edge = hull.edges.at(*it);
				if (edge.j == p1) std::swap(edge.i, edge.j);
				connectedEdges.push_back(*it);
			}

			// Extracting loop orientation.
			const auto& e0 = hull.edges.at(connectedEdges[0]);
			const auto& e1 = hull.edges.at(connectedEdges[1]);
			const auto t0 = Tools::makeVec(pnts[e0.i], pnts[e0.j]);
			const auto t1 = Tools::makeVec(pnts[e1.i], pnts[e1.j]);
			const auto n0 = Tools::makeVec(pnts[e0.i], pnts[pntIdx]);
			const bool zplus = Tools::dot(n0, Tools::cross(t0, t1)) > GEOM_EPSILON;

			// Creating set of new edges.
			std::vector<EdgeId> newEdges;
			newEdges.reserve(boundary.size());
			for (const auto& e : connectedEdges)
			{
				const auto pi = hull.edges.at(e).i;
				if(zplus) newEdges.push_back(addEdge(hull, pi, pntIdx));
				else newEdges.push_back(addEdge(hull, pntIdx, pi));
			}

			// Creating new faces.
			std::vector<FaceId> newFaces;
			newFaces.reserve(boundarySz);
			for (size_t i = 0; i < boundarySz-1; ++i)
				newFaces.push_back(addFace(pnts, hull, connectedEdges[i], newEdges[i], newEdges[i + 1]));
			newFaces.push_back(addFace(pnts, hull, connectedEdges.back(), newEdges.back(), newEdges.front()));

			// Updating neighbouring relations.
			for (size_t i = 0; i < boundarySz; ++i)
			{
				// Adding new faces to neihbouring faces.
				for (const auto& face : hull.edgeIncidentFaces[connectedEdges[i]])
					hull.faceNeighbours[face].push_back(newFaces[i]);

				// Updating new faces neighbouring tables by adding incident faces of boundary edges and
				hull.faceNeighbours[newFaces[i]] = hull.edgeIncidentFaces[connectedEdges[i]];

				// Adding new faces in incident faces tables of connected edges.
				hull.edgeIncidentFaces[connectedEdges[i]].push_back(newFaces[i]);
			}

			// Adding new faces to neighbour tables.
			{ 
				// Special treatment for first face.
				auto& faceNeighbours0 = hull.faceNeighbours[newFaces[0]];
				faceNeighbours0.push_back(newFaces.back());
				faceNeighbours0.push_back(newFaces[1]);
			}
			for (size_t i = 1; i < boundarySz - 1; ++i)
			{
				auto& faceNeighbours = hull.faceNeighbours[newFaces[i]];
				faceNeighbours.push_back(newFaces[i + 1]);
				faceNeighbours.push_back(newFaces[i - 1]);
			}
			{
				// Special treatment for last face.
				auto& faceNeighboursN = hull.faceNeighbours[newFaces.back()];
				faceNeighboursN.push_back(newFaces[boundarySz - 2]);
				faceNeighboursN.push_back(newFaces[0]);
			}

			// Creating incident faces tables of new edges.
			hull.edgeIncidentFaces[newEdges[0]] = { newFaces[0], newFaces.back() };
			for (size_t i = 1; i < boundarySz; ++i)
				hull.edgeIncidentFaces[newEdges[i]] = { newFaces[i], newFaces[i - 1] };

			return newFaces;
		}

		/*! \brief Extracting boundary of a connected set of faces.
		*/
		std::vector<EdgeId> getBoundary(const ConvexHull3D& hull, const std::unordered_set<FaceId>& facesId)
		{
			std::unordered_map<EdgeId, unsigned int> count;
			for (const auto& faceId : facesId)
			{
				const auto& face = hull.faces.at(faceId);

				// first edge.
				{
					const auto it = count.find(face.e0);
					if (it == count.end()) count[face.e0] = 1;
					else ++(it->second);
				}

				// second edge.
				{
					const auto it = count.find(face.e1);
					if (it == count.end()) count[face.e1] = 1;
					else ++(it->second);
				}

				// third edge.
				{
					const auto it = count.find(face.e2);
					if (it == count.end()) count[face.e2] = 1;
					else ++(it->second);
				}
			}

			// Computing boundary from occurence counts.
			std::vector<EdgeId> boundary;
			boundary.reserve(3 * facesId.size());
			for (const auto& c : count)
				if (c.second == 1)
					boundary.push_back(c.first);
			boundary.shrink_to_fit();
			
			return boundary;
		}

		/*! \brief Adding a new edge in hull from two point indexes.
		*/
		EdgeId addEdge(ConvexHull3D& hull, unsigned int p0, unsigned int p1)
		{
			EdgeId id = maxEdgeId_;
			hull.edges[id] = { p0, p1 };
			++maxEdgeId_;
			return id;
		}

		/*! \brief Adding a new face in hull from three edges.
		*/
		FaceId addFace(const std::vector<float3>& pnts, ConvexHull3D& hull, EdgeId e0, EdgeId e1, EdgeId e2)
		{
			// Creating id.
			FaceId faceId = maxFaceId_;
			++maxFaceId_;

			// Creating face
			hull.faces[faceId] = Face3D{};
			auto& face = hull.faces.at(faceId);

			// Setting edge ids.
			face.e0 = e0; face.e1 = e1; face.e2 = e2;

			// Extracting pnt index.
			const auto& edge0 = hull.edges[e0];
			const auto& edge1 = hull.edges[e1];
			const auto p2 = (edge1.i == edge0.i || edge1.i == edge0.j) ? edge1.j : edge1.i;

			// Computing center.
			const float coef = float(1.0 / 3.0);
			face.center.x = face.center.y = face.center.z = 0.;
			Tools::addIn(face.center, coef, pnts[edge0.i]);
			Tools::addIn(face.center, coef, pnts[edge0.j]);
			Tools::addIn(face.center, coef, pnts[p2]);

			// Computing normal.
			const auto du = Tools::makeVec(pnts[edge0.i], pnts[edge0.j]);
			const auto dv = Tools::makeVec(pnts[edge1.i], pnts[edge1.j]);
			face.normal = Tools::cross(du, dv);

			return faceId;
		}

		/*! \brief Deleting a face from its face id while maintaining the adjacency table.
		*/
		void deleteFace(ConvexHull3D& hull, FaceId faceId)
		{
			// Removing face from neighbours.
			for(const auto& neighbourFaceId : hull.faceNeighbours[faceId])
			{
				auto& neighbours = hull.faceNeighbours[neighbourFaceId];
				neighbours.erase(std::find(neighbours.begin(), neighbours.end(), faceId));
			}

			// Removing associated neighbour table.
			hull.faceNeighbours.erase(hull.faceNeighbours.find(faceId));

			// Locating face in hull faces.
			const auto it = hull.faces.find(faceId);
			const auto& face = it->second;

			// Removing face from first edge incident table.
			{
				const auto edgeIt = hull.edgeIncidentFaces.find(face.e0);
				auto& incFaces = edgeIt->second;
				incFaces.erase(std::find(incFaces.begin(), incFaces.end(), faceId));
				if (incFaces.empty())
				{
					hull.edgeIncidentFaces.erase(edgeIt);
					hull.edges.erase(hull.edges.find(face.e0));
				}
			}

			// Removing face from second edge incident table.
			{
				const auto edgeIt = hull.edgeIncidentFaces.find(face.e1);
				auto& incFaces = edgeIt->second;
				incFaces.erase(std::find(incFaces.begin(), incFaces.end(), faceId));
				if (incFaces.empty())
				{
					hull.edgeIncidentFaces.erase(edgeIt);
					hull.edges.erase(hull.edges.find(face.e1));
				}
			}

			// Removing face from third edge incident table.
			{
				const auto edgeIt = hull.edgeIncidentFaces.find(face.e2);
				auto& incFaces = edgeIt->second;
				incFaces.erase(std::find(incFaces.begin(), incFaces.end(), faceId));
				if (incFaces.empty())
				{
					hull.edgeIncidentFaces.erase(edgeIt);
					hull.edges.erase(hull.edges.find(face.e2));
				}
			}

			// Removing face.
			hull.faces.erase(it);
		}
	};
}

#endif
