#ifndef _CONVEX_HULL_BUILDER_3D_
#define _CONVEX_HULL_BUILDER_3D_

#include <vector>
#include <array>
#include <map>

#include "Face3D.h"
#include "OutsideSet.h"
#include "Tools.h"

namespace Delaunay2D {

	/*! \brief Definition of Face and edge Ids.
	*/
	using FaceId = uint;
	using EdgeId = uint;

	/*! \brief Definition of a convex hull builder.
	*/
	class ConvexHullBuilder3D {

	public:

		/*! \brief Initializing first simplex.
		*/
		std::array<FaceId, 4> init(const std::vector<float3>& pnts);

		/*! \brief Accessing a face from its Id.
		*/
		const Face3D& getFace(FaceId faceId) const
		{
			return faces_.at(faceId);
		}

		/*! brief Accessing neighbours of a face.
		*/
		const std::vector<FaceId>& getNeighbours(FaceId faceId) const
		{
			return faceNeighbours_.at(faceId);
		}

		/*! \brief Building new faces from every edges in the input set of faces and a point.
		*/
		std::vector<FaceId> addFacesFromBoundary(const std::vector<float3>& /*pnts*/, std::vector<FaceId> /*faces*/, uint /*pntIdx*/)
		{
			// To Do !
			return {};
		}

		/*! \brief Deleting a face from its face id while maintaining the adjacency table.
		*/
		void deleteFace(FaceId faceId)
		{
			// Removing face from neighbours.
			for(const auto& neighbourFace : getNeighbours(faceId))
			{
				auto& neighbours = faceNeighbours_[neighbourFace];
				neighbours.erase(std::find(neighbours.begin(), neighbours.end(), faceId));
			}

			// Removing associated neighbours.
			faceNeighbours_.erase(faceNeighbours_.find(faceId));

			// Removing associated outside set.
			faceOutsideSets_.erase(faceOutsideSets_.find(faceId));

			// Removing edges only if they did not have any incident faces.
			const auto& face = getFace(faceId);
			if(edgeIncidentFaces_[face.i].empty()) deleteEdge(face.i);
			if(edgeIncidentFaces_[face.j].empty()) deleteEdge(face.j);
			if(edgeIncidentFaces_[face.k].empty()) deleteEdge(face.k);

			// Removing face.
			faces_.erase(faces_.find(faceId));
		}

		/*! \brief Deleting a edge from its id.
		*/
		void deleteEdge(EdgeId /*edgeId*/)
		{
			// To Do !
		}

		/*! \brief Setting set of outside points index.
		*/
		void setFaceOutsideSet(FaceId faceId, OutsideSet&& outsideSet)
		{
			faceOutsideSets_[faceId] = std::move(outsideSet);
		}

		/*! \brief Accessing outside points associated to a face id.
		*/
		const OutsideSet& getFaceOutsideSet(FaceId faceId) const
		{
			return faceOutsideSets_.at(faceId);
		}

		/*! \brief Accessing built convex hull.
		*/
		std::vector<Face3D> get()
		{
			// To Do !
			return {};
		}

	private:

		/*! \brief Associated faces.
		*/
		std::map<FaceId, Face3D> faces_;

		/*! \brief Mapping from face id to face neighbours.
		*/
		std::map<FaceId, std::vector<FaceId>> faceNeighbours_;

		/*! \brief Mapping from face id to outside set of points.
		*/
		std::map<FaceId, OutsideSet> faceOutsideSets_;

		/*! \brief Associated edges.
		*/
		std::map<EdgeId, uint2> edges_;

		/*! \brief Mapping from edge id to incident faces.
		*/
		std::map<EdgeId, std::vector<FaceId>> edgeIncidentFaces_;

		/*! \brief Mapping from edge id to edge neighbours.
		*/
		std::map<EdgeId, std::vector<EdgeId>> edgeNeighbours_;
	};
}

#endif
