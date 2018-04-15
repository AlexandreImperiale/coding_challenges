#ifndef _CONVEX_HULL_3D_
#define _CONVEX_HULL_3D_

#include <vector>
#include <array>
#include <map>
#include <algorithm>

struct float3 { float x; float y; float z; };

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
	static float3 makeVec(const float3& p, const float3& q)
	{
		return { q.x - p.x,  q.y - p.y,  q.z - p.z };
	}

	/*! \brief Computing cross product from two input vectors.
	*/
	static float3 cross(const float3& u, const float3& v)
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
	static Face3D makeFace(const std::vector<float3>& pnts, uint p0, uint p1, uint p2)
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
	static std::array<uint2, 3> getEdges(const Face3D& face)
	{
		return {{ { face.p0, face.p1}, { face.p1, face.p2}, { face.p2, face.p0} }};
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

		if( dot > float(0.0) ) return face.normal.z > -ConvexHull3DTools::GEOM_EPSILON;
		else return face.normal.z < ConvexHull3DTools::GEOM_EPSILON;
	}

	/*! \brief Definition of Face and edge Ids.
	*/
	using FaceId = uint;
	using EdgeId = uint;

	/*! \brief Definition face outside set.
	*/
	struct OutsideSet {

		/*! \brief Adding a point in the outside set & updating the index of the
			furthest point.
		*/
		void addPointIndex(const std::vector<float3>& pnts, uint pntIdx);

		/*! \brief Associated index of outside points.
		*/
		std::vector<uint> pntIdxs;

		/*! \brief Index of furthest point.
		*/
		uint furthestPntIdx;
	};

	/*! \brief Definition of a convex hull builder.
	*/
	class ConvexHullBuilder3D {

	public:

		/*! \brief Initializing first simplex.
		*/
		std::array<FaceId, 4> init(const std::vector<float3>& pnts);

		/*! \brief Accessing a face from its Id.
		*/
		const Face& getFace(FaceId faceId) const;

		/*! brief Accessing neighbours of a face.
		*/
		const std::vector<FaceId>& getNeighbours(FaceId faceId) const;

		/*! \brief Building a new face from an edge id and a point index while updating
			the adjacency table.
		*/
		FaceId addFace(const std::vector<flot3>& pnts, EdgeId edgeId, uint pntIdx);

		/*! \brief Deleting a face from its face id, while maintaining the adjacency table.
		*/
		void deleteFace(FaceId faceId);

		/*! \brief Extracting boundaries associated to a set of faces as a set of edge ids.
		*/
		std::vector<EdgeId> getBoundary(const std::vector<FaceId>& faces) const;

		/*! \brief Setting set of outside points index.
		*/
		void setFaceOutsideSet(FaceId faceId, OutsideSet&& outsideSet);

		/*! \brief Accessing outside points associated to a face id.
		*/
		const OutsideSet& getFaceOutsideSet(FaceId faceId) const;

		/*! \brief Accessing built convex hull.
		*/
		std::vector<Face3D> get();
	};

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
				if(isAbove(pnts[i], builder.getFace(newFaces[j])))
					initOutsideSets[j].addPointIndex(pnts, i);

		// Initializing stack of faces to be processed
		std::stack<FaceId> facesInProcess;
		for(size_t i = 0; i < 4; ++i)
		{
			if(!initOutsideSets[i].empty()) facesInProcess.push(initFaces[i]);
			builder.setFaceOutsideSet(initFaces[i], std::move(initOutsideSets[i]));
		}

		// Core loop of quickhull algo.
		while(!facesInProcess.empty())
		{
			// Poping face id.
			const auto currentFace = facesInProcess.pop();

			// Accessing outside set.
			const auto& outsideSet = builder.getFaceOutsideSet(currentFace);

			// Extracting furthest points.
			const auto& furthestPnt = pnts[outsideSet.furthestPntIdx];

			// Building set of visible faces from the furthest point.
			std::vector<FaceID> visibleFaces;
			visibleFaces.push(currentFace);

			std::set<FaceID> visitedCandidates;
			visitedCandidates.insert(currentFace);

			std::stack<FaceID> candidates;
			for(const auto& neighbour : builder.getNeighbours(currentFace))
				candidates.push(neighbour);

			while(!candidates.empty())
			{
				const auto candidateFace = candidates.pop();

				visitedCandidates.push(candidateFace);
				if(isAbove(furthestPnt, builder.getFace(candidateFace)))
					visibleFaces.push(candidateFace);

				for(const auto& neighbour : builder.getNeighbours(candidateFace))
					if(visitedCandidates.find(neighbour) == visitedCandidates.end())
						candidates.push(neighbour);
			}

			// Extracting boundary of the visible faces.
			const auto boundary = builder.getBoundary(visibleFaces);

			// Creating new faces in the hull builder from the boundary of the visible faces.
			std::vector<FaceID> newFaces;
			newFaces.reserve(boundary.size());
			for(const auto& edgeId : boundary)
				newFaces.push_back(builder.addFace(pnts, edgeId, outsideSet.furthestPntIdx));

			// Creating outside sets associated to new faces.
			const size_t nNewFaces = newFaces.size();
			std::vector<OutsideSet> newOutsideSets(nNewFaces);
			for(const auto& visibleFace : visibleFaces)
				for(const auto& outsidePntIdx : builder.getFaceOutsideSet(visibleFace).pntIdxs)
					for(size_t i = 0; i < nNewFaces; ++i)
						if(isAbove(pnts[outsidePntIdx], builder.getFace(newFaces[i])))
							newOutsideSets[i].addPointIndex(pnts, outsidePntIdx);

			// Setting outside set of new faces and updating faces to be processed.
			for(size_t i = 0, ni = newFaces.size(); i < ni; ++i)
			{
				if(!newOutsideSets[i].empty()) facesInProcess.push(newFaces[i]);
				builder.setFaceOutsideSet(newFaces[i], std::move(newOutsideSets[i]));
			}

			// Deleting visible faces.
			for(const auto& face : visibleFaces) builder.deleteFace(face);
		};

		return builder.get();
	}
}

#endif
