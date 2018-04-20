#ifndef _CONVEX_HULL_3D_
#define _CONVEX_HULL_3D_

#include <fstream>
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
		std::unordered_map<FaceId, std::vector<FaceId>> faceNeighbours;

		/*! \brief Mapping from edge id to incident faces.
		*/
		std::unordered_map<EdgeId, std::vector<FaceId>> edgeIncidentFaces;
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

		/*!
			\brief Writing hull in .vtk file.
		*/
		static void write(const std::vector<float3>& pnts, const ConvexHull3D& hull, const std::string& fileName)
		{
			std::ofstream ofs(fileName);

			// Writing header.
			ofs << "# vtk DataFile Version 2.0\n";
			ofs << "Echo 3D scene geometry\n";
			ofs << "ASCII\n";
			ofs << "DATASET UNSTRUCTURED_GRID\n";

			// Writing points.
			ofs << "POINTS " << pnts.size() << " double \n";
			for (const auto& p : pnts)
					ofs << p.x << " " << p.y << " " << p.z << std::endl;

			// Writing triangles.
			ofs << "CELLS " << hull.faces.size() << " " << 4 * hull.faces.size() << std::endl;
			for (const auto& face : hull.faces)
			{
				const auto& edge0 = hull.edges.at(face.second.e0);
				const auto& edge1 = hull.edges.at(face.second.e1);
				const auto p2 = (edge1.i == edge0.i || edge1.i == edge0.j) ? edge1.j : edge1.i;
				ofs << 3 << " " << edge0.i << " " << edge0.j << " " << p2 << std::endl;
			}

			// Writing cell types.
			ofs << "CELL_TYPES " << hull.faces.size() << std::endl;
			for (size_t i = 0, ni = hull.faces.size(); i < ni; ++i) ofs << 5 << std::endl;

			// Writing normals.
			ofs << "CELL_DATA " << hull.faces.size() << std::endl;
			ofs << "VECTORS NORMALS double" << std::endl;
			for (const auto& face : hull.faces)
				ofs << face.second.normal.x << " " << face.second.normal.y << " " << face.second.normal.z << std::endl;

			// Closing file.
			ofs.close();
		}
	};
}

#endif
