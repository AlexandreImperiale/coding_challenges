#ifndef _TOOLS_
#define _TOOLS_

#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

struct float2 { float x; float y; };
struct float3 { float x; float y; float z; };

using uint = size_t;
struct uint2 { uint i; uint j; };

namespace Tools {

	/*! \brief Definition of Face and edge Ids.
	*/
	using FaceId = unsigned int;
	using EdgeId = unsigned int;

	/*! \brief Floating point geometric epsilon.
	*/
	static const float GEOM_EPSILON = std::numeric_limits<double>::epsilon();


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

	/*! \brief Computing dot product from two input vectors.
	*/
	static float dot(const float3& u, const float3& v)
	{
		return u.x * v.x + u.y * v.y + u.z * v.z;
	}

	/*! \brief Computing square distance between two points.
	*/
	static float squareDistance(const float2& p, const float2& q)
	{
		const float dx = q.x - p.x;
		const float dy = q.y - p.y;
		return dx * dx + dy * dy;
	}

	/*!
		\brief Generating points randomly in unit square.
		\param nPts is the number of points to be generated.
	*/
	std::vector<float2> generateSquarePoints(unsigned int nPts, float radius = 0.5)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> dis(-radius, radius);

		std::vector<float2> pnts;
		pnts.reserve(nPts);
		for (size_t i = 0; i < nPts; ++i)
			pnts.push_back({ dis(gen), dis(gen) });

		return pnts;
	}


	/*!
		\brief Writing edges into VTK file.
		\params points are the points coordinates in the graph.
		\params edges are the set of edges linking points in graph.
	*/
	static void write(const std::vector<float2>& pnts, const std::vector<uint2>& edges, const std::string& fileName)
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
			ofs << p.x << " " << p.y << " " << 0. << std::endl;

		// Writing triangles.
		ofs << "CELLS " << edges.size() << " " << 3 * edges.size() << std::endl;
		for (const auto& edge : edges)
			ofs << 2 << " " << edge.i << " " << edge.j << std::endl;

		// Writing cell types.
		ofs << "CELL_TYPES " << edges.size() << std::endl;
		for (size_t i = 0, ni = edges.size(); i < ni; ++i) ofs << 3 << std::endl;

		// Closing file.
		ofs.close();
	}
}

#endif
