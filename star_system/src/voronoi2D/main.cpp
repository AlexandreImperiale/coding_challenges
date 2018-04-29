#include <memory>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <stack>
#include <unordered_set>

using uint = size_t;
struct uint2 { uint i; uint j; };
struct float2 { float x; float y; };

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

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

struct VEdge {
	float2 p0, p1;
	float2 tangent;
	bool isHalfLine;
};

struct VCell {

	uint generator, ne0, ne1;
	
	std::vector<uint> neighbors0, neighbors1;
	
	std::vector<VEdge> edges0, edges1;

};

struct VDiagram {

	std::vector<VCell> cells;
};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static bool intersect(const float2& p, const float2& normal, const VEdge& edge)
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static uint getNearestNeighbor(const VDiagram& diagram, const std::vector<float2>& generators, uint idx, uint guess = 0)
{
	const auto& p = generators[idx];

	bool found;
	uint i = guess;

	do
	{
		uint k = i;
		found = true;
		auto di = squareDistance(p, generators[i]);

		for (auto j : diagram.cells[i].neighbors)
		{
			auto dj = squareDistance(p, generators[j]);
			if (dj < di)
			{
				k = j;
				di = dj;
				found = false;
			}
		}
		i = k;

	} while (!found);

	return i;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void main() {


}
