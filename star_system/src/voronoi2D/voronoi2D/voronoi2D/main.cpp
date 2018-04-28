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
 \params pnts are the points coordinates in the graph.
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
	float p0x, p0y, p1x, p1y, tx, ty;
	bool isHalfLine;
	
	VEdge(float a, float b, float c, float d)
	: p0x(a), p0y(b), p1x(c), p1y(d), tx(c - a), ty(d - b), isHalfLine(false) {}
	
	VEdge(float a, float b, float c, float d, bool)
	: p0x(a), p0y(b), p1x(0.f), p1y(0.f), tx(c), ty(d), isHalfLine(true) {}
	
};

struct VCell {
	
	std::vector<uint> neighbors;
	
	std::vector<VEdge> edges;
	
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

static const uint AVERAGE_NEDGE_PER_CELL = 6;

static VDiagram makeDiagram(const std::vector<float2>& pnts)
{
	const uint np = pnts.size();
	
	// Initializing diagram.
	VDiagram diagram;
	diagram.cells.reserve(np);
	
	// Initializing first three cells.
	// To Do !
	
	// Declaration of algo variables.
	uint c0, c0T, c1T, c0B, c1B, ie, ne;
	bool A, B, C;
	float bpx, bpy, bnx, bny, rx, ry, a, b, topx, topy, botx, boty;
	
	// Loop on points.
	for(uint ip = 3; ip < np; ++ip)
	{
		VCell newCell;
		newCell.edges.reserve(AVERAGE_NEDGE_PER_CELL);
		newCell.neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
		
		// INIT /////////////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		c0 = getNearestNeighbor(diagram, pnts, ip);
		
		const VCell& cell0 = diagram.cells[c0];
		
		VCell newCell0;
		newCell0.edges.reserve(AVERAGE_NEDGE_PER_CELL);
		newCell0.neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
		
		const float2& p = pnts[ip];
		const float2& q = pnts[c0];
		
		bpx = 0.5 * (q.x + p.x);
		bpy = 0.5 * (q.y + p.y);
		bnx = q.x - p.x;
		bny = q.y - p.y;
		
		for(ie = 0, ne = cell0.edges.size(); ie < ne; ++ie)
		{
			const VEdge& e = cell0.edges[ie];
			
			if(e.isHalfLine)
			{
				a = bnx * (e.p0x - bpx) + bny * (e.p0y - bpy);
				b = e.tx * bnx + e.ty * bny;
				
				A = a > 0.f;
				B = b > 0.f;
				
				if(A && B)
				{
					newCell0.edges.push_back(e);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if( (A && !B) || (!A && B) )
				{
					a /= b;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( bnx * (ry - bpy) - bny * (rx - bpx) > 0.f)
					{
						topx = rx; topy = ry;
						c0T = cell0.neighbors[ie];
					}
					else
					{
						botx = rx; boty = ry;
						c0B = cell0.neighbors[ie];
					}
					
					newCell0.edges.emplace_back(e.p0x, e.p0y, rx, ry);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
			}
			else
			{
				a = bnx * (e.p0x - bpx) + bny * (e.p0y - bpy);
				b = bnx * (e.p1x - bpx) + bny * (e.p1y - bpy);
				
				A = a > 0.f;
				B = b > 0.f;
				
				if (A && B)
				{
					newCell0.edges.push_back(e);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(A && !B)
				{
					a /= e.tx * bnx + e.ty * bny;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( bnx * (ry - bpy) - bny * (rx - bpx) > 0.f)
					{
						topx = rx; topy = ry;
						c0T = cell0.neighbors[ie];
					}
					else
					{
						botx = rx; boty = ry;
						c0B = cell0.neighbors[ie];
					}
					
					newCell0.edges.emplace_back(e.p0x, e.p0y, rx, ry);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(!A && B)
				{
					a /= e.tx * bnx + e.ty * bny;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( bnx * (ry - bpy) - bny * (rx - bpx) > 0.f)
					{
						topx = rx; topy = ry;
						c0T = cell0.neighbors[ie];
					}
					else
					{
						botx = rx; boty = ry;
						c0B = cell0.neighbors[ie];
					}
					
					newCell0.edges.emplace_back(rx, ry, e.p1x, e.p1y);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
			}
		}
		
		newCell0.edges.emplace_back(botx, boty, topx, topy);
		newCell0.neighbors.push_back(ip);
		diagram.cells[c0] = newCell0;
		
		newCell.edges.emplace_back(botx, boty, topx, topy);
		newCell.neighbors.push_back(c0);
		
		C = true;
		
		// TOP LOOP /////////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		
		while (C)
		{
			const VCell& cell = diagram.cells[c0T];
			
			for(ie = 0, ne = cell.edges.size(); ie < ne; ++ie)
			{
				
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int main() {
	
	return 0;
}

