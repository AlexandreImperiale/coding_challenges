#include <memory>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <stack>
#include <unordered_set>

// using uint = size_t;
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

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

struct VEdge {
	double p0x, p0y, p1x, p1y, tx, ty;
	bool isHalfLine;
	
	VEdge(double a, double b, double c, double d)
	: p0x(a), p0y(b), p1x(c), p1y(d), tx(c - a), ty(d - b), isHalfLine(false) {}
	
	VEdge(double a, double b, double c, double d, bool)
	: p0x(a), p0y(b), p1x(0.), p1y(0.), tx(c), ty(d), isHalfLine(true) {}
	
};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static const uint AVERAGE_NEDGE_PER_CELL = 6;

struct VCell {
	
	std::vector<uint> neighbors;
	
	std::vector<VEdge> edges;
	
	VCell()
	{
		neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
		edges.reserve(AVERAGE_NEDGE_PER_CELL);
	}
	
	VCell(VCell&& cell) : neighbors(std::move(cell.neighbors)), edges(std::move(cell.edges)) {}
	VCell& operator=(VCell&& cell)
	{
		neighbors = std::move(cell.neighbors);
		edges = std::move(cell.edges);
		return *this;
	}
};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static uint getNearestNeighbor(const std::vector<VCell>& diagram, const std::vector<float2>& generators, uint idx, uint guess = 0)
{
	const auto& p = generators[idx];
	
	bool found;
	uint i = guess;
	
	do
	{
		uint k = i;
		found = true;
		auto di = squareDistance(p, generators[i]);
		
		for (auto j : diagram[i].neighbors)
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

static void initDiagram(const std::vector<float2>& pnts, std::vector<VCell>& diagram)
{
	const float2& p0 = pnts[0];
	const float2& p1 = pnts[1];
	const float2& p2 = pnts[2];
	
	double p0p1x = p1.x - p0.x;
	double p0p1y = p1.y - p0.y;
	double p0p2x = p2.x - p0.x;
	double p0p2y = p2.y - p0.y;
	double p1p2x = p2.x - p1.x;
	double p1p2y = p2.y - p1.y;
	
	double b01x = 0.5 * (p0.x + p1.x);
	double b01y = 0.5 * (p0.y + p1.y);
	double b02x = 0.5 * (p0.x + p2.x);
	double b02y = 0.5 * (p0.y + p2.y);
	
	double d = p0p1y * p0p2x - p0p1x * p0p2y;
	double a = p0p2x * (b02x - b01x) + p0p2y * (b02y - b01y);
	a /= d;
	double cx = b01x + a * p0p1y;
	double cy = b01y - a * p0p1x;
	double o = d > 0. ? 1.0 : -1.0;
	
	auto& c0 = diagram[0];
	auto& c1 = diagram[1];
	auto& c2 = diagram[2];
	
	c0.edges.emplace_back(cx, cy, -o * p0p1y, o * p0p1x, true);
	c0.edges.emplace_back(cx, cy, o * p0p2y, -o * p0p2x, true);
	c0.neighbors.push_back(1);
	c0.neighbors.push_back(2);
	
	c1.edges.emplace_back(cx, cy, -o * p0p1y, o * p0p1x, true);
	c1.edges.emplace_back(cx, cy, -o * p1p2y, o * p1p2x, true);
	c1.neighbors.push_back(0);
	c1.neighbors.push_back(2);
	
	c2.edges.emplace_back(cx, cy, o * p0p2y, -o * p0p2x, true);
	c2.edges.emplace_back(cx, cy, -o * p1p2y, o * p1p2x, true);
	c2.neighbors.push_back(0);
	c2.neighbors.push_back(1);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static const double GEOM_EPS = 1e-12;

struct CutResult {
	double rx0, ry0, rx1, ry1;
	uint nei0, nei1;
	bool found0, found1;
};

static void cutCell(
		double bpx, double bpy,
		uint cellIdx, 			VCell& cell, 		const float2& q,
		uint cuttingCellIdx, 	VCell& cuttingCell, const float2& p,
		CutResult& cutResult)
{
	VCell newCell;
	
	double bnx = q.x - p.x;
	double bny = q.y - p.y;
	double nbn2 = bnx * bnx + bny * bny;
	
	double dx, dy, a, b, c, na, nb, nc;
	bool A, B;
	size_t ie0;
	
	cutResult.found0 = cutResult.found1 = false;
	
	for(size_t ie = 0, ne = cell.edges.size(); ie < ne; ++ie)
	{
		const VEdge& e = cell.edges[ie];
		
		if(e.isHalfLine)
		{
			dx = e.p0x - bpx;
			dy = e.p0y - bpy;
			a = dx * bnx + dy * bny;
			na = std::sqrt(nbn2 * (dx * dx + dy * dy));
			
			b = e.tx * bnx + e.ty * bny;
			nb = std::sqrt(nbn2 * (e.tx * e.tx + e.ty * e.ty));
			
			A = a / na > GEOM_EPS;
			B = b / nb > GEOM_EPS;
			
			if(fabs(b) / nb < GEOM_EPS)
			{
				if(A)
				{
					newCell.edges.push_back(e);
					newCell.neighbors.push_back(cell.neighbors[ie]);
				}
			}
			else
			{
				if(A && B)
				{
					newCell.edges.push_back(e);
					newCell.neighbors.push_back(cell.neighbors[ie]);
				}
				else if(A && !B)
				{
					a /= b;
					
					if(cutResult.found0)
					{
						cutResult.rx1 = e.p0x - a * e.tx;
						cutResult.ry1 = e.p0y - a * e.ty;
						cutResult.nei1 = cell.neighbors[ie];
						cutResult.found1 = true;
						
						newCell.edges.emplace_back(e.p0x, e.p0y, cutResult.rx1, cutResult.ry1);
						newCell.neighbors.push_back(cutResult.nei1);
					}
					else
					{
						cutResult.rx0 = e.p0x - a * e.tx;
						cutResult.ry0 = e.p0y - a * e.ty;
						cutResult.nei0 = cell.neighbors[ie];
						cutResult.found0 = true;
						
						newCell.edges.emplace_back(e.p0x, e.p0y, cutResult.rx0, cutResult.ry0);
						newCell.neighbors.push_back(cutResult.nei0);
						ie0 = ie;
					}
				}
				else if(!A && B)
				{
					a /= b;
					
					if(cutResult.found0)
					{
						cutResult.rx1 = e.p0x - a * e.tx;
						cutResult.ry1 = e.p0y - a * e.ty;
						cutResult.nei1 = cell.neighbors[ie];
						cutResult.found1 = true;
						
						newCell.edges.emplace_back(cutResult.rx1, cutResult.ry1, e.tx, e.ty, true);
						newCell.neighbors.push_back(cutResult.nei1);
					}
					else
					{
						cutResult.rx0 = e.p0x - a * e.tx;
						cutResult.ry0 = e.p0y - a * e.ty;
						cutResult.nei0 = cell.neighbors[ie];
						cutResult.found0 = true;
						
						newCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, e.tx, e.ty, true);
						newCell.neighbors.push_back(cutResult.nei0);
						ie0 = ie;
					}
				}
			}
		}
		else
		{
			dx = e.p0x - bpx;
			dy = e.p0y - bpy;
			a = dx * bnx + dy * bny;
			na = std::sqrt(nbn2 * (dx * dx + dy * dy));
			
			dx = e.p1x - bpx;
			dy = e.p1y - bpy;
			b = dx * bnx + dy * bny;
			nb = std::sqrt(nbn2 * (dx * dx + dy * dy));
			
			c = e.tx * bnx + e.ty * bny;
			nc = std::sqrt(nbn2 * (bnx * bnx + bny * bny));
			
			A = a / na > GEOM_EPS;
			B = b / nb > GEOM_EPS;
			
			if(fabs(c) / nc < GEOM_EPS)
			{
				if(A)
				{
					newCell.edges.push_back(e);
					newCell.neighbors.push_back(cell.neighbors[ie]);
				}
			}
			else
			{
				if(A && B)
				{
					newCell.edges.push_back(e);
					newCell.neighbors.push_back(cell.neighbors[ie]);
				}
				else if(A && !B)
				{
					a /= c;
					
					if(cutResult.found0)
					{
						cutResult.rx1 = e.p0x - a * e.tx;
						cutResult.ry1 = e.p0y - a * e.ty;
						cutResult.nei1 = cell.neighbors[ie];
						cutResult.found1 = true;
						
						newCell.edges.emplace_back(e.p0x, e.p0y, cutResult.rx1, cutResult.ry1);
						newCell.neighbors.push_back(cutResult.nei1);
					}
					else
					{
						cutResult.rx0 = e.p0x - a * e.tx;
						cutResult.ry0 = e.p0y - a * e.ty;
						cutResult.nei0 = cell.neighbors[ie];
						cutResult.found0 = true;
						
						newCell.edges.emplace_back(e.p0x, e.p0y, cutResult.rx0, cutResult.ry0);
						newCell.neighbors.push_back(cutResult.nei0);
						ie0 = ie;
					}
				}
				else if(!A && B)
				{
					a /= c;
					
					if(cutResult.found0)
					{
						cutResult.rx1 = e.p0x - a * e.tx;
						cutResult.ry1 = e.p0y - a * e.ty;
						cutResult.nei1 = cell.neighbors[ie];
						cutResult.found1 = true;
						
						newCell.edges.emplace_back(cutResult.rx1, cutResult.ry1, e.p1x, e.p1y);
						newCell.neighbors.push_back(cutResult.nei1);
					}
					else
					{
						cutResult.rx0 = e.p0x - a * e.tx;
						cutResult.ry0 = e.p0y - a * e.ty;
						cutResult.nei0 = cell.neighbors[ie];
						cutResult.found0 = true;
						
						newCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, e.p1x, e.p1y);
						newCell.neighbors.push_back(cutResult.nei0);
						ie0 = ie;
					}
				}
			}
		}
	}

	if(cutResult.found1)
	{
		cuttingCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, cutResult.rx1, cutResult.ry1);
		newCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, cutResult.rx1, cutResult.ry1);
	}
	else
	{
		const auto& edge0 = cell.edges[ie0];
		
		nc = edge0.ty * edge0.ty + edge0.tx * edge0.tx;
		
		dx = q.x - edge0.p0x;
		dy = q.y - edge0.p0y;
		a = edge0.ty * dx - edge0.tx * dy;
		na = std::sqrt(nc * (dx * dx + dy * dy));
		
		b = edge0.ty * bny + edge0.tx * bnx;
		nb = std::sqrt(nc * nbn2);
		
		A = a / na > GEOM_EPS;
		B = b / nb > GEOM_EPS;
		
		if(A == B)
		{
			cuttingCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, bny, -bnx, true);
			newCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, bny, -bnx, true);
		}
		else
		{
			cuttingCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, -bny, bnx, true);
			newCell.edges.emplace_back(cutResult.rx0, cutResult.ry0, -bny, bnx, true);
		}
	}
	
	newCell.neighbors.push_back(cuttingCellIdx);
	cuttingCell.neighbors.push_back(cellIdx);
	
	cell = std::move(newCell);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static void analyzeCutResult(const CutResult& cutRes1, uint& currentCell, uint& nextCell, double& rx, double& ry)
{
	if(cutRes1.found0)
	{
		if(cutRes1.found1)
		{
			if(cutRes1.nei1 == currentCell)
			{
				currentCell = nextCell;
				nextCell = cutRes1.nei0;
				rx = cutRes1.rx0;
				ry = cutRes1.ry0;
			}
			else
			{
				currentCell = nextCell;
				nextCell = cutRes1.nei1;
				rx = cutRes1.rx1;
				ry = cutRes1.ry1;
			}
		}
		else
			if (cutRes1.nei0 != currentCell)
				throw std::exception(); // Irrelevant cut result !
	}
	else throw std::exception(); // Can't cut cell !
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static std::vector<VCell> makeDiagram(const std::vector<float2>& pnts)
{
	// Allocating cells in diagram.
	const size_t np = pnts.size();
	std::vector<VCell> diagram;
	diagram.resize(np);
	
	// Initializing diagram.
	initDiagram(pnts, diagram);
	
	// Declartion of cut results.
	CutResult cutRes0, cutRes1;
	uint currentCell, nextCell;
	double rx, ry;
	
	// Loop on points.
	for(uint ip = 3; ip < np; ++ip)
	{
		// Reference to current cutting cell.
		auto& cuttingCell = diagram[ip];
		
		// Extracting nearest cell.
		auto c0 = getNearestNeighbor(diagram, pnts, ip);
		
		// Cutting first cell.
		const float2& p = pnts[ip];
		const float2& q = pnts[c0];
		cutCell(
			0.5 * (q.x + p.x), 	// bpx
			0.5 * (q.y + p.y), 	// bpy
			c0, diagram[c0], q,
			ip, cuttingCell, p,
			cutRes0);
		
		currentCell = c0;
		nextCell = cutRes0.nei0;
		
		// Propagating through first cutting point.
		if(cutRes0.found0)
		{
			rx = cutRes0.rx0;
			ry = cutRes0.ry0;
			
			do
			{
				cutCell(
					rx, ry,
					nextCell, diagram[nextCell], pnts[nextCell],
					ip, cuttingCell, p,
					cutRes1);
				
				analyzeCutResult(cutRes1, currentCell, nextCell, rx, ry);
				
			} while(cutRes1.found1 && nextCell != c0);
		}
		
		// Propagating through second cutting point.
		if(cutRes0.found1 && nextCell != c0)
		{
			currentCell = c0;
			nextCell = cutRes0.nei1;
			rx = cutRes0.rx1;
			ry = cutRes0.ry1;
			
			do
			{
				cutCell(
					rx, ry,
					nextCell, diagram[nextCell], pnts[nextCell],
					ip, cuttingCell, p,
					cutRes1);
				
				analyzeCutResult(cutRes1, currentCell, nextCell, rx, ry);
				
			} while(cutRes1.found1 && nextCell != c0);
		}
	}
	
	return diagram;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

/*!
 \brief Writing edges into VTK file.
 \params pnts are the points coordinates in the graph.
*/
static void write(const std::vector<VCell>& diagram, const std::vector<float2>& pnts, const std::string& fileName, float HalfLineDistance = 2.0)
{
	std::ofstream ofs(fileName);
	
	uint nedge = 0;
	for(const auto& cell : diagram)
		nedge += cell.edges.size();
	
	// Writing header.
	ofs << "# vtk DataFile Version 2.0\n";
	ofs << "Voronoi 2D\n";
	ofs << "ASCII\n";
	ofs << "DATASET UNSTRUCTURED_GRID\n";
	
	// Writing points.
	ofs << "POINTS " << 2 * nedge << " double \n";
	for(const auto& cell : diagram)
		for(const auto& edge : cell.edges)
		{
			if(edge.isHalfLine)
			{
				float c = HalfLineDistance / std::sqrt(edge.tx * edge.tx + edge.ty * edge.ty);
				ofs << edge.p0x << " " << edge.p0y << " " << 0. << std::endl;
				ofs << edge.p0x + c * edge.tx << " " << edge.p0y + c * edge.ty << " " << 0. << std::endl;
			}
			else
			{
				ofs << edge.p0x << " " << edge.p0y << " " << 0. << std::endl;
				ofs << edge.p1x << " " << edge.p1y << " " << 0. << std::endl;
			}
		}
	
	// Writing edges.
	ofs << "CELLS " << nedge << " " << 3 * nedge << std::endl;
	uint ip = 0;
	for(const auto& cell : diagram)
		for(const auto& edge : cell.edges)
		{
			ofs << 2 << " " << ip << " " << ip + 1 << std::endl;
			ip += 2;
		}
	
	// Writing cell types.
	ofs << "CELL_TYPES " << nedge << std::endl;
	for (size_t i = 0; i < nedge; ++i) ofs << 3 << std::endl;
	
	// Closing file.
	ofs.close();
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

/*!
 \brief Writing edges into VTK file.
 \params pnts are the points coordinates in the graph.
*/
static void writeDual(const std::vector<VCell>& diagram, const std::vector<float2>& pnts, const std::string& fileName)
{
	std::ofstream ofs(fileName);
	
	// Writing header.
	ofs << "# vtk DataFile Version 2.0\n";
	ofs << "Voronoi 2D\n";
	ofs << "ASCII\n";
	ofs << "DATASET UNSTRUCTURED_GRID\n";
	
	// Writing points.
	ofs << "POINTS " << pnts.size() << " double \n";
	for(const auto& p : pnts)
		ofs << p.x << " " << p.y << " " << 0. << std::endl;
	
	uint nedge = 0;
	for(const auto& cell : diagram)
		nedge += cell.edges.size();
	
	// Writing edges.
	ofs << "CELLS " << nedge << " " << 3 * nedge << std::endl;
	for (size_t i = 0, ni = diagram.size(); i < ni; ++i)
		for (size_t j = 0, nj = diagram[i].neighbors.size(); j < nj; ++j)
			ofs << 2 << " " << i << " " << diagram[i].neighbors[j] << std::endl;
	
	// Writing cell types.
	ofs << "CELL_TYPES " << nedge << std::endl;
	for (size_t i = 0; i < nedge; ++i) ofs << 3 << std::endl;
	
	// Closing file.
	ofs.close();
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int main() {
	
	// const auto pnts = generateSquarePoints(5000);
	/*const std::vector<float2> pnts = {
		{0.0f, 0.0f},
		{1.0f, 0.0f},
		{0.0f, 1.0f},
		{1.0f, -0.5f},
		{1.0f, 0.25f},
	};*/
	
	/*const auto diagram = makeDiagram(pnts);
	write(diagram, pnts, "voronoi_dbg.vtk");
	writeDual(diagram, pnts, "delaunay_dbg.vtk");*/
	
	for(size_t i = 0; i < 100; ++i)
	{
		const auto pnts = generateSquarePoints(50000);
		// const auto pnts = generateSquarePoints(50);
		const auto diagram = makeDiagram(pnts);
		std::cout << i << std::endl;
		// write(diagram, pnts, "VTK/voronoi_" + std::to_string(i) + ".vtk");
	}
	
	return 0;
}

