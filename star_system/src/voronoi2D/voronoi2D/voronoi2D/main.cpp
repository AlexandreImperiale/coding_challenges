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
	
	VCell() = default;
	VCell(VCell&& cell) : neighbors(std::move(cell.neighbors)), edges(std::move(cell.edges)) {}
	VCell& operator=(VCell&& cell)
	{
		neighbors = std::move(cell.neighbors);
		edges = std::move(cell.edges);
		return *this;
	}
};

struct VDiagram {
	
	std::vector<VCell> cells;
};

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
static const float REL_GEOM_EPS = 1e-6;

static VDiagram makeDiagram(const std::vector<float2>& pnts)
{
	const uint np = pnts.size();
	
	// Initializing diagram.
	VDiagram diagram;
	diagram.cells.reserve(np);
	
	// Initializing first three cells.
	{
		const float2& p0 = pnts[0];
		const float2& p1 = pnts[1];
		const float2& p2 = pnts[2];
		
		float p0p1x = p1.x - p0.x;
		float p0p1y = p1.y - p0.y;
		float p0p2x = p2.x - p0.x;
		float p0p2y = p2.y - p0.y;
		float p1p2x = p2.x - p1.x;
		float p1p2y = p2.y - p1.y;
		
		float b01x = 0.5 * (p0.x + p1.x);
		float b01y = 0.5 * (p0.y + p1.y);
		float b02x = 0.5 * (p0.x + p2.x);
		float b02y = 0.5 * (p0.y + p2.y);
		float b12x = 0.5 * (p1.x + p2.x);
		float b12y = 0.5 * (p1.y + p2.y);
		
		float d = p0p1y * p0p2x - p0p1x * p0p2y;
		float a = p0p2x * (b02x - b01x) + p0p2y * (b02y - b01y);
		a /= d;
		float cx = b01x + a * p0p1y;
		float cy = b01y - a * p0p1x;
		float o = d > 0.f ? -1.0f : 1.0f;
		
		VCell c0, c1, c2;

		c0.edges.emplace_back(cx, cy, -d * p0p1y, d * p0p1x, true);
		c0.edges.emplace_back(cx, cy, d * p0p2y, -d * p0p2x, true);
		c0.neighbors.push_back(1);
		c0.neighbors.push_back(2);
		
		c1.edges.emplace_back(cx, cy, -d * p0p1y, d * p0p1x, true);
		c1.edges.emplace_back(cx, cy, -d * p1p2y, d * p1p2x, true);
		c1.neighbors.push_back(0);
		c1.neighbors.push_back(2);
		
		c2.edges.emplace_back(cx, cy, d * p0p2y, -d * p0p2x, true);
		c2.edges.emplace_back(cx, cy, -d * p1p2y, d * p1p2x, true);
		c2.neighbors.push_back(0);
		c2.neighbors.push_back(1);
		
		diagram.cells.push_back(std::move(c0));
		diagram.cells.push_back(std::move(c1));
		diagram.cells.push_back(std::move(c2));
	}

	// Declaration of algo variables.
	uint c0, c0T, c1T, c0B, c1B, ie, ne;
	bool A, B, foundT, foundB, reachedInitCell;
	float bpx, bpy, bnx, bny, rx, ry, a, b, c;
	float top0x, top0y, bot0x, bot0y, top1x, top1y, bot1x, bot1y;
	
	// Loop on points.
	for(uint ip = 3; ip < np; ++ip)
	{
		VCell newCell;
		newCell.edges.reserve(AVERAGE_NEDGE_PER_CELL);
		newCell.neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
		
		// Extracting nearest cell.
		c0 = getNearestNeighbor(diagram, pnts, ip);
		
		// INIT /////////////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////

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
		foundT = foundB = reachedInitCell = false;
		
		for(ie = 0, ne = cell0.edges.size(); ie < ne; ++ie)
		{
			const VEdge& e = cell0.edges[ie];
			
			if(e.isHalfLine)
			{
				a = bnx * (e.p0x - bpx) + bny * (e.p0y - bpy);
				b = e.tx * bnx + e.ty * bny;
				c = 1.0f / (bnx * bnx + bny * bny);
				
				A = a * c > REL_GEOM_EPS;
				B = b * c > REL_GEOM_EPS;
				
				if(A && B)
				{
					newCell0.edges.push_back(e);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(A && !B)
				{
					a /= b;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( (bnx * (ry - bpy) - bny * (rx - bpx)) * c > REL_GEOM_EPS)
					{
						top0x = rx; top0y = ry;
						c0T = cell0.neighbors[ie];
						foundT = true;
					}
					else
					{
						bot0x = rx; bot0y = ry;
						c0B = cell0.neighbors[ie];
						foundB = true;
					}
					
					newCell0.edges.emplace_back(e.p0x, e.p0y, rx, ry);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(!A && B)
				{
					a /= b;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( (bnx * (ry - bpy) - bny * (rx - bpx)) * c > REL_GEOM_EPS)
					{
						top0x = rx; top0y = ry;
						c0T = cell0.neighbors[ie];
						foundT = true;
					}
					else
					{
						bot0x = rx; bot0y = ry;
						c0B = cell0.neighbors[ie];
						foundB = true;
					}
					
					newCell0.edges.emplace_back(rx, ry, e.tx, e.ty, true);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
			}
			else
			{
				a = bnx * (e.p0x - bpx) + bny * (e.p0y - bpy);
				b = bnx * (e.p1x - bpx) + bny * (e.p1y - bpy);
				c = 1.0f / (bnx * bnx + bny * bny);
				
				A = a * c > REL_GEOM_EPS;
				B = b * c > REL_GEOM_EPS;
				
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
					
					if( (bnx * (ry - bpy) - bny * (rx - bpx)) * c > REL_GEOM_EPS)
					{
						top0x = rx; top0y = ry;
						c0T = cell0.neighbors[ie];
						foundT = true;
					}
					else
					{
						bot0x = rx; bot0y = ry;
						c0B = cell0.neighbors[ie];
						foundB = true;
					}
					
					newCell0.edges.emplace_back(e.p0x, e.p0y, rx, ry);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(!A && B)
				{
					a /= e.tx * bnx + e.ty * bny;
					rx = e.p0x - a * e.tx;
					ry = e.p0y - a * e.ty;
					
					if( (bnx * (ry - bpy) - bny * (rx - bpx)) * c > REL_GEOM_EPS)
					{
						top0x = rx; top0y = ry;
						c0T = cell0.neighbors[ie];
						foundT = true;
					}
					else
					{
						bot0x = rx; bot0y = ry;
						c0B = cell0.neighbors[ie];
						foundB = true;
					}
					
					newCell0.edges.emplace_back(rx, ry, e.p1x, e.p1y);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
			}
		}
		
		if(!foundT)
		{
			newCell0.edges.emplace_back(bot0x, bot0y, -bny, bnx, true);
			newCell.edges.emplace_back(bot0x, bot0y, -bny, bnx, true);
		}
		else if(!foundB)
		{
			newCell0.edges.emplace_back(top0x, top0y, bny, -bnx, true);
			newCell.edges.emplace_back(top0x, top0y, bny, -bnx, true);
		}
		else
		{
			newCell0.edges.emplace_back(bot0x, bot0y, top0x, top0y);
			newCell.edges.emplace_back(bot0x, bot0y, top0x, top0y);
		}
		
		newCell.neighbors.push_back(c0);
		newCell0.neighbors.push_back(ip);
		diagram.cells[c0] = std::move(newCell0);
		
		// TOP LOOP /////////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		
		while (foundT && !reachedInitCell)
		{
			const VCell& cellT = diagram.cells[c0T];
			
			VCell newCellT;
			newCellT.edges.reserve(AVERAGE_NEDGE_PER_CELL);
			newCellT.neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
			
			const float2& q = pnts[c0T];
			
			bnx = q.x - p.x;
			bny = q.y - p.y;
			foundT = false;
			
			for(ie = 0, ne = cellT.edges.size(); ie < ne; ++ie)
			{
				const VEdge& e = cellT.edges[ie];
				
				if(e.isHalfLine)
				{
					a = bnx * (e.p0x - top0x) + bny * (e.p0y - top0y);
					b = e.tx * bnx + e.ty * bny;
					c = 1.0f / (bnx * bnx + bny * bny);
					
					A = a * c > REL_GEOM_EPS;
					B = b * c > REL_GEOM_EPS;
					
					if (A && B)
					{
						newCellT.edges.push_back(e);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(A && !B)
					{
						a /= b;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bnx * (ry - top0y) - bny * (rx - top0x)) * c > REL_GEOM_EPS)
						{
							top1x = rx; top1y = ry;
							c1T = cellT.neighbors[ie];
							foundT = true;
						}

						newCellT.edges.emplace_back(e.p0x, e.p0y, rx, ry);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(!A && B)
					{
						a /= b;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bnx * (ry - top0y) - bny * (rx - top0x)) * c > REL_GEOM_EPS)
						{
							top1x = rx; top1y = ry;
							c1T = cellT.neighbors[ie];
							foundT = true;
						}
						
						newCellT.edges.emplace_back(rx, ry, e.tx, e.ty, true);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
				}
				else
				{
					a = bnx * (e.p0x - top0x) + bny * (e.p0y - top0y);
					b = bnx * (e.p1x - top0x) + bny * (e.p1y - top0y);
					c = 1.0f / (bnx * bnx + bny * bny);
					
					A = a * c > REL_GEOM_EPS;
					B = b * c > REL_GEOM_EPS;
					
					if (A && B)
					{
						newCellT.edges.push_back(e);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(A && !B)
					{
						a /= e.tx * bnx + e.ty * bny;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bnx * (ry - top0y) - bny * (rx - top0x)) * c > REL_GEOM_EPS)
						{
							top1x = rx; top1y = ry;
							c1T = cellT.neighbors[ie];
							foundT = true;
						}
						
						newCellT.edges.emplace_back(e.p0x, e.p0y, rx, ry);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(!A && B)
					{
						a /= e.tx * bnx + e.ty * bny;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bnx * (ry - top0y) - bny * (rx - top0x)) * c > REL_GEOM_EPS)
						{
							top1x = rx; top1y = ry;
							c1T = cellT.neighbors[ie];
							foundT = true;
						}
						
						newCellT.edges.emplace_back(rx, ry, e.p1x, e.p1y);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
				}
			}
			
			if(foundT)
			{
				newCellT.edges.emplace_back(top0x, top0y, top1x, top1y);
				newCellT.neighbors.push_back(ip);
				diagram.cells[c0T] = std::move(newCellT);
				
				newCell.edges.emplace_back(top0x, top0y, top1x, top1y);
				newCell.neighbors.push_back(c0T);
				
				reachedInitCell = c1T == c0;
				c0T = c1T;
				top0x = top1x;
				top0y = top1y;
				
			}
			else
			{
				newCellT.edges.emplace_back(top0x, top0y, -bny, bnx, true);
				newCellT.neighbors.push_back(ip);
				diagram.cells[c0T] = std::move(newCellT);
				
				newCell.edges.emplace_back(top0x, top0y, -bny, bnx, true);
				newCell.neighbors.push_back(c0T);
			}
		}
		
		// BOTTOM LOOP //////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		
		while (foundB && !reachedInitCell)
		{
			const VCell& cellB = diagram.cells[c0B];
			
			VCell newCellB;
			newCellB.edges.reserve(AVERAGE_NEDGE_PER_CELL);
			newCellB.neighbors.reserve(AVERAGE_NEDGE_PER_CELL);
			
			const float2& q = pnts[c0B];
			
			bnx = q.x - p.x;
			bny = q.y - p.y;
			foundB = false;
			
			for(ie = 0, ne = cellB.edges.size(); ie < ne; ++ie)
			{
				const VEdge& e = cellB.edges[ie];
				
				if(e.isHalfLine)
				{
					a = bnx * (e.p0x - bot0x) + bny * (e.p0y - bot0y);
					b = e.tx * bnx + e.ty * bny;
					c = 1.0f / (bnx * bnx + bny * bny);
					
					A = a * c > REL_GEOM_EPS;
					B = b * c > REL_GEOM_EPS;
					
					if (A && B)
					{
						newCellB.edges.push_back(e);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(A && !B)
					{
						a /= b;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bny * (rx - bot0x) - bnx * (ry - bot0y)) * c > REL_GEOM_EPS)
						{
							bot1x = rx; bot1y = ry;
							c1B = cellB.neighbors[ie];
							foundB = true;
						}
						
						newCellB.edges.emplace_back(e.p0x, e.p0y, rx, ry);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(!A && B)
					{
						a /= b;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bny * (rx - bot0x) - bnx * (ry - bot0y)) * c > REL_GEOM_EPS)
						{
							bot1x = rx; bot1y = ry;
							c1B = cellB.neighbors[ie];
							foundB = true;
						}
						
						newCellB.edges.emplace_back(rx, ry, e.tx, e.ty, true);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
				}
				else
				{
					a = bnx * (e.p0x - bot0x) + bny * (e.p0y - bot0y);
					b = bnx * (e.p1x - bot0x) + bny * (e.p1y - bot0y);
					c = 1.0f / (bnx * bnx + bny * bny);
					
					A = a * c > REL_GEOM_EPS;
					B = b * c > REL_GEOM_EPS;
					
					if (A && B)
					{
						newCellB.edges.push_back(e);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(A && !B)
					{
						a /= e.tx * bnx + e.ty * bny;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bny * (rx - bot0x) - bnx * (ry - bot0y)) * c > REL_GEOM_EPS)
						{
							bot1x = rx; bot1y = ry;
							c1B = cellB.neighbors[ie];
							foundB = true;
						}
						
						newCellB.edges.emplace_back(e.p0x, e.p0y, rx, ry);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(!A && B)
					{
						a /= e.tx * bnx + e.ty * bny;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						if( (bny * (rx - bot0x) - bnx * (ry - bot0y)) * c > REL_GEOM_EPS)
						{
							bot1x = rx; bot1y = ry;
							c1B = cellB.neighbors[ie];
							foundB = true;
						}
						
						newCellB.edges.emplace_back(rx, ry, e.p1x, e.p1y);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
				}
			}
			
			if(foundB)
			{
				newCellB.edges.emplace_back(bot0x, bot0y, bot1x, bot1y);
				newCellB.neighbors.push_back(ip);
				diagram.cells[c0B] = std::move(newCellB);
				
				newCell.edges.emplace_back(bot0x, bot0y, bot1x, bot1y);
				newCell.neighbors.push_back(c0B);
				
				reachedInitCell = c1B == c0;
				c0B = c1B;
				bot0x = bot1x;
				bot0y = bot1y;
			}
			else
			{
				newCellB.edges.emplace_back(bot0x, bot0y, bny, -bnx, true);
				newCellB.neighbors.push_back(ip);
				diagram.cells[c0B] = std::move(newCellB);
				
				newCell.edges.emplace_back(bot0x, bot0y, bny, -bnx, true);
				newCell.neighbors.push_back(c0B);
			}
		}
		
		// FIN INCREMENTAL ALGO ////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////
		
		diagram.cells.push_back(std::move(newCell));
	}
	
	return diagram;
}

/*!
 \brief Writing edges into VTK file.
 \params pnts are the points coordinates in the graph.
 \params edges are the set of edges linking points in graph.
 */
static void write(const VDiagram& diagram, const std::vector<float2>& pnts, const std::string& fileName, float HalfLineDistance = 2.0)
{
	std::ofstream ofs(fileName);
	
	uint nedge = 0;
	for(const auto& cell : diagram.cells)
		nedge += cell.edges.size();
	
	// Writing header.
	ofs << "# vtk DataFile Version 2.0\n";
	ofs << "Voronoi 2D\n";
	ofs << "ASCII\n";
	ofs << "DATASET UNSTRUCTURED_GRID\n";
	
	// Writing points.
	ofs << "POINTS " << 2 * nedge + pnts.size() << " double \n";
	for(const auto& cell : diagram.cells)
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
	
	for(const auto& p : pnts)
		ofs << p.x << " " << p.y << " " << 0. << std::endl;
	
	// Writing triangles.
	ofs << "CELLS " << nedge << " " << 3 * nedge << std::endl;
	uint ip = 0;
	for(const auto& cell : diagram.cells)
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

int main() {
	
	const auto pnts = generateSquarePoints(200);
	/*const std::vector<float2> pnts = {
		{0.0f, 0.0f},
		{1.0f, 0.0f},
		{0.0f, 1.0f},
		{1.0f, -0.5f},
		{0.9f, 0.25f},
	};*/
	const auto diagram = makeDiagram(pnts);
	write(diagram, pnts, "voronoi_dbg.vtk");

	return 0;
}

