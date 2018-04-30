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
	const size_t np = pnts.size();
	
	/// Initializing diagram.
	diagram.reserve(np);
	
	// Initializing first three cells.
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
	
	VCell c0, c1, c2;
	
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
	
	diagram.push_back(std::move(c0));
	diagram.push_back(std::move(c1));
	diagram.push_back(std::move(c2));
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

static const uint AVERAGE_NEDGE_PER_CELL = 6;
static const double REL_GEOM_EPS = 1e-12;

static std::vector<VCell> makeDiagram(const std::vector<float2>& pnts)
{
	const uint np = pnts.size();
	
	// Initializing diagram.
	std::vector<VCell> diagram;
	initDiagram(pnts, diagram);
	
	// Declaration of algo variables.
	uint c0, c0T, c1T, c0B, c1B, ie, ne;
	bool A, B, C, foundT, foundB, reachedInitCell;
	double bpx, bpy, bnx, bny, rx, ry, a, b, c, d, na, nb, nd, n0, dx0, dy0, dx1, dy1;
	double top0x, top0y, bot0x, bot0y, top1x, top1y, bot1x, bot1y;
	
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

		const VCell& cell0 = diagram[c0];
		
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
				dx0 = e.p0x - bpx;
				dy0 = e.p0y - bpy;
				n0 = bnx * bnx + bny * bny;
				
				a = bnx * dx0 + bny * dy0;
				na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
				
				b = e.tx * bnx + e.ty * bny;
				nb = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
				
				A = (a / na) > REL_GEOM_EPS;
				B = (b / nb) > REL_GEOM_EPS;
				C = (fabs(b) / nb) > REL_GEOM_EPS;
				
				if((A && B) || (A && !C))
				{
					newCell0.edges.push_back(e);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(C)
				{
					if(A && !B)
					{
						a /= b;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						dx0 = (rx - bpx);
						dy0 = (ry - bpy);
						a = bnx * dy0 - bny * dx0;
						na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // *
						
						if((a / na) > REL_GEOM_EPS)
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
						
						dx0 = (rx - bpx);
						dy0 = (ry - bpy);
						a = bnx * dy0 - bny * dx0;
						na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // **
						
						if((a / na) > REL_GEOM_EPS)
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
			}
			else
			{
				dx0 = e.p0x - bpx;
				dy0 = e.p0y - bpy;
				dx1 = e.p1x - bpx;
				dy1 = e.p1y - bpy;
				n0 = bnx * bnx + bny * bny;
				
				a = bnx * dx0 + bny * dy0;
				na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
				
				b = bnx * dx1 + bny * dy1;
				nb = std::sqrt(n0 * (dx1 * dx1 + dy1 * dy1));
				
				d = e.tx * bnx + e.ty * bny;
				nd = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
				
				A = (a / na) > REL_GEOM_EPS;
				B = (b / nb) > REL_GEOM_EPS;
				C = (fabs(d) / nd) > REL_GEOM_EPS;
				
				if ((A && B) || (A && !C))
				{
					newCell0.edges.push_back(e);
					newCell0.neighbors.push_back(cell0.neighbors[ie]);
				}
				else if(C)
				{
					if(A && !B)
					{
						a /= d;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						dx0 = (rx - bpx);
						dy0 = (ry - bpy);
						a = bnx * dy0 - bny * dx0;
						na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // ***
						
						if((a / na) > REL_GEOM_EPS)
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
						a /= d;
						rx = e.p0x - a * e.tx;
						ry = e.p0y - a * e.ty;
						
						dx0 = (rx - bpx);
						dy0 = (ry - bpy);
						a = bnx * dy0 - bny * dx0;
						na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // ****
						
						if((a / na) > REL_GEOM_EPS)
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
		diagram[c0] = std::move(newCell0);
		
		// TOP LOOP /////////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		
		while (foundT && !reachedInitCell)
		{
			const VCell& cellT = diagram[c0T];
			
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
					dx0 = e.p0x - top0x;
					dy0 = e.p0y - top0y;
					n0 = bnx * bnx + bny * bny;
					
					a = bnx * dx0 + bny * dy0;
					na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
					
					b = e.tx * bnx + e.ty * bny;
					nb = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
					
					A = (a / na) > REL_GEOM_EPS;
					B = (b / nb) > REL_GEOM_EPS;
					C = (fabs(b) / nb) > REL_GEOM_EPS;

					if ((A && B) || (A && !C))
					{
						newCellT.edges.push_back(e);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(C)
					{
						if(A && !B)
						{
							a /= b;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - top0x);
							dy0 = (ry - top0y);
							
							// a = bnx * dy0 - bny * dx0;
							// na = 1.0; // std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; // std::sqrt(top0x * top0x + top0y * top0y);
							
							if((a / na) > REL_GEOM_EPS)
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
							
							dx0 = (rx - top0x);
							dy0 = (ry - top0y);
							
							// a = bnx * dy0 - bny * dx0;
							// na = 1.0; // std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; // std::sqrt(top0x * top0x + top0y * top0y);
							
							if((a / na) > REL_GEOM_EPS)
							{
								top1x = rx; top1y = ry;
								c1T = cellT.neighbors[ie];
								foundT = true;
							}
							
							newCellT.edges.emplace_back(rx, ry, e.tx, e.ty, true);
							newCellT.neighbors.push_back(cellT.neighbors[ie]);
						}
					}
				}
				else
				{
					dx0 = e.p0x - top0x;
					dy0 = e.p0y - top0y;
					dx1 = e.p1x - top0x;
					dy1 = e.p1y - top0y;
					n0 = bnx * bnx + bny * bny;
					
					a = bnx * dx0 + bny * dy0;
					na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
					
					b = bnx * dx1 + bny * dy1;
					nb = std::sqrt(n0 * (dx1 * dx1 + dy1 * dy1));
					
					d = e.tx * bnx + e.ty * bny;
					nd = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
					
					A = (a / na) > REL_GEOM_EPS;
					B = (b / nb) > REL_GEOM_EPS;
					C = (fabs(d) / nd) > REL_GEOM_EPS;

					if ((A && B) || (A && !C))
					{
						newCellT.edges.push_back(e);
						newCellT.neighbors.push_back(cellT.neighbors[ie]);
					}
					else if(C)
					{
						if(A && !B)
						{
							a /= d;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - top0x);
							dy0 = (ry - top0y);
							
							//a = bnx * dy0 - bny * dx0;
							//na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(top0x * top0x + top0y * top0y);
							
							if((a / na) > REL_GEOM_EPS)
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
							a /= d;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - top0x);
							dy0 = (ry - top0y);
							
							//a = bnx * dy0 - bny * dx0;
							//na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(top0x * top0x + top0y * top0y);
							
							if((a / na) > REL_GEOM_EPS)
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
			}
			
			if(foundT)
			{
				newCellT.edges.emplace_back(top0x, top0y, top1x, top1y);
				newCellT.neighbors.push_back(ip);
				diagram[c0T] = std::move(newCellT);
				
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
				diagram[c0T] = std::move(newCellT);
				
				newCell.edges.emplace_back(top0x, top0y, -bny, bnx, true);
				newCell.neighbors.push_back(c0T);
			}
		}
		
		// BOTTOM LOOP //////////////////////////////
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		
		while (foundB && !reachedInitCell)
		{
			const VCell& cellB = diagram[c0B];
			
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
					dx0 = e.p0x - bot0x;
					dy0 = e.p0y - bot0y;
					n0 = bnx * bnx + bny * bny;
					
					a = bnx * dx0 + bny * dy0;
					na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
					
					b = e.tx * bnx + e.ty * bny;
					nb = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
					
					A = (a / na) > REL_GEOM_EPS;
					B = (b / nb) > REL_GEOM_EPS;
					C = (fabs(b) / nb) > REL_GEOM_EPS;

					if ((A && B) || (A && !C))
					{
						newCellB.edges.push_back(e);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(C)
					{
						if(A && !B)
						{
							a /= b;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - bot0x);
							dy0 = (ry - bot0y);
							
							//a = bny * dx0 - bnx * dy0;
							//na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(bot0x * bot0x + bot0y * bot0y);
							
							if((a / na) > REL_GEOM_EPS)
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
							
							dx0 = (rx - bot0x);
							dy0 = (ry - bot0y);
							
							// a = bny * dx0 - bnx * dy0;
							// na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(bot0x * bot0x + bot0y * bot0y);
							
							if((a / na) > REL_GEOM_EPS)
							{
								bot1x = rx; bot1y = ry;
								c1B = cellB.neighbors[ie];
								foundB = true;
							}
							
							newCellB.edges.emplace_back(rx, ry, e.tx, e.ty, true);
							newCellB.neighbors.push_back(cellB.neighbors[ie]);
						}
					}
				}
				else
				{
					dx0 = e.p0x - bot0x;
					dy0 = e.p0y - bot0y;
					dx1 = e.p1x - bot0x;
					dy1 = e.p1y - bot0y;
					n0 = bnx * bnx + bny * bny;
					
					a = bnx * dx0 + bny * dy0;
					na = std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0));
					
					b = bnx * dx1 + bny * dy1;
					nb = std::sqrt(n0 * (dx1 * dx1 + dy1 * dy1));
					
					d = e.tx * bnx + e.ty * bny;
					nd = std::sqrt(n0 * (e.tx * e.tx + e.ty * e.ty));
					
					A = (a / na) > REL_GEOM_EPS;
					B = (b / nb) > REL_GEOM_EPS;
					C = (fabs(d) / nd) > REL_GEOM_EPS;
					
					if ((A && B) || (A && !C))
					{
						newCellB.edges.push_back(e);
						newCellB.neighbors.push_back(cellB.neighbors[ie]);
					}
					else if(C)
					{
						if(A && !B)
						{
							a /= d;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - bot0x);
							dy0 = (ry - bot0y);
							
							//a = bny * dx0 - bnx * dy0;
							//na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(bot0x * bot0x + bot0y * bot0y);
							
							if((a / na) > REL_GEOM_EPS)
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
							a /= d;
							rx = e.p0x - a * e.tx;
							ry = e.p0y - a * e.ty;
							
							dx0 = (rx - bot0x);
							dy0 = (ry - bot0y);
							
							//a = bny * dx0 - bnx * dy0;
							//na = 1.0; //std::sqrt(n0 * (dx0 * dx0 + dy0 * dy0)); // !
							a = std::sqrt(dx0 * dx0 + dy0 * dy0);
							na = 1.0; //std::sqrt(bot0x * bot0x + bot0y * bot0y);
							
							if((a / na) > REL_GEOM_EPS)
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
			}
			
			if(foundB)
			{
				newCellB.edges.emplace_back(bot0x, bot0y, bot1x, bot1y);
				newCellB.neighbors.push_back(ip);
				diagram[c0B] = std::move(newCellB);
				
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
				diagram[c0B] = std::move(newCellB);
				
				newCell.edges.emplace_back(bot0x, bot0y, bny, -bnx, true);
				newCell.neighbors.push_back(c0B);
			}
		}
		
		// FIN INCREMENTAL ALGO ////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////
		
		diagram.push_back(std::move(newCell));
	}
	
	return diagram;
}

/*!
 \brief Writing edges into VTK file.
 \params pnts are the points coordinates in the graph.
 \params edges are the set of edges linking points in graph.
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

/*!
 \brief Writing edges into VTK file.
 \params pnts are the points coordinates in the graph.
 \params edges are the set of edges linking points in graph.
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
	
	// const auto pnts = generateSquarePoints(500);
	/*const std::vector<float2> pnts = {
		{0.0f, 0.0f},
		{1.0f, 0.0f},
		{0.0f, 1.0f},
		{1.0f, -0.5f},
		{1.0f, 0.25f},
	};
	const auto diagram = makeDiagram(pnts);
	write(diagram, pnts, "voronoi_dbg.vtk");
	writeDual(diagram, pnts, "delaunay_dbg.vtk");*/
	
	for(size_t i = 0; i < 100; ++i)
	{
		const auto pnts = generateSquarePoints(50000);
		// const auto pnts = generateSquarePoints(50);
		const auto diagram = makeDiagram(pnts);
		std::cout << i << std::endl;
		//write(diagram, pnts, "VTK/voronoi_" + std::to_string(i) + ".vtk");
	}
	
	return 0;
}

