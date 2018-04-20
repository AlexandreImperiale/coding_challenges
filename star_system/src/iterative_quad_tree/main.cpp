#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <string>
#include "Tools.h"

static const uint MIN_STENCIL_RADIUS = 2;
static const uint MAX_STENCIL_RADIUS = 10;

struct QuadTree {

	QuadTree() : hasData(false) {}

	float2 center, delta;
	std::vector<QuadTree> leaves;

	bool hasData;
	float2 dataPnt;
	uint dataIdx;
};

enum class PointPositionX { LEFT, RIGHT };
enum class PointPositionY { BOTTOM, TOP };

static QuadTree makeQuadTree(const std::vector<float2>& pnts, const float2& center, const float2& delta, const std::vector<uint>& idx)
{
	QuadTree tree;

	tree.center = center;
	tree.delta = delta;

	if (!pnts.empty())
	{
		if (pnts.size() == 1)
		{
			tree.hasData = true;
			tree.dataPnt = pnts.front();
			tree.dataIdx = idx.front();
		}
		else
		{
			std::vector<float2> pntsLB, pntsRB, pntsLT, pntsRT;
			std::vector<uint> idxLB, idxRB, idxLT, idxRT;

			for (size_t i = 0, ni = pnts.size(); i < ni; ++i)
			{
				const auto& p = pnts[i];
				const auto posX = p.x <= center.x ? PointPositionX::LEFT : PointPositionX::RIGHT;
				const auto posY = p.y <= center.y ? PointPositionY::BOTTOM : PointPositionY::TOP;

				if (posX == PointPositionX::LEFT && posY == PointPositionY::BOTTOM)
				{
					pntsLB.push_back(p);
					idxLB.push_back(idx[i]);
				}
				else if (posX == PointPositionX::RIGHT && posY == PointPositionY::BOTTOM)
				{
					pntsRB.push_back(p);
					idxRB.push_back(idx[i]);
				}
				else if (posX == PointPositionX::LEFT && posY == PointPositionY::TOP)
				{
					pntsLT.push_back(p);
					idxLT.push_back(idx[i]);
				}
				else if (posX == PointPositionX::RIGHT && posY == PointPositionY::TOP)
				{
					pntsRT.push_back(p);
					idxRT.push_back(idx[i]);
				}
			}

			const float2 leafDelta = { 0.5 * delta.x, 0.5 * delta.y };

			tree.leaves.push_back(makeQuadTree(pntsLB, { center.x - leafDelta.x, center.y - leafDelta.y }, leafDelta, idxLB));
			tree.leaves.push_back(makeQuadTree(pntsRB, { center.x + leafDelta.x, center.y - leafDelta.y }, leafDelta, idxRB));
			tree.leaves.push_back(makeQuadTree(pntsLT, { center.x - leafDelta.x, center.y + leafDelta.y }, leafDelta, idxLT));
			tree.leaves.push_back(makeQuadTree(pntsRT, { center.x + leafDelta.x, center.y + leafDelta.y }, leafDelta, idxRT));
		}
	}

	return tree;
}

static void getNeibhours(const QuadTree& tree, const float2& center, const float2& delta, std::vector<uint>& neighbours)
{
	if (tree.hasData)
	{
		const auto vx = tree.dataPnt.x - center.x;
		const auto vy = tree.dataPnt.y - center.y;
		if( fabs(vx) - delta.x <= 0. && fabs(vy) - delta.y <= 0.)
			neighbours.push_back(tree.dataIdx);
	}
	else
	{
		if (!tree.leaves.empty())
		{
			for (const auto& leaf : tree.leaves)
			{
				const auto vx = center.x - leaf.center.x;
				const auto vy = center.y - leaf.center.y;
				const auto dx = leaf.delta.x + delta.x;
				const auto dy = leaf.delta.y + delta.y;

				if (fabs(vx) <= dx && fabs(vy) <= dy)
					getNeibhours(leaf, center, delta, neighbours);
			}
		}
	}
}

static std::array<float2, 2> getBox(const std::vector<float2>& pnts, float enlarge = 0.1)
{
	const auto minmaxX = std::minmax_element(pnts.begin(), pnts.end(), [](const float2& p, const float2& q) { return p.x < q.x; });
	const auto minmaxY = std::minmax_element(pnts.begin(), pnts.end(), [](const float2& p, const float2& q) { return p.y < q.y; });
	const auto step = std::max(minmaxX.second->x - minmaxX.first->x, minmaxX.second->y - minmaxX.first->y);

	const float2 center = { 0.5 * (minmaxX.first->x + minmaxX.second->x), 0.5 * (minmaxY.first->y + minmaxY.second->y) };
	const float2 delta = { 0.5 * step + enlarge,  0.5 * step + enlarge };

	return { { center, delta } };
}

static float2 getSamplingDelta(const float2& dataBoxDelta, uint nPnts, uint stencilRadius = MIN_STENCIL_RADIUS)
{
	const auto approxStep = static_cast<float>((2.0 * dataBoxDelta.x) / sqrt(nPnts));
	return { stencilRadius * approxStep, stencilRadius * approxStep };
}

static std::vector<uint2> getSamplingEdges(const QuadTree& tree, const std::vector<float2>& pnts, uint pntIdx, const float2& samplingDelta)
{
	std::vector<uint> samplingNeighbours;
	getNeibhours(tree, pnts[pntIdx], samplingDelta, samplingNeighbours);

	std::vector<uint2> samplingEdges;
	samplingEdges.reserve(samplingNeighbours.size() - 1);
	for (const auto& idx : samplingNeighbours)
		if (idx != pntIdx)
			samplingEdges.push_back({ pntIdx, idx });

	return samplingEdges;
}

static std::vector<uint2> getSamplingGraph(const QuadTree& tree, const std::vector<float2>& pnts, const float2& samplingDelta, size_t nNeighbors = (2 * MIN_STENCIL_RADIUS + 1) * (2 * MIN_STENCIL_RADIUS + 1))
{
	std::vector<uint2> samplingGraph;
	samplingGraph.reserve(pnts.size() * nNeighbors);

	for (size_t i = 0, ni = pnts.size(); i < ni; ++i)
	{
		std::vector<uint> samplingNeighbours;
		getNeibhours(tree, pnts[i], samplingDelta, samplingNeighbours);

		for (const auto& idx : samplingNeighbours)
			if (idx != i)
				samplingGraph.push_back({ i, idx });
	}

	return samplingGraph;
}

static std::vector<uint> getInitialIndexes(uint nPts)
{
	std::vector<uint> idx(nPts);
	std::iota(idx.begin(), idx.end(), 0);
	return idx;
}

static void removeDoublons(std::vector<uint2>& edges, uint maxPntIdx)
{
	for (auto& e : edges)
		if (e.i > e.j)
			std::swap(e.i, e.j);

	std::sort(edges.begin(), edges.end(), [maxPntIdx](const uint2& e0, const uint2& e1)
	{
		return (e0.j * maxPntIdx + e0.i) < (e1.j * maxPntIdx + e1.i);
	});

	edges.erase(std::unique(edges.begin(), edges.end(), [](const uint2& e0, const uint2& e1)
	{
		return e0.i == e1.i && e0.j == e1.j;
	}), edges.end());
}

struct UnionFind
{
	std::vector<uint> parent, rank;
};

static UnionFind makeUnionFind(size_t nPts)
{
	UnionFind uf;
	uf.rank = std::vector<uint>(nPts, 0);
	uf.parent = std::vector<uint>(nPts);
	std::iota(uf.parent.begin(), uf.parent.end(), 0);
	return uf;
}

static uint findParent(UnionFind& uf, uint i)
{
	if (i != uf.parent[i]) uf.parent[i] = findParent(uf, uf.parent[i]);
	return uf.parent[i];
}

static void mergeByRank(UnionFind& uf, uint i, uint j)
{
	i = findParent(uf, i);
	j = findParent(uf, j);

	const auto rki = uf.rank[i];
	const auto rkj = uf.rank[j];

	if (rki > rkj) uf.parent[j] = i;
	else uf.parent[i] = j;

	if(rki == rkj) ++uf.rank[j];
}

static bool applyKruskal(const std::vector<float2>& pnts, std::vector<uint2>& in, std::vector<uint2>& out)
{
	std::sort(in.begin(), in.end(), [&pnts](const uint2& e0, const uint2& e1) {
		return Tools::squareDistance(pnts[e0.i], pnts[e0.j]) < Tools::squareDistance(pnts[e1.i], pnts[e1.j]);
	});

	auto uf = makeUnionFind(pnts.size());

	size_t ie = 0, je = 0, nie = in.size(), nje = pnts.size() - 1;
	out.reserve(nje);
	while (ie < nie && je < nje)
	{
		const auto& e = in[ie];
		const auto set0 = findParent(uf, e.i);
		const auto set1 = findParent(uf, e.j);
		if (set0 != set1)
		{
			mergeByRank(uf, set0, set1);
			out.push_back(e);
			++je;
		}
		++ie;
	}
	out.shrink_to_fit();

	return out.size() == nje;
}

static std::vector<uint2> makeEMST(const std::vector<float2>& pnts)
{
	std::vector<uint2> emst;
	
	const auto box = getBox(pnts);

	const auto tree = makeQuadTree(pnts, box[0], box[1], getInitialIndexes(pnts.size()));

	bool success = false;
	uint stencilRadius = MIN_STENCIL_RADIUS;
	while (!success & stencilRadius < MAX_STENCIL_RADIUS)
	{
		emst.clear();

		const auto samplingDelta = getSamplingDelta(box[1], pnts.size(), stencilRadius);

		auto samplingGraph = getSamplingGraph(tree, pnts, samplingDelta, (2 * stencilRadius + 1) * (2 * stencilRadius + 1));

		removeDoublons(samplingGraph, pnts.size());

		success = applyKruskal(pnts, samplingGraph, emst);

		++stencilRadius;
	}

	if (!success)
		std::cout << "FAILED TO APPLY KRUSKAL" << std::endl;

	return emst;
}

void test0()
{
	const auto pnts = std::vector<float2>
	{
		{ 0., 0. },
		{ 1., 1. },
		{ 6., 0. },
		{ 5., 1. },
		{ 5., 6. },
		{ 7., 8. },
	};

	const auto box = getBox(pnts);
	const auto tree = makeQuadTree(pnts, box[0], box[1], getInitialIndexes(pnts.size()));
	const auto samplingDelta = getSamplingDelta(box[1], pnts.size(), 3);

	auto samplingGraph = getSamplingGraph(tree, pnts, samplingDelta);
	removeDoublons(samplingGraph, pnts.size());
	Tools::write(pnts, samplingGraph, "Blank/sampling_graph_0.vtk");

	std::vector<uint2> emst;
	if (!applyKruskal(pnts, samplingGraph, emst))
		std::cout << "FAILED TO APPLY KRUSKAL" << std::endl;
	Tools::write(pnts, emst, "Blank/emst_0.vtk");
}

void test1()
{
	const auto pnts = Tools::generateSquarePoints(50000);

	const auto box = getBox(pnts);
	const auto tree = makeQuadTree(pnts, box[0], box[1], getInitialIndexes(pnts.size()));
	const auto samplingDelta = getSamplingDelta(box[1], pnts.size());

	auto samplingGraph = getSamplingGraph(tree, pnts, samplingDelta);
	removeDoublons(samplingGraph, pnts.size());

	std::vector<uint2> emst;
	if (!applyKruskal(pnts, samplingGraph, emst))
		std::cout << "FAILED TO APPLY KRUSKAL" << std::endl;
	Tools::write(pnts, emst, "Blank/emst_1.vtk");
}

void test2()
{
	const auto pnts = Tools::generateSquarePoints(50000, 10);
	Tools::write(pnts, makeEMST(pnts), "Blank/emst_2.vtk");
}

void test3()
{
	for (size_t i = 0; i < 50; ++i)
	{
		const auto pnts = Tools::generateSquarePoints(50000, 10);
		makeEMST(pnts);
		std::cout << i << std::endl;
	}
}

int main()
{
	test3();
	return 0;
}
