
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <math.h>

// Type definition.
using uint = unsigned int;
struct uint2 { uint i; uint j; };
struct float2 { float x; float y; };

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Geometric tools
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct Geom2DTools {

  /*! \brief Computing square distance between two points.
  */
  static float getSquareDistance(const float2& p, const float2& q);

  /*! \brief Removing doublons in a set of edges.
  */
  static void removeDoublons(std::vector<uint2>& edges);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Definition of 2D AABBs.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct Box2D {

  /*! \brief Associated center of the box.
  */
	float2 center;

  /*! \brief Associated half lengths of the box in X & Y directions.
  */
  float2 delta;
};

struct Box2DTools {

  /*! \brief Computing box containing every input points plus enlargement.
  */
  static Box2D getBox(const std::vector<float2>& pnts, float enlarge);

  /*! \brief Check if a box contains a point.
  */
  static bool contains(const Box2D& box, const float2& p);

  /*! \brief Check if two boxes are intersecting.
  */
  static bool intersects(const Box2D& box0, const Box2D& box1);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Definition of quad trees.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct QuadTree {

  /*! \brief Associated box.
  */
  Box2D box;

  /*! \brief Associated leaves.
  */
  std::vector<QuadTree> leaves;

  /*! \brief Potentially Associated points coordinates and index. Note that this
    could be more efficiently implemented using utilities such as boost::optional
    but for simplicity we restrict ourselves to a garding bool "hasData" which is
    set to true only if some data point are associated to a quad tree.
  */
  bool hasData;
  uint dataIdx;
};

struct QuadTreeTools {

  /*! \brief Enum used to define position of a data point w.r.t the center of a box.
  */
  enum class PointPositionX { LEFT, RIGHT };
  enum class PointPositionY { BOTTOM, TOP };

  /*! \brief Building data tree from a set of points and a enclosing box.
  */
  static QuadTree make(const Box2D& box, const std::vector<float2>& pnts, std::vector<uint>&& idxs);

  /*! \brief Creating initial set of points indexes used when computed quad tree.
  */
  static std::vector<uint> getInitialIndexes(uint nPts);

  /*! \brief Extracting indexes of points contained within a box.
  */
  static void getPointIndexes(const QuadTree& tree, const std::vector<float2>& pnts, const Box2D& box, std::vector<uint>& neighbours);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Definition of a Union - Find structure.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct UnionFind {

  /*! \brief Associated node parents and ranks.
  */
	std::vector<uint> parent, rank;
};

struct UnionFindTools {

  /*! \brief Creating union find structure.
  */
  static UnionFind make(uint nPts);

  /*! \brief Finding parent of an associated node.
  */
  static uint findParent(UnionFind& uf, uint i);

  /*! \brief merging to set by rank.
  */
  static void mergeByRank(UnionFind& uf, uint i, uint j);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterative algorithm details.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct AlgoTools {

  /*! \brief Definition of the minimal & maximal stencil radius use in the iterative process.
  */
  static const uint MIN_STENCIL_RADIUS = 2;
  static const uint MAX_STENCIL_RADIUS = 10;

  /* \brief Estimating average number of neighbors from an input stencil radius.
  */
  static uint getAverageNNeighbor(uint stencilRadius);

  /* \brief Computing sampling deltas used to build the sampling graph.
  */
  static float2 getSamplingDelta(const float2& dataBoxDelta, uint nPnts, uint stencilRadius = MIN_STENCIL_RADIUS);

  /*
    \brief Computing the sampling graph associated to a set of points.
    \params tree is a quad tree used to accelerate search of points contained in sampling boxes.
    \params pnts are the set of data points.
    \params samplingDelta are the delta of the sampling box used to compute neighors of each points.
    \params nNeighbors is an estimation of the average number of neighbors of a point, typically used to avoid large number of reallocations.
    \returns The sampling graph associated to a set of points and computed using a specific sampling delta.
  */
  static std::vector<uint2> getSamplingGraph(const QuadTree& tree, const std::vector<float2>& pnts, const float2& samplingDelta, uint nNeighbors);

  /*
    \brief Attempt to apply Kruskal algorithm from an input graph.
    \param in is the set of input edges that will be sorted in the process, typically a sampling graph.
    \param out is the set of edges forming the EMST computed using Kruskal's algorithm.
    \returns true if the algorithm was successfull.
  */
  static bool applyKruskal(const std::vector<float2>& pnts, std::vector<uint2>& in, std::vector<uint2>& out);

  /*
    \brief Computing EMST from an input set of points. This iterative algorithm
    is based upon the assumption that the points are uniformly distributed within
    a 2D box.
  */
  static std::vector<uint2> makeEMST(const std::vector<float2>& pnts);
};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Implementations of iterative algorithm details.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------//
uint AlgoTools::getAverageNNeighbor(uint stencilRadius)
{
  return (2 * stencilRadius + 1) * (2 * stencilRadius + 1);
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
float2 AlgoTools::getSamplingDelta(const float2& dataBoxDelta, uint nPnts, uint stencilRadius)
{
  const auto approxStep = (2.0f * dataBoxDelta.x) / sqrtf(float(nPnts));
  return { stencilRadius * approxStep, stencilRadius * approxStep };
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
std::vector<uint2> AlgoTools::getSamplingGraph(const QuadTree& tree, const std::vector<float2>& pnts, const float2& samplingDelta, uint nNeighbors)
{
  // Reserving memory from estimation of the number of neighbors per point.
  std::vector<uint2> samplingGraph;
  samplingGraph.reserve(pnts.size() * nNeighbors);

  // Creating sampling box centered on every point in loop.
  Box2D samplingBox;
  samplingBox.delta = samplingDelta;

  // Loop on every points to extract neighbors from sampling box.
  for (uint i = 0, ni = pnts.size(); i < ni; ++i)
  {
    // Centering sampling box on current point.
    samplingBox.center = pnts[i];

    // Extracting indexes of points inside the sampling box.
    std::vector<uint> samplingNeighbours;
	samplingNeighbours.reserve(nNeighbors);
    QuadTreeTools::getPointIndexes(tree, pnts, samplingBox, samplingNeighbours);

    // Creating edges from current point to every point in the sampling box.
    for (const auto& idx : samplingNeighbours)
      if (idx != i)
        samplingGraph.push_back({ i, idx });
  }

  return samplingGraph;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
bool AlgoTools::applyKruskal(const std::vector<float2>& pnts, std::vector<uint2>& in, std::vector<uint2>& out)
{
  // Sorting edges by square length.
	std::sort(in.begin(), in.end(), [&pnts](const uint2& e0, const uint2& e1) {
		return Geom2DTools::getSquareDistance(pnts[e0.i], pnts[e0.j]) < Geom2DTools::getSquareDistance(pnts[e1.i], pnts[e1.j]);
	});

  // Initializing Union-Find structure.
	auto uf = UnionFindTools::make(pnts.size());

  // Initializing loop variabales.
	uint ie = 0, je = 0, nie = in.size(), nje = pnts.size() - 1;

  // Applying Kruskal's algorithm.
	out.reserve(nje);
	while (ie < nie && je < nje)
	{

    // Extracting edge and associated sets.
		const auto& e = in[ie];
		const auto set0 = UnionFindTools::findParent(uf, e.i);
		const auto set1 = UnionFindTools::findParent(uf, e.j);

    // Pushing edge into otuput graph only of it is not forming a cycle.
		if (set0 != set1)
		{
      out.push_back(e);
      ++je;

      // Merging sets.
			UnionFindTools::mergeByRank(uf, set0, set1);
		}
		++ie;
	}
	out.shrink_to_fit();

  // Kruskal's algorithm formed a complete EMST only if every points are connected.
	return out.size() == nje;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
std::vector<uint2> AlgoTools::makeEMST(const std::vector<float2>& pnts)
{
  // Declaration of EMSt to be built.
  std::vector<uint2> emst;

  // Computing box containing every points plus enlargement.
	const auto box = Box2DTools::getBox(pnts, 0.1);

  // Building quad tree.
	const auto tree = QuadTreeTools::make(box, pnts, QuadTreeTools::getInitialIndexes(pnts.size()));

  // Applying iterative algorithm.
	bool success = false;
	uint stencilRadius = MIN_STENCIL_RADIUS;
	while (!success && stencilRadius < MAX_STENCIL_RADIUS)
	{
    // Clearing previously built & incomplete EMST.
		emst.clear();

    // Computing current sampling delta.
		const auto samplingDelta = AlgoTools::getSamplingDelta(box.delta, pnts.size(), stencilRadius);

    // Computing associated sampling graph.
		auto samplingGraph = AlgoTools::getSamplingGraph(tree, pnts, samplingDelta, AlgoTools::getAverageNNeighbor(stencilRadius));

    // Remoubing doublons in graph.
		Geom2DTools::removeDoublons(samplingGraph);

    // Applying Kruskal's algorithm on sampling graph.
		success = AlgoTools::applyKruskal(pnts, samplingGraph, emst);

    // Increasing sampling stencil radius.
		++stencilRadius;
	}

	return emst;
}
// ---------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Implementations of Union - Find structure.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------//
UnionFind UnionFindTools::make(uint nPts)
{
  UnionFind uf;
  uf.rank = std::vector<uint>(nPts, 0);
  uf.parent = std::vector<uint>(nPts);
  std::iota(uf.parent.begin(), uf.parent.end(), 0);
  return uf;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
uint UnionFindTools::findParent(UnionFind& uf, uint i)
{
  if (i != uf.parent[i]) uf.parent[i] = findParent(uf, uf.parent[i]);
	return uf.parent[i];
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
void UnionFindTools::mergeByRank(UnionFind& uf, uint i, uint j)
{
  i = findParent(uf, i);
  j = findParent(uf, j);

  const auto rki = uf.rank[i];
  const auto rkj = uf.rank[j];

  if (rki > rkj) uf.parent[j] = i;
  else uf.parent[i] = j;

  if(rki == rkj) ++uf.rank[j];
}
// ---------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Implementations of Quad tree tools.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------//
QuadTree QuadTreeTools::make(const Box2D& box, const std::vector<float2>& pnts, std::vector<uint>&& idxs)
{
  QuadTree tree;

  // Setting box associated to current leaf.
  tree.box = box;

  // By default, tree has no associated data point.
  tree.hasData = false;

	if (!idxs.empty())
	{
		if (idxs.size() == 1)
		{
			tree.hasData = true;
			tree.dataIdx = idxs.front();
		}
		else
		{
			std::vector<uint> idxLB, idxRB, idxLT, idxRT;

			for (uint i = 0, ni = idxs.size(); i < ni; ++i)
			{
        // Exrtracting position of current point.
				const auto& p = pnts[idxs[i]];
				const auto posX = p.x <= box.center.x ? PointPositionX::LEFT : PointPositionX::RIGHT;
				const auto posY = p.y <= box.center.y ? PointPositionY::BOTTOM : PointPositionY::TOP;

        // Assigning current index to one of the leaves in tree.
				if (posX == PointPositionX::LEFT && posY == PointPositionY::BOTTOM) idxLB.push_back(idxs[i]);
				else if (posX == PointPositionX::RIGHT && posY == PointPositionY::BOTTOM) idxRB.push_back(idxs[i]);
				else if (posX == PointPositionX::LEFT && posY == PointPositionY::TOP) idxLT.push_back(idxs[i]);
				else if (posX == PointPositionX::RIGHT && posY == PointPositionY::TOP) idxRT.push_back(idxs[i]);
			}

      // Clearing input index vector.
      idxs.clear();

      // Creating leaves boxes.
      const float2 leafDelta = { 0.5f * box.delta.x, 0.5f * box.delta.y };
      const auto boxLB = Box2D{ { box.center.x - leafDelta.x, box.center.y - leafDelta.y }, leafDelta };
      const auto boxRB = Box2D{ { box.center.x + leafDelta.x, box.center.y - leafDelta.y }, leafDelta };
      const auto boxLT = Box2D{ { box.center.x - leafDelta.x, box.center.y + leafDelta.y }, leafDelta };
      const auto boxRT = Box2D{ { box.center.x + leafDelta.x, box.center.y + leafDelta.y }, leafDelta };

      // Creating leaves.
			tree.leaves.push_back(make(boxLB, pnts, std::move(idxLB)));
			tree.leaves.push_back(make(boxRB, pnts, std::move(idxRB)));
			tree.leaves.push_back(make(boxLT, pnts, std::move(idxLT)));
			tree.leaves.push_back(make(boxRT, pnts, std::move(idxRT)));
		}
	}

	return tree;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
std::vector<uint> QuadTreeTools::getInitialIndexes(uint nPts)
{
  std::vector<uint> idx(nPts);
  std::iota(idx.begin(), idx.end(), 0);
  return idx;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
void QuadTreeTools::getPointIndexes(const QuadTree& tree, const std::vector<float2>& pnts, const Box2D& box, std::vector<uint>& neighbours)
{
  if (tree.hasData)
  {
    if(Box2DTools::contains(box, pnts[tree.dataIdx]))
      neighbours.push_back(tree.dataIdx);
  }
  else
  {
    if (!tree.leaves.empty())
      for (const auto& leaf : tree.leaves)
        if(Box2DTools::intersects(box, leaf.box))
          getPointIndexes(leaf, pnts, box, neighbours);
  }
}
// ---------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Implementations of box tools.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------//
Box2D Box2DTools::getBox(const std::vector<float2>& pnts, float enlarge)
{
  const auto minmaxX = std::minmax_element(pnts.begin(), pnts.end(), [](const float2& p, const float2& q) { return p.x < q.x; });
  const auto minmaxY = std::minmax_element(pnts.begin(), pnts.end(), [](const float2& p, const float2& q) { return p.y < q.y; });
  const auto step = std::max(minmaxX.second->x - minmaxX.first->x, minmaxX.second->y - minmaxX.first->y);

  const float2 center = { 0.5f * (minmaxX.first->x + minmaxX.second->x), 0.5f * (minmaxY.first->y + minmaxY.second->y) };
  const float2 delta = { 0.5f * step + enlarge,  0.5f * step + enlarge };

  return { center, delta };
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
bool Box2DTools::contains(const Box2D& box, const float2& p)
{
  const auto vx = p.x - box.center.x;
  const auto vy = p.y - box.center.y;
  return fabs(vx) <= box.delta.x && fabs(vy) <= box.delta.y;
}
// ---------------------------------------------------------------------------//

// ---------------------------------------------------------------------------//
bool Box2DTools::intersects(const Box2D& box0, const Box2D& box1)
{
  const auto vx = box0.center.x - box1.center.x;
  const auto vy = box0.center.y - box1.center.y;
  const auto dx = box0.delta.x + box1.delta.x;
  const auto dy = box0.delta.y + box1.delta.y;
  return fabs(vx) <= dx && fabs(vy) <= dy;
}
// ---------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Implementations of geometric tools.
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------//
float Geom2DTools::getSquareDistance(const float2& p, const float2& q)
{
  const float dx = q.x - p.x;
  const float dy = q.y - p.y;
  return dx * dx + dy * dy;
}

void Geom2DTools::removeDoublons(std::vector<uint2>& edges)
{
  for (auto& e : edges)
    if (e.i > e.j)
      std::swap(e.i, e.j);

  std::sort(edges.begin(), edges.end(), [](const uint2& e0, const uint2& e1)
  {
    return e0.i < e1.i && e0.j < e1.j;
  });

  edges.erase(std::unique(edges.begin(), edges.end(), [](const uint2& e0, const uint2& e1)
  {
    return e0.i == e1.i && e0.j == e1.j;
  }), edges.end());
}
// ---------------------------------------------------------------------------//
