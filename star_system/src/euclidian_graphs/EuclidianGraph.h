#include <vector>
#include <algorithm>

struct float2 { float x; float y; };
struct uint2 { size_t i; size_t j; };

namespace EuclidianGraph {

  /*!
    \brief Returns true if two edges are equivalent, i.e. they link the same points.
  */
  static auto areEquivalent(const uint2& e0, const uint2& e1) -> bool
  {
    return (e0.i == e1.i && e0.j == e1.j) || (e0.i == e1.j && e0.j == e1.i);
  }

  /*!
    \brief Extracting square length of a specific edge linking two points.
  */
  static auto getSquareLength(const std::vector<float2>& points, const uint2& e) -> float
  {
    const auto& p0 = points[e.i];
    const auto& p1 = points[e.j];
    const float dx = p1.x - p0.x;
    const float dy = p1.y - p0.y;
    return dx * dx + dy * dy;
  }

  /*!
    \brief Sorting edges in euclidian graph by length.
    \params points are the points coordinates in the graph.
    \params edges are a set of edges linking a sub set of points.
  */
  static auto sortByLength(const std::vector<float2>& points, std::vector<uint2>& edges) -> void
  {
    std::sort(edges.begin(), edges.end(), [&points](const uint2& e0, const uint2& e1)
    {
      return getSquareLength(points, e0) < getSquareLength(points, e1);
    });
  }

  /*!
    \brief Returns true if the union of input external edge with input graph will lead to a cycling graph.
    \params edges are a set of edges in an euclidian graph.
    \params e is the external edge.
  */
  static auto unionsIsCycling(const std::vector<uint2>& edges, const uint2& e) -> bool
  {
	  bool foundI = false, foundJ = false;
	  size_t i = 0;
	  const size_t nEdge = edges.size();

	  while (i < nEdge && (!foundI || !foundJ))
	  {
		  const auto& ei = edges[i];
		  if (!foundI) foundI = ei.i == e.i || ei.j == e.i;
		  if (!foundJ) foundJ = ei.i == e.j || ei.j == e.j;
		  ++i;
	  }

	  return foundI && foundJ;
  }

  /*!
    \brief Add new edge in input euclidan graph if it leads to a non cylcyling graph.
    \params edges are a set of edges linking a sub set of points.
    \params e0 first point index in new edge.
    \params e1 second point index in new edge.
    \returns True if the new edge was indeed added, which means that the graph is not cycling.
  */
  static auto notCyclingAdd(std::vector<uint2>& edges, const uint2& e) -> bool
  {
    if(unionsIsCycling(edges, e)) return false;
    else
    {
      edges.push_back(e);
      return true;
    }
  }

  /*!
    \brief Computing a Euclidian Minimal Spanning Tree (EMST) in the euclidian graph. This procedure
    is implemented using Kruskal's algorithm.
    \params points are the points coordinates in the graph.
    \params edges are the set of edges linking points in graph.
    \returns a sub set of edges representing a (MST) subgraph.
  */
  static auto makeEMST(const std::vector<float2>& points, std::vector<uint2>& edges) -> std::vector<uint2>
  {
    // Sorting edges using their lengths.
    sortByLength(points, edges);

    // Extracting size of EMST.
    const size_t nPts = points.size();
    const size_t nEmstEdge = nPts - 1;

    // Allocating EMST edges.
    std::vector<uint2> emst;
    emst.reserve(nEmstEdge);

    // Applying Krushkal's algorithm.
    size_t i = 0, j = 0;
    const size_t nEdge = edges.size();
    while(i < nEmstEdge && j < nEdge)
    {
      if(notCyclingAdd(emst, edges[j])) ++i;
      ++j;
    }
    emst.shrink_to_fit();

    // Returning built EMST.
    return emst;
  }

  /*!
	\brief Extracting minimal length enabling full browsing of input graph. This implementation is based upon Kruskal's algorithm.
	\params points are the points coordinates in the graph.
	\params edges are the set of edges linking points in graph.
	\returns the minimal browsing length.
  */
  static auto getBrowsingLength(const std::vector<float2>& points, std::vector<uint2>& edges) -> float
  {
	  return std::sqrt(getSquareLength(points, makeEMST(points, edges).back()));
  }
}
