#include <iostream>
#include <string>
#include "EuclidianGraph.h"

/*! Testing if two floating points values are equal.
*/
static bool areEqual(float a, float b, float eps = 1e-6)
{
	return std::abs(a - b) < eps;
}

/* Helper class for running test functions.
*/
template<typename Test> static void runTest(Test test, const std::string& testName)
{
	if(!test()) std::cout << testName << ": NOT SUCCESSFUL" << std::endl;
	else std::cout << testName << ": SUCCESSFUL" << std::endl;
}

/* Testing equivalent relations between edges.
*/
static bool testAreEquivalent0() { return EuclidianGraph::areEquivalent({0, 1}, {0, 1}); }
static bool testAreEquivalent1() { return EuclidianGraph::areEquivalent({0, 1}, {1, 0}); }
static bool testAreEquivalent2() { return !EuclidianGraph::areEquivalent({0, 1}, {2, 1}); }

/*! Testing computation of square length.
*/
static bool testGetSquareLength0()
{
	const std::vector<float2> points = {
		{float(0.2), float(-0.1)},
		{float(1.0), float(2.3)}
	};

	return areEqual(EuclidianGraph::getSquareLength(points, {0, 1}), float(6.4));
}

/*! Testing sorting set of edges using their length.
*/
static bool testSortByLength0()
{
	const std::vector<float2> points = {
		{ float(0.0), float(0.0) },
		{ float(0.0), float(2.0) },
		{ float(0.0), float(3.0) },
		{ float(0.0), float(3.5) }
	};

	std::vector<uint2> edges = { { 0, 1 }, { 1, 2 }, { 2, 3 } };
	EuclidianGraph::sortByLength(points, edges);

	return
		EuclidianGraph::areEquivalent(edges[0], { 2, 3 }) &&
		EuclidianGraph::areEquivalent(edges[1], { 1, 2 }) &&
		EuclidianGraph::areEquivalent(edges[2], { 0, 1 });
}

/*! Testing recognition of cycling graphs.
*/
static bool testUnionIsCycling0() { return EuclidianGraph::unionsIsCycling({ { 0, 1 },{ 1, 2 } }, { 2, 0 }); }
static bool testUnionIsCycling1() { return EuclidianGraph::unionsIsCycling({ { 0, 1 },{ 0, 2 } }, { 1, 2 }); }
static bool testUnionIsCycling2() { return !EuclidianGraph::unionsIsCycling({ { 0, 1 },{ 1, 2 } }, { 2, 3 }); }

/*! Testing to add edges so that the result is acycling.
*/
static bool testNotCyclingAdd0()
{
	std::vector<uint2> edges = { { 0, 1 },{ 1, 2 } };
	EuclidianGraph::notCyclingAdd(edges, { 2, 0 });
	return edges.size() == 2;
}

static bool testNotCyclingAdd1()
{
	std::vector<uint2> edges = { { 0, 1 },{ 0, 2 } };
	EuclidianGraph::notCyclingAdd(edges, { 2, 1 });
	return edges.size() == 2;
}

static bool testNotCyclingAdd2()
{
	std::vector<uint2> edges = { { 0, 1 },{ 1, 2 } };
	EuclidianGraph::notCyclingAdd(edges, { 2, 3 });
	return edges.size() == 3;
}

/*! Testing construction of EMST.
*/
static bool testMakeEMST0()
{
	const std::vector<float2> points = {
		{ float(0.0), float(0.0) },
		{ float(0.0), float(1.0) },
		{ float(1.0), float(0.0) },
	};

	std::vector<uint2> edges = { { 0, 1 },{ 1, 2 },{ 2, 0 } };

	const auto emst = EuclidianGraph::makeEMST(points, edges);

	return
		EuclidianGraph::areEquivalent(emst[0], { 0, 1 }) &&
		EuclidianGraph::areEquivalent(emst[1], { 2, 0 });
}

/*! Testing computation of minimal browsing length.
*/
static bool testGetBrowingLength0()
{
	const std::vector<float2> points = {
		{ float(0.0), float(0.0) },
		{ float(0.0), float(3.2) },
		{ float(1.0), float(0.0) },
	};

	std::vector<uint2> edges = { { 0, 1 },{ 1, 2 },{ 2, 0 } };

	return areEqual(EuclidianGraph::getBrowsingLength(points, edges), float(3.2));
}

/*! Launching tests for EuclidianGraph tools.
*/
int main()
{
	std::cout << std::endl;

	// Testing equivalent relations between edges.
	runTest(testAreEquivalent0, "areEquivalent0");
	runTest(testAreEquivalent1, "areEquivalent1");
	runTest(testAreEquivalent2, "areEquivalent2");
	std::cout << std::endl;

	// Testing square length computation.
	runTest(testGetSquareLength0, "getSquareLength0");
	std::cout << std::endl;

	// Testing sorting edges by length.
	runTest(testSortByLength0, "sortByLength0");
	std::cout << std::endl;

	// Testing cycle recognition in graphs.
	runTest(testUnionIsCycling0, "UnionIsCycling0");
	runTest(testUnionIsCycling1, "UnionIsCycling1");
	runTest(testUnionIsCycling2, "UnionIsCycling2");
	std::cout << std::endl;

	// Testing not cycling add.
	runTest(testNotCyclingAdd0, "NotCyclingAdd0");
	runTest(testNotCyclingAdd1, "NotCyclingAdd1");
	runTest(testNotCyclingAdd2, "NotCyclingAdd2");
	std::cout << std::endl;

	// Testing building EMST.
	runTest(testMakeEMST0, "MakeEMST0");
	std::cout << std::endl;

	// Testing browsing length computation.
	runTest(testGetBrowingLength0, "GetBrowsingLength0");
	std::cout << std::endl;

	return 0;
}
