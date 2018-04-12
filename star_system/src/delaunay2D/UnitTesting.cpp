#include <iostream>
#include <string>
#include "Delaunay2D.h"

/*! Testing if two floating points values are equal.
*/
static bool areEqual(float a, float b, float eps = 1e-5)
{
	return std::abs(a - b) < eps;
}

/* Helper class for running test functions.
*/
template<typename Test> static void runTest(Test test, const std::string& testName)
{
	if (!test()) std::cout << testName << ": NOT SUCCESSFUL" << std::endl;
	else std::cout << testName << ": SUCCESSFUL" << std::endl;
}

/*! \brief Testing lift operations.
*/
static bool testLift0()
{
	const auto pnts3D = Delaunay2D::lift({
		{ float(0.00), float(1.2) }
	});

	return areEqual(pnts3D[0].x, float(0.0)) && areEqual(pnts3D[0].y, float(1.2)) && areEqual(pnts3D[0].z, float(1.44));
}

static bool testLift1()
{
	const auto pnts3D = Delaunay2D::lift({
		{ float(0.00), float(1.2) },
		{ float(-2.1), float(3.6) },
	});

	return areEqual(pnts3D[0].z, float(1.44)) && areEqual(pnts3D[1].z, float(17.37));
}

/*! \brief Testing projection of underside faces of a convex hull.
*/
static bool testProjectUnderside0()
{
	const auto pnts3D = std::vector<float3>{
		{ float(0.0), float(0.0), float(1.0) },
		{ float(1.0), float(0.0), float(2.0) },
		{ float(0.0), float(1.0), float(2.0) },
		{ float(-1.), float(0.0), float(2.0) }
	};

	// Creating faces.
	const auto faces = std::vector<ConvexHull3D::Face3D>{
		ConvexHull3D::makeFace(pnts3D, 0, 1, 2),
		ConvexHull3D::makeFace(pnts3D, 0, 2, 3),
		ConvexHull3D::makeFace(pnts3D, 3, 1, 2),
		ConvexHull3D::makeFace(pnts3D, 0, 1, 3),
	};

	// Projecting underside.
	const auto edges = Delaunay2D::projectUnderside(pnts3D, faces);

	return edges.size() == 5;
}

int main()
{
	std::cout << std::endl;

	// Testing lift operations.
	runTest(testLift0, "Lift0");
	runTest(testLift1, "Lift1");
	std::cout << std::endl;

	// Testing projection of underside.
	runTest(testProjectUnderside0, "ProjectUnderside0");

	return 0;
}
