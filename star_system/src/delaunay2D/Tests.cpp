#include <iostream>
#include <string>
#include<random>
#include<cmath>
#include<chrono>

#include "ConvexHullBuilder3D.h"
#include "DelaunayBuilder2D.h"

void test0()
{
	std::vector<float3> pnts = {
		{ float(0), float(0), float(0) },
		{ float(1), float(0), float(0) },
		{ float(0), float(1), float(0) },
		{ float(0), float(0), float(1) },
		{ float(0.3), float(0.3), float(0.3) },
	};

	const auto hull = Delaunay2D::ConvexHullBuilder3D{}.make(pnts);
	Delaunay2D::ConvexHull3DTools::write(pnts, hull, "Blank/test0.vtk");
}

void test1()
{
	std::vector<float3> pnts = {
		{ float(0), float(0), float(0) },
		{ float(1), float(0), float(0) },
		{ float(0), float(1), float(0) },
		{ float(0), float(0), float(1) },
		{ float(0.3), float(2.0), float(0.3) },
	};

	const auto hull = Delaunay2D::ConvexHullBuilder3D{}.make(pnts);
	Delaunay2D::ConvexHull3DTools::write(pnts, hull, "Blank/test1.vtk");
}

std::vector<float3> generateSpherePnts(unsigned int nPts)
{
	const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);
	std::uniform_real_distribution<float> uniform(0.0, 1.0);

	std::vector<float3> pnts;
	pnts.reserve(nPts);

	for (size_t i = 0; i < nPts; ++i) 
	{
		const float theta = 2.0 * M_PI * uniform(generator);
		const float phi = acos(1.0 - 2.0 * uniform(generator));
		pnts.push_back({ sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi) });
	}
	
	return pnts;
}

void test2()
{
	const auto pnts = generateSpherePnts(1000);
	const auto hull = Delaunay2D::ConvexHullBuilder3D{}.make(pnts);
	Delaunay2D::ConvexHull3DTools::write(pnts, hull, "Blank/test2.vtk");
}

std::vector<float3> generateParaboloid(unsigned int nPts)
{
	const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);
	std::uniform_real_distribution<float> uniform(0.0, 1.0);

	std::vector<float2> pnts;
	pnts.reserve(nPts);

	for (size_t i = 0; i < nPts; ++i)
		pnts.push_back({ uniform(generator), uniform(generator) });

	return Delaunay2D::DelaunayBuilder2D::lift(pnts);
}

void test3()
{
	const auto pnts = generateParaboloid(100);
	const auto hull = Delaunay2D::ConvexHullBuilder3D{}.make(pnts);
	Delaunay2D::ConvexHull3DTools::write(pnts, hull, "Blank/test3.vtk");
}

int main()
{
	test3();
	return 0;
}
