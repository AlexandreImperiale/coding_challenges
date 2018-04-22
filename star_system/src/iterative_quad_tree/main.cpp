#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <string>
#include "Solution.h"
#include "Tools.h"

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

	Tools::write(pnts, AlgoTools::makeEMST(pnts), "emst_0.vtk");
}

void test1()
{
	const auto pnts = Tools::generateSquarePoints(10000, 10);
	Tools::write(pnts, AlgoTools::makeEMST(pnts), "emst_1.vtk");
}

void test2()
{
	for (size_t i = 0; i < 50; ++i)
	{
		const auto pnts = Tools::generateSquarePoints(50000, 10);
		AlgoTools::makeEMST(pnts);
		std::cout << i << std::endl;
	}
}

int main()
{
	test1();
	return 0;
}
