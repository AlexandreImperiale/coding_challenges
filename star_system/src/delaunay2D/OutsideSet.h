#ifndef _OUTSIDESET_
#define _OUTSIDESET_

#include <vector>
#include "Tools.h"
#include "Face3D.h"

namespace Delaunay2D {

	/*! \brief Definition face outside set.
	*/
	struct OutsideSet {

		/*! \brief Default constructors & assignement operator.
		*/
		OutsideSet() = default;
		OutsideSet& operator=(const OutsideSet&) = default;

		/*! \brief Move constructor and move assignement definition.
		*/
		OutsideSet(OutsideSet&& outsideSet)
			: pntIdxs(std::move(outsideSet.pntIdxs)), furthestPntIdx(outsideSet.furthestPntIdx) {}
		OutsideSet& operator=(OutsideSet&& outsideSet)
		{
			pntIdxs = std::move(outsideSet.pntIdxs);
			furthestPntIdx = outsideSet.furthestPntIdx;
			return *this;
		}

		/*! \brief Associated index of outside points.
		*/
		std::vector<unsigned int> pntIdxs;

		/*! \brief Index of furthest point.
		*/
		unsigned int furthestPntIdx;
	};

	struct OutsideSetTools {

		/*! \brief Adding a point in the outside set & updating the index of the furthest point.
		*/
		static void addPointIndex(OutsideSet& outsideSet, const Face3D& face, const std::vector<float3>& pnts, unsigned int pntIdx)
		{
			outsideSet.pntIdxs.push_back(pntIdx);
			if(outsideSet.pntIdxs.size() == 1) outsideSet.furthestPntIdx = pntIdx;
			else
			{
				const float d0 = Tools::squareDistance(pnts[outsideSet.furthestPntIdx], face.center);
				const float d1 = Tools::squareDistance(pnts[pntIdx], face.center);
				if( d0 < d1 ) outsideSet.furthestPntIdx = pntIdx;
			}
		}

		/*! \brief Testing is outside set is empty.
		*/
		static bool isEmpty(const OutsideSet& outsideSet)
		{
			return outsideSet.pntIdxs.empty();
		}
	};
}

#endif
