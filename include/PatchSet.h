#pragma once

#include "common.h"
#include "FencedRegion.h"

#include <map>
#include <vector>

struct EnteringEdge
{
	size_t patchIdx;
	int enteringViaEdgeColor;

	EnteringEdge() { }

	EnteringEdge(size_t patchIdx, int enteringViaEdgeColor)
		: patchIdx(patchIdx), enteringViaEdgeColor(enteringViaEdgeColor)
	{ }
};

struct PatchSet
{
	std::map<HEMesh::HalfedgeHandle, EnteringEdge> enteringEdgeToPatch;

	std::vector<FencedRegion> patches;

	void clear()
	{
		patches.clear();
		enteringEdgeToPatch.clear();
	}
};