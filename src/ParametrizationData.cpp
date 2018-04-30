#include "ParametrizationData.h"

FaceConstraintInfo::FaceConstraintInfo(ParametrizationData* optData, int face, ConstraintDirection direction)
	: optData(optData)
{
	int patchSides[2] = { direction == LeftRight ? 0 : 1, direction == LeftRight ? 2 : 3 };

	auto& f = optData->patches[face];
	if (f.PatchSides().size() != 4)
		return;
	if (f.PatchSides()[patchSides[0]].size() == 0 || f.PatchSides()[patchSides[1]].size() == 0)
		return;

	for (int k = 0; k < 2; ++k)
	{
		auto currentSide = patchSides[k];
		canGrow[k] = f.CanSideGrow(currentSide);
		for (int i = 0; i < f.PatchSides()[currentSide].size(); ++i)
		{
			auto arcIdx = f.PatchSides()[currentSide][i];

			auto& arc = optData->graph->Halfarcs()[arcIdx];
			auto arcEdge = optData->graph->MotorcycleHalfedge(*arc.begin());

			bool isBoundary = optData->mesh->is_boundary(arcEdge) || optData->mesh->is_boundary(optData->mesh->opposite_halfedge_handle(arcEdge));			
			
			if(!isBoundary)				
				arcs[k].push_back(arcIdx);
		}		
	}
}