#include "ParametrizationStrategyLSCMForInteriorScaffoldForBoundary.h"

#include "ParametrizationHelper.h"

#include <nsessentials/util/TimedBlock.h>
#include <nsessentials/math/LeastSquaresSystem.h>
#include <nsessentials/util/MathematicaFormatter.h>
#include <fstream>

ParametrizationStrategyLSCMForInteriorScaffoldForBoundary::ParametrizationStrategyLSCMForInteriorScaffoldForBoundary(float factor)
	: factor(factor)
{ }

void ParametrizationStrategyLSCMForInteriorScaffoldForBoundary::CalculateParameterization(ParametrizationData& optData, float parametrizationErrorThreshold)
{	
	nse::util::TimedBlock b("Parametrizing patches ..");
//#pragma omp parallel for
	for (int i = 0; i < optData.patches.size(); ++i)
	{
		auto& patch = optData.patches[i];
		auto& texCoords = optData.texCoords[i];

		if (patch.PatchSides().size() != 4)
			continue;

		float patchSize[2] = { 0.0f, 0.0f };
		for (int i = 0; i < 2; ++i)
		{
			auto side = (patch.PatchSides()[i].empty() ? i + 2 : i);
			for (int iArc = 0; iArc < patch.PatchSides()[side].size(); ++iArc)
			{
				auto arcIdx = patch.PatchSides()[side][iArc];
				patchSize[i] += optData.parametricHalfarcLengths[arcIdx];
			}
			if (patchSize[i] == 0)
				patchSize[i] = 1;
		}
		/*if (patchSize[0] == 0.0f && patchSize[1] == 0.0f)
		{
		for (int side = 0; side < 4; ++side)
		{
		bool constraintOnSide = false;
		for (int iArc = 0; iArc < patch.PatchSides()[side].size(); ++iArc)
		{
		auto arcIdx = patch.PatchSides()[side][iArc];
		auto& arc = graph->Halfarcs()[arcIdx];
		patchSize[side % 2] += 0.5f * arc.length;
		}
		}
		}*/

		bool hasOpenBoundary = false;
		int fixedSides = 0;
		for (int side = 0; side < 4; ++side)
		{
			bool sideIsFixed = false;
			for (auto& arcIdx : patch.PatchSides()[side])
			{
				auto& arc = optData.graph->Halfarcs()[arcIdx];
				auto h = optData.graph->MotorcycleHalfedge(arc.segments.front().location);
				if (optData.mesh->is_boundary(optData.mesh->opposite_halfedge_handle(h)))
					hasOpenBoundary = true;
				else
					sideIsFixed = true;
			}
			if (sideIsFixed)
				++fixedSides;
		}

		texCoords.ResetTextureCoordinates();		

		//set up boundary constraints
		auto useSourceLineConstraint = [&](int side, size_t iArc) { return !optData.arcConstraints[patch.PatchSides()[side][iArc] / 2].broken || (fixedSides <= 2 && iArc == 0); };
		auto useTargetLineConstraint = [&](int side, size_t iArc) { return !optData.arcConstraints[patch.PatchSides()[side][iArc] / 2].broken || (fixedSides <= 2 && iArc == patch.PatchSides()[side].size() - 1); };
		auto useInteriorLineConstraint = [&](int side, size_t iArc) { return !optData.arcConstraints[patch.PatchSides()[side][iArc] / 2].broken; };
		SetupBoundaryConstraints(texCoords, patchSize, optData.parametricHalfarcLengths, *optData.mesh, *optData.graph, useSourceLineConstraint, useTargetLineConstraint, useInteriorLineConstraint);
		
		float error = 0;
		if (!hasOpenBoundary)
			error = CalculateLSCM(texCoords, patchSize, *optData.mesh, *optData.graph);
		else
			CalculateScaffoldMap(texCoords, patchSize, *optData.mesh, *optData.graph, factor);		

		if (error > parametrizationErrorThreshold)
		{
			bool deactivated = false;
			int maxSide = 0;
			for (int i = 0; i < 4; ++i)
				if (patch.PatchSides()[i].size() > maxSide)
					maxSide = i;
			for (auto arcIdx : patch.PatchSides()[maxSide])
			{
				if (!optData.arcConstraints[arcIdx / 2].broken)
					deactivated = true;
				optData.arcConstraints[arcIdx / 2].broken = true;
			}
			if(deactivated)
				--i; //calculate this patch again
		}
	}
}