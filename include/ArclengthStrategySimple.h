#pragma once

#include "Parametrization.h"
#include "ParametrizationInterfaces.h"

#include <nsessentials/util/TimedBlock.h>

class ArclengthStrategySimple : public IArclengthStrategy
{
public:
	void CalculateParametricLengths(ParametrizationData& optData)
	{
		nse::util::TimedBlock b("Using simple strategy to set parametric arc lengths ..");

#pragma omp parallel for
		for (int i = 0; i < optData.arcConstraints.size(); ++i)
			optData.arcConstraints[i].broken = true;

#pragma omp parallel for
		for (int i = 0; i < optData.patches.size(); ++i)
		{
			auto& patch = optData.patches[i];
			if (patch.PatchSides().size() != 4)
				continue;
			float targetSideLengths[4] = { 0 };
			for (int side = 0; side < 4; ++side)
			{
				for (auto arcIdx : patch.PatchSides()[side])
					targetSideLengths[side] += optData.parametricHalfarcTargetLengths[arcIdx];
			}
			for (int side = 0; side < 4; ++side)
			{
				float sideLength = 0.5f * (targetSideLengths[side] + targetSideLengths[(side + 2) % 2]);
				for (auto arcIdx : patch.PatchSides()[side])
					optData.parametricHalfarcLengths[arcIdx] = optData.parametricHalfarcTargetLengths[arcIdx] * sideLength / targetSideLengths[side];
			}
		}
	}
};