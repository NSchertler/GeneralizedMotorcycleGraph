#pragma once

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"

//Strategy that sets the target lengths of all halfarcs to their respective number of comprising edges.
class TargetLengthStrategySimple : public ITargetlengthStrategy
{
public:
	TargetLengthStrategySimple()
	{ }

	void SetParametricTargetLength(ParametrizationData& optData)
	{
#pragma omp parallel for
		for (int i = 0; i < optData.halfarcs.size(); ++i)
		{
			optData.parametricHalfarcTargetLengths[i] = optData.halfarcs[i].Length();
		}
	}
};