#pragma once

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"

#include <nsessentials/util/TimedBlock.h>

class MultiplierStrategyAllInvalid : public IMultiplierStrategy
{
public:
	void CalculateMultipliers(ParametrizationData& optData)
	{
		nse::util::TimedBlock b("Invalidating all arc constraints ..");

#pragma omp parallel for
		for (int i = 0; i < optData.arcConstraints.size(); ++i)
			optData.arcConstraints[i].broken = true;
	}
};