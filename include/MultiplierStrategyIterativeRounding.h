#pragma once

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"

//Deprecated
class MultiplierStrategyIterativeRounding : public IMultiplierStrategy
{
public:
	void CalculateMultipliers(ParametrizationData& optData);
};