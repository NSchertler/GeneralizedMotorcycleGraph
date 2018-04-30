#pragma once

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"

//Parametrization strategy that uses LSCM for patches without an open boundary and
//scaffold map for all other patches.
class ParametrizationStrategyLSCMForInteriorScaffoldForBoundary : public IParameterizationStrategy
{
public:
	ParametrizationStrategyLSCMForInteriorScaffoldForBoundary(float factor);
	void CalculateParameterization(ParametrizationData& data, float parametrizationErrorThreshold);

private:
	float factor;
};