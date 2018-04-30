#pragma once

#include "ParametrizationData.h"

class ITargetlengthStrategy
{
public:
	virtual void SetParametricTargetLength(ParametrizationData& optData) = 0;
};


class IMultiplierStrategy
{
public:
	virtual void CalculateMultipliers(ParametrizationData& optData) = 0;
};


class IArclengthStrategy
{
public:
	virtual void CalculateParametricLengths(ParametrizationData& optData) = 0;
};


class IParameterizationStrategy
{
public:
	//
	// parametrizationErrorThreshold - the patch energy above which more seams are made visible in order to reduce distortion
	virtual void CalculateParameterization(ParametrizationData& data, float parametrizationErrorThreshold) = 0;
};