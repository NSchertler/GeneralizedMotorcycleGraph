#pragma once

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"
#include <chrono>

//Represents a strategy that calculates integer arc lengths using a global constrained quadratic energy and Gurobi.
class ArclengthStrategyGurobi : public IArclengthStrategy
{
public:
	//Instantiates the strategy.
	//discreteOptimizationTimeLimit - the time limit after which optimization is stopped; there may or may not be a result after that
	ArclengthStrategyGurobi(const std::chrono::steady_clock::duration& discreteOptimizationTimeLimit);
	
	void CalculateParametricLengths(ParametrizationData& optData);

private:
	//The time limit after which optimization is stopped.
	const std::chrono::steady_clock::duration& discreteOptimizationTimeLimit;
};