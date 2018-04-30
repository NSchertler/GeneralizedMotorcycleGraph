#include "ParametrizationInterfaces.h"

//Represents a strategy that calculates integer arc lengths using separate quadratic energies per patch and Gurobi.
class ArclengthStrategyGurobiSeparate : public IArclengthStrategy
{
public:
	void CalculateParametricLengths(ParametrizationData& optData);
};