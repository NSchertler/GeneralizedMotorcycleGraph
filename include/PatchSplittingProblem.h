#pragma once

#include <vector>
#include <cassert>

//Represents the optimization problem of splitting a patch at potential cuts. The objective
//minimizes isometric distortion in the patch while allowing multipliers between the segments.
//Formally, the objective is:
//  D = (l - l1)^2 + (l - l2)^2 + Sum_i (m_i * l - l_i+1)^2 + (m_i * l - l_i+2)^2
//  minimize over l and m_i
class PatchSplittingProblem
{
public:
	//Initializes the patch splitting problem with the given potential cuts (the entries
	//correspond to the cuts' lengths, including the boundaries)
	PatchSplittingProblem(const std::vector<float>& potentialCuts)
		: potentialCuts(potentialCuts)
	{ 
		assert(potentialCuts.size() >= 2);

		lengthSolutionNumerator = potentialCuts[0] + potentialCuts[1];
		lengthSolutionDenominator = 2;
		multiplierFixed.resize(potentialCuts.size() - 2, false);
		multipliers.resize(potentialCuts.size() - 2);
		fixedMultipliersCount = 0;
	}

	float Length() const { return lengthSolutionNumerator / lengthSolutionDenominator; }

	//Optimizes for the first length and all unfixed multipliers. Returns
	//the energy of the solution.
	float Optimize()
	{
		float length = Length();

		float residual1 = length - potentialCuts[0];
		float residual2 = length - potentialCuts[1];
		float energy = residual1 * residual1 + residual2 * residual2;

		for (int i = 0; i < multiplierFixed.size(); ++i)
		{
			if (!multiplierFixed[i])
				multipliers[i] = (potentialCuts[i + 1] + potentialCuts[i + 2]) / (2 * length);
			residual1 = multipliers[i] * length - potentialCuts[i + 1];
			residual2 = multipliers[i] * length - potentialCuts[i + 2];
			energy += residual1 * residual1 + residual2 * residual2;
		}

		return energy;
	}

	void UnfixMultiplier(int i)
	{
		if (!multiplierFixed[i])
			return;
		--fixedMultipliersCount;
		multiplierFixed[i] = false;
		lengthSolutionDenominator -= 2 * multipliers[i] * multipliers[i];
		lengthSolutionNumerator -= multipliers[i] * (potentialCuts[i + 1] + potentialCuts[i + 2]);
	}

	void FixMultiplier(int i, float value)
	{
		if (multiplierFixed[i])
			UnfixMultiplier(i);
		++fixedMultipliersCount;		
		multiplierFixed[i] = true;
		multipliers[i] = value;
		lengthSolutionDenominator += 2 * multipliers[i] * multipliers[i];
		lengthSolutionNumerator += multipliers[i] * (potentialCuts[i + 1] + potentialCuts[i + 2]);
	}

	size_t NumberOfMultipliers() const { return multipliers.size(); }
	float GetMultiplier(int i) const { return multipliers[i]; }
	bool IsMultiplierFixed(int i) const { return multiplierFixed[i]; }
	bool HasFixedMultipliers() const { return fixedMultipliersCount != 0; }
	bool HasUnfixedMultipliers() const { return fixedMultipliersCount != multipliers.size(); }

private:	
	std::vector<bool> multiplierFixed;
	size_t fixedMultipliersCount;

	std::vector<float> multipliers;
	const std::vector<float>& potentialCuts;

	float lengthSolutionNumerator;
	float lengthSolutionDenominator;
};