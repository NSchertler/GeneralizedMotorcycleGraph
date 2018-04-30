#include "MultiplierStrategyIterativeRounding.h"

#ifdef WITH_NLOPT
#include <nlopt.hpp>
#endif
#include <nsessentials/util/TimedBlock.h>

#include <iostream>

#define LOG_SPACE_MULTIPLIERS

double MinimizationObjective(unsigned n, const double *x, double *grad, void *_data)
{	
	ParametrizationData* data = static_cast<ParametrizationData*>(_data);

	double objective = 0.0;

	for (int i = 0; i < data->halfarcs.size(); ++i)
	{
		double weight = data->geometricArcFitWeight[i] * data->weightLengthFit;

		double lengthResidual = x[i] - data->parametricHalfarcTargetLengths[i];
		objective += weight * lengthResidual * lengthResidual;
		if (grad)
			grad[i] = 2 * weight * lengthResidual;
	}

	for (int i = data->halfarcs.size(); i < n; ++i)
	{
#ifdef LOG_SPACE_MULTIPLIERS
		double multiplierResidual = x[i];
#else
		double multiplierResidual = x[i] - 1;
#endif
		objective += data->weightArcMultipliers * multiplierResidual * multiplierResidual;
		if (grad)
			grad[i] = 2 * data->weightArcMultipliers * multiplierResidual;

		//#ifdef LOG_SPACE_MULTIPLIERS
		//		double multiplier = std::exp(x[i]);
		//		double multiplierResidual = multiplier - 1;
		//		double multiplierResidualDeriv = multiplier;		
		//#else
		//		double multiplier = x[i];
		//		double multiplierResidual = multiplier - 1;
		//		double multiplierResidualDeriv = 1.0;
		//#endif		
		//		objective += data->weightArcMultipliers * multiplierResidual * multiplierResidual;
		//		if (grad)
		//			grad[i] = 2 * data->weightArcMultipliers * multiplierResidual * multiplierResidualDeriv;
	}

	return objective;
}

double FaceConstraint(unsigned n, const double* x, double* grad, void* data)
{
	if (grad)
		memset(grad, 0, sizeof(double) * n);
	double objective = 0;

	FaceConstraintInfo* info = static_cast<FaceConstraintInfo*>(data);

	for (int k = 0; k < 2; ++k)
	{
		int sign = (k == 0 ? 1 : -1);
		for (auto var : info->arcs[k])
		{
			objective += sign * x[var];
			if (grad)
				grad[var] = sign;
		}
	}

	return objective;
}

double ArcConstraint(unsigned n, const double* x, double* grad, void* data)
{
	if (grad)
		memset(grad, 0, sizeof(double) * n);

	ArcConstraintInfo* info = static_cast<ArcConstraintInfo*>(data);
	if (info->broken)
		return 0.0;

	auto lengthLhs = info->arc;
	auto lengthRhs = info->arc + 1;

#ifdef LOG_SPACE_MULTIPLIERS
	auto multiplierIdxLhs = info->arc / 2 + info->optData->halfarcs.size();
	double multiplierLhs = std::exp(x[multiplierIdxLhs]);
	double multiplierRhs = 1.0;
	double multiplierLhsDeriv = multiplierLhs;
#else
	auto multiplierIdxLhs = info->arc + info->optData->halfarcs.size();
	auto multiplierIdxRhs = info->arc + 1 + info->optData->halfarcs.size();
	double multiplierLhs = x[multiplierIdxLhs];
	double multiplierRhs = x[multiplierIdxRhs];
	double multiplierLhsDeriv = 1.0;
	double multiplierRhsDeriv = 1.0;
#endif

	double objective = multiplierLhs * x[lengthLhs] - multiplierRhs * x[lengthRhs];
	if (grad)
	{
		grad[multiplierIdxLhs] = multiplierLhsDeriv * x[lengthLhs];
		grad[lengthLhs] = multiplierLhs;
		grad[lengthRhs] = -multiplierRhs;
#ifndef LOG_SPACE_MULTIPLIERS
		grad[multiplierIdxRhs] = -multiplierRhsDeriv * x[lengthRhs];
#endif
	}

	return objective;
}

double logRound(double logValue, double& error)
{
	int sign = 1;
	if (logValue < 0)
	{
		sign = -1;
		logValue = -logValue;
	}
	double value = std::exp(logValue);
	double rounded = std::round(value);
	error = rounded - value;
	return sign * std::log(rounded);
}

void MultiplierStrategyIterativeRounding::CalculateMultipliers(ParametrizationData& optData)
{
	for (int i = 0; i < optData.arcConstraints.size(); ++i)
	{
		optData.arcConstraints[i].multiplier = 1;
	}
	return;

#ifdef WITH_NLOPT
	nse::util::TimedBlock b("Calculating arc multipliers using iterative rounding ..");

	bool verbose = true;

	//parameters, for which the rounding error is below this number are rounded immediately
	const double roundErrorThreshold = 0.01;

	//the relative increase of the energy that is allowed after fixing a multiplier (if the 
	//threshold is exceeded, the corresponding arc constraint is disabled)
	const double fixEnergyIncreaseThreshold = 0.01;

#ifdef LOG_SPACE_MULTIPLIERS
	int nVariables = optData.halfarcs.size() + optData.halfarcs.size() / 2;
#else
	int nVariables = optData.halfarcs.size() * 2;
#endif

	std::vector<bool> fixed(nVariables, false);

#ifdef LOG_SPACE_MULTIPLIERS
	optData.weightArcMultipliers = 2.0;
#else
	optData.weightArcMultipliers = 1.0;
#endif
	optData.weightLengthFit = 1.0;

	nlopt::opt opt(nlopt::AUGLAG, nVariables);
	nlopt::opt innerOpt(nlopt::LD_MMA, nVariables);
	innerOpt.set_xtol_rel(1e-4);
	opt.set_local_optimizer(innerOpt);

	opt.set_min_objective(&MinimizationObjective, &optData);
	opt.set_xtol_rel(1e-4);

	for (int i = 0; i < optData.faceConstraints.size(); ++i)
		opt.add_equality_constraint(&FaceConstraint, &optData.faceConstraints[i], 0.01);

	for (int i = 0; i < optData.arcConstraints.size(); ++i)
		opt.add_equality_constraint(&ArcConstraint, &optData.arcConstraints[i], 0.01);

	std::vector<double> lowerBounds(nVariables);
	std::vector<double> upperBounds(nVariables);

	std::vector<double> parameters(nVariables);
#pragma omp parallel for
	for (int i = 0; i < optData.halfarcs.size(); ++i)
	{
		lowerBounds[i] = 0.1;
		upperBounds[i] = 1000;
		parameters[i] = optData.parametricHalfarcTargetLengths[i] + 0.01; //Starting a bit next to the optimal solution helps to avoid singularities
	}
#ifdef LOG_SPACE_MULTIPLIERS
	for (int i = 0; i < optData.halfarcs.size() / 2; ++i)
	{
		lowerBounds[optData.halfarcs.size() + i] = -1.4; //corresponds to multiplier 1/4
		upperBounds[optData.halfarcs.size() + i] = 1.4; //corresponds to multiplier 4
		parameters[optData.halfarcs.size() + i] = 0.0;
	}
#else
	for (int i = 0; i < optData.halfarcs.size(); ++i)
	{
		lowerBounds[optData.halfarcs.size() + i] = 1.0;
		upperBounds[optData.halfarcs.size() + i] = 4.0;
		parameters[optData.halfarcs.size() + i] = 1.0;
	}
#endif

	std::cout << "Optimization problem with " << nVariables << " variables, " << optData.faceConstraints.size() << " face constraints, and " << optData.arcConstraints.size() << " arc constraints." << std::endl;

	int arcLengthsFixed = 0;
	int multipliersFixed = 0;
	int iterations = 0;
	int totalMultipliers = nVariables - optData.halfarcs.size();

	auto& printResult = [&]()
	{
		std::cout << "Parametric lengths: " << std::endl;
		for (int i = 0; i < optData.halfarcs.size(); ++i)
		{
			std::cout << parameters[i] << "  (optimal length: " << optData.parametricHalfarcTargetLengths[i] << ")";
			if (fixed[i])
				std::cout << " x";
			std::cout << std::endl;
		}
		std::cout << "------------------" << std::endl;
		std::cout << "Arc multipliers:" << std::endl;
		for (int i = 0; i < optData.halfarcs.size() / 2; ++i)
		{
			auto idx = optData.halfarcs.size() + i;
			std::cout << parameters[idx];
#ifdef LOG_SPACE_MULTIPLIERS
			std::cout << " --> " << std::exp(parameters[idx]);
#endif
			if (fixed[idx])
				std::cout << " x";
#ifdef LOG_SPACE_MULTIPLIERS
			auto& correspondingArcConstraint = optData.arcConstraints[i];
#else
			auto& correspondingArcConstraint = optData.arcConstraints[i / 2];
#endif
			if (correspondingArcConstraint.broken)
				std::cout << " (deactivated)";
			std::cout << std::endl;
		}		
	};

	opt.set_lower_bounds(lowerBounds);
	opt.set_upper_bounds(upperBounds);

	double oldMinimum;

	try
	{
		if (opt.optimize(parameters, oldMinimum) < 0)
		{
			std::cout << "Optimization failed." << std::endl;
			return;
		}
	}
	catch (std::exception& e)
	{
		std::cout << "Error optimizing: " << e.what() << std::endl;
		return;
	}

	while (multipliersFixed < totalMultipliers)
	{
		++iterations;

		std::vector<int> fixedVariablesInThisIteration;
		int minErrorVariable = 0;
		double minError = std::numeric_limits<double>::infinity();
		for (int i = optData.halfarcs.size(); i < nVariables; ++i)
		{
			if (fixed[i])
				continue;
#ifdef LOG_SPACE_MULTIPLIERS
			double roundError;
			auto rounded = logRound(parameters[i], roundError);
#else
			double rounded = std::round(parameters[i]);
			double roundError = rounded - parameters[i];
#endif			
			if (std::abs(roundError) < minError)
			{
				minError = std::abs(roundError);
				minErrorVariable = i;
			}

			if (std::abs(roundError) < roundErrorThreshold)
			{
				fixed[i] = true;
				lowerBounds[i] = rounded;
				upperBounds[i] = rounded;
				parameters[i] = rounded;
				fixedVariablesInThisIteration.push_back(i);
				if (i < optData.halfarcs.size())
					++arcLengthsFixed;
				else
					++multipliersFixed;
			}
		}

		if (fixedVariablesInThisIteration.empty())
		{
#ifdef LOG_SPACE_MULTIPLIERS
			double error;
			auto rounded = logRound(parameters[minErrorVariable], error);
#else
			auto rounded = std::round(parameters[minErrorVariable]);
			double error = rounded - parameters[minErrorVariable];
#endif
			if(verbose)
				std::cout << "Rounding with an error of " << error << std::endl;

			fixed[minErrorVariable] = true;
			lowerBounds[minErrorVariable] = rounded;
			upperBounds[minErrorVariable] = rounded;
			parameters[minErrorVariable] = rounded;
			fixedVariablesInThisIteration.push_back(minErrorVariable);
			if (minErrorVariable < optData.halfarcs.size())
				++arcLengthsFixed;
			else
				++multipliersFixed;
		}

		opt.set_lower_bounds(lowerBounds);
		opt.set_upper_bounds(upperBounds);

		//std::cout << "Fixed " << fixedVariablesInThisIteration.size() << " variables." << std::endl;
		//std::cout << arcLengthsFixed << " of " << optData.halfarcs.size() << " arc lengths fixed" << std::endl;
		std::cout << multipliersFixed << " of " << totalMultipliers << " multipliers fixed" << std::endl;

		{
			try
			{
				double newMinimum;
				if (opt.optimize(parameters, newMinimum) < 0)
				{
					std::cout << "Optimization failed." << std::endl;
					break;
				}

				auto relativeChange = (newMinimum - oldMinimum) / oldMinimum;
				if(verbose)
					std::cout << "Optimal objective changed from " << oldMinimum << " to " << newMinimum << "(" << relativeChange * 100.0 << " %)" << std::endl;

				if (relativeChange > fixEnergyIncreaseThreshold)
				{
					if(verbose)
						std::cout << "Relative change is too big. Removing arc constraint .." << std::endl;
					for (auto i : fixedVariablesInThisIteration)
						optData.arcConstraints[i - optData.halfarcs.size()].broken = true;
					{
						if (opt.optimize(parameters, newMinimum) < 0)
						{
							std::cout << "Optimization failed." << std::endl;
							break;
						}
					}
					if(verbose)
						std::cout << "New energy after deactivating constraints: " << newMinimum << std::endl;
				}

				oldMinimum = newMinimum;
			}
			catch (std::exception& e)
			{
				std::cout << "Error optimizing: " << e.what() << std::endl;
				break;
			}
		}

		bool writeResult = (iterations - 1) % 50 == 0;
		if(verbose)
			std::cout << "Optimization successful." << std::endl;
		if (writeResult && verbose)
			printResult();
	}

	for (int i = 0; i < optData.arcConstraints.size(); ++i)
	{
#ifdef LOG_SPACE_MULTIPLIERS
		optData.arcConstraints[i].multiplier = std::exp(parameters[optData.halfarcs.size() + i]);
#else
		optData.arcConstraints[i].multiplier = parameters[optData.halfarcs.size() + 2 * i] / parameters[optData.halfarcs.size() + 2 * i + 1];
#endif
	}

	int deactivatedArcConstraints = 0;
	for (int i = 0; i < optData.arcConstraints.size(); ++i)
	{
		if (optData.arcConstraints[i].broken)
			++deactivatedArcConstraints;
		}
	std::cout << deactivatedArcConstraints << " deactivated arc constraints." << std::endl;

	if(verbose)
		printResult();
#endif
}