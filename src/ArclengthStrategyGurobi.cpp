#include "ArclengthStrategyGurobi.h"

#include <nsessentials/util/TimedBlock.h>
#ifdef WITH_GUROBI
#include <gurobi_c++.h>
#endif

#include <random>

ArclengthStrategyGurobi::ArclengthStrategyGurobi(const std::chrono::steady_clock::duration& discreteOptimizationTimeLimit)
	: discreteOptimizationTimeLimit(discreteOptimizationTimeLimit)
{ }

void ArclengthStrategyGurobi::CalculateParametricLengths(ParametrizationData& optData)
{
#ifdef WITH_GUROBI
	//Ratio threshold between actual halfarc length and target length. If the threshold is
	//exceeded, the arc is made a visible seam. An allowableLengthRatioError of t means that
	//the actual ratio is allowed in the range [1/t, t].
	const float allowableLengthRatioError = 2.0f;

	nse::util::TimedBlock b("Finding parametric sizes of patches using Gurobi ..");

	bool verbose = false;

	//We have one variable per halfarc, variable indices correspond to halfarc indices.
	int nVariables = optData.halfarcs.size();

	std::vector<GRBVar> variables(nVariables);
	try
	{
		//Set up the Gurobi environment
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);		

		model.set(GRB_IntParam_OutputFlag, verbose ? 1 : 0);
		model.set(GRB_DoubleParam_MIPGap, 0.0025);

		auto finishBy = std::chrono::steady_clock::now() + discreteOptimizationTimeLimit;

		//Construct the objective function
		GRBQuadExpr objective = 0.0;
		
		for (int i = 0; i < optData.halfarcs.size(); ++i)
		{
			//variable for the arc length
			variables[i] = model.addVar(1.0, 5000.0, 0.0, GRB_INTEGER);
			//number of faces included in the patch
			size_t patchSize = 0;
			if (optData.halfarcs[i].face != (size_t)-1)
				patchSize = optData.patches[optData.halfarcs[i].face].Faces().size();

			//add a quadratic fitting term to the energy
			objective = objective + 
				patchSize * optData.geometricArcFitWeight[i] * (variables[i] - optData.parametricHalfarcTargetLengths[i]) * (variables[i] - optData.parametricHalfarcTargetLengths[i]);
		}		
		model.setObjective(objective, GRB_MINIMIZE);

		//constraints that keep both sides of an arc equal
		std::vector<GRBConstr> gurobiArcConstraints(optData.arcConstraints.size());	

		//Determines if an arc constraint is active
		std::vector<bool> useConstraint(optData.arcConstraints.size());

		//pentalties for breaking an arc constraint
		std::vector<double> arcConstraintPenalties(optData.arcConstraints.size());
				
		//Set up arc constraints
		for (int i = 0; i < optData.arcConstraints.size(); ++i)
		{
			auto& c = optData.arcConstraints[i];
			if (c.broken)
				continue;			
			gurobiArcConstraints[i] = model.addConstr(c.multiplier * variables[c.arc] - variables[c.arc + 1] == 0, "Arc" + std::to_string(i));
			arcConstraintPenalties[i] = optData.halfarcs[2 * i].Length();
			useConstraint[i] = true;
		}
		//Set up face constraints that keep both sides of a face equal
		std::vector<GRBConstr> gurobiFaceConstraints(optData.faceConstraints.size());
		for (int i = 0; i < optData.faceConstraints.size(); ++i)
		{
			auto& c = optData.faceConstraints[i];
			GRBLinExpr expr;
			for (int k = 0; k < 2; ++k)
			{
				double sign = (k == 0 ? +1 : -1);
				for (auto arc : c.arcs[k])
					expr = expr + sign * variables[arc];
			}
			//if both sides can grow, the face does not impose any constraints
			if (c.canGrow[0] && c.canGrow[1])
				continue;
			char sense;
			if (!c.canGrow[0] && !c.canGrow[1])
				sense = GRB_EQUAL;
			else if (c.canGrow[0])
				sense = GRB_LESS_EQUAL;
			else
				sense = GRB_GREATER_EQUAL;
			gurobiFaceConstraints[i] = model.addConstr(expr, sense, 0.0);
		}

		bool ratiosAreGood = false;
		int iterations = 0;
		while (!ratiosAreGood && iterations++ < 20)
		{
			//solve the current model
			model.set(GRB_DoubleParam_TimeLimit, std::chrono::duration_cast<std::chrono::milliseconds>(finishBy - std::chrono::steady_clock::now()).count() / 1000.0);
			model.optimize();
			auto status = model.get(GRB_IntAttr_Status);
			if (status == GRB_TIME_LIMIT)
				break;
			if (status == GRB_INFEASIBLE)
			{
				//if the model is infeasible, solve a feasibility relaxation
				std::cout << "Feasibility relaxation" << std::endl;
				GRBModel modelCopy(model);
				std::vector<GRBConstr> breakableConstraints; breakableConstraints.reserve(gurobiArcConstraints.size());
				std::vector<double> breakableConstraintPenalties; breakableConstraintPenalties.reserve(gurobiArcConstraints.size());
				for (int i = 0; i < gurobiArcConstraints.size(); ++i)
				{
					if (!useConstraint[i])
						continue;
					breakableConstraints.push_back(gurobiArcConstraints[i]);
					breakableConstraintPenalties.push_back(arcConstraintPenalties[i]);
				}
				modelCopy.feasRelax(2, false, 0, nullptr, nullptr, nullptr, breakableConstraints.size(), breakableConstraints.data(), breakableConstraintPenalties.data());
				modelCopy.optimize();
				if (modelCopy.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
					break;

				//determine what arc constraints have been broken by the feasibility relaxation and
				//disable them
				auto feasibilityTolerance = model.get(GRB_DoubleParam_FeasibilityTol);
				for (int i = 0; i < breakableConstraints.size(); ++i)
				{
					auto& cname = breakableConstraints[i].get(GRB_StringAttr_ConstrName);
					auto ArtP = modelCopy.getVarByName("ArtP_" + cname).get(GRB_DoubleAttr_X);
					auto ArtN = modelCopy.getVarByName("ArtN_" + cname).get(GRB_DoubleAttr_X);
					if (ArtP > feasibilityTolerance || ArtN > feasibilityTolerance)
					{
						std::cout << "Broken arc " << i << std::endl;
						auto arcConstraintId = std::stoi(cname.substr(3));
						model.remove(breakableConstraints[i]);
						useConstraint[arcConstraintId] = false;
						optData.arcConstraints[arcConstraintId].broken = true;
					}
				}

				//solve again with the disabled arc constraints
				model.set(GRB_DoubleParam_TimeLimit, std::chrono::duration_cast<std::chrono::seconds>(finishBy - std::chrono::steady_clock::now()).count());
				model.optimize();
				if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
					break;
			}

			//check if actual / target length ratios are in the desired range
			float minRatio = 1000000;
			float maxRatio = 0;
			std::vector<std::pair<float, size_t>> failedMinRatioArcs;
			for (int i = 0; i < optData.halfarcs.size(); ++i)
			{
				auto ratio = (variables[i].get(GRB_DoubleAttr_X) / optData.parametricHalfarcTargetLengths[i]);
				if (ratio < minRatio)
				{
					minRatio = ratio;
					if (ratio < 1.0f / allowableLengthRatioError)
						failedMinRatioArcs.emplace_back(ratio, i);
				}
				if (ratio > maxRatio)
					maxRatio = ratio;
				//if (verbose)
				//	std::cout << "Length: " << variables[i].get(GRB_DoubleAttr_X) << " (optimal length: " << optData.halfarcs[i].parametricTargetLength << ", ratio: " << ratio << ", vertex: " << optData.mesh->from_vertex_handle(optData.graph->MotorcycleHalfedge(*optData.halfarcs[i].begin())) << ")" << std::endl;
			}
			std::cout << "Min ratio: " << minRatio << std::endl;
			std::cout << "Max ratio: " << maxRatio << std::endl;

			if (minRatio > 1.0f / allowableLengthRatioError)
			{
				ratiosAreGood = true;
				break;
			}

			//try to disable arc constraints in order to achieve better ratios
			std::sort(failedMinRatioArcs.begin(), failedMinRatioArcs.end());
			bool removedArc = false;
			for (auto& failed : failedMinRatioArcs)
			{
				if (removedArc)
					break;
				auto arcIdx = failed.second;
				auto fIdx = optData.halfarcs[arcIdx].face;
				if (fIdx == (size_t)-1)
					continue;
				auto& f = optData.patches[fIdx];
				int sideOfPatch = 0;
				for (int side = 0; side < 4; ++side)
				{
					for (auto arc : f.PatchSides()[side])
						if (arc == arcIdx)
							sideOfPatch = side;
				}
				for (auto arc : f.PatchSides()[sideOfPatch])
				{
					auto cIdx = arc / 2;
					if (!useConstraint[cIdx])
						std::cout << "Arc not constrained." << std::endl;
					else
					{
						std::cout << "Removing constraint " << cIdx << " for arc " << arc << std::endl;
						optData.arcConstraints[arc / 2].broken = true;
						model.remove(gurobiArcConstraints[cIdx]);
						removedArc = true;
						useConstraint[cIdx] = false;
					}
				}
			}
			if (!removedArc)
			{
				std::cout << "None of the arcs with bad length ratio are constrained." << std::endl;
				break;
			}
		}		

		//store the solution
		for (int i = 0; i < optData.halfarcs.size(); ++i)		
			optData.parametricHalfarcLengths[i] = variables[i].get(GRB_DoubleAttr_X);

		if (verbose)
		{
			int violatedArcConstraints = 0, violatedFaceConstraints = 0;
			for (int i = 0; i < gurobiArcConstraints.size(); ++i)
			{			
				if (!useConstraint[i])
					++violatedArcConstraints;
			}
		
			std::cout << "Violated " << violatedArcConstraints << " arc constraints." << std::endl;
			std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
		}
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;
	}
#else
	throw std::runtime_error("Gurobi is not available.");
#endif
}