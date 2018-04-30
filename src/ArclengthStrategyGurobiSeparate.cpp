#include "ArclengthStrategyGurobiSeparate.h"

#include <nsessentials/util/TimedBlock.h>
#ifdef WITH_GUROBI
#include <gurobi_c++.h>
#endif

void ArclengthStrategyGurobiSeparate::CalculateParametricLengths(ParametrizationData& optData)
{
#ifdef WITH_GUROBI
	nse::util::TimedBlock b("Finding parametric sizes of patches using separate Gurobi ..");

	bool verbose = false;

#pragma omp parallel for
	for (int i = 0; i < optData.arcConstraints.size(); ++i)
		optData.arcConstraints[i].broken = true;	

	try
	{
		GRBEnv env = GRBEnv();

		for (int iPatch = 0; iPatch < optData.patches.size(); ++iPatch)
		{
			auto& patch = optData.patches[iPatch];
			if (patch.PatchSides().size() != 4)
				continue;

			GRBModel model = GRBModel(env);

			model.set(GRB_IntParam_OutputFlag, verbose ? 1 : 0);

			GRBQuadExpr objective = 0.0;

			std::map<size_t, GRBVar> variables;
			for (auto& side : patch.PatchSides())
			{
				for (auto arcIdx : side)
				{
					auto var = model.addVar(1.0, 5000.0, 0, GRB_INTEGER);
					objective = objective + optData.geometricArcFitWeight[arcIdx] * (var - optData.parametricHalfarcTargetLengths[arcIdx]) * (var - optData.parametricHalfarcTargetLengths[arcIdx]);
					variables[arcIdx] = var;
				}
			}

			for (int k = 0; k < 2; ++k)
			{
				GRBLinExpr expr;
				auto& c = optData.faceConstraints[2 * iPatch + k];
				for (int i = 0; i < 2; ++i)
				{
					double sign = (i == 0 ? +1 : -1);
					for (auto arc : c.arcs[i])
						expr = expr + sign * variables.at(arc);
				}
				if (!(c.canGrow[0] && c.canGrow[1]))
				{

					char sense;
					if (!c.canGrow[0] && !c.canGrow[1])
						sense = GRB_EQUAL;
					else if (c.canGrow[0])
						sense = GRB_LESS_EQUAL;
					else
						sense = GRB_GREATER_EQUAL;
					model.addConstr(expr, sense, 0.0);
				}
			}
						
			model.setObjective(objective, GRB_MINIMIZE);

			model.optimize();

			for(auto& entry : variables)
			{
				if (verbose)
					std::cout << "Length: " << entry.second.get(GRB_DoubleAttr_X) << " (optimal length: " << optData.parametricHalfarcTargetLengths[entry.first] << ")" << std::endl;
				optData.parametricHalfarcLengths[entry.first] = entry.second.get(GRB_DoubleAttr_X);
			}
			
			if (verbose)
				std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
		}
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << std::endl;
	}
#endif
}