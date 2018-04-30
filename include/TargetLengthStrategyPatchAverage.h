#pragma once

#include "common.h"
#include "MotorcycleGraph.h"
#include "Parametrization.h"
#include "ParametrizationInterfaces.h"

class TargetLengthStrategyPatchAverage : public ITargetlengthStrategy
{
public:
	TargetLengthStrategyPatchAverage(LengthMeasure lengthMeasure, const HEMesh& mesh, const MotorcycleGraph& graph, float multiplier = 1)
		: lengthMeasure(lengthMeasure), mesh(mesh), graph(graph), multiplier(multiplier)
	{ }

	void SetParametricTargetLength(ParametrizationData& optData)
	{
#pragma omp parallel for
		for (int i = 0; i < optData.halfarcs.size(); ++i)
			optData.parametricHalfarcTargetLengths[i] = std::numeric_limits<float>::quiet_NaN();

#pragma omp parallel for
		for (int i = 0; i < optData.patches.size(); ++i)
		{
			auto& patch = optData.patches[i];
			if (patch.PatchSides().size() != 4)
				continue;
			for (int side = 0; side < 4; ++side)
			{
				int weight;
				float average = patch.GetAverageInteriorPatchLength((side + 1) % 2, weight);
				if (weight == 0)
					continue;

				int totalSideLength = 0;
				for (auto arcIdx : patch.PatchSides()[side])
					totalSideLength += optData.halfarcs[arcIdx].Length();
				for (auto arcIdx : patch.PatchSides()[side])
					optData.parametricHalfarcTargetLengths[arcIdx] = multiplier * (average * optData.halfarcs[arcIdx].Length()) / totalSideLength;
			}
		}

#pragma omp parallel for
		for (int i = 0; i < optData.halfarcs.size(); ++i)
		{		
			if (std::isnan(optData.parametricHalfarcTargetLengths[i]))
			{
				optData.parametricHalfarcTargetLengths[i] = 0;
				for (auto p : optData.halfarcs[i])
				{
					switch (lengthMeasure)
					{
					case Topological:
						optData.parametricHalfarcTargetLengths[i] += 1;
						break;
					case Geometrical:
						optData.parametricHalfarcTargetLengths[i] += multiplier * mesh.calc_edge_length(graph.MotorcycleHalfedge(p));
						break;
					}
				}
			}
		}
	}

private:
	LengthMeasure lengthMeasure;
	const HEMesh& mesh;
	const MotorcycleGraph& graph;
	float multiplier;
};