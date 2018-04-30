#include "Parametrization.h"

#include <nsessentials/util/TimedBlock.h>

Parametrization::Parametrization(const HEMesh& mesh, const MotorcycleGraph& graph, const std::vector<TexturePatch>& patches, 
	std::vector<TextureCoordinatesStorage>& texCoords, float parametrizationErrorThreshold, ITargetlengthStrategy* targetlengthStrategy, 
	IArclengthStrategy* arclengthStrategy, IParameterizationStrategy* parameterizationStrategy)
	: optData(graph.Halfarcs(), patches, texCoords, mesh, graph), parametrizationErrorThreshold(parametrizationErrorThreshold),
	targetlengthStrategy(targetlengthStrategy), arclengthStrategy(arclengthStrategy), parameterizationStrategy(parameterizationStrategy)
{ }

void Parametrization::run()
{
	nse::util::TimedBlock b("Calculating parametrization ..");

	CalculateConstraints();
	targetlengthStrategy->SetParametricTargetLength(optData);

	optData.geometricArcFitWeight.resize(optData.halfarcs.size());
#pragma omp parallel for
	for (int i = 0; i < optData.halfarcs.size(); ++i)
	{
		//Calculate fitting weight such that the energy of a residual of 50% of the target length equals 1
		double residualForEnergy1 = 0.5 * optData.parametricHalfarcTargetLengths[i];
		optData.geometricArcFitWeight[i] = 1.0 / (residualForEnergy1 * residualForEnergy1);
	}

	arclengthStrategy->CalculateParametricLengths(optData);
	parameterizationStrategy->CalculateParameterization(optData, parametrizationErrorThreshold);
}

void Parametrization::CalculateConstraints()
{
	if (optData.patches.size() == 0)
		return;

	nse::util::TimedBlock b("Calculating parametrization constraints ..");

	optData.faceConstraints.resize(2 * optData.patches.size());
#pragma omp parallel for
	for (int i = 0; i < optData.patches.size(); ++i)
	{
		optData.faceConstraints[2 * i + 0] = FaceConstraintInfo(&optData, i, FaceConstraintInfo::LeftRight);
		optData.faceConstraints[2 * i + 1] = FaceConstraintInfo(&optData, i, FaceConstraintInfo::TopBottom);
	}

	optData.arcConstraints.resize(optData.halfarcs.size() / 2);
#pragma omp parallel for
	for (int i = 0; i < optData.halfarcs.size(); i += 2)
	{
		optData.arcConstraints[i / 2] = ArcConstraintInfo(&optData, i);
		optData.arcConstraints[i / 2].multiplier = optData.halfarcs[i].multiplier /*/ optData.halfarcs[i + 1].multiplier*/;
		//std::cout << optData.arcConstraints[i / 2].multiplier << std::endl;
	}
}


std::vector<size_t> Parametrization::BrokenHalfarcs() const
{
	std::vector<size_t> broken;
	for (auto& c : optData.arcConstraints)
		if (c.broken)
			broken.push_back(c.arc);
	return broken;
}