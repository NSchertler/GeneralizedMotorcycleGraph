#pragma once

#include "common.h"
#include "MotorcycleGraph.h"
#include "TexturePatch.h"

#include "ParametrizationData.h"
#include "ParametrizationInterfaces.h"

//Performs all necessary computations to calculate a parametrization for a given set of patches.
class Parametrization
{
public:
	//Instantiates the parametrization calculation.
	//mesh - the underlying mesh
	//graph - the underlying motorcycle graph
	//patches - the set of extracted patches
	//texCoords - a list of texture coordinates storage objects for each patch in which to save the result	
	//parametrizationErrorThreshold - determines the upper threshold of the parametrization energy above which 
	//                                parametrization is re-attempted after relaxing the system. Details depend 
	//                                on the chosen parametrization strategy
	//targetlengthStrategy - the strategy to use to determine parametric target lengths
	//arclengthStrategy - the strategy to use to determine parametric arc lengths
	//parametrizationStrategy - the strategy to use to calculate the actual parametrization
	Parametrization(const HEMesh& mesh, const MotorcycleGraph& graph, const std::vector<TexturePatch>& patches, std::vector<TextureCoordinatesStorage>& texCoords,
		float parametrizationErrorThreshold, ITargetlengthStrategy* targetlengthStrategy, IArclengthStrategy* arclengthStrategy, IParameterizationStrategy* parameterizationStrategy);

	//Calculates the parametrization
	void run();
	
	//Returns a list of halfarc ids that are broken (i.e. visible seams)
	std::vector<size_t> BrokenHalfarcs() const;

private:		
	//Sets up optimization constraints for the parametric lengths
	void CalculateConstraints();

	ParametrizationData optData;
	
	ITargetlengthStrategy* targetlengthStrategy;
	IArclengthStrategy* arclengthStrategy;
	IParameterizationStrategy* parameterizationStrategy;

	float parametrizationErrorThreshold;
};