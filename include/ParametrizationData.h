#pragma once

#include <vector>

#include "common.h"
#include "MotorcycleGraph.h"

struct ParametrizationData;

//Constraint of the form  multiplier * arcLength - oppositeArcLength == 0
struct ArcConstraintInfo
{
	ParametrizationData* optData;

	//The arc id that generated this constraint
	int arc;

	//Determines if this arc constraint is broken, i.e. it is a visible seam
	bool broken = false;

	//Determines parametric size multiplier between the two sides of the arc
	double multiplier;

	ArcConstraintInfo() { }

	ArcConstraintInfo(ParametrizationData* optData, int arc)
		: optData(optData), arc(arc)
	{ }
};

//Constraint of the form  Sum of left/top arc lengths == Sum of right/bottom arc lengths
struct FaceConstraintInfo
{
	enum ConstraintDirection
	{
		//Uses the sides 0 and 2 of the face
		LeftRight,

		//Uses the sides 1 and 3 of the face
		TopBottom,
	};

	ParametrizationData* optData;

	//The two sets of arc indices that are relevant for this constraint
	std::vector<int> arcs[2];

	//Determines for the two face sides if the side can grow, i.e. if there is an open boundary on the side
	bool canGrow[2];

	FaceConstraintInfo() { }

	//Instantiates the face constraint
	//face - specify the index of the face/patch for this constraint
	//direction - specify the side pair for this constraint
	FaceConstraintInfo(ParametrizationData* optData, int face, ConstraintDirection direction);
};

//Holds all relevant data for parametrization
struct ParametrizationData
{
	ParametrizationData(const std::vector<MotorcycleGraph::HalfArc>& halfarcs, const std::vector<TexturePatch>& patches, 
		std::vector<TextureCoordinatesStorage>& texCoords, const HEMesh& mesh, const MotorcycleGraph& graph)
		: halfarcs(halfarcs), patches(patches), texCoords(texCoords), mesh(&mesh), graph(&graph), 
		parametricHalfarcLengths(halfarcs.size(), -1), parametricHalfarcTargetLengths(halfarcs.size(), -1)
	{ }

	const HEMesh* mesh;
	const MotorcycleGraph* graph;

	const std::vector<MotorcycleGraph::HalfArc>& halfarcs;
	const std::vector<TexturePatch>& patches;
	std::vector<TextureCoordinatesStorage>& texCoords;

	//A weight for the geometric fitting energy for every arc
	std::vector<double> geometricArcFitWeight;

	std::vector<ArcConstraintInfo> arcConstraints;
	std::vector<FaceConstraintInfo> faceConstraints;

	std::vector<float> parametricHalfarcLengths, parametricHalfarcTargetLengths;

	//Deprecated - weight for the length fitting term for arc lengths
	double weightLengthFit;

	//Deprecated - weight for the unit multiplier term
	double weightArcMultipliers;
};