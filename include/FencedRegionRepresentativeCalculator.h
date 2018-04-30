#pragma once

#include "common.h"
#include "MotorcycleGraph.h"
#include "FencedRegion.h"
#include "MotorcycleOptionsInPatch.h"

//Calculates the representative singular point for a fenced region along with the initial paths.
class FencedRegionRepresentativeCalculator
{
	//Run motorcycles from the patch boundary to the inside of the patch and find a good point where 
	//all motorcycles meet.

public:
	FencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp);

	//Calculates the representative for a fenced region and returns the corresponding cost.
	float CalculateRepresentative();

	//Fills a vector with the path of the motorcycle with the given orientation until it leaves the fenced region.
	//Returns if the resulting motorcycle will be active.
	bool FillWithSolution(int orientation, std::vector<HEMesh::HalfedgeHandle>& result) const;

protected:

	struct CycleAdditionalData
	{
		//The orientation through which the motorcycle entered/leaves the fenced region
		int orientation;

		//The first edge outside of the fenced region, oriented outwards
		HEMesh::HalfedgeHandle firstEdgeOutside;

		//the index of the motorcycle's entry point on the patch's boundary
		size_t entryPoint; 

		CycleAdditionalData() = default;
		CycleAdditionalData(int orientation, HEMesh::HalfedgeHandle firstEdgeOutside, size_t entryPoint)
			: orientation(orientation), firstEdgeOutside(firstEdgeOutside), entryPoint(entryPoint)
		{ }
	};

	//Record of a motorcycle passing through a vertex
	struct PassedThroughMotorcycle
	{
		//The index of the motorcycle
		size_t motorcycleIdx;

		//The length of the motorcycle path at the time of pass-through
		size_t lengthofMotorcyclePath;

		//The cost of the path up to the pass-through
		float cost;

		PassedThroughMotorcycle() { }

		PassedThroughMotorcycle(size_t motorcycleIdx, size_t lengthOfMotorcyclePath, float cost)
			: motorcycleIdx(motorcycleIdx), lengthofMotorcyclePath(lengthOfMotorcyclePath), cost(cost)
		{ }

		bool operator<(const PassedThroughMotorcycle& other) const { return cost < other.cost; }
	};

	//Information stored for a vertex in the fenced region
	struct VertexInfo
	{
		//a vector with an entry for each possible color
		//entry: a list of motorcycles that passed through this vertex sorted by their cost
		std::vector<std::set<PassedThroughMotorcycle>> possiblePaths;

		//Number of edges of the shortest path from the region's fence to the vertex
		int distanceFromBoundary;
	};

	const HEMesh& mesh;
	const FencedRegion& patch;
	const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp;

	//degree of the fenced region
	size_t distinctOrientations;

	//simulated motorcycle within the fenced region
	MotorcycleOptionsInPatch<CycleAdditionalData> cycles;

	//key: vertex that lies within the patch
	std::map<HEMesh::VertexHandle, VertexInfo> vertexInfo;	

	//For every orientation, references the motorcycle that had the best emanating path.
	std::vector<PassedThroughMotorcycle> bestPaths;

	//Initializes motorcycles on the boundary of the fenced region and calculates vertex distances
	//from the boundary.
	void InitializeMotorcyclesAndDistances();

	//Advances the motorcycles until they leave the fenced region.
	void FindAllPossibleMotorcyclePaths();

	//Calculates the best set of emanating motorcycles for the entire fenced region. Returns
	//the respective cost.
	float CalculateBestEmanatingPaths();

	//Calculates the centeredness energy for a given distance of a vertex from the region fence.
	float GetCenterednessEnergy(int distanceFromBoundary) const;

	//Calculates the perpendicularity energy between two motorcycles
	float GetAngleEnergy(HEMesh::HalfedgeHandle h1, HEMesh::HalfedgeHandle h2) const;

	//Calculates the the set of emanating motorcycles for a given vertex by updating bestPaths.
	//Returns the respective cost.
	virtual float CalculateBestEmanatingPathsForVertex(HEMesh::VertexHandle v, const VertexInfo&, float bestCostSoFar) = 0;

	//Returns the first halfedge on the path of an emanating motorcycle (starting at the representative, going outwards)
	HEMesh::HalfedgeHandle GetFirstPathEdge(const PassedThroughMotorcycle& motorcycle) const;
};

//Simple brute force calculation of the representative. May yield invalid results.
class SimpleFencedRegionRepresentativeCalculator : public FencedRegionRepresentativeCalculator
{
public:
	SimpleFencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp);
protected:
	float CalculateBestEmanatingPathsForVertex(HEMesh::VertexHandle v, const VertexInfo&, float bestCostSoFar);
};

//Dynamic Programming calculation of the representative as presented in the paper.
class DPFencedRegionRepresentativeCalculator : public FencedRegionRepresentativeCalculator
{
public:
	DPFencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp);
protected:
	float CalculateBestEmanatingPathsForVertex(HEMesh::VertexHandle v, const VertexInfo&, float bestCostSoFar);
};
