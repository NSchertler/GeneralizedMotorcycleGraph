#pragma once

#include <vector>
#include <queue>

#include "common.h"
#include "FencedRegion.h"

template <typename TAdditionalDataForMotorcycles>
class MotorcycleOptionsInPatch;

//Represents a motorcycle during multidirectional Dijkstra
//TAdditionalData - The data type with which to augment the internal motorcycles.
template <typename TAdditionalData>
class PatchMotorcycle
{
public:
	//Initializes the motorcycle starting at the given edge.
	PatchMotorcycle(HEMesh::HalfedgeHandle edge)
		: costSoFar(0)
	{
		edges.push_back(edge);
	}

	TAdditionalData additionalData;
	
	//Returns the cost accumulated by this motorcycle so far.
	float CostSoFar() const { return costSoFar; }

	//Returns the edges traversed by this motorcycle so far.
	const std::vector<HEMesh::HalfedgeHandle>& Edges() const { return edges; }

private:
	//All edges traversed by this motorcycle so far.
	std::vector<HEMesh::HalfedgeHandle> edges;	

	//The cost accumulated by this motorcycle so far.
	float costSoFar;

	//Determines if the motorcycle has traversed any non-original edge.
	bool hasMovedAlongNonOriginalEdge = false;
	
	template <typename TAdditionalDataForMotorcycles>
	friend class MotorcycleOptionsInPatch;
};

namespace MotorcycleOptionsHelper
{
	static bool VoidVertexAction(size_t, HEMesh::VertexHandle) { return true; }
	static void VoidHalfedgeAction(size_t, HEMesh::HalfedgeHandle) { }
	static float DefaultMaxCostQuery() { return std::numeric_limits<float>::infinity(); };

	struct VoidStruct {};
}

//Performs multidirectional Dijkstra on the conjugate graph of a fenced region by virtually tracing motorcycles.
//Motorcycles that travel straight are prioritized.
//TAdditionalDataForMotorcycles - The data type with which to augment the internal motorcycles.
template <typename TAdditionalDataForMotorcycles = MotorcycleOptionsHelper::VoidStruct>
class MotorcycleOptionsInPatch
{	
	//The additional cost that is added to a motorcycle if it uses a non-original edge
	const float NonOriginalCost = 100;

public:
	MotorcycleOptionsInPatch()
		: activeMotorcycles(PriorityComparer(motorcycles))
	{ }

	//Adds a new motorcycle to the initialization.
	//Returns the id of the new motorcycle.
	//edge - The location of the motorcycle
	//mesh - The underlying mesh
	//isOriginalEdgeProp - mesh property
	//additionalDataArgs - arguments that are used for calling the constructor of the additional data of motorcycles
	template <typename ... TArgs>
	size_t AddMotorcycle(HEMesh::HalfedgeHandle edge, const HEMesh& mesh, const OpenMesh::EPropHandleT<bool> isOriginalEdgeProp, TArgs... additionalDataArgs)
	{
		//Add a new motorcycle
		motorcycles.emplace_back(edge);

		//Initialize the additional data with the provided arguments
		motorcycles.back().additionalData = TAdditionalDataForMotorcycles(additionalDataArgs...);

		//Check if the initial edge is non-original
		if (!mesh.property(isOriginalEdgeProp, mesh.edge_handle(edge)))
		{
			motorcycles.back().costSoFar += NonOriginalCost;
			motorcycles.back().hasMovedAlongNonOriginalEdge = true;
		}

		//Add the new motorcycle to the active set.
		activeMotorcycles.push(motorcycles.size() - 1);
		return motorcycles.size() - 1;
	}

	//Returns a list of all (active and inactive) motorcycles.
	const std::vector<PatchMotorcycle<TAdditionalDataForMotorcycles>>& Motorcycles() const { return motorcycles; }

	/*
	 * Runs Dijkstra based on the current initialization.
	 * CyclePassThroughVertexAction: bool operator()(size_t cycleIdx, HEMesh::VertexHandle vertex)
	 * CycleLeavePatchAction: void operator()(size_t cycleIdx, HEMesh::HalfedgeHandle firstEdgeOutsidePatch)
	 * MaxCostQuery: float operator()()
	 * mesh - the underlying mesh
	 * isOriginalEdgeProp - mesh property
	 * patch - the fenced region for which to perform tracing
	 * cyclePassThroughVertexAction - called each time a motorcycle passes through a vertex. Must return if the motorcycle is allowed to continue.
	 * cycleLeavePatchAction - called each time a motorcycle leaves the fenced region.
	 * maxCostQuery - tracing is stopped as soon as all motorcycles have a cost that is greater than maxCostQuery().
	 * minCosDeviationAngle - threshold that determines how much the motorcycles are allowed to deviate from straight continuations.
	 * allowTurnsOnRegularVertices - determines if motorcycles are allowed to deviate from the topologically straight path at irregular vertices.
	 */
	template <typename CyclePassThroughVertexAction, typename CycleLeavePatchAction, typename MaxCostQuery>
	void SimulateCyclesInPatch(const HEMesh& mesh, const OpenMesh::EPropHandleT<bool> isOriginalEdgeProp, const FencedRegion& patch, CyclePassThroughVertexAction&& cyclePassThroughVertexAction, CycleLeavePatchAction&& cycleLeavePatchAction, MaxCostQuery&& maxCostQuery, const float minCosDeviationAngle = 0.5f, bool allowTurnsOnRegularVertices = false)
	{
		//Holds possible continuations for the currently examined motorcycle
		std::vector<PathContinuation> continuations;

		//run the motorcycles as long as their are active ones
		while (!activeMotorcycles.empty())
		{
			//Take the motorcycle with the least cost so far
			auto cycleIdx = activeMotorcycles.top();
			activeMotorcycles.pop();

			//current location of the current motorcycle
			auto motorcycleHalfedge = motorcycles[cycleIdx].edges.back();
			auto target = mesh.to_vertex_handle(motorcycleHalfedge);

			//if the accumulated cost is too much, terminate the process
			if (motorcycles[cycleIdx].costSoFar > std::forward<MaxCostQuery>(maxCostQuery)())
				return;

			//check the callback if the motorcycle is allowed to continue from its current location
			bool canContinue = std::forward<CyclePassThroughVertexAction>(cyclePassThroughVertexAction)(cycleIdx, target);
			if (!canContinue)
				continue;

			//If the motorcycle passes through a vertex on a boundary, it must leave the region for consistency reasons.
			bool mustExitPatch = patch.IsToVertexOnBoundary(motorcycleHalfedge);

			//Let the motorcycle continue
			float costBeforeContinuation = motorcycles[cycleIdx].costSoFar;			
			FindPossibleContinuations(mesh, motorcycleHalfedge, continuations, minCosDeviationAngle, allowTurnsOnRegularVertices);
			bool cycleContinued = false;
			for (auto& continuation : continuations)
			{
				if (!continuation.edge.is_valid() && !patch.IsToVertexOnBoundary(motorcycleHalfedge, false))
					continue; //inner hole in patch; TODO: continue on the other side
				
				if (continuation.edge.is_valid())
				{
					//check if the motor cycle forms a cyclic path
					bool isCycle = false;
					auto continuationToVertex = mesh.to_vertex_handle(continuation.edge);
					if (mesh.from_vertex_handle(motorcycles[cycleIdx].edges[0]) == continuationToVertex)
						isCycle = true;
					else
						for (auto e : motorcycles[cycleIdx].edges)
							if (mesh.to_vertex_handle(e) == continuationToVertex)
							{
								isCycle = true;
								break;
							}
					if (isCycle)
						continue;
				}

				bool edgeIsOutsideOfPatch = true;
				if (continuation.edge.is_valid())
				{
					//check if the edge is outside of the patch
					auto face1 = mesh.face_handle(continuation.edge);
					auto face2 = mesh.opposite_face_handle(continuation.edge);
					edgeIsOutsideOfPatch = (face1.is_valid() && !patch.IsFaceIncluded(face1)) || (face2.is_valid() && !patch.IsFaceIncluded(face2));
				}
				if (edgeIsOutsideOfPatch)
				{
					if (continuation.edge.is_valid())
					{
						//check if motorcycle made a turn on the boundary
						if (mesh.face_handle(motorcycleHalfedge).is_valid() && continuation.edge == mesh.next_halfedge_handle(motorcycleHalfedge))
							continue; // this is an invalid turn
						if (mesh.opposite_face_handle(motorcycleHalfedge).is_valid() && continuation.edge == mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(motorcycleHalfedge))))
							continue; //this is an invalid turn
					}

					//record this cycle
					motorcycles.push_back(motorcycles[cycleIdx]);
					motorcycles.back().costSoFar = costBeforeContinuation;
					if (cycleContinued) //we already re-used this cycle - delete the last path segment
						motorcycles.back().edges.erase(motorcycles.back().edges.end() - 1);

					std::forward<CycleLeavePatchAction>(cycleLeavePatchAction)(motorcycles.size() - 1, continuation.edge);
					continue;
				}
			
				//at this point, we have a continuation that does not leave the patch

				if (mustExitPatch)
					continue;

				float continuationCost = std::acos(OpenMesh::sane_aarg(continuation.score));
				if (!mesh.property(isOriginalEdgeProp, mesh.edge_handle(continuation.edge)) && !motorcycles[cycleIdx].hasMovedAlongNonOriginalEdge)
					continuationCost += NonOriginalCost;				

				if (!cycleContinued)
				{
					//Re-use the existing motorcycle and just append a new segment to the path
					motorcycles[cycleIdx].edges.push_back(continuation.edge);
					motorcycles[cycleIdx].costSoFar += continuationCost;
					activeMotorcycles.push(cycleIdx);
					cycleContinued = true;
				}
				else
				{
					//create a new motorcycle
					motorcycles.push_back(motorcycles[cycleIdx]); //copy the current motorcycle
					motorcycles.back().edges.back() = continuation.edge;
					motorcycles.back().costSoFar = costBeforeContinuation + continuationCost;
					activeMotorcycles.push(motorcycles.size() - 1);
				}
			}
		}
	}

private:
	std::vector<PatchMotorcycle<TAdditionalDataForMotorcycles>> motorcycles;

	class PriorityComparer
	{
	public:
		PriorityComparer(std::vector<PatchMotorcycle<TAdditionalDataForMotorcycles>>& motorcycles)
			: motorcycles(motorcycles) { }

		bool operator()(size_t first, size_t second) const { return motorcycles[first].costSoFar > motorcycles[second].costSoFar; }
	private:
		std::vector<PatchMotorcycle<TAdditionalDataForMotorcycles>>& motorcycles;
	};

	std::priority_queue<size_t, std::vector<size_t>, PriorityComparer> activeMotorcycles;
};