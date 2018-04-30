#include "FencedRegionRepresentativeCalculator.h"

FencedRegionRepresentativeCalculator::FencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp)
	: mesh(mesh), patch(patch), isOriginalEdgeProp(isOriginalEdgeProp)
{
	auto& loop = patch.Loops()[0];
	distinctOrientations = loop.Degree();
	bestPaths.resize(distinctOrientations);
}

float FencedRegionRepresentativeCalculator::CalculateRepresentative()
{
	InitializeMotorcyclesAndDistances();

	FindAllPossibleMotorcyclePaths();

	return CalculateBestEmanatingPaths();
}

//Returns if the resulting motorcycle will be active
bool FencedRegionRepresentativeCalculator::FillWithSolution(int orientation, std::vector<HEMesh::HalfedgeHandle>& result) const
{
	auto& pathDescriptor = bestPaths[orientation];
	result.clear();
	for (int i = pathDescriptor.lengthofMotorcyclePath - 1; i >= 0; --i)
		result.push_back(mesh.opposite_halfedge_handle(cycles.Motorcycles().at(pathDescriptor.motorcycleIdx).Edges().at(i)));
	auto edgeOutside = cycles.Motorcycles().at(pathDescriptor.motorcycleIdx).additionalData.firstEdgeOutside;
	if (edgeOutside.is_valid())
	{
		result.push_back(edgeOutside);
		return true;
	}
	else
		return false; //this motorcycle terminates at the boundary
}

void FencedRegionRepresentativeCalculator::InitializeMotorcyclesAndDistances()
{
	std::queue<std::pair<HEMesh::VertexHandle, int>> distanceBFS;

	//Initialize motorcycles from the boundary
	auto& loop = patch.Loops()[0];
	size_t nextEntryPoint = 0;
	for (int i = 0; i < loop.Edges().size(); ++i)
	{
		auto perpendicularEdge = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(loop.Edges()[i].edge));
		auto v = mesh.from_vertex_handle(perpendicularEdge);
		auto& sourceInfo = vertexInfo[v];
		sourceInfo.distanceFromBoundary = 0;
		sourceInfo.possiblePaths.resize(distinctOrientations);
		auto toVertex = mesh.to_vertex_handle(perpendicularEdge);
		if (!patch.IsToVertexOnBoundary(perpendicularEdge, false))
		{
			auto inserted = vertexInfo.emplace(toVertex, VertexInfo());
			if (inserted.second)
			{
				inserted.first->second.distanceFromBoundary = 1;
				inserted.first->second.possiblePaths.resize(distinctOrientations);
				distanceBFS.emplace(toVertex, 1);
			}
		}

		HEMesh::HalfedgeHandle continuationOutside(-1);
		if (!mesh.is_boundary(mesh.opposite_halfedge_handle(loop.Edges()[i].edge)))
		{
			RegularContinuation(mesh, mesh.opposite_halfedge_handle(perpendicularEdge), continuationOutside);
		}
		auto cycleIdx = cycles.AddMotorcycle(perpendicularEdge, mesh, isOriginalEdgeProp, loop.Edges()[i].orientation, continuationOutside, nextEntryPoint++);
		if(mesh.is_boundary(mesh.opposite_halfedge_handle(loop.Edges()[i].edge)) && mesh.is_boundary(mesh.opposite_halfedge_handle(loop.Edges()[(i - 1 + loop.Edges().size()) % loop.Edges().size()].edge)))
			sourceInfo.possiblePaths[loop.Edges()[i].orientation].emplace(cycleIdx, 0, 0);

		auto nextEdge = (i + 1) % loop.Edges().size();
		if (loop.Edges()[i].orientation != loop.Edges()[nextEdge].orientation)
		{
			perpendicularEdge = mesh.next_halfedge_handle(loop.Edges()[i].edge);
			v = mesh.from_vertex_handle(perpendicularEdge);
			toVertex = mesh.to_vertex_handle(perpendicularEdge);
			if (!patch.IsToVertexOnBoundary(perpendicularEdge, false))
			{
				auto inserted = vertexInfo.emplace(toVertex, VertexInfo());
				if (inserted.second)
				{
					inserted.first->second.distanceFromBoundary = 1;
					inserted.first->second.possiblePaths.resize(distinctOrientations);
					distanceBFS.emplace(toVertex, 1);
				}
			}
			continuationOutside = HEMesh::HalfedgeHandle(-1);
			if (!mesh.is_boundary(loop.Edges()[i].edge))
			{
				RegularContinuation(mesh, mesh.opposite_halfedge_handle(perpendicularEdge), continuationOutside);
			}
			cycleIdx = cycles.AddMotorcycle(perpendicularEdge, mesh, isOriginalEdgeProp, loop.Edges()[i].orientation, continuationOutside, nextEntryPoint++);
			//sourceInfo.possiblePaths[loop.Edges()[i].color].emplace(cycleIdx, 0, 0);
		}
	}

	//Calculate vertex distances
	while (!distanceBFS.empty())
	{
		auto entry = distanceBFS.front();
		distanceBFS.pop();
		for (auto v : mesh.vv_range(entry.first))
		{
			auto inserted = vertexInfo.emplace(v, VertexInfo());
			if (inserted.second)
			{
				inserted.first->second.distanceFromBoundary = entry.second + 1;
				inserted.first->second.possiblePaths.resize(distinctOrientations);
				distanceBFS.emplace(v, entry.second + 1);
			}
		}
	}
}

void FencedRegionRepresentativeCalculator::FindAllPossibleMotorcyclePaths()
{
	auto cyclePassThroughVertexAction = [&](size_t cycleIdx, HEMesh::VertexHandle vertex)
	{
		if (patch.IsToVertexOnBoundary(cycles.Motorcycles()[cycleIdx].Edges().back()))
			return true; //don't allow representative on the fence

		//Remember that this motorcycle passed through the current vertex
		auto& info = vertexInfo.at(vertex);
		//if there are too many options, consider an incremental approach like the following:
		//Make the dynamic programming problem incremental (i.e. allow to add nodes). Then simulate the motorcycles
		//until we can find a solution (recalculating the solution incrementally). After that, simulate the motorcycles
		//to the end as long as their cost does not exceed the cost that could make any vertex' cost smaller than the
		//initially calculated cost.
		//std::cout << "Recording pass-through at vertex " << vertex.idx() << " (color: " << cycles.Motorcycles()[cycleIdx].additionalData.color << ", colors: " << distinctOrientations << ")" <<  std::endl;

		info.possiblePaths[cycles.Motorcycles()[cycleIdx].additionalData.orientation].emplace(cycleIdx, cycles.Motorcycles()[cycleIdx].Edges().size(), cycles.Motorcycles()[cycleIdx].CostSoFar());
		return true;
	};

	cycles.SimulateCyclesInPatch(mesh, isOriginalEdgeProp, patch, cyclePassThroughVertexAction, &MotorcycleOptionsHelper::VoidHalfedgeAction, &MotorcycleOptionsHelper::DefaultMaxCostQuery);
}

float FencedRegionRepresentativeCalculator::GetCenterednessEnergy(int distanceFromBoundary) const
{
	const float centerednessFactor = 0.2f;
	return -centerednessFactor * distanceFromBoundary * distinctOrientations;
}

float FencedRegionRepresentativeCalculator::GetAngleEnergy(HEMesh::HalfedgeHandle h1, HEMesh::HalfedgeHandle h2) const
{
	if (!h1.is_valid() || !h2.is_valid())
		return 0;
	auto dir1 = mesh.calc_edge_vector(h1); dir1.normalize();
	auto dir2 = mesh.calc_edge_vector(h2); dir2.normalize();
	auto dot = OpenMesh::dot(dir1, dir2);
	//this does not consider the angle's orientation but that should be fine for mostly regular meshes
	auto angleDeviation = std::acos(OpenMesh::sane_aarg(dot)) - M_PI / 2;
	return angleDeviation * angleDeviation;
}

float FencedRegionRepresentativeCalculator::CalculateBestEmanatingPaths()
{
	float bestCost = std::numeric_limits<float>::infinity();
	for (auto& ventry : vertexInfo)
	{
		bestCost = CalculateBestEmanatingPathsForVertex(ventry.first, ventry.second, bestCost);
	}
	return bestCost;
}

HEMesh::HalfedgeHandle FencedRegionRepresentativeCalculator::GetFirstPathEdge(const PassedThroughMotorcycle& motorcycle) const
{
	if (motorcycle.lengthofMotorcyclePath > 0)
	{
		auto h = cycles.Motorcycles()[motorcycle.motorcycleIdx].Edges()[motorcycle.lengthofMotorcyclePath - 1];
		if (h.is_valid())
			h = mesh.opposite_halfedge_handle(h);
		return h;
	}
	else
		return cycles.Motorcycles()[motorcycle.motorcycleIdx].additionalData.firstEdgeOutside;
};

SimpleFencedRegionRepresentativeCalculator::SimpleFencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp)
	: FencedRegionRepresentativeCalculator(mesh, patch, isOriginalEdgeProp)
{ }


float SimpleFencedRegionRepresentativeCalculator::CalculateBestEmanatingPathsForVertex(HEMesh::VertexHandle vertex, const VertexInfo& info, float bestCostSoFar)
{	
	//find the best meeting point	
	auto& perColorPaths = info.possiblePaths;

	if (mesh.valence(vertex) + (mesh.is_boundary(vertex) ? 1 : 0) < distinctOrientations)
		bestCostSoFar;

	//check if we have paths for all of the colors
	bool hasPathForAllColors = true;
	for (int orientation = 0; orientation < distinctOrientations; ++orientation)
	{
		if (perColorPaths[orientation].size() == 0)
		{
			hasPathForAllColors = false;
			break;
		}
	}

	if (!hasPathForAllColors)
		return bestCostSoFar;

	//the path that is taken
	std::vector<std::set<PassedThroughMotorcycle>::iterator> pathPerColor(distinctOrientations);
	//std::cout << "Possible paths: ";
	int possiblePathsCount = 1;
	static int maxPathsCount = 0;
	for (int orientation = 0; orientation < distinctOrientations; ++orientation)
	{
		//if (color != 0)
		//	std::cout << " * ";
		//std::cout << perColorPaths[color].size();
		//possiblePathsCount *= perColorPaths[color].size();
		pathPerColor[orientation] = perColorPaths[orientation].begin();
	}
	//maxPathsCount = std::max(possiblePathsCount, maxPathsCount);
	//std::cout << " = " << possiblePathsCount << " (max seen so far: " << maxPathsCount << ")" << std::endl;		
	//TODO: investigate other paths

	//check if the paths do not cross
	std::set<HEMesh::VertexHandle> crossedVertices;
	bool pathIsValid = true;
	for (int orientation = 0; orientation < distinctOrientations; ++orientation)
	{
		auto& pathDescriptor = *pathPerColor[orientation];
		for (int i = 0; i < pathDescriptor.lengthofMotorcyclePath; ++i)
		{
			auto vertex = mesh.from_vertex_handle(cycles.Motorcycles()[pathDescriptor.motorcycleIdx].Edges()[i]);
			auto inserted = crossedVertices.insert(vertex);
			if (!inserted.second)
			{
				//The paths cross, do not consider this configuration
				pathIsValid = false;
				break;
			}
		}
		if (!pathIsValid)
			break;
	}

	if (!pathIsValid)
		return bestCostSoFar;

	float cost = GetCenterednessEnergy(info.distanceFromBoundary);
	int minPathLength = std::numeric_limits<int>::max();
	for (int orientation = 0; orientation < distinctOrientations; ++orientation)
	{
		//add cost for path
		cost += pathPerColor[orientation]->cost;

		if (pathPerColor[orientation]->lengthofMotorcyclePath < minPathLength)
			minPathLength = pathPerColor[orientation]->lengthofMotorcyclePath;

		//add cost for perpendicularity at meeting vertex
		int next = (orientation + 1) % distinctOrientations;
		auto e1 = GetFirstPathEdge(*pathPerColor[orientation]);
		auto e2 = GetFirstPathEdge(*pathPerColor[next]);

		cost += GetAngleEnergy(e1, e2);
	}

	if (cost < bestCostSoFar)
	{
		for (int orientation = 0; orientation < distinctOrientations; ++orientation)
			bestPaths[orientation] = *pathPerColor[orientation];
		return cost;
	}

	return bestCostSoFar;
}

DPFencedRegionRepresentativeCalculator::DPFencedRegionRepresentativeCalculator(const HEMesh& mesh, const FencedRegion& patch, const OpenMesh::EPropHandleT<bool>& isOriginalEdgeProp)
	: FencedRegionRepresentativeCalculator(mesh, patch, isOriginalEdgeProp)
{ }

float DPFencedRegionRepresentativeCalculator::CalculateBestEmanatingPathsForVertex(HEMesh::VertexHandle v, const VertexInfo& info, float bestCostSoFar)
{
	//For each color, a list of motorcycles (.second) that arrived at this vertex, ordered by their entry point (.first)
	std::vector<std::vector<std::pair<size_t, PassedThroughMotorcycle>>> motorcyclesPerColor(distinctOrientations);

	int maxEntryPoint = 0;
	for (int orientation = 0; orientation < distinctOrientations; ++orientation)
	{
		for (auto& path : info.possiblePaths[orientation])
		{
			if (motorcyclesPerColor[orientation].size() > 20)
				break; //only use the 20 best motorcycles per color to improve performance.
			auto entryPoint = cycles.Motorcycles()[path.motorcycleIdx].additionalData.entryPoint;
			if (entryPoint > maxEntryPoint)
				maxEntryPoint = entryPoint;
			motorcyclesPerColor[orientation].push_back(std::make_pair(entryPoint, path));
		}
		if (motorcyclesPerColor[orientation].size() == 0)
			return bestCostSoFar;
		std::sort(motorcyclesPerColor[orientation].begin(), motorcyclesPerColor[orientation].end());
	}	

	float centerednessEnergy = GetCenterednessEnergy(info.distanceFromBoundary);

	//Incrementally calculate the function f(color, lastMotorcycle, colorZeroMotorcycle) that is
	//the least possible cost to select motorcycles for all edge colors up to color, where the last
	//used motorcycle is lastMotorcycle, and the first motorcycle is colorZeroMotorcycle.
	//The best choice is then the argmin_(lastMotorcycle, colorZeroMotorcycle) f(lastColor, lastMotorcycle, colorZeroMotorcycle)

	struct DPEntry
	{
		float cost = std::numeric_limits<float>::infinity();
		size_t predecessorNode;
	};

	//dpTable[i, j] = f(color - 1, nodeIndexOfLastMotorcycle, currentColorZeroMotorcycle)
	//(no entries for first color as this is fixed by colorZeroMotorcycle)
	std::vector<std::vector<DPEntry>> dpTable(distinctOrientations - 1); 

	for (int colorZeroMotorcycle = 0; colorZeroMotorcycle < motorcyclesPerColor[0].size(); ++colorZeroMotorcycle)
	{
		auto& firstChoice = motorcyclesPerColor[0][colorZeroMotorcycle];
		auto firstChoiceEntryPoint = firstChoice.first;
		auto& firstMotorcycleDesc = firstChoice.second;
		auto& firstMotorcycle = cycles.Motorcycles()[firstMotorcycleDesc.motorcycleIdx];
		float firstCost = firstMotorcycleDesc.cost;
		auto firstMotorcycleEdge = GetFirstPathEdge(firstMotorcycleDesc);

		//resort the other motorcycles and ensure that all entrypoints are larger than the first
		for (int orientation = 1; orientation < distinctOrientations; ++orientation)
		{
			while (motorcyclesPerColor[orientation].front().first < firstChoiceEntryPoint)
			{
				auto entry = motorcyclesPerColor[orientation].front();
				entry.first += maxEntryPoint + 1;
				motorcyclesPerColor[orientation].erase(motorcyclesPerColor[orientation].begin());
				motorcyclesPerColor[orientation].push_back(entry);
			}
		}

		for (int orientation = 1; orientation < distinctOrientations; ++orientation)
		{
			auto& currentColumn = dpTable[orientation - 1];
			currentColumn.clear();
			currentColumn.resize(motorcyclesPerColor[orientation].size());
			
			int minCurrentChoiceIndex;
			float minCurrentChoiceCost = std::numeric_limits<float>::infinity();
			for (int currentChoiceIndex = 0; currentChoiceIndex < motorcyclesPerColor[orientation].size(); ++currentChoiceIndex)
			{
				auto& currentChoice = motorcyclesPerColor[orientation][currentChoiceIndex];
				auto currentChoiceEntryPoint = currentChoice.first;
				auto& currentMotorcycleDesc = currentChoice.second;
				auto& currentMotorcycle = cycles.Motorcycles()[currentMotorcycleDesc.motorcycleIdx];
				float currentCost = currentMotorcycleDesc.cost;
				auto currentMotorcycleEdge = GetFirstPathEdge(currentMotorcycleDesc);

				//record all vertices visited by the current motorcycle to evaluate crossings
				std::set<HEMesh::VertexHandle> coveredVertices;
				for (int i = 0; i < currentMotorcycleDesc.lengthofMotorcyclePath; ++i)
					coveredVertices.insert(mesh.from_vertex_handle(currentMotorcycle.Edges().at(i)));

				if (orientation == distinctOrientations - 1)
				{
					//if this is the last color, go around to the beginning

					//check if the current and first motorcycles collide
					bool collision = false;
					for (int i = 0; i < firstMotorcycleDesc.lengthofMotorcyclePath; ++i)
						if (coveredVertices.find(mesh.from_vertex_handle(firstMotorcycle.Edges().at(i))) != coveredVertices.end())
						{
							collision = true;
							break;
						}
					if (collision)
						continue;

					currentCost += GetAngleEnergy(currentMotorcycleEdge, firstMotorcycleEdge);
				}

				//Find the best predecessor
				int prevChoicesEnd = motorcyclesPerColor[orientation - 1].size();
				if (orientation == 1)
					prevChoicesEnd = 1;
				int minPrevChoice = -1;
				float minPrevChoiceCost = std::numeric_limits<float>::infinity();
				for (int prevChoiceIndex = 0; prevChoiceIndex < prevChoicesEnd; ++prevChoiceIndex)
				{
					auto& prevChoice = (orientation == 1 ? firstChoice : motorcyclesPerColor[orientation - 1][prevChoiceIndex]);
					auto prevChoiceEntryPoint = prevChoice.first;
					if (prevChoiceEntryPoint >= currentChoiceEntryPoint)
						break; //wrong order
					auto& prevMotorcycleDesc = prevChoice.second;
					auto& prevMotorcycle = cycles.Motorcycles()[prevMotorcycleDesc.motorcycleIdx];
					float prevCost = (orientation == 1 ? firstCost : dpTable[orientation - 1 - 1][prevChoiceIndex].cost);
					if (std::isinf(prevCost))
						continue;
					auto prevMotorcycleEdge = GetFirstPathEdge(prevMotorcycleDesc);

					//check if the current and previous motorcycles collide
					bool collision = false;
					for (int i = 0; i < prevMotorcycleDesc.lengthofMotorcyclePath; ++i)
						if (coveredVertices.find(mesh.from_vertex_handle(prevMotorcycle.Edges().at(i))) != coveredVertices.end())
						{
							collision = true;
							break;
						}
					if (collision)
						continue;

					prevCost += GetAngleEnergy(prevMotorcycleEdge, currentMotorcycleEdge);

					if (prevCost < minPrevChoiceCost)
					{
						minPrevChoiceCost = prevCost;
						minPrevChoice = prevChoiceIndex;
					}
				}

				if (std::isinf(minPrevChoiceCost))
					continue;

				auto& dpEntry = dpTable[orientation - 1][currentChoiceIndex];
				dpEntry.cost = currentCost + minPrevChoiceCost;
				dpEntry.predecessorNode = minPrevChoice;

				if (dpEntry.cost < minCurrentChoiceCost)
				{
					minCurrentChoiceCost = dpEntry.cost;
					minCurrentChoiceIndex = currentChoiceIndex;
				}
			}

			if (std::isinf(minCurrentChoiceCost))
				break; //no solution possible, go to next colorZeroMotorcycle

			if (orientation == distinctOrientations - 1)
			{
				//we finished the propagation of the function
				if (minCurrentChoiceCost + centerednessEnergy < bestCostSoFar)
				{
					bestCostSoFar = minCurrentChoiceCost + centerednessEnergy;
					//track back
					auto choice = minCurrentChoiceIndex;
					for (int orientation = distinctOrientations - 1; orientation > 0; --orientation)
					{
						bestPaths[orientation] = motorcyclesPerColor[orientation][choice].second;
						choice = dpTable[orientation - 1][choice].predecessorNode;
					}
					bestPaths[0] = motorcyclesPerColor[0][colorZeroMotorcycle].second;
				}
			}
		}
	}

	return bestCostSoFar;
}