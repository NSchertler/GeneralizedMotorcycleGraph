#include "FencedRegion.h"

#include <iostream>
#include <queue>

// ------------  Loop  ------------

FencedRegion::Loop::Loop(int degree, std::vector<LoopEdge>&& edges, std::vector<HEMesh::HalfedgeHandle>&& concaveCornersOutgoing)
	: degree(degree), edges(edges), concaveCornersOutgoing(concaveCornersOutgoing)
{
	for(int i = 0; i < edges.size(); ++i)
		edgeToIndex[edges[i].edge] = i;
}

const std::vector<HEMesh::HalfedgeHandle>& FencedRegion::Loop::ConcaveCornersOutgoing() const
{
	return concaveCornersOutgoing;
}

size_t FencedRegion::Loop::FindEdge(HEMesh::HalfedgeHandle e) const
{
	auto it = edgeToIndex.find(e);
	if (it == edgeToIndex.end())
		return (size_t)-1;
	else
		return it->second;
}

bool FencedRegion::Loop::AddFaceGeneratesHandle(HEMesh::FaceHandle f, const HEMesh& mesh) const
{
	//adding a face generates a handle (or a hole) if not all consecutive face
	//edges are consecutive in the loop
	std::vector<size_t> edgeIndices;
	int firstMatchedEdge = -1;
	for (auto h : mesh.fh_range(f))
	{
		auto opp = mesh.opposite_halfedge_handle(h);
		edgeIndices.push_back(FindEdge(opp));
		if (firstMatchedEdge == -1 && edgeIndices.back() != (size_t)-1)
			firstMatchedEdge = edgeIndices.size() - 1;
	}
	if (firstMatchedEdge == -1)
		return true;

	bool hadGap = false;
	for (int i = firstMatchedEdge; i < firstMatchedEdge + edgeIndices.size(); ++i)
	{
		auto current = i % edgeIndices.size();
		auto next = (i + 1) % edgeIndices.size();
		if (edgeIndices[current] == (size_t)-1)
			hadGap = true;
		else
		{
			if (hadGap)
				return true;
			auto currentEdgeIndex = edgeIndices[current];
			auto nextEdgeIndex = edgeIndices[next];
			if (nextEdgeIndex == (size_t)-1)
				continue;
			if (currentEdgeIndex < nextEdgeIndex)
				currentEdgeIndex += edges.size();
			if (currentEdgeIndex != nextEdgeIndex + 1)
				return true;
		}
	}
	return false;
}

void FencedRegion::Loop::NormalizeEdgeOrientations()
{
	if (degree != 0)
	{
		//normalize edge orientations (positive and mod degree)
		auto distinctOrientations = std::abs(degree);
		for (auto& e : edges)
		{
			while (e.orientation < 0)
				e.orientation += distinctOrientations;
			e.orientation = e.orientation % distinctOrientations;
		}
	}
}

// ------------  FencedRegion  ------------

FencedRegion::FencedRegion(const HEMesh& mesh, const ManifoldnessAwareVertexMap<size_t>& vertexToSingularityIndex, bool verbose)
	: mesh(&mesh), verbose(verbose), vertexToSingularityIndex(&vertexToSingularityIndex)
{ }

const std::set<size_t>& FencedRegion::CoveredSingularities() const { return coveredSingularities; }

bool FencedRegion::HasSingularitiesOnBoundary() const { return !singularitiesOnBoundary.empty(); }

const size_t FencedRegion::NextSingularity()
{
	return *singularitiesOnBoundary.begin();
}

void FencedRegion::UpdateSingularitiesFromChanges()
{
	for (auto h : edgesToCheck)
	{
		const size_t* singPtr;		
		if (vertexToSingularityIndex->TryAccessAtToVertex(h, singPtr))
		{
			//There is a singularity incident to the halfedge

			size_t sing = *singPtr;

			coveredSingularities.insert(sing);
			if (IsToVertexOnBoundary(h))
				singularitiesOnBoundary.insert(sing);
			else
				singularitiesOnBoundary.erase(sing);

		}
	}
	edgesToCheck.clear();	
}

void FencedRegion::Clear()
{
	faces.clear();
	loops.clear();
	outgoing.clear();
	totalAdditionalEdges = 0;
	coveredSingularities.clear();
	singularitiesOnBoundary.clear();
}

const FencedRegion::OutgoingEdgeInfo& FencedRegion::OutgoingInfo(OpenMesh::VertexHandle v) const
{
	return outgoing.at(v);
}

size_t FencedRegion::BoundaryLength() const { return outgoing.size(); }
const std::set<HEMesh::FaceHandle>& FencedRegion::IncludedFaces() const { return faces; }

void FencedRegion::AddFace(HEMesh::FaceHandle f)
{
	if (faces.find(f) != faces.end())
		return; //face is already part of the patch
	faces.insert(f);
	for (auto h : mesh->fh_range(f))
	{
		auto hInverse = mesh->opposite_halfedge_handle(h);
		auto from_vertex = mesh->from_vertex_handle(h);
		auto to_vertex = mesh->to_vertex_handle(h);

		//check if the inverse edge exists
		bool inverseExists = DeleteEdge(hInverse);

		if (!inverseExists)
		{
			//The inverse of h does not exists, add h
			auto inserted = outgoing.insert(std::make_pair(from_vertex, OutgoingEdgeInfo(h)));
			if (!inserted.second)
			{
				//from_vertex already has outgoing edges
				inserted.first->second.additionalEdges.push_back(h);
				++totalAdditionalEdges;
			}
			edgesToCheck.push_back(h);
		}
	}
	UpdateSingularitiesFromChanges();
}

bool FencedRegion::FillMeshHole(HEMesh::HalfedgeHandle h, size_t maxHoleEdges)
{
	//add the faces incident to the hole	
	auto currentHalfedge = h;
	std::vector<HEMesh::HalfedgeHandle> holeHalfedges;
	std::vector<HEMesh::FaceHandle> facesToAdd;
	do
	{
		holeHalfedges.push_back(currentHalfedge);
		if (holeHalfedges.size() > maxHoleEdges)
			return false;

		CirculateBackwardUntil<true>(currentHalfedge, *mesh, [&](HEMesh::HalfedgeHandle h)
		{
			if (!mesh->is_boundary(h))
				facesToAdd.push_back(mesh->face_handle(h));
			return false;
		});
	} while (currentHalfedge != h);

	for (auto f : facesToAdd)
		AddFace(f);

	//remove the hole loop
	for (auto h : holeHalfedges)
	{
		DeleteEdge(mesh->opposite_halfedge_handle(h));
	}

	UpdateSingularitiesFromChanges();

	return true;
}

Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> FencedRegion::GenerateBoundaryIndices() const
{
	Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> result(2, outgoing.size() + totalAdditionalEdges);
	int i = 0;
	for (auto& edge : outgoing)
	{
		result(0, i) = edge.first.idx();
		result(1, i) = mesh->to_vertex_handle(edge.second.edge).idx();
		++i;
		for (auto& additional : edge.second.additionalEdges)
		{
			result(0, i) = edge.first.idx();
			result(1, i) = mesh->to_vertex_handle(additional).idx();
			++i;
		}
	}
	return result;
}

std::vector<unsigned int> FencedRegion::GenerateFaceIndices() const
{
	std::vector<unsigned int> indices;
	indices.reserve(4 * faces.size());
	for (auto& f : faces) //for each face
	{
		//span triangle fan
		OpenMesh::VertexHandle base;
		for (auto h : mesh->fh_range(f))
		{
			if (base.idx() == -1)
				base = mesh->from_vertex_handle(h);
			else if (mesh->to_vertex_handle(h) == base)
				break;
			else
			{
				indices.push_back(base.idx());
				indices.push_back(mesh->from_vertex_handle(h).idx());
				indices.push_back(mesh->to_vertex_handle(h).idx());
			}
		}
	}
	return indices;
}

bool FencedRegion::IsRectilinear() const
{
	for (auto& l : loops)
		if (l.IsSingular())
			return false;
	return true;
}

bool FencedRegion::IsFaceIncluded(HEMesh::FaceHandle face) const
{
	return faces.find(face) != faces.end();
}

void FencedRegion::MergeWith(FencedRegion& other)
{
	//TODO: there is probably a more efficient way...

	//std::set<HEMesh::HalfedgeHandle> edgesOnBoundary;
	//for (auto& o : outgoing)
	//{
	//	edgesOnBoundary.insert(o.second.edge);
	//	for (auto e : o.second.additionalEdges)
	//		edgesOnBoundary.insert(e);
	//}
	//for (auto& o : other.outgoing)
	//{
	//	edgesOnBoundary.insert(o.second.edge);
	//	for (auto e : o.second.additionalEdges)
	//		edgesOnBoundary.insert(e);
	//}


	for (auto f : other.faces)
		AddFace(f);

	//remove edges that had not been boundaries before (e.g. filled holes)
	//std::vector<HEMesh::HalfedgeHandle> edgesToDelete;
	//for (auto& o : outgoing)
	//{
	//	if (edgesOnBoundary.find(o.second.edge) == edgesOnBoundary.end())
	//		edgesToDelete.push_back(o.second.edge);
	//	for(auto e : o.second.additionalEdges)
	//		if (edgesOnBoundary.find(e) == edgesOnBoundary.end())
	//			edgesToDelete.push_back(e);
	//}
	//for (auto e : edgesToDelete)
	//	DeleteEdge(e);	
}

OpenMesh::HalfedgeHandle FencedRegion::NextOnBoundary(OpenMesh::HalfedgeHandle h) const
{
	auto out = outgoing.at(mesh->to_vertex_handle(h));
	if (out.additionalEdges.size() == 0)
		return out.edge;
	else
	{
		//non-manifold vertex; find the correct outgoing edge
		h = mesh->next_halfedge_handle(h);
		while (h != out.edge && std::find(out.additionalEdges.begin(), out.additionalEdges.end(), h) == out.additionalEdges.end())
			h = mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(h));
		return h;
	}
}

void FencedRegion::EstablishLoops()
{
	loops.clear();
	std::set<OpenMesh::HalfedgeHandle> handled;
	//iterate all edges on the boundary of this patch
	for (auto edge : outgoing)
	{
		auto currentH = edge.second.edge;
		//if the edge is not already part of a loop, trace a new loop starting from there
		if (handled.find(currentH) == handled.end())
			loops.push_back(CalculateLoop(currentH, handled));

		//do the same for the additional edges
		for (auto e : edge.second.additionalEdges)
		{
			if (handled.find(e) != handled.end())
				continue;

			loops.push_back(CalculateLoop(e, handled));
		}
	}

	MakeSimple();
}

void FencedRegion::MakeSimple()
{
	const int maxHoleSize = 200; //number of faces

	if (loops.size() <= 1)
		return;
	std::cout << "Non-simple patch at vertex " << mesh->from_vertex_handle(loops[0].Edges().front().edge).idx() << std::endl;

	//Check how many faces need to be added for all of the loops (outside of the patch) and keep the loop 
	//where this number is largest (and will consequentially not be filled).

	struct PerLoopInfo
	{
		std::set<HEMesh::FaceHandle> holeFaces;
		std::set<HEMesh::FaceHandle> facesOutsideOfHole;

		std::queue<HEMesh::FaceHandle> bfsQueue;
	};
	std::vector<PerLoopInfo> loopInfo(loops.size());

	//Perform simultaneous BFS on all loops until we found the largest loop

	//Initialize BFS
	for (int i = 0; i < loops.size(); ++i)
	{
		for (auto& e : loops[i].Edges())
		{
			auto f = mesh->opposite_face_handle(e.edge);
			if (f.is_valid() && loopInfo[i].holeFaces.find(f) == loopInfo[i].holeFaces.end())
			{
				loopInfo[i].holeFaces.insert(f);
				loopInfo[i].bfsQueue.push(f);
			}
			loopInfo[i].facesOutsideOfHole.insert(mesh->face_handle(e.edge));
		}
	}

	//do the actual BFS
	int largestHole = 0;
	while (true)
	{
		int loopsWithNonEmptyQueue = 0;
		int holesAboveMaxSize = 0;
		for (int i = 0; i < loops.size(); ++i)
		{
			if (!loopInfo[i].bfsQueue.empty())
			{
				auto f = loopInfo[i].bfsQueue.front();
				loopInfo[i].bfsQueue.pop();
				for (auto neighborH : mesh->fh_range(f))
				{
					auto neighborF = mesh->opposite_face_handle(neighborH);
					if (loopInfo[i].facesOutsideOfHole.find(neighborF) != loopInfo[i].facesOutsideOfHole.end())
						continue;
					if (loopInfo[i].holeFaces.find(neighborF) == loopInfo[i].holeFaces.end())
					{
						loopInfo[i].holeFaces.insert(neighborF);
						loopInfo[i].bfsQueue.push(neighborF);
					}
				}
			}
			if (!loopInfo[i].bfsQueue.empty())
				++loopsWithNonEmptyQueue;
			if (loopInfo[i].holeFaces.size() > loopInfo[largestHole].holeFaces.size())
				largestHole = i;
			if (loopInfo[i].holeFaces.size() > maxHoleSize)
				++holesAboveMaxSize;
		}
		if (holesAboveMaxSize >= 2)
		{
			std::cout << "Cannot make patch simple. The holes are too big." << std::endl;
			return;
		}
		if (loopsWithNonEmptyQueue == 0)
			break;

		if ((loopsWithNonEmptyQueue == 1 && !loopInfo[largestHole].bfsQueue.empty()) || loopsWithNonEmptyQueue == 0)
			break; //if the currently largest hole is the only growing one, we found the largest hole
	}

	//fill all loops except the one with the largest hole
	for (int i = 0; i < loops.size(); ++i)
	{
		if (i == largestHole)
			continue;

		//TODO: simple workaround...
		if (loopInfo[i].holeFaces.size() > maxHoleSize)
			return;
		//std::cout << "Closing hole with " << loopInfo[i].holeFaces.size() << " faces at " << mesh->point(mesh->from_vertex_handle(loops[i].Edges().front().edge)) << "." << std::endl;
		//continue;

		for (auto& e : loops[i].Edges())
			DeleteEdge(e.edge);

		for (auto f : loopInfo[i].holeFaces)
			for (auto h : mesh->fh_range(f))
			{
				const size_t* singPtr;
				if (vertexToSingularityIndex->TryAccessAtToVertex(h, singPtr))
					coveredSingularities.insert(*singPtr);
			}

		this->faces.insert(loopInfo[i].holeFaces.begin(), loopInfo[i].holeFaces.end());
	}
	if (largestHole < loops.size() - 1)
		loops.erase(loops.begin() + (largestHole + 1), loops.end());
	if (largestHole > 0)
		loops.erase(loops.begin(), loops.begin() + largestHole);
}

bool FencedRegion::IsToVertexOnBoundary(HEMesh::HalfedgeHandle h, bool considerPaddedBoundary) const
{
	auto v = mesh->to_vertex_handle(h);
	auto outgoingIt = outgoing.find(v);
	if (outgoingIt == outgoing.end())
		return false; //There is no record of the vertex in the boundary

	if (!considerPaddedBoundary)
		return true;

	if (!mesh->is_boundary(mesh->opposite_halfedge_handle(outgoingIt->second.edge)))
		return true;

	//check if the previous boundary edge is also in the loop (then we add padding to the boundary and the vertex is not on the patch's boundary).

	HEMesh::HalfedgeHandle outgoingEdge;
	if (mesh->is_manifold(v))
	{
		outgoingEdge = outgoingIt->second.edge;
		if (outgoingIt->second.additionalEdges.size() > 0)
			return true;
	}
	else
	{
		std::set<HEMesh::HalfedgeHandle> validEdges; //The edges in this surface patch
		auto action = [&](HEMesh::HalfedgeHandle halfedge)
		{
			validEdges.insert(halfedge);
			return false;
		};
		auto hCopy = h;
		CirculateForwardUntil<false>(hCopy, *mesh, action);

		if (validEdges.find(outgoingIt->second.edge) != validEdges.end())
			outgoingEdge = outgoingIt->second.edge;
		for (auto& e : outgoingIt->second.additionalEdges)
		{
			if (validEdges.find(e) != validEdges.end())
			{
				if (outgoingEdge.is_valid())
					return true;
				else
					outgoingEdge = e;
			}
		}
	}

	auto prevEdge = mesh->opposite_halfedge_handle(outgoingEdge);
	CirculateBackwardUntil<true>(prevEdge, *mesh, [](HEMesh::HalfedgeHandle) { return false; });
	prevEdge = mesh->opposite_halfedge_handle(prevEdge);
	auto prevOutgoingIt = outgoing.find(mesh->from_vertex_handle(prevEdge));
	if (prevOutgoingIt == outgoing.end())
		return true;
	if (prevOutgoingIt->second.edge == prevEdge)
		return false;
	for (auto h : prevOutgoingIt->second.additionalEdges)
		if (h == prevEdge)
			return false;
	return true;
}

const std::vector<FencedRegion::Loop>& FencedRegion::Loops() const { return loops; }

//Tries to remove an edge from the boundary. Returns if the edge existed.
bool FencedRegion::DeleteEdge(HEMesh::HalfedgeHandle h)
{
	auto fromVertex = mesh->from_vertex_handle(h);
	auto outgoingIt = outgoing.find(fromVertex);
	if (outgoingIt != outgoing.end())
	{
		//fromVertex has outgoing edges
		//check if one of the outgoing edges is h
		if (outgoingIt->second.edge == h)
		{
			//remove this edge
			if (outgoingIt->second.additionalEdges.size() > 0)
			{
				//There are more edges, move the last edge to the primary spot
				outgoingIt->second.edge = std::move(outgoingIt->second.additionalEdges.back());
				outgoingIt->second.additionalEdges.pop_back();
				--totalAdditionalEdges;
			}
			else
			{
				outgoing.erase(outgoingIt);
			}
			edgesToCheck.push_back(mesh->opposite_halfedge_handle(h));
			return true;
		}
		else
		{
			auto findHIt = std::find(outgoingIt->second.additionalEdges.begin(), outgoingIt->second.additionalEdges.end(), h);
			if (findHIt != outgoingIt->second.additionalEdges.end())
			{
				//h exists, remove it
				*findHIt = std::move(outgoingIt->second.additionalEdges.back());
				outgoingIt->second.additionalEdges.pop_back();
				--totalAdditionalEdges;
				edgesToCheck.push_back(mesh->opposite_halfedge_handle(h));
				return true;
			}
		}
	}
	return false;
}

FencedRegion::Loop FencedRegion::CalculateLoop(OpenMesh::HalfedgeHandle h, std::set<OpenMesh::HalfedgeHandle>& handled) const
{
	//accumulate the patch degree
	int degree = 0;

	std::vector<LoopEdge> edges;
	std::vector<HEMesh::HalfedgeHandle> concaveCornersOutgoing;
	
	while (true)
	{
		handled.insert(h);

		edges.emplace_back(h, degree);
		auto nextH = NextOnBoundary(h);

#if __DEBUG
		assert(!mesh->is_boundary(h));
#endif

		auto turns = TurnsBetweenEdges(h, nextH, *mesh, true);
		if (turns < 0)
			concaveCornersOutgoing.push_back(nextH);
		degree += turns;
		h = nextH;
		if (handled.find(nextH) != handled.end())
		{
			//loop is closed
			break;
		}
	}
	Loop result(degree, std::move(edges), std::move(concaveCornersOutgoing));
	result.NormalizeEdgeOrientations();
	return result;
}