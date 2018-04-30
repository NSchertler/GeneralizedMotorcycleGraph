#include "TexturePatch.h"
#include "MotorcycleGraph.h"

#include <nsessentials/math/LeastSquaresSystem.h>
#include <fstream>

#include <nsessentials/util/MathematicaFormatter.h>

#include "PolygonTriangulation.h"

bool IsGapBetween(const MotorcycleGraph::HalfArc& arc1, const MotorcycleGraph::HalfArc& arc2, const MotorcycleGraph& graph, const HEMesh& mesh)
{
	auto lastOf1 = graph.MotorcycleHalfedge(arc1.LastPathSegment());
	auto firstOf2 = graph.MotorcycleHalfedge(*arc2.begin());
	return mesh.to_vertex_handle(lastOf1) != mesh.from_vertex_handle(firstOf2);
}

TexturePatch::TexturePatch(MotorcycleGraph& graph)
	: graph(&graph)
{ }

void TexturePatch::PrepareBuild(std::vector<std::vector<size_t>>&& patchSides, std::set<HEMesh::FaceHandle>&& faces, std::vector<HEMesh::HalfedgeHandle>&& openBoundaryEdges)
{
	this->patchSides = patchSides;
	this->openBoundaryEdges = openBoundaryEdges;
	this->faces = faces;

	contextToIndex.clear();
	vertices.clear();
	filledHoles.clear();
}

HEMesh::HalfedgeHandle TexturePatch::GetContextEdgeOfToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	auto v = mesh.to_vertex_handle(h);
	if (verticesOnArcs.find(v) == verticesOnArcs.end())
		return mesh.halfedge_handle(v);
	else
	{
		CirculateBackwardUntil<true>(h, mesh, [&](HEMesh::HalfedgeHandle checkH)
		{
			return halfedgesOnArcs.find(mesh.opposite_halfedge_handle(checkH)) != halfedgesOnArcs.end();
		});
		return h;
	}
}

HEMesh::HalfedgeHandle TexturePatch::GetContextEdgeOfFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	return GetContextEdgeOfToVertex(mesh.prev_halfedge_handle(h), mesh);	
}

void TexturePatch::Build(const HEMesh& mesh)
{	
	//Record the vertices and halfedges that lie on the comprising halfarcs.
	for (auto& side : patchSides)
	{
		for (auto arcIdx : side)
		{
			auto& arc = graph->Halfarcs()[arcIdx];
			for (auto p : arc)
			{
				auto h = graph->MotorcycleHalfedge(p);
				verticesOnArcs.insert(mesh.to_vertex_handle(h));
				halfedgesOnArcs.insert(h);
			}
		}
	}

	//Assign contiguous vertex indices
	int nextVertexIndex = 0;
	for (auto f : faces)
	{
		for (auto h : mesh.fh_range(f))
		{
			auto context = GetContextEdgeOfToVertex(h, mesh);
			auto inserted = contextToIndex.emplace(context, nextVertexIndex);
			if (inserted.second)
			{
				vertices.push_back(mesh.to_vertex_handle(h));
				++nextVertexIndex;
			}
		}
	}

	//find which of the sides can grow (i.e. if they have open boundaries or the incident sides do not constrain them

	if (patchSides.size() == 4)
	{
		for (int k = 0; k < 4; ++k)
		{
			sidesCanGrow[k] = false;
			bool hasNonBoundaryArc = false;
			for (int i = 0; i < patchSides[k].size(); ++i)
			{
				auto arcIdx = patchSides[k][i];
				auto& arc = graph->Halfarcs()[arcIdx];
				auto arcEdge = graph->MotorcycleHalfedge(*arc.begin());

				bool isBoundary = mesh.is_boundary(arcEdge) || mesh.is_boundary(mesh.opposite_halfedge_handle(arcEdge));

				if (isBoundary)
				{
					sidesCanGrow[k] = true;
					break;
				}
				else
					hasNonBoundaryArc = true;

				//check if there is an open boundary between this arc and the previous one
				if (i != 0 && IsGapBetween(graph->Halfarcs()[patchSides[k][i - 1]], arc, *graph, mesh))
					sidesCanGrow[k] = true;
			}

			if (!hasNonBoundaryArc)
			{
				sidesCanGrow[k] = true;
				continue;
			}

			//check if the previous and next sides are fixed
			auto nextSide = (k + 1) % 4;
			auto prevSide = (k + 3) % 4;

			if (patchSides[nextSide].empty())
				sidesCanGrow[k] = true;
			else
			{
				auto& myLastArc = graph->Halfarcs()[patchSides[k].back()];
				auto& firstArcOfNext = graph->Halfarcs()[patchSides[nextSide].front()];
				if (IsGapBetween(myLastArc, firstArcOfNext, *graph, mesh))
					sidesCanGrow[k] = true;
			}

			if (patchSides[prevSide].empty())
				sidesCanGrow[k] = true;
			else
			{
				auto& myFirstArc = graph->Halfarcs()[patchSides[k].front()];
				auto& lastArcOfPrev = graph->Halfarcs()[patchSides[prevSide].back()];
				if (IsGapBetween(lastArcOfPrev, myFirstArc, *graph, mesh))
					sidesCanGrow[k] = true;
			}
		}
	}

	FillHoles(mesh);
}

void TexturePatch::FillHoles(const HEMesh& mesh)
{
	//try to fill holes	
	std::set<HEMesh::HalfedgeHandle> handledBoundary;
	std::set<HEMesh::HalfedgeHandle> removeFromBoundary;
	for (auto h : openBoundaryEdges)
	{
		if (handledBoundary.find(h) != handledBoundary.end())
			continue;

		h = mesh.opposite_halfedge_handle(h);
		std::vector<HEMesh::VertexHandle> holeVertices;
		std::vector<Eigen::Vector3f> holeCoordinates;
		std::vector<HEMesh::HalfedgeHandle> holeHalfedges;
		auto startVertex = mesh.from_vertex_handle(h);
		auto p = mesh.point(startVertex);
		holeCoordinates.emplace_back(p[0], p[1], p[2]);
		holeVertices.push_back(startVertex);
		auto v = mesh.to_vertex_handle(h);
		bool isInsidePatch;
		do
		{
			auto opp = mesh.opposite_halfedge_handle(h);
			handledBoundary.insert(opp);
			holeHalfedges.push_back(opp);

			holeVertices.push_back(v);
			p = mesh.point(v);
			holeCoordinates.emplace_back(p[0], p[1], p[2]);

			CirculateBackwardUntil<true>(h, mesh, [](HEMesh::HalfedgeHandle) { return false; });
			//h = mesh.next_halfedge_handle(h);
			v = mesh.to_vertex_handle(h);
			auto context = GetContextEdgeOfToVertex(h, mesh);
			isInsidePatch = contextToIndex.find(context) != contextToIndex.end();
		} while (v != startVertex && isInsidePatch);

		if (v == startVertex)
		{
			//we have found a hole
			for (auto h : holeHalfedges)
				removeFromBoundary.insert(h);

			auto triangulation = TriangulatePolygon(holeCoordinates);
			for (auto& tri : triangulation)
				filledHoles.push_back({ holeVertices[tri[0]], holeVertices[tri[1]], holeVertices[tri[2]] });
		}
	}

	openBoundaryEdges.erase(std::remove_if(openBoundaryEdges.begin(), openBoundaryEdges.end(),
		[&](HEMesh::HalfedgeHandle h) { return removeFromBoundary.find(h) != removeFromBoundary.end(); }), openBoundaryEdges.end());	
}

size_t TexturePatch::IdOfToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	auto context = GetContextEdgeOfToVertex(h, mesh);
	auto id = contextToIndex.find(context);
	if (id == contextToIndex.end())
		throw std::runtime_error("The specified vertex does not belong to this patch");
	return id->second;
}

size_t TexturePatch::IdOfFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	auto context = GetContextEdgeOfFromVertex(h, mesh);
	auto id = contextToIndex.find(context);
	if (id == contextToIndex.end())
		throw std::runtime_error("The specified vertex does not belong to this patch");
	return id->second;
}

size_t TexturePatch::IdOfInnerVertex(HEMesh::VertexHandle v, const HEMesh& mesh) const
{
	auto context = mesh.halfedge_handle(v);
	auto id = contextToIndex.find(context);
	if (id == contextToIndex.end())
		throw std::runtime_error("The specified vertex does not belong to this patch");

	return id->second;
}

HEMesh::VertexHandle TexturePatch::Vertex(size_t id) const
{
	return vertices[id];
}

void TexturePatch::InsertArcAfter(size_t arcToInsert, size_t afterArc)
{
	//Find afterArc
	bool found = false;
	for (auto& side : patchSides)
	{
		for (auto it = side.begin(); it != side.end(); ++it)
		{
			if (*it == afterArc)
			{
				//We found afterArc; insert arcToInsert
				++it;
				side.insert(it, arcToInsert);
				found = true;
				break;
			}
		}
		if (found)
			break;
	}
}

void TexturePatch::InsertArcBefore(size_t arcToInsert, size_t beforeArc)
{
	//Find beforeArc
	bool found = false;
	for (auto& side : patchSides)
	{
		for (auto it = side.begin(); it != side.end(); ++it)
		{
			if (*it == beforeArc)
			{
				//We found beforeArc; insert arcToInsert
				side.insert(it, arcToInsert);
				found = true;
				break;
			}
		}
		if (found)
			break;
	}
}

void TexturePatch::AddInteriorPatchLength(int dimension, float length)
{
	assert(dimension == 0 || dimension == 1);
	sumInteriorPatchLengths[dimension] += length;
	++countInteriorEdgeLengths[dimension];
}

float TexturePatch::GetAverageInteriorPatchLength(int dimension, int & outWeight) const
{
	assert(dimension == 0 || dimension == 1);
	outWeight = countInteriorEdgeLengths[dimension];
	return sumInteriorPatchLengths[dimension] / outWeight;
}

size_t TexturePatch::NumberOfVertices() const
{
	return contextToIndex.size();
}

// ----------  TextureCoordinatesStorage  ----------

TextureCoordinatesStorage::TextureCoordinatesStorage(const TexturePatch& patch)
	: patch(&patch)
{
	ResetTextureCoordinates();
}

const TexturePatch& TextureCoordinatesStorage::Patch() const
{
	return *patch;
}

const Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	return textureCoordinates.at(patch->IdOfToVertex(h, mesh));
}

const Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const
{
	return textureCoordinates.at(patch->IdOfFromVertex(h, mesh));
}

const Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtInnerVertex(HEMesh::VertexHandle v, const HEMesh& mesh) const
{
	return textureCoordinates.at(patch->IdOfInnerVertex(v, mesh));
}

const Eigen::Vector2f& TextureCoordinatesStorage::TexCoord(size_t i) const
{
	return textureCoordinates[i];
}

Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) { return const_cast<Eigen::Vector2f&>(static_cast<const TextureCoordinatesStorage*>(this)->TexCoordAtToVertex(h, mesh)); }
Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) { return const_cast<Eigen::Vector2f&>(static_cast<const TextureCoordinatesStorage*>(this)->TexCoordAtFromVertex(h, mesh)); }
Eigen::Vector2f& TextureCoordinatesStorage::TexCoordAtInnerVertex(HEMesh::VertexHandle v, const HEMesh& mesh) { return const_cast<Eigen::Vector2f&>(static_cast<const TextureCoordinatesStorage*>(this)->TexCoordAtInnerVertex(v, mesh)); }
Eigen::Vector2f& TextureCoordinatesStorage::TexCoord(size_t i) { return const_cast<Eigen::Vector2f&>(static_cast<const TextureCoordinatesStorage*>(this)->TexCoord(i)); }

size_t TextureCoordinatesStorage::TexCoordCount() const
{
	return textureCoordinates.size();
}

void TextureCoordinatesStorage::ResetTextureCoordinates()
{
	textureCoordinates.clear();
	textureCoordinates.resize(patch->NumberOfVertices(), Eigen::Vector2f::Constant(std::numeric_limits<float>::quiet_NaN()));
}

void TextureCoordinatesStorage::GetMinMaxTextureCoordinates(Eigen::Vector2f& min, Eigen::Vector2f& max) const
{
	min.setConstant(std::numeric_limits<float>::infinity());
	max.setConstant(-std::numeric_limits<float>::infinity());
	for (auto& entry : textureCoordinates)
	{
		for (int i = 0; i < 2; ++i)
		{
			if (entry[i] < min[i])
				min[i] = entry[i];
			if (entry[i] > max[i])
				max[i] = entry[i];
		}
	}

	if (std::isinf(min.x()))
	{
		min.setZero();
		max.setZero();
	}
}