#pragma once

#include "common.h"
#include <set>

class MotorcycleGraph;

//Represents a rectangular domain for parametrization.
class TexturePatch
{
public:
	//Instantiates the texture patch.
	TexturePatch(MotorcycleGraph& graph);

	//Prepares the build of the patch by saving the provided data in the object.
	void PrepareBuild(std::vector<std::vector<size_t>>&& patchSides, std::set<HEMesh::FaceHandle>&& faces, std::vector<HEMesh::HalfedgeHandle>&& openBoundaryEdges);

	//Builds the patch doing the following steps:
	// * Gathers verticesOnArcs and halfedgesOnArcs
	// * Assigns contiguous vertex indices to every vertex included in this patch.
	// * Determines which of the sides can grow (i.e. that have an open boundary)
	// * Fills holes in the patch ( == holes in the original mesh)
	void Build(const HEMesh& mesh);

	//Return the number of corners for this patch. Should always be 4.
	size_t Corners() const { return patchSides.size(); }

	//Returns the patch-internal index of the specified vertex.
	size_t IdOfToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;
	//Returns the patch-internal index of the specified vertex.
	size_t IdOfFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;
	//Returns the patch-internal index of the specified vertex.
	size_t IdOfInnerVertex(HEMesh::VertexHandle h, const HEMesh& mesh) const;
	
	//Returns the vertex with the given patch-internal index.
	HEMesh::VertexHandle Vertex(size_t id) const;
	
	//Inserts a new arc after another arc on the patch boundary. If the other
	//arc is not part of this patch, nothing is done.
	void InsertArcAfter(size_t arcToInsert, size_t afterArc);

	//Inserts a new arc before another arc on the patch boundary. If the other
	//arc is not part of this patch, nothing is done.
	void InsertArcBefore(size_t arcToInsert, size_t beforeArc);

	//Returns the number of vertices in this patch. If a vertex is present twice (e.g. on the left edge and the right
	//edge), it is counted twice.
	size_t NumberOfVertices() const;

	//Records an interior patch length (i.e. length from one side to the opposite).
	//If dimension=0, the length is from side0 to side2. If dimension=1, the length is from side1 to side3.
	void AddInteriorPatchLength(int dimension, float length);

	//Returns the average of all recorded interior patch lengths for the given dimension.
	float GetAverageInteriorPatchLength(int dimension, int& outWeight) const;

	//Returns a list of (4) patch sides. Each patch side contains a list of the comprising halfarc indices.
	const std::vector<std::vector<size_t>>& PatchSides() const { return patchSides; }

	//Returns if a given side can grow (i.e. if it has an open boundary)
	bool CanSideGrow(size_t side) const { return sidesCanGrow[side]; }

	//Returns a list of all faces covered by this patch.
	const std::set<OpenMesh::FaceHandle>& Faces() const { return faces; }

	//Returns a lost of all open boundary edges of this patch (i.e. edges that are not part of any arc).
	const std::vector<HEMesh::HalfedgeHandle>& OpenBoundaryEdges() const { return openBoundaryEdges; }

	//Returns triangles that have been generated to fill holes.
	const std::vector<std::array<HEMesh::VertexHandle, 3>>& FilledHoles() const { return filledHoles; }
private:

	//Fills holes in the patch by triangulation.
	void FillHoles(const HEMesh& mesh);

	MotorcycleGraph* graph;
	
	//Sum of motorcycle lengths traversing the patch in the two possible directions
	float sumInteriorPatchLengths[2] = { 0 };
	int countInteriorEdgeLengths[2] = { 0 };

	std::set<HEMesh::VertexHandle> verticesOnArcs;
	std::set<HEMesh::HalfedgeHandle> halfedgesOnArcs;

	//Maps context edges to vertex indices
	std::map<HEMesh::HalfedgeHandle, size_t> contextToIndex;

	HEMesh::HalfedgeHandle GetContextEdgeOfToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;
	HEMesh::HalfedgeHandle GetContextEdgeOfFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;	

	std::vector<HEMesh::VertexHandle> vertices;
	std::set<OpenMesh::FaceHandle> faces;
	std::vector<std::vector<size_t>> patchSides; //indices to motorcycle graph arcs
	std::vector<HEMesh::HalfedgeHandle> openBoundaryEdges;

	bool sidesCanGrow[4];

	std::vector<std::array<HEMesh::VertexHandle, 3>> filledHoles;
};

//Represents storage space for the texture coordinates of a patch.
class TextureCoordinatesStorage
{
public:
	TextureCoordinatesStorage(const TexturePatch& patch);

	//Returns the patch associated with this storage object.
	const TexturePatch& Patch() const;

	Eigen::Vector2f& TexCoordAtToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh);
	Eigen::Vector2f& TexCoordAtFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh);
	Eigen::Vector2f& TexCoordAtInnerVertex(HEMesh::VertexHandle v, const HEMesh& mesh);
	Eigen::Vector2f& TexCoord(size_t i);

	const Eigen::Vector2f& TexCoordAtToVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;
	const Eigen::Vector2f& TexCoordAtFromVertex(HEMesh::HalfedgeHandle h, const HEMesh& mesh) const;
	const Eigen::Vector2f& TexCoordAtInnerVertex(HEMesh::VertexHandle v, const HEMesh& mesh) const;
	const Eigen::Vector2f& TexCoord(size_t i) const;

	//Returns the number of texture coordinates stored in this object.
	size_t TexCoordCount() const;

	//Resizes the internal container to hold enough texture coordinates for the associated set and
	//sets all to nan.
	void ResetTextureCoordinates();

	//Calculates the minimum and maximum texture coordinate (i.e. the bounding box).
	void GetMinMaxTextureCoordinates(Eigen::Vector2f& min, Eigen::Vector2f& max) const;
private:
	const TexturePatch* patch;

	std::vector<Eigen::Vector2f> textureCoordinates;
};