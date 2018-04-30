#pragma once

#include "common.h"
#include "ManifoldnessAwareVertexMap.h"
#include <set>

//Represents a fenced region that covers one or more singularities.
class FencedRegion
{
public:
	FencedRegion(const HEMesh& mesh, const ManifoldnessAwareVertexMap<size_t>& vertexToSingularityIndex, bool verbose = true);

	FencedRegion(FencedRegion&& movedFrom) = default;
	FencedRegion& operator=(FencedRegion&& movedFrom) = default;

	//Returns a set of singularity indices covered by this fenced region.
	const std::set<size_t>& CoveredSingularities() const;

	//Returns if the fenced region has singularities on its boundary and is thus invalid.
	bool HasSingularitiesOnBoundary() const;

	//Returns an arbitrary singularity on the boundary of this fenced region.
	const size_t NextSingularity();

	//Returns if this fenced region has not been deactivated
	bool IsActive() const { return isActive; }

	//Lets the user deactivate a fenced region.
	void SetActive(bool active) { isActive = active; }

	//Represents incidence relations between a vertex and all outgoing
	//edges on the fence.
	struct OutgoingEdgeInfo
	{
		OutgoingEdgeInfo(OpenMesh::HalfedgeHandle edge)
			: edge(edge)
		{ }

		//An arbitrary outgoing boundary edge
		OpenMesh::HalfedgeHandle edge;

		//All other outgoing edges if there are more than one (in the case of
		//a non-manifold boundary)
		std::vector<OpenMesh::HalfedgeHandle> additionalEdges;
	};

	//Represents an edge in a boundary loop of a fenced region.
	struct LoopEdge
	{
		HEMesh::HalfedgeHandle edge;

		//Integer orientation of this edge in [0, n), where
		//n is the loop's degree.
		int orientation;

		LoopEdge() { }
		LoopEdge(HEMesh::HalfedgeHandle e, unsigned int orientation)
			: edge(e), orientation(orientation)
		{ }
	};

	//Represents a boundary loop of a fenced region.
	class Loop
	{
	public:
		Loop(int degree, std::vector<LoopEdge>&& edges, std::vector<HEMesh::HalfedgeHandle>&& concaveCornersOutgoing);

		int Degree() const { return degree; }

		//Returns the ordered set of edges of this loop
		const std::vector<LoopEdge>& Edges() const { return edges; }		

		//Returns a list of concave corners of this loop (represented as outgoing halfedges)
		const std::vector<HEMesh::HalfedgeHandle>& ConcaveCornersOutgoing() const;

		//Returns if this loop has a degree other than 4
		bool IsSingular() const { return degree != 4 && degree != -4; }

		//Returns the index of the given halfedge within this loop or (size_t)-1 if the
		//halfedge is not part of this loop.
		size_t FindEdge(HEMesh::HalfedgeHandle e) const;

		//Returns if adding a face inside this loop generates a handle on the boundary.
		//This happens if the face merges two parts of the enclosed region.
		bool AddFaceGeneratesHandle(HEMesh::FaceHandle f, const HEMesh& mesh) const;		

		//Reduces all edge modulo Degree.
		void NormalizeEdgeOrientations();

	private:
		//number of corners of this loop
		int degree = 0;

		std::vector<LoopEdge> edges;
		std::map<HEMesh::HalfedgeHandle, size_t> edgeToIndex;

		//All edges outgoing from concave vertices of this loop
		std::vector<HEMesh::HalfedgeHandle> concaveCornersOutgoing;
	};

	//Returns a list of outgoing boundary edges for a vertex on the fence
	const OutgoingEdgeInfo& OutgoingInfo(OpenMesh::VertexHandle v) const;

	//Returns the number of edges comprising the fence
	size_t BoundaryLength() const;

	//Returns a list of faces that are covered by this region
	const std::set<HEMesh::FaceHandle>& IncludedFaces() const;

	//Adds a face to the fenced region
	void AddFace(HEMesh::FaceHandle f);

	//Adds the faces surrounding a hole to the current patch and removes the hole loop.
	//h is a boundary halfedge that specifies the hole to fill.
	//Returns true if the hole has been added. Returns false if the hole was bigger than the given threshold.
	bool FillMeshHole(HEMesh::HalfedgeHandle h, size_t maxHoleEdges);

	//Returns an index list (of vertices) for rendering the boundary,
	Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> GenerateBoundaryIndices() const;

	//Returns an index list (of vertices) for rendering the faces as triangles.
	std::vector<unsigned int> GenerateFaceIndices() const;

	//Returns the halfedge following a given halfedge on the boundary.
	OpenMesh::HalfedgeHandle NextOnBoundary(OpenMesh::HalfedgeHandle h) const;

	//Calculates the loops for this fenced region and calls MakeSimple afterwards.
	void EstablishLoops();

	//Fills any holes in the patch and adds the according faces to the patch.
	void MakeSimple();

	//Returns a list of loops that border this fenced region. A valid fenced region has a
	//single loop.
	const std::vector<Loop>& Loops() const;

	//Returns if all loops are of degree 4.
	bool IsRectilinear() const;

	//Returns if a face is covered by this region.
	bool IsFaceIncluded(HEMesh::FaceHandle face) const;

	//Merges this region with another one by copying all faces of other into
	//this region. Invalidates existing loops.
	void MergeWith(FencedRegion& other);

	//Returns if the to vertex of h is on the fence of this region. If considerPaddedBoundary
	//is true, a vertex is only considered on the boundary if the surrounding edges are not virtually 
	//padded.
	bool IsToVertexOnBoundary(HEMesh::HalfedgeHandle h, bool considerPaddedBoundary = true) const;

	//Resets the region.
	void Clear();

private:
	//Updates the set of covered singularities from the changes recorded in edgesToCheck.
	void UpdateSingularitiesFromChanges();

	//Modified edges that need to be checked for covered singularities
	std::vector<HEMesh::HalfedgeHandle> edgesToCheck;

	//Tries to remove an edge from the boundary. Returns if the edge existed.
	bool DeleteEdge(HEMesh::HalfedgeHandle h);

	//Calculates a loop starting at a given halfedge. Adds all encountered halfedges to handled.
	Loop CalculateLoop(OpenMesh::HalfedgeHandle h, std::set<OpenMesh::HalfedgeHandle>& handled) const;

	//Determines if functions of this object should output debug information
	bool verbose;

	//represent boundary loop as a linked list
	std::map<OpenMesh::VertexHandle, OutgoingEdgeInfo> outgoing;

	//Total number of all additional edges in the outgoing structure.
	size_t totalAdditionalEdges = 0;

	//Set of all covered faces
	std::set<OpenMesh::FaceHandle> faces;

	//Set of generated boundary loops for this fenced region
	std::vector<Loop> loops;

	//The underlying mesh
	const HEMesh* mesh;	

	//indices into singularities
	std::set<size_t> singularitiesOnBoundary; 

	//indices into singularities	
	std::set<size_t> coveredSingularities;	

	const ManifoldnessAwareVertexMap<size_t>* vertexToSingularityIndex;

	bool isActive = true;
};