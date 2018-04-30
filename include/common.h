#pragma once

#include <Eigen/Core>
#include <vector>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

typedef Eigen::Matrix<float, 2, Eigen::Dynamic>					Matrix2Xf;
typedef Eigen::Matrix<float, 3, Eigen::Dynamic>					Matrix3Xf;
typedef Eigen::Matrix<float, 4, Eigen::Dynamic>					Matrix4Xf;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;
typedef Eigen::Matrix<unsigned short, 3, Eigen::Dynamic>		Matrix3Xus;

typedef std::vector<std::vector<unsigned int>> FaceList;

typedef OpenMesh::PolyMesh_ArrayKernelT<>  HEMesh;

//Converts an OpenMesh vector to an Eigen vector
static Eigen::Vector3f ToEigenVector(const OpenMesh::Vec3f& v)
{
	return Eigen::Vector3f(v[0], v[1], v[2]);
}

enum LengthMeasure
{
	Geometrical,
	Topological,
};

//Returns the number of faces incident to the given vertex
extern int FaceValence(OpenMesh::VertexHandle v, const HEMesh& mesh);

//Returns if a manifold vertex is a singularity
extern bool IsSingularityManifoldVertex(OpenMesh::VertexHandle v, const HEMesh& mesh, int& valenceDefect);
//Returns if a manifold vertex is a singularity
extern bool IsSingularityManifoldVertex(OpenMesh::VertexHandle v, const HEMesh& mesh);

//Returns if a non-manifold vertex is a singularity. Non-manifold vertices are split into several manifold ones per
//incident boundary loop. The context edge determines what representative to process.
//contextOutgoingBoundary - A boundary edge (outgoing from the examined vertex) that represents the examined boundary loop
extern bool IsSingularityNonManifoldVertex(HEMesh::HalfedgeHandle contextOutgoingBoundary, const HEMesh& mesh, int& valenceDefect);

//Returns if a non-manifold vertex is a singularity. Non-manifold vertices are split into several manifold ones per
//incident boundary loop. The context edge determines what representative to process.
//contextOutgoingBoundary - A boundary edge (outgoing from the examined vertex) that represents the examined boundary loop
extern bool IsSingularityNonManifoldVertex(HEMesh::HalfedgeHandle contextOutgoingBoundary, const HEMesh& mesh);

//Returns if the to vertex of h is a singularity
extern bool IsToVertexSingularity(HEMesh::HalfedgeHandle h, const HEMesh& mesh);

//Performs topological Catmull Clark subdivision of the input mesh.
//subdivisionInfo - Output variable. Stores for every input face the indices of all subdivided output faces.
extern void CatmullClarkSubdivide(FaceList& F, Matrix3Xf& V, std::vector<std::vector<unsigned int>>& subdivisionInfo);

//Merges neighboring triangles in the given mesh if their union becomes a rectangle.
//cosAngleThreshold - specifies the angle tolerance for the final quads (deviation from 90°)
extern void MergeTriangulatedQuads(FaceList& F, const Matrix3Xf& V, float cosAngleThreshold);

//Represents a possible continuation for a motorcycle path
struct PathContinuation
{
	//The next edge to traverse
	HEMesh::HalfedgeHandle edge;

	//A score for this continuation that measures straightness deviation in [-1, 1]
	float score;

	PathContinuation() { }

	PathContinuation(HEMesh::HalfedgeHandle e)
		: edge(e), score(1)
	{ }

	PathContinuation(HEMesh::HalfedgeHandle e, float score)
		: edge(e), score(score)
	{ }
};

//Finds all possible continuations for a motorcycle sitting at the given edge.
extern void FindPossibleContinuations(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, std::vector<PathContinuation>& outContinuations, const float minCosDeviationAngle = 0.5f, bool allowTurnsOnRegularVertices = false);

//Finds the best (in terms of straightness) continuation for a motorcycle sitting at the given edge.
//Returns true if a valid continuation has been found.
extern bool FindBestContinuation(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, PathContinuation& outContinuation, const float minCosDeviationAngle = 0.5f, bool allowTurnsOnRegularVertices = false);

//Finds the continuation of a motorcycle sitting at the given edge, assuming that the target vertex is regular.
//Returns true if a valid continuation has been found.
extern bool RegularContinuation(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, HEMesh::HalfedgeHandle& outContinuation);

//Projects an arbitrary vector onto a plane defined by the given normal.
extern void ProjectToTangentSpace(OpenMesh::Vec3f& vec, const OpenMesh::Vec3f& normalizedNormal);

//Prints a halfedge representation to the given stream in the form "(from vertex) -- halfedge index --> "
extern void PrintHalfedge(std::ostream& stream, HEMesh::HalfedgeHandle h, const HEMesh& mesh);

//Prints a full halfedge representation to the given stream in the form "(from vertex) -- halfedge index --> (to vertex)"
extern void PrintFullHalfedge(std::ostream& stream, HEMesh::HalfedgeHandle h, const HEMesh& mesh);

//Calculates the number of turns to go from h1 to h2 (to_vertex(h1) == from_vertex(h2)).
//Both h1 and h2 must be non-boundary halfedges.
extern int TurnsBetweenEdges(HEMesh::HalfedgeHandle h1, HEMesh::HalfedgeHandle h2, const HEMesh& mesh, bool useMargin);

//Calculates the vertex normal as a weighted average of face normals (using angle weights)
extern void CalcVertexNormalAngleWeights(const HEMesh& mesh, OpenMesh::VertexHandle _vh, OpenMesh::Vec3f& _n);

//Describes the result of the circulation functions
enum CirculationResult
{
	//The circulation has been stopped because it has reached its starting point.
	ReachedStartingPoint,

	//The circulation has been stopped because it has reached a boundary.
	StoppedByBoundary,

	//The circulation has been stopped by the stopping condition.
	StoppedByCondition,
};

//Takes a halfedge (incoming to vertex v) and circulates the outgoing halfedge in positive direction while a condition is met.
//The input halfedge will be the last in the order (if the circulation is not stopped earlier).
//Condition: bool(HEMesh::HalfedgeHandle)
template <bool StopAtBoundary, typename Condition>
CirculationResult CirculateForwardUntil(HEMesh::HalfedgeHandle& h, const HEMesh& mesh, Condition&& stopCondition);

//Takes a halfedge (incoming to vertex v) and circulates the outgoing halfedge in negative direction while a condition is met.
//The input halfedge will be the first in the order (if the circulation is not stopped earlier).
//Condition: bool(HEMesh::HalfedgeHandle)
template <bool StopAtBoundary, typename Condition>
CirculationResult CirculateBackwardUntil(HEMesh::HalfedgeHandle& h, const HEMesh& mesh, Condition&& stopCondition);

//Behavior of the forward circulation
template <bool StopAtBoundary>
struct CirculateForwardTraits
{ };

template <>
struct CirculateForwardTraits <true>
{
	static void ResolveBoundary(const HEMesh& mesh, HEMesh::HalfedgeHandle& h) 
	{
		h = mesh.opposite_halfedge_handle(h);
	}
	static bool ReturnAtBoundary() { return true; }
};

template <>
struct CirculateForwardTraits <false>
{
	static void ResolveBoundary(const HEMesh& mesh, HEMesh::HalfedgeHandle& h)
	{
		CirculateBackwardUntil<true>(h, mesh, [](HEMesh::HalfedgeHandle) { return false; });
	}
	static bool ReturnAtBoundary() { return false; }
};

//Takes a halfedge (incoming to vertex v) and circulates the outgoing halfedge in positive direction while a condition is met.
//The input halfedge will be the last in the order (if the circulation is not stopped earlier).
//Condition: bool(HEMesh::HalfedgeHandle)
template <bool StopAtBoundary, typename Condition>
CirculationResult CirculateForwardUntil(HEMesh::HalfedgeHandle& h, const HEMesh& mesh, Condition&& stopCondition)
{
	auto startEdge = mesh.opposite_halfedge_handle(h);
	h = startEdge;
	do
	{
		h = mesh.opposite_halfedge_handle(h);
		if (mesh.is_boundary(h))
		{			
			CirculateForwardTraits<StopAtBoundary>::ResolveBoundary(mesh, h);
			if (CirculateForwardTraits<StopAtBoundary>::ReturnAtBoundary())
				return StoppedByBoundary;
		}
		else
			h = mesh.next_halfedge_handle(h);
		if (std::forward<Condition>(stopCondition)(h))
			return StoppedByCondition;
	} while (h != startEdge);
	return ReachedStartingPoint;
}

template <bool stopAtBoundary>
struct CirculateBackwardTraits
{ };

template <>
struct CirculateBackwardTraits <true>
{ 
	static void ResolveBoundary(const HEMesh& mesh, HEMesh::HalfedgeHandle& h)	{ }
	static bool ReturnAtBoundary() { return true; }
};

template <>
struct CirculateBackwardTraits <false>
{
	static void ResolveBoundary(const HEMesh& mesh, HEMesh::HalfedgeHandle& h) 
	{
		h = mesh.opposite_halfedge_handle(h);
		CirculateForwardUntil<true>(h, mesh, [](HEMesh::HalfedgeHandle) { return false; });
	}
	static bool ReturnAtBoundary() { return false; }
};

//Takes a halfedge (incoming to vertex v) and circulates the outgoing halfedge in negative direction while a condition is met.
//The input halfedge will be the first in the order (if the circulation is not stopped earlier).
//Condition: bool(HEMesh::HalfedgeHandle)
template <bool stopAtBoundary, typename Condition>
CirculationResult CirculateBackwardUntil(HEMesh::HalfedgeHandle& h, const HEMesh& mesh, Condition&& stopCondition)
{
	auto startEdge = mesh.opposite_halfedge_handle(h);
	h = startEdge;
	do
	{
		if (std::forward<Condition>(stopCondition)(h))
			return StoppedByCondition;

		if (mesh.is_boundary(h))
		{
			CirculateBackwardTraits<stopAtBoundary>::ResolveBoundary(mesh, h);
			if(CirculateBackwardTraits<stopAtBoundary>::ReturnAtBoundary())
				return StoppedByBoundary;
		}
		else
			h = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(h));		
		
	} while (h != startEdge);
	return ReachedStartingPoint;
}