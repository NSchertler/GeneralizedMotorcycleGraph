#pragma once

#include "common.h"
#include "PatchSet.h"
#include "TexturePatch.h"
#include "ManifoldnessAwareVertexMap.h"
#include "Singularity.h"

#include <vector>
#include <set>

#include <nsessentials/util/Observer.h>

class MotorcycleGraph
{
public:

	//Represents a motorcycle that travels across the mesh together with its trail.
	struct Motorcycle
	{
		Motorcycle() { }

		//Instantiates a motorcycle
		//h - The current position of the motorcycle
		Motorcycle(HEMesh::HalfedgeHandle h, bool highPriority)
			: currentPosition(h), highPriority(highPriority)
		{ }

		//The current position of the motorcycle. The motorcycle is assumed to sit
		//at the beginning of this halfedge and is going to traverse it.
		HEMesh::HalfedgeHandle currentPosition;		

		//saves the locations on the path where other motorcycles crashed into this motorcycle
		std::set<size_t> crashesOnPath; 

		//Index of a motorcycle that continues in a 180° angle (not a 90° angle as usual) at
		//the start of this motorcycle. This happens when two motorcycles collide head-on or 
		//if a valence-2 singularity is introduced. A value of (size_t)-1 specifies no continuation.
		size_t straightContinuationStart = (size_t)-1;

		//Index of a motorcycle that continues in a 180° angle (not a 90° angle as usual) at
		//the end of this motorcycle. This happens when two motorcycles collide head-on or 
		//if a valence-2 singularity is introduced. A value of (size_t)-1 specifies no continuation.
		size_t straightContinuationEnd = (size_t)-1;

		//Determines if this motorcycle has already terminated and will not move further.
		bool terminated = false;

		//Determines if this motorcycle collided with the padded boundary.
		bool collideWithPaddedBoundary = false;

		//Determines if this motorcycle start from the padded boundary.
		bool startInPaddedBoundary = false;

		bool highPriority = false;

		//Determines if this motorcycle is ignored for patch extraction.
		bool isDeactivated = false; 

		//temporary data, only used for sorting
		int temporaryContextEdge; 

		//The multiplier of the parametric sizes between the two sides of the motorcycle.
		float multiplier = 1; 

		//Returns a list of all previously traversed halfedges.
		const std::vector<HEMesh::HalfedgeHandle>& Path() const { return path; }

		//Adds a halfedge to the path of this motorcycle.
		void AddToPath(HEMesh::HalfedgeHandle h, const HEMesh& mesh);

		//Removes the last halfedge from the path of this motorcycle.
		void RemoveLastPathSegment(const HEMesh& mesh);

		//Returns the number of non-boundary edges in the path of this motorcycle.
		size_t NonBoundaryEdges() const { return nonBoundaryEdges; }

	private:
		std::vector<HEMesh::HalfedgeHandle> path;
		size_t nonBoundaryEdges = 0;
	};

	//Represents a point on the path of a motorcycle
	struct MotorcycleStop
	{
		MotorcycleStop() = default;
		MotorcycleStop(size_t motorcycleId, size_t locationInMotorcyclePath);

		//The id of the motorcycle
		size_t motorcycle;

		//The index of the location on the motorcycle's path in [0, n], where n is the 
		//number of path halfedges. 0 specifies the start of the path. n specifies the end.
		size_t locationInMotorcyclePath;
	};

	//Represents a directed location on the path of a motorcycle.
	struct LocationOnPath
	{
		//The id of the motorcycle
		size_t motorcycle;

		//The number of the halfedge within the motorcycle's path
		size_t pathSegment;

		//Determines if this location is directed in the same direction as the
		//original halfedge. If so, the location specifies the beginning of the
		//halfedge (going forward). Otherwise, the location specifies the end
		//of the halfedge (going backward).
		bool goForward;

		LocationOnPath();
		LocationOnPath(size_t motorcycle, size_t pathSegment, bool goForward);

		bool operator==(const LocationOnPath& other) const;
	};

	//Holds data gathered during patch extraction
	struct ExtractionStatistics
	{
		//The total number of faces covered by patches
		size_t totalFaces = 0;

		//The total number of open patches (i.e. patches whose
		//extraction failed).
		size_t openPatches = 0;

		//The maximum degree of all extracted patches.
		int maxPatchDegree = 0;

		//The number of patches, grouped by their degree
		std::vector<uint32_t> patchCountPerNumberOfCorners;

		//The total number of faces covered by patches of the respective degree
		std::vector<uint32_t> patchSizesPerNumberOfCorners;

		//A list of all patch sizes (number of covered faces)
		std::vector<size_t> patchSizes;		

		void Clear();
	};

	//Represents a half arc of the patch graph. The arc itself is made up of partial
	//motorcycle trails.
	struct HalfArc
	{
		//Iterates the halfedges of a halfarc
		class PathSegmentIterator : public std::iterator<std::forward_iterator_tag, LocationOnPath>
		{			
		public:
			PathSegmentIterator() = default;
			PathSegmentIterator(int iMotorcycle, int iPathSegment, const HalfArc* arc);

			bool operator==(const PathSegmentIterator& other) const;
			bool operator!=(const PathSegmentIterator& other) const;
			LocationOnPath operator*() const;
			PathSegmentIterator& operator++();

		private:
			int iMotorcycle;
			int iPathSegmentOnMotorcycle;
			const HalfArc* arc;

			friend struct HalfArc;
		};

		//Represents a segment of a halfarc, i.e. a contiguous set of halfedges on a single motorcycle
		struct Segment
		{
			LocationOnPath location;
			size_t length;

			Segment() = default;
			Segment(const LocationOnPath& location, size_t length);
		};

		//The segments that constitute this halfarc
		std::vector<Segment> segments;

		//The id of the patch incident to this halfarc
		size_t face = (size_t)-1;

		//The parametric size multiplier for this halfarc. All multipliers are 1 except the ones that
		//have been introduced by patch splitting.
		float multiplier;

		
		HalfArc() = default;
		HalfArc(const PathSegmentIterator& first, const PathSegmentIterator& last);

		//Returns the length (number of halfedges) of this halfarc.
		size_t Length() const;

		//Merges two consecutive halfarcs into one by moving the segments of the provided
		//argument to this halfarc.
		void MergeWith(HalfArc&& merged);

		//Returns the last halfedge of this halfarc
		LocationOnPath LastPathSegment() const;				

		//Trims the halfarc by removing everything before start. start must be
		//part of the halfarc.
		void SetStartInclusive(const PathSegmentIterator& start);

		//Trims the halfarc by removing end and everything after end. end must
		//be part of the halfarc.
		void SetEndExclusive(const PathSegmentIterator& end);

		PathSegmentIterator begin() const;
		PathSegmentIterator end() const;

		friend class PathSegmentIterator;
	};

	//Represents the various states that a motorcycle graph might be in.
	enum TracingStatus
	{
		//The graph is in the process of being calculated.
		InProcess,

		//Graph construction is finished.
		Finished,

		//There was an internal error while calculating the graph.
		InternalFailure,

		//Graph construction failed because a motorcycle could not be traced through
		//a certain fenced region. Try again after removing MotorcycleGraph::ProblematicFencedRegion().
		CannotTraceThroughRegion,
	};	

	//Constructs a new empty motorcycle graph. After construction, add motorcycles with AddMotorcycle() or
	//AddMotorcycleWithHistory() and call AdvanceToEnd() to calculate the graph.
	//mesh - The mesh on which to build the graph
	//patchSet - A set of all fenced regions
	//metaSingularities - A set of all meta singularities
	//vertexToMetaSingularityIndex - A map from vertices in the mesh to the corresponding meta singularity index
	//isOriginalEdgeProp - A mesh property that determines if a given edge is original (and not introduced by Catmull-Clark subdivision)
	//lengthMeasure - The measure that determines if patches should be split
	MotorcycleGraph(const HEMesh& mesh, PatchSet& patchSet, const std::vector<MetaSingularity>& metaSingularities,
	const ManifoldnessAwareVertexMap<size_t>& vertexToMetaSingularityIndex, const OpenMesh::EPropHandleT<bool> isOriginalEdgeProp, LengthMeasure lengthMeasure);

	//Returns the index of the halfarc opposite to the provided halfarc.
	size_t OppositeHalfarc(size_t arc) const;

	//Adds a new active motorcycle starting at the given halfedge to the graph. If the halfedge has already been used
	//in the same direction, an exception is thrown. If the halfedge has already been used in the opposite direction
	//and the other motorcycle terminated just after that, the new motorcycle is merged with the other one and made
	//inactive.
	//Returns the index of the added motorcycle (can be (size_t)-1 if the motorcycle could not be added).
	size_t AddMotorcycle(HEMesh::HalfedgeHandle h, bool highPriority);

	//Adds a new pair of motorcycles starting at the given halfedges to the graph and marks them as straight continuations
	//of on another. If the halfedge has already been used in the same direction, an exception is thrown. If the halfedge
	//has already been used in the opposite direction and the other motorcycle terminated just after that, the new motorcycle 
	//is merged with the other one and made inactive. This function is used for artificial valence-2 singularities.
	//Returns the index of the first added motorcycle. The index of the second added motorcycle is the return value + 1.
	size_t AddMotorcyclePair(const std::array<HEMesh::HalfedgeHandle, 2>&, bool highPriority);

	//Adds a new motorcycle with an existing path to the graph. The path must not be occupied by any other motorcycle. The
	//motorcycle is marked active if collideWithBoundary is false and if it does not immediately collide with another motorcycle.
	void AddMotorcycleWithHistory(const std::vector<HEMesh::HalfedgeHandle>& path, bool highPriority, bool collideWithBoundary);

	//Returns if there are active motorcycles in the graph.
	bool HasActiveMotorcycles() const;

	//Continues motorcycle tracing for a single step for every motorcycle.
	void AdvanceAllMotorcycles();	

	//Continues motorcycle tracing until all motorcycles have terminated.
	void AdvanceToEnd();

	//Returns a list of all motorcycles.
	const std::vector<Motorcycle>& Motorcycles() const;

	//Returns a list of all extracted patches. These are only available after a call to ExtractPatches().
	const std::vector<TexturePatch>& Patches() const;

	//Tries to extract the texture patches from the finalized motorcycle graph. The graph must be in
	//the Finished state. The provided statistics object is filled during extraction. The extracted
	//patches are available via Patches()
	void ExtractPatches(ExtractionStatistics& statistics);	

	//Returns a list of all halfarcs of the patch graph.
	const std::vector<HalfArc>& Halfarcs() const { return halfArcs; }

	//Returns the halfedge corresponding to a location on a motorcycle's path
	HEMesh::HalfedgeHandle MotorcycleHalfedge(const LocationOnPath&) const;

	//Returns the current status of this motorcycle graph.
	TracingStatus Status() const { return status; }

	//Returns the index of the region that caused the CannotTraceThroughRegion status.
	size_t ProblematicFencedRegion() const { return problematicFencedRegion; }

	//Deactivates all unnecessary motorcycles, i.e. those that have no collisions on their path.
	//These motorcycles will be ignored during patch extraction.
	void DeactivateUnnecessaryMotorcycles();

	//Splits highly non-rectangular patches into multiple patches in order to reduce parametric distortion.
	void SplitNonRectangularPatches();

	//Loads an external graph into this data structure. No processing functions can be used on this graph.
	//It is only available for visualization.
	void LoadExternal(std::vector<TexturePatch>&& patches, std::vector<HalfArc>&& halfarcs);

	//This event is raised whenever the set of patches of this motorcycle graph changes.
	nse::util::Observer<> PatchesChanged;

private:
	//Describes the result of a routing process (through a fenced region)
	enum RouteResultType
	{
		//Could not find a route through the region
		NoRoute,

		//Could find a valid route where the motorcycle leaves the region on the other side
		RouteThrough,

		//Could find a valid route where the motorcycle collides inside the region
		RouteWithCollision,

		//The motorcycle enters the region at an edge where another motorcycle is just leaving it
		HeadOnCollision,
	};

	//Describes the result of a routing process (through a fenced region)
	struct RouteResult
	{
		//The type of this result
		RouteResultType type;

		//The calculated route through the region
		std::vector<HEMesh::HalfedgeHandle> route;

		//If the route ends with a collision with another motorcycle, this
		//is its motorcycle index
		size_t collisionWith;

		RouteResult() { }
		RouteResult(RouteResultType type) : type(type) { }
		RouteResult(RouteResultType type, std::vector<HEMesh::HalfedgeHandle>&& route) : type(type), route(route) { }
		RouteResult(RouteResultType type, std::vector<HEMesh::HalfedgeHandle>&& route, size_t collisionWith) : type(type), route(route), collisionWith(collisionWith) { }
	};

	//Stopping condition for CirculateForwardUntil that continues as long as the current halfedge
	//has not been visited by a motorcycle. Objects of this struct will save the first motorcycle
	//location that they encounter during circulation.
	struct CirculateUntilMotorcycleCondition
	{
		CirculateUntilMotorcycleCondition(const MotorcycleGraph& graph)
			: graph(graph)
		{ }

		//A list of motorcycles that passed through the target vertex of the query halfedge.
		//This must be set by the caller before circulation in order to be able to find the
		//required motorcycles.
		const std::vector<MotorcycleStop>* motorcyclesAtQueryTargetVertex = nullptr;

		//This member holds the first encountered motorcycle after circulation.
		LocationOnPath result;

		const MotorcycleGraph& graph;

		bool operator()(HEMesh::HalfedgeHandle h);
	};

	//Represents a pair of stopping motorcycles and their validity. They are used as a 
	//fallback strategy for tracing motorcycles through fenced regions. If the pair is
	//spawned at the mesh boundary, one of the two motorcycles can be left away.
	struct StoppingMotorcycles
	{
		//Instantiates a pair of stopping motorcycles
		//mesh - the underlying mesh
		//stopBeforeEdge - the halfedge whose traversal to prevent
		StoppingMotorcycles(const HEMesh& mesh, HEMesh::HalfedgeHandle stopBeforeEdge, OpenMesh::EPropHandleT<bool> isOriginalEdgeProp);

		//Tries to find a route for the motorcycle pair until both are inside a regular region.
		//The motorcycles might be in a regular region right from the beginning.
		//graph - the underlying motorcycle graph
		//tryAlternative - determines if a route with stopBeforeEdge as the motorcycle direction should be attempted
		void TryToRoute(MotorcycleGraph& graph, bool tryAlternative);

		//Realizes the routes for the pair of stopping motorcycles.
		void Realize(MotorcycleGraph& graph, bool highPriority);

		//The edge whose traversal to prevent
		HEMesh::HalfedgeHandle stopBeforeEdge;

		//The two halfedges along which the motorcycle pair starts
		std::array<HEMesh::HalfedgeHandle, 2> newSingularityH;

		//The calculated routes for the motorcycle pair
		std::array<RouteResult, 2> newSingularityHRoutes;

		//Set in the constructor; determines if the motorcycle pair uses non-original edges.
		bool useNonOriginalEdge;

		//Set in TryToRoute(); determines if all motorcycles have a valid route.
		bool hasValidRoute = false;

		//If the pair has an invalid route, determines the index of the motorcycle that failed routing (either 0 or 1)
		int failedStoppingMotorcycle = -1;
	};	

	//Extracts all arcs from the current motorcycle graph and stores them in halfarcs.
	void ExtractArcs();

	//Partitions the path of a given motorcycle into its arcs.
	//iMotorcycle - the motorcycle index to process
	//arcMerges - Output variable. If the motorcycle has continuations at the start or end, stores pairs of halfarc indices that can be merged
	//checkCollisionPoins - If set to true, checks all recorded collisions on the path and creates a new node only if there exist non-deactivated motorcycles.
	//                      Otherwise, all recorded collisions become nodes.
	void ExtractArcsForMotorcycle(size_t iMotorcycle, std::vector<std::pair<size_t, size_t>>& arcMerges, bool checkCollisionPoints = true);

	//Splits the halfarc and its opposite halfarc after the specified segment. Reassigns the faces and replaces the arc with two split arcs in the referenced faces.
	void SplitArcAfterSegment(size_t iArc, HalfArc::PathSegmentIterator it);

	//Determines if a given location on the path of a motorcycle is about to collide with the bounary (either directly or 
	//through a straightly continuing motorcycle with an empty path)
	bool CollidesWithBoundary(const LocationOnPath& location) const;

	//Returns the index of the next unterminated straight continuation of the given motorcycle in the given direction.
	//If there is no straight continuation, motorcycleIdx is returned. This process skips any motorcycle with an empty
	//path.
	//motorcycleIdx - the index of the initial motorcycle to check
	//forward - the direction in which to continue on the motorcycle's path
	size_t FindStraightContinuation(size_t motorcycleIdx, bool forward) const;

	//Returns the index of the next unterminated straight continuation of the given motorcycle in the given direction.
	//If there is no straight continuation, motorcycleIdx is returned. This process skips any motorcycle with an empty
	//path.
	//motorcycleIdx - the index of the initial motorcycle to check
	//forward - the direction in which to continue on the motorcycle's path
	//outForward - Output variable. Determines if the returned motorcycle path is traversed in its forward direction.
	size_t FindStraightContinuation(size_t motorcycleIdx, bool inForward, bool& outForward) const;

	//Calculates the possible exits for a motorcycle entering a patch at a given location.
	void FindPossibleExits(const FencedRegion& patch, const std::pair<HEMesh::HalfedgeHandle, EnteringEdge>& enteringInfo, std::set<HEMesh::HalfedgeHandle>& possibleExits, bool verbose = false);

	//Calculates a route through a fenced region. The last path segment will be inside the region with its target vertex on the fence.
	RouteResult RouteThroughPatch(HEMesh::HalfedgeHandle edge, const FencedRegion& patch, const std::pair<HEMesh::HalfedgeHandle, EnteringEdge>& enteringInfo, bool desperate = false, bool verbose = false);
	RouteResult RouteThroughPatch(HEMesh::HalfedgeHandle edge, const FencedRegion& patch, const float minCosDeviationAngle, bool allowTurnsOnRegularVertices, const std::set<HEMesh::HalfedgeHandle>& possibleExits);

	//Performs all necessary changes for a motorcycle to follow a calculated route, i.e. marks the corresponding vertices and halfedges as visited,
	//records straight continuations for head-on collisions, adds the edges to the motorcycle's path, and terminates the motorcycle if the route
	//ends with a collision.
	void RealizeRoute(size_t cycleIdx, const RouteResult& route);

	//Finds the motorcycle whose path contains the given halfedge (either in forward or backward direction).
	//h - the edge to check
	//motorcyclesAtSourceVertex - a list of motorcycle stops for the from vertex of h
	//outForward - Output variable. Determines if the motorcycle traversed this edge in its forward direction.
	std::vector<MotorcycleGraph::MotorcycleStop>::const_iterator FindMotorcycleAtEdge(HEMesh::HalfedgeHandle h, const std::vector<MotorcycleGraph::MotorcycleStop>& motorcyclesAtSourceVertex, bool& outForward) const;

	//Converts an edge that has been traversed by a motorcycle to its corresponding LocationOnPath.
	//h - the edge to convert
	//motorcyclesAtSourceVertex - a list of motorcycle stops for the from vertex of h
	MotorcycleGraph::LocationOnPath EdgeToMotorcycleLocation(HEMesh::HalfedgeHandle h, const std::vector<MotorcycleGraph::MotorcycleStop>& motorcyclesAtSourceVertex) const;

	//Returns an iterator to the next recorded collision on a motorcycle's path starting at (but not including) the given location.
	std::set<size_t>::iterator NextCrashSite(const LocationOnPath&) const;	

	//Returns the next edge on the boundary of a patch in counter-clockwise direction. The search can skip open boundaries that
	//have no motorcycles.
	//current - the current location on the patch boundary
	//outOrientationChange - Output variable. Determines how many topological turns have been made in order to go to the next edge.
	//skippedMotorcycles - Output variable. A list of deactivated motorcycles that have been skipped.
	LocationOnPath NextEdgeInPatch(const LocationOnPath& current, int& outColorChange, std::vector<size_t>& skippedMotorcycles, bool verbose) const;
	LocationOnPath NextEdgeInPatchInternal(const LocationOnPath& current, int& outColorChange, bool verbose) const;

	//Extracts all patches enclosed by the given set of halfarcs. Separate connected components are emitted as separate patches.
	//PatchOutput: void(std::vector<std::vector<size_t>>&& patchSides, std::set<OpenMesh::FaceHandle>&& faces, std::vector<HEMesh::HalfedgeHandle>&& openBoundary)
	//halfarcs - the halfarcs to process. This list must be ordered by their orientation (and usually resembles a rectangle). There may be halfarcs with a 
	//           zero orientation at the end of the list
	//halfarcOrientations - The integer arc orientations for the final patches.	
	//firstArcOnFirstSide - Index of the first arc on the side with a zero orientation (usually 0 except if there are halfarcs with a zero orientation at the end of halfarcs)	
	//maxOrientation - maximum of halfarcOrientations
	//additionalBoundaries - additional boundary edges that border the final patches (edges that are not part of halfarcs)
	//patchOutput - callback that is called once for every extracted patch
	template <typename PatchOutput>
	void BuildPatchesFromArcs(const std::vector<size_t>& halfarcs, const std::vector<int>& halfarcOrientations, int firstArcOnFirstSide, int maxOrientation, const std::vector<HEMesh::HalfedgeHandle>& additionalBoundaries, PatchOutput&& patchOutput) const;

	//Assigns the given patch index to all halfarcs comprising the given patch.
	void AssignFaceIndexToHalfarcs(size_t patchIdx);	

	//Calculates if the transition from a given motorcycle to a given edge on the boundary is valid. A transition
	//is considered valid if one of the following applies:
	//  - This motorcycle collided with a singularity or starts there and the singularity has an outgoing edge to the boundary
	//  - The motorcycle collided with the boundary directly
	//  - The new motorcycle collided with another motorcycle that immediately collided with the boundary
	//location - the incoming motorcycle to check
	//edgeOnBoundary - the prospectice edge on the boundary
	//motorCyclesAtTransition - a list of all motorcycles that passed through the transitioning point (i.e. the source vertex of edgeOnBoundary)
	//forward - determines if the motorcycle is traversed in its forward direction.
	//outOrientationChange - Output variable. Determines the number of turns that are required to make the transition.
	bool IsValidTransitionFromMotorcycleToBoundary(const LocationOnPath& location, HEMesh::HalfedgeHandle edgeOnBoundary, const std::vector<MotorcycleStop>& motorCyclesAtTransition, bool forward, int& outOrientationChange, bool verbose) const;

	//Marks a given edge as occupied by a given motorcycle
	void OccupyEdge(HEMesh::EdgeHandle, size_t motorcycle);

	//The internal status of this graph
	TracingStatus status = InProcess;

	//The problematic fenced region. Only set if the according state is set.
	size_t problematicFencedRegion;

	//The patches extracted from the graph.
	std::vector<TexturePatch> patches;

	//Determines if there are motorcycles with high priority in the graph
	bool hasHighPriorityMotorcycles = false;

	//Reference to the underlying mesh
	const HEMesh* mesh;

	//Reference to the underlying set of fenced regions
	PatchSet* patchSet;

	//Property that determines if an edge is original (i.e. not introduced by subdivision)
	const OpenMesh::EPropHandleT<bool> isOriginalEdgeProp;

	//Set of all motorcycles, both active and terminated
	std::vector<Motorcycle> motorcycles;

	//Indices of active motorcycles
	std::vector<size_t> activeMotorcycles;

	typedef std::vector<MotorcycleStop> VisitedEntry;
	//Stores for every vertex the motorcycles that passed through it and when
	ManifoldnessAwareVertexMap<VisitedEntry> visitedVertices;

	//Stores all halfedges that have been visited by any motorcycle
	std::set<HEMesh::HalfedgeHandle> visitedHalfEdges;

	//Set of all halfarcs extracted from this graph
	std::vector<HalfArc> halfArcs;

	//Maps the first halfedge of every halfarc to the halfarc index.
	std::map<HEMesh::HalfedgeHandle, size_t> firstHalfedgeToHalfArc;

	//Maps occupied edges to the according motorcycle's id
	std::map<HEMesh::EdgeHandle, size_t> occupiedEdges;

	//reference to the underlying list of meta singularities
	const std::vector<MetaSingularity>& metaSingularities;
	//reference to the underlying structure that maps vertices to meta singularities
	const ManifoldnessAwareVertexMap<size_t>& vertexToMetaSingularityIndex;		

	//The length measure to use in order to determine patch splits
	LengthMeasure lengthMeasure;	
};

extern std::ostream& operator<<(std::ostream& stream, const MotorcycleGraph::ExtractionStatistics& s);