#include "MotorcycleGraph.h"

#include <iostream>
#include <bitset>
#include <random>
#include <algorithm>
#include <cmath>

#include <nsessentials/util/TimedBlock.h>
#include <nsessentials/data/Parallelization.h>

#include "MotorcycleOptionsInPatch.h"
#include "ParametrizationHelper.h"

#include "PatchSplittingProblem.h"

// ------------  Motorcycle  ------------

void MotorcycleGraph::Motorcycle::AddToPath(HEMesh::HalfedgeHandle h, const HEMesh& mesh)
{
	path.push_back(h);
	if (!mesh.is_boundary(h) && !mesh.is_boundary(mesh.opposite_halfedge_handle(h)))
		++nonBoundaryEdges;
}

void MotorcycleGraph::Motorcycle::RemoveLastPathSegment(const HEMesh& mesh)
{
	if (!mesh.is_boundary(path.back()) && !mesh.is_boundary(mesh.opposite_halfedge_handle(path.back())))
		--nonBoundaryEdges;
	path.pop_back();
}

// ------------  MotorcycleStop  ------------

MotorcycleGraph::MotorcycleStop::MotorcycleStop(size_t motorcycleId, size_t locationInMotorcyclePath)
	: motorcycle(motorcycleId), locationInMotorcyclePath(locationInMotorcyclePath)
{ }

// ------------  LocationOnPath  ------------

MotorcycleGraph::LocationOnPath::LocationOnPath() { }
MotorcycleGraph::LocationOnPath::LocationOnPath(size_t motorcycle, size_t pathSegment, bool goForward)
	: motorcycle(motorcycle), pathSegment(pathSegment), goForward(goForward)
{ }

bool MotorcycleGraph::LocationOnPath::operator==(const LocationOnPath& other) const
{
	return motorcycle == other.motorcycle && pathSegment == other.pathSegment && goForward == other.goForward;
}

// ------------  ExtractionStatistics  ------------

void MotorcycleGraph::ExtractionStatistics::Clear()
{
	totalFaces = 0;
	openPatches = 0;
	maxPatchDegree = 0;
	patchCountPerNumberOfCorners.clear();
	patchSizesPerNumberOfCorners.clear();
	patchSizes.clear();
}

std::ostream& operator<<(std::ostream& stream, const MotorcycleGraph::ExtractionStatistics& s)
{
	std::cout << s.openPatches << " open patches." << std::endl;
	for (int corners = 0; corners < s.patchCountPerNumberOfCorners.size(); ++corners)
		std::cout << corners << " corners: " << (100.0f * s.patchSizesPerNumberOfCorners[corners] / s.totalFaces) << "% (" << s.patchCountPerNumberOfCorners[corners] << " patches)" << std::endl;

	size_t steps[] = { 0, 1, 10, 50, 100, 500, 1000, 5000, 10000 };
	size_t stepCount = sizeof(steps) / sizeof(steps[0]);
	for (int i = 0; i <= stepCount; ++i)
	{
		auto lower = (i == 0 ? s.patchSizes.begin() : std::upper_bound(s.patchSizes.begin(), s.patchSizes.end(), steps[i - 1]));
		auto upper = (i == stepCount ? s.patchSizes.end() : std::upper_bound(s.patchSizes.begin(), s.patchSizes.end(), steps[i]));
		size_t faces = 0;
		for (auto it = lower; it != upper; ++it)
			faces += *it;
		if (i < stepCount)
			std::cout << "Faces in patches with <= " << steps[i];
		else
			std::cout << "Faces in patches with > " << steps[i - 1];
		std::cout << " faces: " << 100.0f * (float)faces / s.totalFaces << "% (" << (upper - lower) << " patches)" << std::endl;
	}

	return stream;
}

// ------------  HalfArc  ------------

MotorcycleGraph::HalfArc::Segment::Segment(const LocationOnPath& location, size_t length)
	: location(location), length(length)
{ }

MotorcycleGraph::HalfArc::HalfArc(const PathSegmentIterator& first, const PathSegmentIterator& last)
{
	segments.emplace_back(*first, 0);
	if (first == last)
		return;

	auto it = first;
	++it;
	int length = 1;
	//Find the segments of equal motorcycles that make up the given range
	while (it != last)
	{
		auto newLocation = *it;
		if (newLocation.motorcycle != segments.back().location.motorcycle)
		{
			segments.back().length = length;
			segments.emplace_back(newLocation, 0);
			length = 1;
		}
		else
			++length;
		++it;
	}
	segments.back().length = length;
}

size_t MotorcycleGraph::HalfArc::Length() const
{
	size_t l = 0;
	for (auto& segment : segments)
		l += segment.length;
	return l;
}

void MotorcycleGraph::HalfArc::MergeWith(HalfArc&& merged)
{
	segments.insert(segments.end(), merged.segments.begin(), merged.segments.end());
	merged.segments.clear();
}

MotorcycleGraph::LocationOnPath MotorcycleGraph::HalfArc::LastPathSegment() const
{
	auto& segment = segments.back();
	if (segment.location.goForward)
		return LocationOnPath(segment.location.motorcycle, segment.location.pathSegment + segment.length - 1, segment.location.goForward);
	else
		return LocationOnPath(segment.location.motorcycle, segment.location.pathSegment - segment.length + 1, segment.location.goForward);
}

void MotorcycleGraph::HalfArc::SetStartInclusive(const PathSegmentIterator& start)
{
	assert(start.arc == this);
	auto newStartLocation = *start;
	segments.erase(segments.begin(), segments.begin() + start.iMotorcycle);
	segments.front().location = newStartLocation;
	segments.front().length -= start.iPathSegmentOnMotorcycle;
}

void MotorcycleGraph::HalfArc::SetEndExclusive(const PathSegmentIterator& end)
{
	assert(end.arc == this);
	if (end.iPathSegmentOnMotorcycle == 0)
	{
		segments.erase(segments.begin() + end.iMotorcycle, segments.end());
		return;
	}

	segments.erase(segments.begin() + end.iMotorcycle + 1, segments.end());
	segments.back().length = end.iPathSegmentOnMotorcycle;
}

MotorcycleGraph::HalfArc::PathSegmentIterator MotorcycleGraph::HalfArc::begin() const { return PathSegmentIterator(0, 0, this); }
MotorcycleGraph::HalfArc::PathSegmentIterator MotorcycleGraph::HalfArc::end() const { return PathSegmentIterator(segments.size(), 0, this); }

// ------------  HalfArc::PathSegmentIterator  ------------

MotorcycleGraph::HalfArc::PathSegmentIterator::PathSegmentIterator(int iMotorcycle, int iPathSegment, const HalfArc* arc)
	: iMotorcycle(iMotorcycle), iPathSegmentOnMotorcycle(iPathSegment), arc(arc)
{ }

bool MotorcycleGraph::HalfArc::PathSegmentIterator::operator==(const PathSegmentIterator& other) const
{
	return iMotorcycle == other.iMotorcycle && iPathSegmentOnMotorcycle == other.iPathSegmentOnMotorcycle;
}

bool MotorcycleGraph::HalfArc::PathSegmentIterator::operator!=(const PathSegmentIterator& other) const 
{
	return !(*this == other);
}

MotorcycleGraph::LocationOnPath MotorcycleGraph::HalfArc::PathSegmentIterator::operator*() const
{
	auto& segment = arc->segments[iMotorcycle];
	if (segment.location.goForward)
		return LocationOnPath(segment.location.motorcycle, segment.location.pathSegment + iPathSegmentOnMotorcycle, segment.location.goForward);
	else
		return LocationOnPath(segment.location.motorcycle, segment.location.pathSegment - iPathSegmentOnMotorcycle, segment.location.goForward);
}

MotorcycleGraph::HalfArc::PathSegmentIterator& MotorcycleGraph::HalfArc::PathSegmentIterator::operator++()
{
	++iPathSegmentOnMotorcycle;
	if (iPathSegmentOnMotorcycle >= arc->segments[iMotorcycle].length)
	{
		iPathSegmentOnMotorcycle = 0;
		++iMotorcycle;
	}
	return *this;
}

// ------------  MotorcycleGraph  ------------

MotorcycleGraph::MotorcycleGraph(const HEMesh& mesh, PatchSet& patchSet, const std::vector<MetaSingularity>& metaSingularities,
	const ManifoldnessAwareVertexMap<size_t>& vertexToMetaSingularityIndex, const OpenMesh::EPropHandleT<bool> isOriginalEdgeProp, 
	LengthMeasure lengthMeasure)
	: mesh(&mesh), patchSet(&patchSet), metaSingularities(metaSingularities), vertexToMetaSingularityIndex(vertexToMetaSingularityIndex),
	visitedVertices(mesh), isOriginalEdgeProp(isOriginalEdgeProp), lengthMeasure(lengthMeasure)
{ }

size_t MotorcycleGraph::OppositeHalfarc(size_t arc) const
{
	return arc ^ 1;
}

size_t MotorcycleGraph::AddMotorcycle(HEMesh::HalfedgeHandle h, bool highPriority)
{
	auto& visitedAtSource = visitedVertices.AccessOrCreateAtToVertex(mesh->opposite_halfedge_handle(h));
	int mergeWith = -1;
	if (visitedHalfEdges.find(h) != visitedHalfEdges.end())
	{
		//The halfedge has already been visited
		auto oppH = mesh->opposite_halfedge_handle(h);
		for (auto& stop : visitedAtSource)
		{
			//if the other motorcycle is terminated here, we can just merge with it
			auto& motorcycle = motorcycles[stop.motorcycle];
			if (motorcycle.terminated && motorcycle.straightContinuationEnd == (size_t)-1 && motorcycle.Path().back() == oppH)
				mergeWith = stop.motorcycle;
		}
		if (mergeWith == -1)
		{
			std::cout << "The halfedge ";
			PrintFullHalfedge(std::cout, h, *mesh);
			std::cout << " has already been used." << std::endl;
			throw std::runtime_error("Error adding a motorcycle");
		}
	}	
	motorcycles.emplace_back(h, highPriority);
	if (highPriority)
		hasHighPriorityMotorcycles = true;

	if (mergeWith != -1)
	{
		motorcycles[motorcycles.size() - 1].straightContinuationEnd = mergeWith;
		motorcycles[mergeWith].straightContinuationEnd = motorcycles.size() - 1;
		motorcycles[motorcycles.size() - 1].terminated = true;
	}
	else
		activeMotorcycles.push_back(motorcycles.size() - 1);		

	//Record the collision for all motorcycles passing through the start vertex
	for (auto& stop : visitedAtSource)
		motorcycles[stop.motorcycle].crashesOnPath.insert(stop.locationInMotorcyclePath);
	visitedAtSource.emplace_back(motorcycles.size() - 1, 0);
	visitedHalfEdges.insert(h);

	return motorcycles.size() - 1;
}

size_t MotorcycleGraph::AddMotorcyclePair(const std::array<HEMesh::HalfedgeHandle, 2>& h, bool highPriority)
{
#if __DEBUG
	assert(mesh->from_vertex_handle(h[0]) == mesh->from_vertex_handle(h[1]));
#endif
	
	auto& visitedAtSource = visitedVertices.AccessOrCreateAtToVertex(mesh->opposite_halfedge_handle(h[0]));
	int mergeH0With = -1, mergeH1With = -1;
	if (visitedHalfEdges.find(h[0]) != visitedHalfEdges.end())
	{
		//The halfedge has already been visited
		auto oppH0 = mesh->opposite_halfedge_handle(h[0]);
		for (auto& stop : visitedAtSource)
		{
			//if the other motorcycle is terminated here, we can just merge with it
			auto& motorcycle = motorcycles[stop.motorcycle];
			if (motorcycle.terminated && motorcycle.straightContinuationEnd == (size_t)-1 && motorcycle.Path().back() == oppH0)
				mergeH0With = stop.motorcycle;
		}
		if(mergeH0With == -1)
			std::cout << "The first halfedge at " << mesh->point(mesh->from_vertex_handle(h[0])) << " of the motorcycle pair has already been used." << std::endl;
	}
	if (visitedHalfEdges.find(h[1]) != visitedHalfEdges.end())
	{
		//The halfedge has already been visited
		auto oppH1 = mesh->opposite_halfedge_handle(h[1]);
		for (auto& stop : visitedAtSource)
		{
			//if the other motorcycle is terminated here, we can just merge with it
			auto& motorcycle = motorcycles[stop.motorcycle];
			if (motorcycle.terminated && motorcycle.straightContinuationEnd == (size_t)-1 && motorcycle.Path().back() == oppH1)
				mergeH1With = stop.motorcycle;
		}
		if (mergeH1With == -1)
			std::cout << "The second halfedge at " << mesh->point(mesh->from_vertex_handle(h[1])) << " of the motorcycle pair has already been used." << std::endl;
	}	

	
	for (auto& stop : visitedAtSource)
		motorcycles[stop.motorcycle].crashesOnPath.insert(stop.locationInMotorcyclePath);

	for (int i = 0; i < 2; ++i)
	{
		motorcycles.emplace_back(h[i], highPriority);
		visitedAtSource.emplace_back(motorcycles.size() - 1, 0);
		visitedHalfEdges.insert(h[i]);
	}
	if (highPriority)
		hasHighPriorityMotorcycles = true;

	motorcycles[motorcycles.size() - 2].straightContinuationStart = motorcycles.size() - 1;
	motorcycles[motorcycles.size() - 1].straightContinuationStart = motorcycles.size() - 2;

	if (mergeH0With != -1)
	{
		motorcycles[motorcycles.size() - 2].straightContinuationEnd = mergeH0With;
		motorcycles[mergeH0With].straightContinuationEnd = motorcycles.size() - 2;
		motorcycles[motorcycles.size() - 2].terminated = true;
	}
	else
		activeMotorcycles.push_back(motorcycles.size() - 2);

	if (mergeH1With != -1)
	{
		motorcycles[motorcycles.size() - 1].straightContinuationEnd = mergeH1With;
		motorcycles[mergeH1With].straightContinuationEnd = motorcycles.size() - 1;
		motorcycles[motorcycles.size() - 1].terminated = true;
	}
	else
		activeMotorcycles.push_back(motorcycles.size() - 1);

	return motorcycles.size() - 2;
}


void MotorcycleGraph::AddMotorcycleWithHistory(const std::vector<HEMesh::HalfedgeHandle>& path, bool highPriority, bool collideWithBoundary)
{
	if (path.size() == 0)
		return;

	if (visitedHalfEdges.find(path[0]) != visitedHalfEdges.end())
		return; //already used the beginning of this path	

	motorcycles.emplace_back();	
	motorcycles.back().highPriority = highPriority;
	if (highPriority)
		hasHighPriorityMotorcycles = true;

	auto& sourceVisited = visitedVertices.AccessOrCreateAtToVertex(mesh->opposite_halfedge_handle(path[0]));
	for (auto& stop : sourceVisited)
		motorcycles[stop.motorcycle].crashesOnPath.insert(stop.locationInMotorcyclePath);
	sourceVisited.emplace_back(motorcycles.size() - 1, motorcycles.back().Path().size());

	bool active = !collideWithBoundary;
	for(auto it = path.begin(); it < path.end(); ++it)
	{
		//The last edge is only used if the motorcycle will immediately terminate. Otherwise,
		//it just sits at the edge.
		bool useEdge = it != path.end() - 1 || !active;

		auto inserted = visitedHalfEdges.insert(*it);
		if (!inserted.second)
		{
			//the motorcycle already collided with another motorcycle
			std::cout << "This should not happen." << std::endl;
		}

		if (useEdge)
		{
			OccupyEdge(mesh->edge_handle(*it), motorcycles.size() - 1);
			visitedHalfEdges.insert(mesh->opposite_halfedge_handle(*it));
			motorcycles.back().AddToPath(*it, *mesh);

			auto& targetVisited = visitedVertices.AccessOrCreateAtToVertex(*it);
			for (auto& stop : targetVisited)
				motorcycles[stop.motorcycle].crashesOnPath.insert(stop.locationInMotorcyclePath);
			targetVisited.emplace_back(motorcycles.size() - 1, motorcycles.back().Path().size());
			if (targetVisited.size() > 1)
			{
				active = false;
				//check for head-on collision
				for (auto& stop : targetVisited)
				{
					if (motorcycles[stop.motorcycle].currentPosition == mesh->opposite_halfedge_handle(*it))
					{
						//this is a head-on collision
						motorcycles.back().straightContinuationEnd = stop.motorcycle;
						motorcycles[stop.motorcycle].straightContinuationEnd = motorcycles.size() - 1;
						motorcycles[stop.motorcycle].terminated = true;
						break;
					}
				}				
			
				break;
			}
		}
		else
			motorcycles.back().currentPosition = *it;
	}	

	if (active)
		activeMotorcycles.push_back(motorcycles.size() - 1);
	else
		motorcycles.back().terminated = true;

	if(collideWithBoundary)
		motorcycles.back().collideWithPaddedBoundary = true;	
}

bool MotorcycleGraph::HasActiveMotorcycles() const
{
	return !activeMotorcycles.empty();
}

const std::vector<TexturePatch>& MotorcycleGraph::Patches() const
{
	return patches;
}

void MotorcycleGraph::FindPossibleExits(const FencedRegion& patch, const std::pair<HEMesh::HalfedgeHandle, EnteringEdge>& enteringInfo, std::set<HEMesh::HalfedgeHandle>& possibleExits, bool verbose)
{
	auto& loop = patch.Loops()[0];

	//Walk the fence and record the edges that result in a zero-turn traversal.
	//Motorcycles that pass through the region are considered boundaries. The edge through which
	//the motorcycles enters is marked with -2 turns.

	//The orientation of the current edge on the fence. Can be -1 if we are
	//currently not on the boundary but on a motorcycle.
	int currentRegionEdgeOrientation = enteringInfo.second.enteringViaEdgeColor;

	//The current motorcycle that we are travelling on. currentMotorcycle.motorcycle can be
	//(size_t)-1 if we are currently not on a motorcycle but on the fence.
	LocationOnPath currentMotorcycle; currentMotorcycle.motorcycle = (size_t)-1;

	//Records the number of turns needed to traverse through the fenced region leaving at the
	//current position.
	int currentTurns = -2;
	HEMesh::HalfedgeHandle currentH;
	
	if (!mesh->is_boundary(enteringInfo.first))
		//find the halfedge that is to the left of the entry point
		currentH = mesh->prev_halfedge_handle(enteringInfo.first);
	else
	{
		//if we are entering through a boundary edge, the halfedge to the left
		//of the entry point does not exist. Instead, we use the inverse of the
		//entering edge, which has a turn count of -3 instead of -2.
		currentH = mesh->opposite_halfedge_handle(enteringInfo.first);
		currentRegionEdgeOrientation = (currentRegionEdgeOrientation + loop.Degree() - 1) % loop.Degree();
		currentTurns = -3;		
	}

	{
		//Check if we are also currently on a motorcycle
		VisitedEntry* visited;
		if(visitedVertices.TryAccessAtToVertex(mesh->opposite_halfedge_handle(currentH), visited))
			currentMotorcycle = EdgeToMotorcycleLocation(currentH, *visited);
	}

	auto startVertex = mesh->from_vertex_handle(enteringInfo.first);

	if (verbose)
	{
		std::cout << "Entering patch at ";
		PrintFullHalfedge(std::cout, enteringInfo.first, *mesh);
		std::cout << ", initial halfedge: ";
		PrintFullHalfedge(std::cout, currentH, *mesh);
		std::cout << ", currentQuadPatchEdgeColor: " << currentTurns << std::endl;
	}

	int outerIterations = 0;
	bool hIsInitialEdge = true;
	//Move along the boundary (or motorcycles)
	do
	{
		if (outerIterations++ > 1000)
		{
			std::cout << "Problem with outer iteration when entering patch at " << mesh->point(startVertex) << std::endl;
			throw std::runtime_error("Unknown error while finding appropriate exit points");
		}

		bool currentIsPaddedBorder = mesh->is_boundary(mesh->opposite_halfedge_handle(currentH)) && currentRegionEdgeOrientation != -1;
		//there will never be a motorcycle on the padded boundary
		//DEBUG:
		if (currentIsPaddedBorder && currentMotorcycle.motorcycle != (size_t)-1)
			throw std::runtime_error("There was a motorcycle on the padded boundary.");

		int iterations = 0;
		int nextRegionEdgeOrientation = -1;
		LocationOnPath nextMotorcycle; nextMotorcycle.motorcycle = (size_t)-1;

		VisitedEntry* visited;
		auto isVisited = visitedVertices.TryAccessAtToVertex(currentH, visited);

		int hopsFromCurrentToNext = 0;
		int boundaryOrientationChange = 0;

		//Find the edge onto which to continue
		auto stopCondition = [&](HEMesh::HalfedgeHandle h)
		{
			if (iterations++ > 10)
			{
				std::cout << "Problem with inner iteration when entering patch at " << mesh->point(mesh->from_vertex_handle(enteringInfo.first)) << std::endl;
				throw std::runtime_error("Unknown error while finding appropriate exit points");
			}

			if (verbose)
			{
				std::cout << "Checking halfedge ";
				PrintFullHalfedge(std::cout, h, *mesh);
				std::cout << std::endl;
			}

			++hopsFromCurrentToNext;
			boundaryOrientationChange = 0;

			//check if we continue on the fence
			bool currentHOnPaddedBoundary = false;
			auto edgeOnFence = patch.Loops()[0].FindEdge(h);
			if (edgeOnFence != (size_t)-1)
			{
				nextRegionEdgeOrientation = patch.Loops()[0].Edges()[edgeOnFence].orientation;
				if (mesh->is_boundary(mesh->opposite_halfedge_handle(h)))
				{
					currentHOnPaddedBoundary = true;
					if (currentMotorcycle.motorcycle != (size_t)-1)
					{
						//check if this is a valid transition
						if (!IsValidTransitionFromMotorcycleToBoundary(currentMotorcycle, currentH, *visited, true, boundaryOrientationChange, verbose) && !patch.IsToVertexOnBoundary(mesh->opposite_halfedge_handle(h), true))
						{
							nextRegionEdgeOrientation = -1;
							currentHOnPaddedBoundary = false;
						}
						else if (boundaryOrientationChange == 0)
							boundaryOrientationChange = 2;
					}
				}
			}

			//check if we continue on a motorcycle			
			if (isVisited)
			{
				nextMotorcycle = EdgeToMotorcycleLocation(h, *visited);
				if (nextMotorcycle.motorcycle != (size_t)-1)
				{
					if (currentIsPaddedBorder)
					{
						//check if this is a valid transition
						auto copy = nextMotorcycle;
						copy.goForward = !copy.goForward;
						if (!IsValidTransitionFromMotorcycleToBoundary(copy, mesh->opposite_halfedge_handle(currentH), *visited, false, boundaryOrientationChange, verbose) && !patch.IsToVertexOnBoundary(mesh->opposite_halfedge_handle(h), true))
							nextMotorcycle.motorcycle = (size_t)-1;
						else if (boundaryOrientationChange == 0)
							boundaryOrientationChange = 2;
					}
					if (currentMotorcycle.motorcycle == nextMotorcycle.motorcycle && currentMotorcycle.goForward != nextMotorcycle.goForward)
						nextMotorcycle.motorcycle = (size_t)-1;
				}
				if (nextMotorcycle.motorcycle != (size_t)-1 && currentHOnPaddedBoundary)
					nextRegionEdgeOrientation = -1;
			}

			return nextRegionEdgeOrientation != -1 || nextMotorcycle.motorcycle != (size_t)-1;
		};

		HEMesh::HalfedgeHandle nextH;
		if (stopCondition(mesh->opposite_halfedge_handle(currentH)))
		{
			nextH = mesh->opposite_halfedge_handle(currentH);
			hopsFromCurrentToNext = 0;
		}
		else
		{
			hopsFromCurrentToNext = 0;
			nextH = currentH;
			CirculateForwardUntil<false>(nextH, *mesh, stopCondition);
		}

		bool straightContinuation = false;

		//Find if the current transition is a straight continuation
		if (!hIsInitialEdge)
		{
			if (currentMotorcycle.motorcycle != (size_t)-1 && nextMotorcycle.motorcycle == (size_t)-1)
			{
				if (!motorcycles[currentMotorcycle.motorcycle].terminated && motorcycles[currentMotorcycle.motorcycle].currentPosition == nextH)
					//Although the motorcycle did not yet go there, this is in fact a straight continuation. If we don't
					//do this, this location will be counted as a transition from motorcycle to edge
					straightContinuation = true;
			}
			if (nextMotorcycle.motorcycle != (size_t)-1 && currentMotorcycle.motorcycle == (size_t)-1)
			{
				if (!motorcycles[nextMotorcycle.motorcycle].terminated && mesh->opposite_halfedge_handle(motorcycles[nextMotorcycle.motorcycle].currentPosition) == currentH)
					//Same as above, just in the other direction
					straightContinuation = true;
			}
		}
		else
		{
			hIsInitialEdge = false;
		}

		if (currentIsPaddedBorder && currentMotorcycle.motorcycle == (size_t)-1 && nextMotorcycle.motorcycle != (size_t)-1)
		{
			bool collisionWithPaddedBoundary =
				(motorcycles[nextMotorcycle.motorcycle].startInPaddedBoundary && mesh->edge_handle(nextH) == mesh->edge_handle(motorcycles[nextMotorcycle.motorcycle].Path()[0]))
				|| (motorcycles[nextMotorcycle.motorcycle].collideWithPaddedBoundary && mesh->edge_handle(nextH) == mesh->edge_handle(motorcycles[nextMotorcycle.motorcycle].Path()[motorcycles[nextMotorcycle.motorcycle].Path().size() - 1]));
			if (!collisionWithPaddedBoundary)
				++currentTurns; //two corners collapsed into one (the other one will be counted by the usual check)
		}

		if (mesh->is_boundary(mesh->opposite_halfedge_handle(nextH)) && nextRegionEdgeOrientation != -1 && currentMotorcycle.motorcycle != (size_t)-1 && nextMotorcycle.motorcycle == (size_t)-1)
		{
			bool collisionWithPaddedBoundary =
				(motorcycles[currentMotorcycle.motorcycle].startInPaddedBoundary && mesh->edge_handle(currentH) == mesh->edge_handle(motorcycles[currentMotorcycle.motorcycle].Path()[0]))
				|| (motorcycles[currentMotorcycle.motorcycle].collideWithPaddedBoundary && mesh->edge_handle(currentH) == mesh->edge_handle(motorcycles[currentMotorcycle.motorcycle].Path()[motorcycles[currentMotorcycle.motorcycle].Path().size() - 1]));
			if (!collisionWithPaddedBoundary)
				++currentTurns; //two corners collapsed into one (the other one will be counted by the usual check)
		}

		if ((nextMotorcycle.motorcycle == currentMotorcycle.motorcycle && nextMotorcycle.motorcycle != (size_t)-1) || (nextRegionEdgeOrientation == currentRegionEdgeOrientation && nextRegionEdgeOrientation != -1))
			straightContinuation = true;

		auto findStraightContinuation = [&](size_t motorcycle1, size_t motorcycle2)
		{
			if (motorcycle1 != -1)
			{
				auto continuation = FindStraightContinuation(motorcycle1, true);
				if (continuation != (size_t)-1 && motorcycle2 == continuation)
					straightContinuation = true;
				continuation = FindStraightContinuation(motorcycle1, false);
				if (continuation != (size_t)-1 && motorcycle2 == continuation)
					straightContinuation = true;
				else if (continuation != (size_t)-1 && motorcycles[continuation].Path().size() == 0 && !motorcycles[continuation].terminated)
				{
					//the continuation motorcycle might still be at its start with no path segments at all
					auto continuationEdge = mesh->edge_handle(motorcycles[continuation].currentPosition);
					if (mesh->edge_handle(nextH) == continuationEdge)
						straightContinuation = true;
					if (mesh->edge_handle(currentH) == continuationEdge)
						straightContinuation = true;
				}
			}
		};

		if (!straightContinuation)
			findStraightContinuation(currentMotorcycle.motorcycle, nextMotorcycle.motorcycle);
		if (!straightContinuation)
			findStraightContinuation(nextMotorcycle.motorcycle, currentMotorcycle.motorcycle);

		//Update the current turn count for the current transition
		if (boundaryOrientationChange != 0)
			currentTurns += boundaryOrientationChange;
		else if (straightContinuation)
			; //straight continuation
		else if (currentRegionEdgeOrientation != -1 && nextRegionEdgeOrientation == -1)
			++currentTurns; //switch from fence to motorcycle
		else if (currentRegionEdgeOrientation == -1 && nextRegionEdgeOrientation != -1)
			++currentTurns; //switch from motorcycle to fence
		else if (currentRegionEdgeOrientation != -1 && nextRegionEdgeOrientation != -1 && currentRegionEdgeOrientation != nextRegionEdgeOrientation)
			currentTurns += 2 - hopsFromCurrentToNext;
		else if (currentMotorcycle.motorcycle != nextMotorcycle.motorcycle)
			++currentTurns; //switch from one motorcycle to another

		if (verbose)
			std::cout << "Continued to vertex " << mesh->to_vertex_handle(nextH) << " with regionEdgeColor = " << nextRegionEdgeOrientation << ", motorcycle = " << nextMotorcycle.motorcycle << ", quadPatchColor = " << currentTurns << std::endl;

		if (currentTurns == 0)
		{
			//The current edge is a valid exit point
			possibleExits.insert(nextH);
			int checkMotorcycle = -1;
			HEMesh::HalfedgeHandle edgeOnMotorcycle;
			if (currentRegionEdgeOrientation != -1 && mesh->is_boundary(mesh->opposite_halfedge_handle(currentH)) && nextMotorcycle.motorcycle != (size_t)-1 && currentMotorcycle.motorcycle == (size_t)-1)
			{
				checkMotorcycle = nextMotorcycle.motorcycle;
				edgeOnMotorcycle = nextH;
			}
			if (nextRegionEdgeOrientation != -1 && mesh->is_boundary(mesh->opposite_halfedge_handle(nextH)) && currentMotorcycle.motorcycle != (size_t)-1 && nextMotorcycle.motorcycle == (size_t)-1)
			{
				checkMotorcycle = currentMotorcycle.motorcycle;
				edgeOnMotorcycle = nextH;
			}
			if (checkMotorcycle != -1)
			{
				//Also add the adjacent edge as target (even if it is technically not reachable). This is done
				//because the tracer checks both sides when colliding with a motorcycle
				for (auto& stop : visitedVertices.AccessAtToVertex(currentH))
				{
					if (stop.motorcycle == checkMotorcycle)
					{
						if (stop.locationInMotorcyclePath > 0 && stop.locationInMotorcyclePath < motorcycles[checkMotorcycle].Path().size() - 1)
						{
							if (motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath - 1] == edgeOnMotorcycle || motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath] == edgeOnMotorcycle)
							{
								possibleExits.insert(motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath - 1]);
								possibleExits.insert(motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath]);
							}
							else
							{
								//use the inverse
								possibleExits.insert(mesh->opposite_halfedge_handle(motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath - 1]));
								possibleExits.insert(mesh->opposite_halfedge_handle(motorcycles[checkMotorcycle].Path()[stop.locationInMotorcyclePath]));
							}
						}
					}
				}
			}
		}

		currentRegionEdgeOrientation = nextRegionEdgeOrientation;
		currentMotorcycle = nextMotorcycle;
		currentH = nextH;

	} while (startVertex != mesh->to_vertex_handle(currentH));
}

MotorcycleGraph::RouteResult MotorcycleGraph::RouteThroughPatch(HEMesh::HalfedgeHandle edge, const FencedRegion& patch, const std::pair<HEMesh::HalfedgeHandle, EnteringEdge>& enteringInfo, bool desperate, bool verbose)
{	
	const std::pair<float, bool> attempts[] =
	{
		std::pair<float, bool>(0.3f, false),
		std::pair<float, bool>(-0.3f, false),
		std::pair<float, bool>(-0.5f, true),
		std::pair<float, bool>(-2.0f, true),
	};

	if (patch.Loops().size() > 1)
	{
		std::cout << "Multi-loop patches not supported." << std::endl;
		throw 1;
	}

	auto& loop = patch.Loops()[0];

	//Find an appropriate point of exit
	std::set<HEMesh::HalfedgeHandle> possibleExits; //those are the edges perpendicular to the motorcycle's exit direction
	FindPossibleExits(patch, enteringInfo, possibleExits, verbose);	

	if (possibleExits.empty())
		return RouteResult(NoRoute);

	for (int i = 0; i < (desperate ? 4 : 3); ++i)
	{
		auto result = RouteThroughPatch(edge, patch, attempts[i].first, attempts[i].second, possibleExits);
		if (result.type != NoRoute)
			return result;
	}
	return RouteResult(NoRoute);
}

MotorcycleGraph::RouteResult MotorcycleGraph::RouteThroughPatch(HEMesh::HalfedgeHandle edge, const FencedRegion& patch, const float minCosDeviationAngle, bool allowTurnsOnRegularVertices, const std::set<HEMesh::HalfedgeHandle>& possibleExits)
{
	//Represents a valid routing option
	struct Option
	{
		//cost of the route
		float cost;		

		//route type
		RouteResultType type;

		//index of the corresponding motorcycle within options
		size_t cycleIdx;

		//if the route ends with a collision, this is the index of the motorcycle
		size_t collisionWith;

		//The first edge outside of the fenced region
		HEMesh::HalfedgeHandle edgeOutside;

		Option(float cost, RouteResultType type, size_t cycleIdx, size_t collisionWith, HEMesh::HalfedgeHandle edgeOutside)
			: cost(cost), type(type), cycleIdx(cycleIdx), collisionWith(collisionWith), edgeOutside(edgeOutside)
		{ }

		bool operator<(const Option& other) const
		{
			if (cost != other.cost)
				return cost < other.cost;
			if (type != other.type)
				return type < other.type;
			if (cycleIdx != other.cycleIdx)
				return cycleIdx < other.cycleIdx;
			return false;
		}
	};

	auto& loop = patch.Loops()[0];
	std::set<Option> optionsWithCost;

	//Perform Dijkstra through the fenced region
	MotorcycleOptionsInPatch<> options;
	options.AddMotorcycle(edge, *mesh, isOriginalEdgeProp);

	auto leavePatchAction = [&](size_t localCycleIdx, HEMesh::HalfedgeHandle edgeOutside)
	{
		//whenever a candidate leaves the region, check if it was through a valid exit point

		bool validExit = false;
		
		if (edgeOutside.is_valid())
		{
			//motorcycles continues on the mesh
			auto prevLoopEdge = (mesh->is_boundary(edgeOutside) ? (size_t)-1 : loop.FindEdge(mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(edgeOutside))));
			if (prevLoopEdge != (size_t)-1 && possibleExits.find(loop.Edges()[prevLoopEdge].edge) != possibleExits.end())
				validExit = true;

			auto h = mesh->opposite_halfedge_handle(edgeOutside);
			auto nextLoopEdge = (mesh->is_boundary(h) ? (size_t)-1 : loop.FindEdge(mesh->opposite_halfedge_handle(mesh->next_halfedge_handle(h))));
			if (nextLoopEdge != (size_t)-1 && possibleExits.find(loop.Edges()[nextLoopEdge].edge) != possibleExits.end())
				validExit = true;

			if (validExit)
				optionsWithCost.emplace(options.Motorcycles()[localCycleIdx].CostSoFar(), RouteThrough, localCycleIdx, -1, edgeOutside);
		}
		else
		{
			//motorcycles leaves through mesh boundary
			auto h = options.Motorcycles()[localCycleIdx].Edges().back();
			auto leftOutgoing = h;
			CirculateForwardUntil<true>(leftOutgoing, *mesh, [](HEMesh::HalfedgeHandle) { return false; });
			auto rightIncoming = h;
			CirculateBackwardUntil<true>(rightIncoming, *mesh, [](HEMesh::HalfedgeHandle) { return false; });
			rightIncoming = mesh->opposite_halfedge_handle(rightIncoming);

			if(possibleExits.find(leftOutgoing) != possibleExits.end() && possibleExits.find(rightIncoming) != possibleExits.end())
			{
				optionsWithCost.emplace(options.Motorcycles()[localCycleIdx].CostSoFar(), RouteWithCollision, localCycleIdx, -1, edgeOutside);
			}
		}		
	};

	auto passThroughVertexAction = [&](size_t localCycleIdx, HEMesh::VertexHandle v)
	{
		//whenever a candidate passes through an interior vertex, check if there is a valid collision
		//with a motorcycle

		auto lastHalfedge = options.Motorcycles()[localCycleIdx].Edges().back();

		VisitedEntry* visited;
		bool isVisited = visitedVertices.TryAccessAtToVertex(lastHalfedge, visited);
		if (!isVisited)
			return true; //no other motorcycle has been here before

		//We collided with another motorcycle. Check if this is a valid collision.			
		
		//check for head-on collisions
		for (auto& stop : *visited)
		{
			auto enteringIt = patchSet->enteringEdgeToPatch.find(motorcycles[stop.motorcycle].currentPosition);
			if (stop.locationInMotorcyclePath == motorcycles[stop.motorcycle].Path().size() //if the motorcycle did not already continue ...
				&& !motorcycles[stop.motorcycle].terminated // ... and is still active ...
				&& enteringIt != patchSet->enteringEdgeToPatch.end()) // ... and is about to enter the patch				
			{
				//possible head-on collision
				bool validCollision = false;
				if (!mesh->is_boundary(motorcycles[stop.motorcycle].currentPosition))
				{
					auto perpendicular = mesh->opposite_halfedge_handle(mesh->prev_halfedge_handle(motorcycles[stop.motorcycle].currentPosition));
					if (possibleExits.find(perpendicular) != possibleExits.end())
						validCollision = true;
				}
				auto h = mesh->opposite_halfedge_handle(motorcycles[stop.motorcycle].currentPosition);
				if (!mesh->is_boundary(h))
				{
					auto perpendicular = mesh->next_halfedge_handle(h);
					if (possibleExits.find(perpendicular) != possibleExits.end())
						validCollision = true;
				}
				if (validCollision)
				{
					motorcycles[stop.motorcycle].currentPosition = mesh->opposite_halfedge_handle(lastHalfedge);
					optionsWithCost.emplace(options.Motorcycles()[localCycleIdx].CostSoFar(), HeadOnCollision, localCycleIdx, stop.motorcycle, HEMesh::HalfedgeHandle());
				}
				
				return false;
			}
		}

		bool collideWithOccupied = false;
		bool leftCollisionValid = false;
		bool rightCollisionValid = false;

		size_t collisionWith = -1;

		auto h = mesh->opposite_halfedge_handle(lastHalfedge);
		int iterations = 0;
		while (true)
		{
			if (iterations++ > 10)
			{
				std::cout << "Cannot find occupied edges to the left around vertex " << mesh->from_vertex_handle(h).idx() << " (" << mesh->point(mesh->from_vertex_handle(h)) << ", coming from vertex " << mesh->to_vertex_handle(h).idx() << ")" << std::endl;
				std::cout << "Path of current motorcycle: ";
				for (auto p : options.Motorcycles()[localCycleIdx].Edges())
					PrintHalfedge(std::cout, p, *mesh);
				std::cout << std::endl;
				throw std::runtime_error("Error routing through patch");
			}
			h = mesh->opposite_halfedge_handle(h);
			//h is an outgoing edge
			if (mesh->is_boundary(h))
			{
				//cycle around
				h = lastHalfedge;
				CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) {return false; });
			}
			else
				h = mesh->next_halfedge_handle(h);
			auto it = occupiedEdges.find(mesh->edge_handle(h));
			if (it != occupiedEdges.end())
			{
				collideWithOccupied = true;
				if (possibleExits.find(h) != possibleExits.end())
				{					
					leftCollisionValid = true;
					collisionWith = it->second;
					auto& collidedMotorcycle = motorcycles[it->second];
					if ((collidedMotorcycle.startInPaddedBoundary && collidedMotorcycle.Path()[0] == h)
						|| (collidedMotorcycle.collideWithPaddedBoundary && mesh->opposite_halfedge_handle(collidedMotorcycle.Path().back()) == h))
						rightCollisionValid = true;					
				}
				break;
			}
			auto idxInLoop = loop.FindEdge(h);
			if (idxInLoop != (size_t)-1)
			{
				if (possibleExits.find(h) != possibleExits.end())
					leftCollisionValid = true;
				break;
			}
		}
		h = lastHalfedge;
		iterations = 0;
		while (true)
		{
			if (iterations++ > 10)
			{
				std::cout << "Cannot find occupied edges to the right around " << mesh->point(mesh->to_vertex_handle(h)) << std::endl;
				throw std::runtime_error("Error routing through patch");
			}
			h = mesh->opposite_halfedge_handle(h);
			if (mesh->is_boundary(h))
			{
				//cycle around
				h = lastHalfedge;
				CirculateForwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) {return false; });
				h = mesh->opposite_halfedge_handle(h);
			}
			else
				h = mesh->prev_halfedge_handle(h);
			//h is an incoming edge
			auto it = occupiedEdges.find(mesh->edge_handle(h));
			if (it != occupiedEdges.end())
			{
				collideWithOccupied = true;
				if (possibleExits.find(h) != possibleExits.end())
				{					
					rightCollisionValid = true;
					collisionWith = it->second;
					auto& collidedMotorcycle = motorcycles[it->second];
					if ((collidedMotorcycle.startInPaddedBoundary && mesh->opposite_halfedge_handle(collidedMotorcycle.Path()[0]) == h)
						|| (collidedMotorcycle.collideWithPaddedBoundary && collidedMotorcycle.Path().back() == h))
						leftCollisionValid = true;
				}
				break;
			}
			if (loop.FindEdge(h) != (size_t)-1)
			{
				if (possibleExits.find(h) != possibleExits.end())
					rightCollisionValid = true;
				break;
			}						
		}
		bool valid = false;
		if (collideWithOccupied)
			valid = leftCollisionValid && rightCollisionValid;
		else
			valid = leftCollisionValid || rightCollisionValid;
		if (valid)
			optionsWithCost.emplace(options.Motorcycles()[localCycleIdx].CostSoFar(), RouteWithCollision, localCycleIdx, collisionWith, HEMesh::HalfedgeHandle());

		return false;
	};

	//Run Dijkstra
	auto maxCostQuery = [&]() { if (optionsWithCost.size() == 0) return std::numeric_limits<float>::infinity(); else return optionsWithCost.begin()->cost; };
	options.SimulateCyclesInPatch(*mesh, isOriginalEdgeProp, patch, passThroughVertexAction, leavePatchAction, maxCostQuery, minCosDeviationAngle, allowTurnsOnRegularVertices);
	if (status != InProcess)
		return RouteResult(NoRoute);

	if (optionsWithCost.size() == 0)
		return RouteResult(NoRoute);
	
	auto& bestOption = *optionsWithCost.begin();
	auto bestPathEdges = options.Motorcycles()[bestOption.cycleIdx].Edges();	

	if (bestOption.type == RouteThrough)
		bestPathEdges.push_back(bestOption.edgeOutside);
	return RouteResult(bestOption.type, std::move(bestPathEdges), bestOption.collisionWith);
}

void MotorcycleGraph::RealizeRoute(size_t cycleIdx, const RouteResult & route)
{
	if (route.type == NoRoute)
		return;

	bool useLastSegment = route.type == RouteWithCollision || route.type == HeadOnCollision;

	auto& cycle = motorcycles[cycleIdx];
	for (int i = 0; i < route.route.size(); ++i)
	{
		OccupyEdge(mesh->edge_handle(route.route[i]), cycleIdx);
		if (i < route.route.size() - 1 || useLastSegment)
		{
			cycle.AddToPath(route.route[i], *mesh);
			visitedVertices.AccessOrCreateAtToVertex(cycle.Path().back()).emplace_back(cycleIdx, cycle.Path().size());
		}
		else
			cycle.currentPosition = route.route[i];
	}

	if (useLastSegment)
	{
		auto& passedThrough = visitedVertices.AccessAtToVertex(cycle.Path().back());
		if (passedThrough.size() == 2)
			motorcycles[passedThrough[0].motorcycle].crashesOnPath.insert(passedThrough[0].locationInMotorcyclePath);		
	}
	if (route.type == HeadOnCollision || route.type == RouteWithCollision)
	{
		motorcycles[cycleIdx].terminated = true;
		if (route.collisionWith == (size_t)-1)
			motorcycles[cycleIdx].collideWithPaddedBoundary = true;
	}
	if (route.type == HeadOnCollision)
	{
		motorcycles[cycleIdx].straightContinuationEnd = route.collisionWith;
		motorcycles[route.collisionWith].straightContinuationEnd = cycleIdx;
		motorcycles[route.collisionWith].terminated = true;
	}
}

size_t MotorcycleGraph::FindStraightContinuation(size_t motorcycleIdx, bool forward) const
{
	bool outForward;
	return FindStraightContinuation(motorcycleIdx, forward, outForward);
}

size_t MotorcycleGraph::FindStraightContinuation(size_t motorcycleIdx, bool inForward, bool& outForward) const
{
	outForward = inForward;
	do
	{
		if (outForward)
		{
			if (!motorcycles[motorcycleIdx].terminated)
				break;
			motorcycleIdx = motorcycles[motorcycleIdx].straightContinuationEnd;
		}
		else
			motorcycleIdx = motorcycles[motorcycleIdx].straightContinuationStart;
		outForward = !outForward;
	} while (motorcycleIdx != (size_t)-1 && motorcycles[motorcycleIdx].Path().size() == 0);
	return motorcycleIdx;
}

void MotorcycleGraph::AdvanceAllMotorcycles()
{
	if (status != InProcess)
		return;

	std::vector<PathContinuation> continuations;

	bool verbose = false;

	//check what motorcycles are at the entry of a fenced region and route them through
	bool processedSomething = false;
	for (int iActive = 0; iActive < activeMotorcycles.size();)
	{
		auto cycleIdx = activeMotorcycles[iActive];
		if (hasHighPriorityMotorcycles && !motorcycles[cycleIdx].highPriority)
		{
			++iActive;
			continue;
		}

		if (motorcycles[cycleIdx].terminated)
		{
			activeMotorcycles.erase(activeMotorcycles.begin() + iActive);
			continue;
		}

		processedSomething = true;

		auto enteringEdgeIt = patchSet->enteringEdgeToPatch.find(motorcycles[cycleIdx].currentPosition);
		if (enteringEdgeIt != patchSet->enteringEdgeToPatch.end())
		{			
			auto pathSizeBefore = motorcycles[cycleIdx].Path().size();
			auto& patch = patchSet->patches[enteringEdgeIt->second.patchIdx];
			//The motor cycle just entered a patch. Route the motorcycle through it.
			auto routeResult = RouteThroughPatch(enteringEdgeIt->first, patchSet->patches[enteringEdgeIt->second.patchIdx], *enteringEdgeIt);
			RealizeRoute(cycleIdx, routeResult);
			if (routeResult.type == NoRoute)
			{
				if (status != InProcess)
					return;
				if(verbose)
					std::cout << "Cannot find path through patch for motorcycle " << cycleIdx << " at " << mesh->point(mesh->from_vertex_handle(motorcycles[cycleIdx].currentPosition)) << std::endl;
				//could not find a path through the patch					

				visitedHalfEdges.erase(motorcycles[cycleIdx].currentPosition);

				//TODO: try to track back

				//Add an artificial valence-2 singularity to prevent the motorcycle from entering the region
				//We have two pairs and try to use the one that does not use non-original edges
			
				HEMesh::HalfedgeHandle secondStopBeforeH;
				if (motorcycles[cycleIdx].Path().size() >= 2)
					secondStopBeforeH = motorcycles[cycleIdx].Path().back();
				StoppingMotorcycles stoppingMotorcycles1(*mesh, motorcycles[cycleIdx].currentPosition, isOriginalEdgeProp);
				StoppingMotorcycles stoppingMotorcycles2(*mesh, secondStopBeforeH, isOriginalEdgeProp);

				bool canUseStopping2 = secondStopBeforeH.is_valid() &&
					patchSet->enteringEdgeToPatch.find(mesh->opposite_halfedge_handle(stoppingMotorcycles2.stopBeforeEdge)) == patchSet->enteringEdgeToPatch.end()
					&& (motorcycles[cycleIdx].crashesOnPath.empty() || *motorcycles[cycleIdx].crashesOnPath.rbegin() != motorcycles[cycleIdx].Path().size());

				//Try to find out which pair to use
				StoppingMotorcycles* useStoppingMotorcycles = &stoppingMotorcycles1;
				if (stoppingMotorcycles1.useNonOriginalEdge && canUseStopping2)
					useStoppingMotorcycles = &stoppingMotorcycles2;

				useStoppingMotorcycles->TryToRoute(*this, useStoppingMotorcycles == &stoppingMotorcycles1);
				if (!useStoppingMotorcycles->hasValidRoute)
				{
					if (useStoppingMotorcycles == &stoppingMotorcycles1)
					{
						if(canUseStopping2)
							useStoppingMotorcycles = &stoppingMotorcycles2;
					}
					else
						useStoppingMotorcycles = &stoppingMotorcycles1;
				}
				useStoppingMotorcycles->TryToRoute(*this, useStoppingMotorcycles == &stoppingMotorcycles1);

				//Realize the stopping pair
				if (useStoppingMotorcycles->hasValidRoute)
				{
					activeMotorcycles.erase(activeMotorcycles.begin() + iActive);
					motorcycles[cycleIdx].terminated = true;

					try {
						useStoppingMotorcycles->Realize(*this, motorcycles[cycleIdx].highPriority);
					}
					catch (...)
					{
						std::cout << mesh->point(mesh->from_vertex_handle(motorcycles[cycleIdx].currentPosition)) << std::endl;
						status = InternalFailure;
						return;
					}
					if (useStoppingMotorcycles == &stoppingMotorcycles2)
					{
						//remove the last path segment
						auto& visitedAtSource = visitedVertices.AccessAtToVertex(mesh->opposite_halfedge_handle(motorcycles[cycleIdx].currentPosition));
						visitedAtSource.erase(std::remove_if(visitedAtSource.begin(), visitedAtSource.end(),
							[&](const MotorcycleStop& stop) {return stop.motorcycle == cycleIdx && stop.locationInMotorcyclePath == motorcycles[cycleIdx].Path().size(); }), visitedAtSource.end());
						if (visitedAtSource.empty())
							visitedVertices.EraseAtToVertex(mesh->opposite_halfedge_handle(motorcycles[cycleIdx].currentPosition));
						visitedHalfEdges.erase(motorcycles[cycleIdx].Path().back());
						visitedHalfEdges.erase(mesh->opposite_halfedge_handle(motorcycles[cycleIdx].Path().back()));
						motorcycles[cycleIdx].RemoveLastPathSegment(*mesh);
					}
				}
				else
				{
					std::cout << "The stopping motorcycles cannot be routed through the patch." << std::endl;
					auto alternativeRouteResult = RouteThroughPatch(stoppingMotorcycles1.newSingularityH[useStoppingMotorcycles->failedStoppingMotorcycle], patchSet->patches[enteringEdgeIt->second.patchIdx], *enteringEdgeIt, true);
					RealizeRoute(cycleIdx, alternativeRouteResult);
					if (alternativeRouteResult.type == NoRoute)
					{
						std::cout << "The alternative route does not work at " << mesh->point(mesh->from_vertex_handle(enteringEdgeIt->first)) << "." << std::endl;
						status = CannotTraceThroughRegion;
						problematicFencedRegion = enteringEdgeIt->second.patchIdx;
						return;
					}
					else if (alternativeRouteResult.type == RouteWithCollision || alternativeRouteResult.type == HeadOnCollision)
					{
						activeMotorcycles.erase(activeMotorcycles.begin() + iActive);
						motorcycles[cycleIdx].terminated = true;
					}
					else if (alternativeRouteResult.type == RouteThrough)
						; //this motorcycle might have entered another patch, check it again (do not change iActive)
				}
			} //end if routeResult.type == NoRoute
			else if (routeResult.type == RouteWithCollision || routeResult.type == HeadOnCollision)
			{
				activeMotorcycles.erase(activeMotorcycles.begin() + iActive);
				motorcycles[cycleIdx].terminated = true;					
			}
			else if (routeResult.type == RouteThrough)
				; //this motorcycle might have entered another patch
		}
		else
			++iActive; //motorcycle does not enter a fenced region

	}

	//erase terminated motorcycles from the list of active motorcycles
	activeMotorcycles.erase(std::remove_if(activeMotorcycles.begin(), activeMotorcycles.end(), [&](size_t idx) { return motorcycles[idx].terminated; }), activeMotorcycles.end());	

	if (!processedSomething)
	{
		hasHighPriorityMotorcycles = false;
		return;
	}


	//check collisions on edges
	//sort motorcycles by edge index (motor cycles on the same edge will be next to each other)
	std::sort(activeMotorcycles.begin(), activeMotorcycles.end(),
		[&](const size_t first, const size_t second) { return mesh->edge_handle(motorcycles[first].currentPosition).idx() < mesh->edge_handle(motorcycles[second].currentPosition).idx(); });

	for (int i = activeMotorcycles.size() - 1; i >= 1; --i)
	{
		auto cycleIdx = activeMotorcycles[i];
		auto nextCycleIdx = activeMotorcycles[i - 1];	
		if (mesh->opposite_halfedge_handle(motorcycles[cycleIdx].currentPosition) == motorcycles[nextCycleIdx].currentPosition)
		{
			//there is a collision
			//merge the motor cycles			
			motorcycles[cycleIdx].AddToPath(motorcycles[cycleIdx].currentPosition, *mesh);
			if (motorcycles[nextCycleIdx].Path().size() == 0 && motorcycles[nextCycleIdx].startInPaddedBoundary)
				motorcycles[cycleIdx].collideWithPaddedBoundary = true;
			visitedVertices.AccessOrCreateAtToVertex(motorcycles[cycleIdx].currentPosition).emplace_back(cycleIdx, motorcycles[cycleIdx].Path().size());
			//motorcycles[nextCycleIdx].path.push_back(motorcycles[nextCycleIdx].currentPosition); //last edge is shared by the two motorcycles
			activeMotorcycles.erase(activeMotorcycles.begin() + (i - 1), activeMotorcycles.begin() + (i + 1));
			motorcycles[cycleIdx].terminated = true;
			motorcycles[nextCycleIdx].terminated = true;
			motorcycles[cycleIdx].straightContinuationEnd = nextCycleIdx;
			motorcycles[nextCycleIdx].straightContinuationEnd = cycleIdx;
			--i;
		}
	}

	//add the current edge to each motor cycle's path
//#pragma omp parallel for
	for (auto it = activeMotorcycles.begin(); it != activeMotorcycles.end(); ++it)
	{		
		auto cycleIdx = *it;

		if (hasHighPriorityMotorcycles && !motorcycles[cycleIdx].highPriority)
			continue;

		motorcycles[cycleIdx].AddToPath(motorcycles[cycleIdx].currentPosition, *mesh);
		visitedHalfEdges.insert(mesh->opposite_halfedge_handle(motorcycles[cycleIdx].currentPosition));
	}

	//Set a temporary context edge to differentiate motorcycles that have the same non-manifold target vertex but are on different
	//surface patches.
	for (int i = 0; i < activeMotorcycles.size(); ++i)
	{
		auto& motorcycle = motorcycles[activeMotorcycles[i]];
		if (mesh->is_manifold(mesh->to_vertex_handle(motorcycle.currentPosition)))
			motorcycle.temporaryContextEdge = 0;
		else
		{
			auto h = motorcycle.currentPosition;
			CirculateBackwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });
			motorcycle.temporaryContextEdge = h.idx();
		}
	}

	//check collisions on vertices
	//sort motorcycles by target vertex and temporary context edge
	std::sort(activeMotorcycles.begin(), activeMotorcycles.end(),
		[&](const size_t first, const size_t second) 
	{
		if (motorcycles[first].highPriority && !motorcycles[second].highPriority)
			return true;
		if (motorcycles[second].highPriority && !motorcycles[first].highPriority)
			return false;
		auto firstV = mesh->to_vertex_handle(motorcycles[first].currentPosition).idx();
		auto secondV = mesh->to_vertex_handle(motorcycles[second].currentPosition).idx();
		if (firstV != secondV)
			return firstV < secondV;
		return motorcycles[first].temporaryContextEdge < motorcycles[second].temporaryContextEdge;
	});

	for (int i = activeMotorcycles.size() - 1; i >= 0; --i)
	{
		auto cycleIdx = activeMotorcycles[i];

		if (hasHighPriorityMotorcycles && !motorcycles[cycleIdx].highPriority)
			continue;

		auto currentVertex = mesh->to_vertex_handle(motorcycles[cycleIdx].currentPosition);

		auto& passedThrough = visitedVertices.AccessOrCreateAtToVertex(motorcycles[cycleIdx].currentPosition);
		auto passedThroughBefore = passedThrough.size();
		passedThrough.emplace_back(cycleIdx, motorcycles[cycleIdx].Path().size());

		//check how many motor cycles are at the current vertex
		int j = i - 1;
		while (j >= 0 && currentVertex == mesh->to_vertex_handle(motorcycles[activeMotorcycles[j]].currentPosition) && motorcycles[cycleIdx].temporaryContextEdge == motorcycles[activeMotorcycles[j]].temporaryContextEdge)
		{
			passedThrough.emplace_back(activeMotorcycles[j], motorcycles[activeMotorcycles[j]].Path().size());
			--j;
		}
		//now, j is the first motor cycle that has a different target vertex
		
		size_t cyclesAtVertex = i - j;

		if (passedThroughBefore > 0)
		{
			//vertex has been visited before; stop all motor cycles
			for (int k = j + 1; k <= i; ++k)
				motorcycles[activeMotorcycles[k]].terminated = true;
			motorcycles[passedThrough[0].motorcycle].crashesOnPath.insert(passedThrough[0].locationInMotorcyclePath);
			i = j + 1;
			if(verbose)
				std::cout << "Stopped " << cyclesAtVertex << " motorcycles because the target vertex " << currentVertex << " has already been visited." << std::endl;
		}
		else
		{								
			//vertex has not been visited before
			if (cyclesAtVertex == 1)
			{
				//nothing special here; single motorcycle just continues
				PathContinuation continuation;				
				bool hasContinuation = FindBestContinuation(*mesh, motorcycles[cycleIdx].currentPosition, continuation);
				if (hasContinuation && continuation.edge.is_valid())
				{
					motorcycles[cycleIdx].currentPosition = continuation.edge;
					visitedHalfEdges.insert(continuation.edge);
				}
				else
				{
					motorcycles[cycleIdx].terminated = true;
					if (!continuation.edge.is_valid())
					{
						motorcycles[cycleIdx].collideWithPaddedBoundary = true;
						if (verbose)
							std::cout << "Stopped motorcycle because it collided with the boundary." << std::endl;
					}
					else
						if (verbose)
							std::cout << "Stopped motorcycle because it had no continuation. This should not happen." << std::endl;
				}
			}
			else
			{
				continuations.clear();
				continuations.resize(cyclesAtVertex);				

				//Find out where each cycle would go
				for (int k = 0; k < cyclesAtVertex; ++k)
				{
					PathContinuation continuation;
					bool hasContinuation = FindBestContinuation(*mesh, motorcycles[activeMotorcycles[j + 1 + k]].currentPosition, continuation);
					if (!hasContinuation)
						continuations[k].score = std::numeric_limits<float>::quiet_NaN();
					else
						continuations[k] = continuation;
				}

				bool handled = false;

				//Check if we can merge two cycles, i.e. their continuations align with each other
				for (int k = 0; k < cyclesAtVertex - 1; ++k)
				{
					if (handled)
						break;
					if (!continuations[k].edge.is_valid())
						continue; //no valid continuation
					for (int l = k + 1; l < cyclesAtVertex; ++l)
					{
						if (!continuations[l].edge.is_valid())
							continue; //no valid continuation
						if (continuations[k].edge.is_valid() && continuations[k].edge == mesh->opposite_halfedge_handle(motorcycles[activeMotorcycles[j + 1 + l]].currentPosition)
							&& continuations[l].edge.is_valid() && continuations[l].edge == mesh->opposite_halfedge_handle(motorcycles[activeMotorcycles[j + 1 + k]].currentPosition))
						{							
							//these two cycles can merge
							//k < l
							motorcycles[activeMotorcycles[j + 1 + k]].straightContinuationEnd = activeMotorcycles[j + 1 + l];
							motorcycles[activeMotorcycles[j + 1 + l]].straightContinuationEnd = activeMotorcycles[j + 1 + k];							

							//stop all cycles
							for (int k = j + 1; k <= i; ++k)
								motorcycles[activeMotorcycles[k]].terminated = true;
							if (verbose)
								std::cout << "Stopped motorcycles after merging 2." << std::endl;

							handled = true;
							break;
						}
					}
				}

				if (!handled)
				{

					size_t cyclesWithPossibleContinuation = 0;
					float bestContinuationScore = -std::numeric_limits<float>::infinity();
					int cycleWithBestContinuationScore;

					//check what motor cycles could not continue as planned
					for (int k = 0; k < cyclesAtVertex; ++k)
					{
						if (std::isnan(continuations[k].score))
							continue; //no valid continuation
						for (int l = 0; l < cyclesAtVertex; ++l)
						{
							if (continuations[k].edge == mesh->opposite_halfedge_handle(motorcycles[activeMotorcycles[j + 1 + l]].currentPosition))
							{
								//There is already a motorcycle at the target position
								continuations[k].score = std::numeric_limits<float>::quiet_NaN();
								break;
							}
						}
						if (!std::isnan(continuations[k].score))
						{
							++cyclesWithPossibleContinuation;
							if (continuations[k].score > bestContinuationScore)
							{
								bestContinuationScore = continuations[k].score;
								cycleWithBestContinuationScore = k;
							}
						}
					}

					if (cyclesWithPossibleContinuation == 0)
					{
						//stop all motor cycles
						for (int k = j + 1; k <= i; ++k)
							motorcycles[activeMotorcycles[k]].terminated = true;
						if (verbose)
							std::cout << "Stopping " << cyclesAtVertex << " motorcycles because no motorcycle can continue." << std::endl;
					}
					else
					{
						//stop all others
						for (int k = j + 1; k <= i; ++k)
						{
							if (k != j + 1 + cycleWithBestContinuationScore)
								motorcycles[activeMotorcycles[k]].terminated = true;
						}

						//continue the cycle with the best continuation
						if (continuations[cycleWithBestContinuationScore].edge.is_valid())
						{
							motorcycles[activeMotorcycles[j + 1 + cycleWithBestContinuationScore]].currentPosition = continuations[cycleWithBestContinuationScore].edge;
							visitedHalfEdges.insert(continuations[cycleWithBestContinuationScore].edge);

							motorcycles[activeMotorcycles[j + 1 + cycleWithBestContinuationScore]].crashesOnPath.insert(motorcycles[activeMotorcycles[j + 1 + cycleWithBestContinuationScore]].Path().size());
						}
						else
						{
							motorcycles[activeMotorcycles[j + 1 + cycleWithBestContinuationScore]].terminated = true;
							motorcycles[activeMotorcycles[j + 1 + cycleWithBestContinuationScore]].collideWithPaddedBoundary = true;
						}						
						

						if (verbose)
							std::cout << "Stopping " << cyclesAtVertex << " motorcycles after continuing one." << std::endl;
					}
				}

				i = j + 1;
			}
		}
	}
}

void MotorcycleGraph::AdvanceToEnd()
{
	if (status == Finished)
		status = InProcess;

	nse::util::TimedBlock b("Advancing motorcycles ..");
	try
	{
		while (HasActiveMotorcycles() && status == InProcess)
			AdvanceAllMotorcycles();

		if (status == InProcess)
			status = Finished;
	}
	catch (std::exception& e)
	{
		std::cout << "Exception during motorcycle tracing: " << e.what() << std::endl;
		status = InternalFailure;
	}
}

const std::vector<MotorcycleGraph::Motorcycle>& MotorcycleGraph::Motorcycles() const
{
	return motorcycles;
}

std::vector<MotorcycleGraph::MotorcycleStop>::const_iterator MotorcycleGraph::FindMotorcycleAtEdge(HEMesh::HalfedgeHandle h, const std::vector<MotorcycleGraph::MotorcycleStop>& motorcyclesAtSourceVertex, bool& outForward) const
{
	//motorcycle that continues forward
	auto motorcycleStopIt = std::find_if(motorcyclesAtSourceVertex.begin(), motorcyclesAtSourceVertex.end(), [&](const MotorcycleGraph::MotorcycleStop& stop)
	{
		auto& motorcycle = motorcycles[stop.motorcycle];
		if (stop.locationInMotorcyclePath >= motorcycle.Path().size())
			return false;
		return h == motorcycle.Path()[stop.locationInMotorcyclePath];
	});
	if (motorcycleStopIt != motorcyclesAtSourceVertex.end())
	{
		outForward = true;
		return motorcycleStopIt;
	}

	//motorcycle that continues backward
	motorcycleStopIt = std::find_if(motorcyclesAtSourceVertex.begin(), motorcyclesAtSourceVertex.end(), [&](const MotorcycleStop& stop)
	{
		auto& motorcycle = motorcycles[stop.motorcycle];
		if (stop.locationInMotorcyclePath <= 0)
			return false;
		return h == mesh->opposite_halfedge_handle(motorcycles[stop.motorcycle].Path()[stop.locationInMotorcyclePath - 1]);
	});
	if (motorcycleStopIt != motorcyclesAtSourceVertex.end())
	{
		outForward = false;
		return motorcycleStopIt;
	}

	return motorcycleStopIt;
}

MotorcycleGraph::LocationOnPath MotorcycleGraph::EdgeToMotorcycleLocation(HEMesh::HalfedgeHandle h, const std::vector<MotorcycleGraph::MotorcycleStop>& motorcyclesAtSourceVertex) const
{
	LocationOnPath result;
	auto stopIt = FindMotorcycleAtEdge(h, motorcyclesAtSourceVertex, result.goForward);
	if (stopIt != motorcyclesAtSourceVertex.end())
	{
		result.motorcycle = stopIt->motorcycle;
		if (result.goForward)
			result.pathSegment = stopIt->locationInMotorcyclePath;
		else
			result.pathSegment = stopIt->locationInMotorcyclePath - 1;
	}
	else
	{
		result.motorcycle = -1;
		result.pathSegment = -1;
	}
	return result;
}

MotorcycleGraph::LocationOnPath MotorcycleGraph::NextEdgeInPatch(const LocationOnPath& current, int& outColorChange, std::vector<size_t>& skippedMotorcycles, bool verbose) const
{
	auto location = NextEdgeInPatchInternal(current, outColorChange, verbose);
	while (motorcycles[location.motorcycle].isDeactivated)
	{
		skippedMotorcycles.push_back(location.motorcycle);
		if (verbose)
			std::cout << "Turning around on a deactivated motorcycle." << std::endl;

		outColorChange -= 2;
		location.goForward = !location.goForward;
		int skipColorChange;
		location = NextEdgeInPatchInternal(location, skipColorChange, verbose);
		outColorChange += skipColorChange;
	}
	return location;
}

bool MotorcycleGraph::CollidesWithBoundary(const LocationOnPath& location) const
{
	auto goDirectlyToBoundary = false;
	//Check if this motorcycle collides directly with the boundary
	goDirectlyToBoundary = goDirectlyToBoundary || (location.goForward && motorcycles[location.motorcycle].collideWithPaddedBoundary && location.pathSegment == motorcycles[location.motorcycle].Path().size() - 1);
	goDirectlyToBoundary = goDirectlyToBoundary || (!location.goForward && motorcycles[location.motorcycle].startInPaddedBoundary && location.pathSegment == 0);
	
	//No check possible straight continuations
	auto motorcycle = location.motorcycle;
	auto forward = location.goForward;
	while (true)
	{
		if (forward)
			motorcycle = motorcycles[motorcycle].straightContinuationEnd;
		else
			motorcycle = motorcycles[motorcycle].straightContinuationStart;
		forward = !forward;
		if (motorcycle == (size_t)-1)
			break;

		goDirectlyToBoundary = goDirectlyToBoundary || (forward && motorcycles[motorcycle].collideWithPaddedBoundary && motorcycles[motorcycle].Path().empty());
		goDirectlyToBoundary = goDirectlyToBoundary || (!forward && motorcycles[motorcycle].startInPaddedBoundary && motorcycles[motorcycle].Path().empty());
	}

	return goDirectlyToBoundary;
}

bool MotorcycleGraph::IsValidTransitionFromMotorcycleToBoundary(const LocationOnPath& location, HEMesh::HalfedgeHandle edgeOnBoundary, const std::vector<MotorcycleStop>& motorCyclesAtTransition, bool forward, int& outColorChange, bool verbose) const
{	
	bool isValidTransition = false;
	auto oldMotorcycleH = MotorcycleHalfedge(location);
	outColorChange = 0;
	//first option: this motorcycle collided with a singularity or starts there and the singularity has an outgoing edge to the boundary	
	{
		const size_t* metaSingularityIdxPtr;
		if (vertexToMetaSingularityIndex.TryAccessAtToVertex(oldMotorcycleH, metaSingularityIdxPtr))
		{
			auto& sing = metaSingularities[*metaSingularityIdxPtr];
			auto oldEdgeInEmanatingList = sing.FindOutgoingMotorcycle(mesh->opposite_halfedge_handle(oldMotorcycleH));
			//go to the next emanating motorcycle
			if (forward)
			{				
				if (oldEdgeInEmanatingList == 0)
					oldEdgeInEmanatingList = sing.Degree();
				--oldEdgeInEmanatingList;
			}
			else
			{
				++oldEdgeInEmanatingList;
				if (oldEdgeInEmanatingList == sing.Degree())
					oldEdgeInEmanatingList = 0;
			}
			if (sing.GetEmanatingMotorcycle(oldEdgeInEmanatingList).empty())
			{
				//valid transition to boundary
				outColorChange += 2;
				isValidTransition = true;
				if (verbose)
					std::cout << "Going to boundary over singularity" << std::endl;
			}
		}
	}

	//second option: this motorcycle collided with the boundary directly
	if (!isValidTransition)
	{
		if (CollidesWithBoundary(location))
		{
			outColorChange += 1;
			isValidTransition = true;
			if (verbose)
				std::cout << "Going to boundary directly." << std::endl;
		}
	}

	//third option: the new motorcycle collided with another motorcycle that immediately collided with the boundary
	if (!isValidTransition)
	{
		auto checkEdge = (forward 
			? mesh->prev_halfedge_handle(mesh->opposite_halfedge_handle(oldMotorcycleH))
			: oldMotorcycleH);
		CirculateUntilMotorcycleCondition stopCondition(*this);
		stopCondition.motorcyclesAtQueryTargetVertex = &motorCyclesAtTransition;
		auto circulationResult = (forward
			? CirculateBackwardUntil<true>(checkEdge, *mesh, stopCondition)
			: CirculateForwardUntil<true>(checkEdge, *mesh, stopCondition));
		if(circulationResult == StoppedByCondition)
		{
			stopCondition.result.goForward = !stopCondition.result.goForward;
			if (CollidesWithBoundary(stopCondition.result))
			{
				outColorChange += 2;
				isValidTransition = true;
				if (verbose)
					std::cout << "Going to boundary via other motorcycle" << std::endl;
			}
		}
	}

	return isValidTransition;
}

MotorcycleGraph::LocationOnPath MotorcycleGraph::NextEdgeInPatchInternal(const LocationOnPath& current, int& outColorChange, bool verbose) const
{
	auto h = MotorcycleHalfedge(current);
	auto oldMotorcycleH = h;	

	outColorChange = 0;

	CirculateUntilMotorcycleCondition stopCondition(*this);
	stopCondition.motorcyclesAtQueryTargetVertex = nullptr;

	bool haveResult = false;
	bool isOnPaddedBoundary = false;
	CirculationResult circulationResult;

	//if we are at the end, check for transition to padded boundary
	if ((current.pathSegment == 0 && !current.goForward) || (current.pathSegment == motorcycles[current.motorcycle].Path().size() - 1 && current.goForward))
	{
		stopCondition.motorcyclesAtQueryTargetVertex = &visitedVertices.AccessAtToVertex(h);
		circulationResult = CirculateForwardUntil<true>(h, *mesh, stopCondition);
		if (circulationResult == StoppedByCondition)
			haveResult = true;
		else 
		{
			//stopped by boundary, check if this is a valid transition
			int colorChange;
			isOnPaddedBoundary = IsValidTransitionFromMotorcycleToBoundary(current, h, *stopCondition.motorcyclesAtQueryTargetVertex, true, colorChange, verbose);
			outColorChange += colorChange;
		}
	}

	if (isOnPaddedBoundary && verbose)
	{
		std::cout << "Going over to boundary at " << mesh->point(mesh->to_vertex_handle(oldMotorcycleH)) << " from ";
		PrintFullHalfedge(std::cout, oldMotorcycleH, *mesh);
		std::cout << " to ";
		PrintFullHalfedge(std::cout, h, *mesh);
		std::cout << std::endl;
	}

	//as long as we are moving on the boundary, check for valid transitions back onto a motorcycle
	while (isOnPaddedBoundary)
	{
		auto oldH = h;
		bool isVisited = visitedVertices.TryAccessAtToVertex(h, stopCondition.motorcyclesAtQueryTargetVertex);
		if (!isVisited)
			stopCondition.motorcyclesAtQueryTargetVertex = nullptr;

		if (stopCondition(mesh->opposite_halfedge_handle(h)))
		{
			h = mesh->opposite_halfedge_handle(h);
			circulationResult = StoppedByCondition;
		}
		else
			circulationResult = CirculateForwardUntil<true>(h, *mesh, stopCondition);
		if (circulationResult == StoppedByCondition)
		{
			haveResult = true;
			//check if this is a valid transition from padded boundary back to a motorcycle
			int colorChange;
			auto oppH = mesh->opposite_halfedge_handle(oldH);
			auto oppSeg = stopCondition.result;
			oppSeg.goForward = !oppSeg.goForward;
			isOnPaddedBoundary = !IsValidTransitionFromMotorcycleToBoundary(oppSeg, oppH, *stopCondition.motorcyclesAtQueryTargetVertex, false, colorChange, verbose);
			outColorChange += colorChange;			

			if (isOnPaddedBoundary)
			{
				//we have been stopped by a motorcycle but the transition was invalid. Go to the boundary.
				h = mesh->opposite_halfedge_handle(h);
				CirculateForwardUntil<true>(h, *mesh, [](HEMesh::HalfedgeHandle) { return false; });
			}
		}		
		
		if (verbose)
		{
			std::cout << "Going over from ";
			PrintFullHalfedge(std::cout, oldH, *mesh);
			std::cout << " to ";
			PrintFullHalfedge(std::cout, h, *mesh);
			if(isOnPaddedBoundary)
				std::cout << " (still on boundary)" << std::endl;
			else
				std::cout << " (not on boundary anymore)" << std::endl;
		} 

		if (!isOnPaddedBoundary)
			return stopCondition.result;
	}
	
	if (!haveResult)
	{		
		h = oldMotorcycleH;
		stopCondition.motorcyclesAtQueryTargetVertex = &visitedVertices.AccessAtToVertex(h);
		circulationResult = CirculateForwardUntil<false>(h, *mesh, stopCondition);
	}

	bool isStraightContinuation = false;
	if (!current.goForward && current.pathSegment == 0 && FindStraightContinuation(current.motorcycle, current.goForward) == stopCondition.result.motorcycle)
		isStraightContinuation = true;
	
	if (current.goForward && current.pathSegment == motorcycles[current.motorcycle].Path().size() - 1 && FindStraightContinuation(current.motorcycle, current.goForward) == stopCondition.result.motorcycle)
		isStraightContinuation = true;

	if (current.motorcycle == stopCondition.result.motorcycle && current.goForward == stopCondition.result.goForward
		&& ((current.goForward && current.pathSegment + 1 == stopCondition.result.pathSegment)
		|| (!current.goForward && current.pathSegment - 1 == stopCondition.result.pathSegment)))
		isStraightContinuation = true;

	if (isStraightContinuation)
		outColorChange = 0;
	else
		outColorChange = 1;

	return stopCondition.result;
}

void MotorcycleGraph::AssignFaceIndexToHalfarcs(size_t patchIdx)
{
	auto& patch = patches[patchIdx];

	for (auto& side : patch.PatchSides())
	{
		for (auto arcIdx : side)
		{
			auto& arc = halfArcs[arcIdx];
			arc.face = patchIdx;
		}
	}
}

std::set<size_t>::iterator MotorcycleGraph::NextCrashSite(const LocationOnPath& location) const
{
	auto& motorcycle = motorcycles[location.motorcycle];
	if (location.goForward)
		return std::upper_bound(motorcycle.crashesOnPath.begin(), motorcycle.crashesOnPath.end(), location.pathSegment);
	else
	{
		auto it = std::upper_bound(motorcycle.crashesOnPath.begin(), motorcycle.crashesOnPath.end(), location.pathSegment);
		if (it != motorcycle.crashesOnPath.begin())
			--it;
		else
			it = motorcycle.crashesOnPath.end();
		return it;
	}
}

HEMesh::HalfedgeHandle MotorcycleGraph::MotorcycleHalfedge(const LocationOnPath& location) const
{
	auto halfedge = motorcycles[location.motorcycle].Path()[location.pathSegment];
	if (!location.goForward)
		halfedge = mesh->opposite_halfedge_handle(halfedge);
	return halfedge;
}

void MotorcycleGraph::ExtractArcsForMotorcycle(size_t iMotorcycle, std::vector<std::pair<size_t, size_t>>& arcMerges, bool checkCollisionPoints)
{
	auto& motorcycle = motorcycles[iMotorcycle];

	if (motorcycle.Path().empty())
		return;
	if (motorcycle.isDeactivated)
		return;

	auto closeHalfArc = [&](size_t pathSegment)
	{
		auto motorcycle = halfArcs.back().segments.front().location.motorcycle;
		auto length = pathSegment - halfArcs.back().segments.front().location.pathSegment;;
		halfArcs.back().segments.front().length = length;
		firstHalfedgeToHalfArc[MotorcycleHalfedge(halfArcs.back().segments.front().location)] = halfArcs.size() - 1;

		//add the opposite halfarc
		halfArcs.emplace_back();
		halfArcs.back().segments.emplace_back(LocationOnPath(motorcycle, pathSegment - 1, false), length);
		firstHalfedgeToHalfArc[MotorcycleHalfedge(halfArcs.back().segments.front().location)] = halfArcs.size() - 1;
		halfArcs.back().multiplier = 1;
	};

	halfArcs.emplace_back();
	halfArcs.back().segments.emplace_back(LocationOnPath(iMotorcycle, 0, true), 0);
	halfArcs.back().multiplier = motorcycles[iMotorcycle].multiplier;

	//check if we have a continuing motorcycle at the beginning and schedule for merge if possible
	auto continuingMotorcycleStart = FindStraightContinuation(iMotorcycle, false);
	if (continuingMotorcycleStart != (size_t)-1)
	{
		//find if the other motorcycle already has an arc
		auto arcIt = firstHalfedgeToHalfArc.find(motorcycles[continuingMotorcycleStart].Path().front());
		if (arcIt != firstHalfedgeToHalfArc.end())
		{
			//check if we have other motorcycles colliding at the start point
			bool hasOtherMotorcycles = false;
			for (auto& stop : visitedVertices.AccessAtToVertex(mesh->opposite_halfedge_handle(motorcycle.Path().front())))
			{
				if (stop.motorcycle != iMotorcycle && stop.motorcycle != continuingMotorcycleStart && !motorcycles[stop.motorcycle].isDeactivated)
					hasOtherMotorcycles = true;
			}
			if (!hasOtherMotorcycles)
				arcMerges.emplace_back(OppositeHalfarc(halfArcs.size() - 1), arcIt->second);
		}
	}

	auto nextCollision = motorcycle.crashesOnPath.begin();
	if (nextCollision != motorcycle.crashesOnPath.end() && *nextCollision == 0)
		++nextCollision;
	bool wasOnBorder = false;
	//go along the path and create new arcs whenever there is an intersection with another motorcycle
	for (size_t pathSegment = 0; pathSegment != motorcycle.Path().size(); ++pathSegment)
	{
		auto h = motorcycle.Path()[pathSegment];

		bool newSegment = false;
		if (nextCollision != motorcycle.crashesOnPath.end() && *nextCollision == pathSegment)
		{
			++nextCollision;
			if (pathSegment != 0)
			{
				if (checkCollisionPoints)
				{
					bool hasOtherMotorcycles = false;
					try {
						for (auto& stop : visitedVertices.AccessAtToVertex(motorcycle.Path()[pathSegment - 1]))
						{
							if (!(stop.motorcycle == iMotorcycle && stop.locationInMotorcyclePath == pathSegment) && !motorcycles[stop.motorcycle].isDeactivated)
								hasOtherMotorcycles = true;
						}
					}
					catch (std::runtime_error&)
					{
						std::cout << "Error for halfedge ";
						PrintFullHalfedge(std::cout, motorcycle.Path()[pathSegment - 1], *mesh);
						std::cout << std::endl;
					}
					if (hasOtherMotorcycles)
						newSegment = true;
				}
				else
					newSegment = true;
			}
		}

		bool isOnBorder = mesh->is_boundary(h) || mesh->is_boundary(mesh->opposite_halfedge_handle(h));
		if (wasOnBorder ^ isOnBorder)
			newSegment = true;
		if (!wasOnBorder && !isOnBorder && mesh->is_boundary(mesh->from_vertex_handle(h)))
			newSegment = true;

		if (pathSegment == 0)
			newSegment = false;

		if (newSegment)
		{
			closeHalfArc(pathSegment);

			halfArcs.emplace_back();
			halfArcs.back().segments.emplace_back(LocationOnPath(iMotorcycle, pathSegment, true), 0);
			halfArcs.back().multiplier = motorcycles[iMotorcycle].multiplier;
		}

		wasOnBorder = isOnBorder;
	}
	closeHalfArc(motorcycle.Path().size());

	//check if we have a continuing motorcycle at the end and schedule for merge if possible
	auto continuingMotorcycleEnd = FindStraightContinuation(iMotorcycle, true);
	if (continuingMotorcycleEnd != (size_t)-1)
	{
		//find if the other motorcycle already has an arc
		auto arcIt = firstHalfedgeToHalfArc.find(mesh->opposite_halfedge_handle(motorcycles[continuingMotorcycleEnd].Path().back()));
		if (arcIt != firstHalfedgeToHalfArc.end())
		{
			//check if we have other motorcycles colliding at the start point
			bool hasOtherMotorcycles = false;
			for (auto& stop : visitedVertices.AccessAtToVertex(motorcycle.Path().back()))
			{
				if (stop.motorcycle != iMotorcycle && stop.motorcycle != continuingMotorcycleEnd && !motorcycles[stop.motorcycle].isDeactivated)
					hasOtherMotorcycles = true;
			}
			if (!hasOtherMotorcycles)
				arcMerges.emplace_back(halfArcs.size() - 2, arcIt->second);
		}
	}
}

void MotorcycleGraph::ExtractArcs()
{
	nse::util::TimedBlock b("Extracting arcs of motorcycle graph ..");
	//Extract the graph arcs
	halfArcs.clear();
	firstHalfedgeToHalfArc.clear();

	std::vector<std::pair<size_t, size_t>> arcMerges;	

	//Extract arcs from all motorcycles
	for (size_t iMotorcycle = 0; iMotorcycle != motorcycles.size(); ++iMotorcycle)
	{		
		ExtractArcsForMotorcycle(iMotorcycle, arcMerges);
	}
	std::cout << halfArcs.size() << " halfarcs in motorcycle graph before merge." << std::endl;	

	//Try to merge arcs of straight continuations
	std::map<size_t, size_t> oldArcToNewArc;
	for (auto& merge : arcMerges)
	{
		auto arc1 = merge.first;		
		std::map<size_t, size_t>::iterator it;
		while ((it = oldArcToNewArc.find(arc1)) != oldArcToNewArc.end())
			arc1 = it->second;
	
		auto arc2 = merge.second;		
		while ((it = oldArcToNewArc.find(arc2)) != oldArcToNewArc.end())
			arc2 = it->second;
	
		firstHalfedgeToHalfArc.erase(MotorcycleHalfedge(halfArcs[arc2].segments.front().location));
		halfArcs[arc1].MergeWith(std::move(halfArcs[arc2]));
		oldArcToNewArc[arc2] = arc1;
	
		auto opp1 = OppositeHalfarc(arc1);
		auto opp2 = OppositeHalfarc(arc2);
		firstHalfedgeToHalfArc.erase(MotorcycleHalfedge(halfArcs[opp1].segments.front().location));
		halfArcs[opp2].MergeWith(std::move(halfArcs[opp1]));
		halfArcs[opp1] = std::move(halfArcs[opp2]);
		oldArcToNewArc[opp2] = opp1;
	}
	
	halfArcs.erase(std::remove_if(halfArcs.begin(), halfArcs.end(), [](const HalfArc& a) { return a.segments.empty(); }), halfArcs.end());

	//Record correspondences between halfedges and arcs
	for (int i = 0; i < halfArcs.size(); ++i)
	{
		auto& arc = halfArcs[i];
		firstHalfedgeToHalfArc[MotorcycleHalfedge(arc.segments.front().location)] = i;
	}

	std::cout << halfArcs.size() << " halfarcs in motorcycle graph after merge." << std::endl;
}

void MotorcycleGraph::DeactivateUnnecessaryMotorcycles()
{
	nse::util::TimedBlock b("Deactivating unnecessary motorcycles ..");

	//randomized order
	const int attempts = 20;

	std::mt19937 g;

	std::vector<size_t> bestAttemptDeactivatedMotorcycles;
	size_t bestAttemptDeactivatedLength = 0;
	for (int attempt = 0; attempt < attempts; ++attempt)
	{
		std::deque<size_t> candidatesToRemove;
		std::vector<size_t> thisAttemptDeactivatedMotorcycles;
		size_t thisAttemptDeactivatedLength = 0;
		//Activate all motorcycles and fill the candidates list
		for (size_t i = 0; i < motorcycles.size(); ++i)
		{
			motorcycles[i].isDeactivated = false;
			candidatesToRemove.push_back(i);
		}
		//randomize the list of candidates
		std::shuffle(candidatesToRemove.begin(), candidatesToRemove.end(), g);

		std::vector<size_t> arcsInStraightContinuation;
		while (!candidatesToRemove.empty())
		{
			//check if we can disable this motorcycle, i.e. if it (and all its straight continuations) do
			//not have colliding motorcycles
			auto candidateIdx = candidatesToRemove.front();
			auto* candidate = &motorcycles[candidateIdx];
			candidatesToRemove.pop_front();

			if (candidate->isDeactivated)
				continue;
			if (candidate->Path().size() == 0)
				continue;

			if (candidate->straightContinuationEnd != (size_t)-1 && candidate->straightContinuationStart != (size_t)-1)
				continue; //this motorcycle is in the middle of some straight continuations

			HEMesh::HalfedgeHandle firstH, lastH;

			bool checkForward = candidate->straightContinuationStart == (size_t)-1;
			bool hasNonTrivialCollisions = false;
			arcsInStraightContinuation.clear();

			firstH = (checkForward ? candidate->Path().front() : mesh->opposite_halfedge_handle(candidate->Path().back()));

			while (candidate != nullptr)
			{
				arcsInStraightContinuation.push_back(candidateIdx);
				lastH = (checkForward ? candidate->Path().back() : mesh->opposite_halfedge_handle(candidate->Path().front()));

				auto checkCollision = [&](size_t collisionLocationInMotorcyclePath)
				{
					if (collisionLocationInMotorcyclePath == 0 && candidate->straightContinuationStart == (size_t)-1 && !candidate->startInPaddedBoundary)
						return;
					if (collisionLocationInMotorcyclePath == candidate->Path().size() && candidate->straightContinuationEnd == (size_t)-1 && !candidate->collideWithPaddedBoundary)
						return;
					auto& visited = visitedVertices.AccessAtToVertex((collisionLocationInMotorcyclePath == 0 ? mesh->opposite_halfedge_handle(candidate->Path().front()) : candidate->Path()[collisionLocationInMotorcyclePath - 1]));
					for (auto& stop : visited)
					{
						if (motorcycles[stop.motorcycle].isDeactivated)
							continue;
						if (stop.motorcycle == candidateIdx && stop.locationInMotorcyclePath == collisionLocationInMotorcyclePath)
							continue;
						if (collisionLocationInMotorcyclePath == 0 && stop.motorcycle == candidate->straightContinuationStart)
							continue;
						if (collisionLocationInMotorcyclePath == candidate->Path().size() && stop.motorcycle == candidate->straightContinuationEnd)
							continue;
						hasNonTrivialCollisions = true;
					}
				};

				checkCollision(0);
				if (hasNonTrivialCollisions)
					break;
				checkCollision(candidate->Path().size());
				if (hasNonTrivialCollisions)
					break;
				for (auto collision : candidate->crashesOnPath)
					checkCollision(collision);
				if (hasNonTrivialCollisions)
					break;

				candidateIdx = FindStraightContinuation(candidateIdx, checkForward);
				if (candidateIdx == (size_t)-1)
					candidate = nullptr;
				else
					candidate = &motorcycles[candidateIdx];
				checkForward = !checkForward;
			}

			if (hasNonTrivialCollisions)
				continue;

			//Determines if deactivating a motorcycle around a singularity does not leave a gap of more than
			//one deactivated singularities.
			auto canDeactivateMotorcycleAtSingularity = [&](HEMesh::HalfedgeHandle outgoingH) -> bool
			{
				const size_t* singularityIdx;
				if (vertexToMetaSingularityIndex.TryAccessAtToVertex(mesh->opposite_halfedge_handle(outgoingH), singularityIdx))
				{
					auto& sing = metaSingularities[*singularityIdx];
					auto edgeInEmanatingList = sing.FindOutgoingMotorcycle(outgoingH);
					auto prev = (edgeInEmanatingList == 0 ? sing.Degree() - 1 : edgeInEmanatingList - 1);
					auto next = (edgeInEmanatingList == sing.Degree() - 1 ? 0 : edgeInEmanatingList + 1);
					bool prevActive = sing.GetEmanatingMotorcycle(prev).size() == 0;
					bool nextActive = sing.GetEmanatingMotorcycle(next).size() == 0;
					HEMesh::HalfedgeHandle prevH, nextH;
					if (!prevActive)
						prevH = sing.GetEmanatingMotorcycle(prev).front();
					if (!nextActive)
						nextH = sing.GetEmanatingMotorcycle(next).front();
					for (auto& stop : visitedVertices.AccessAtToVertex(mesh->opposite_halfedge_handle(outgoingH)))
					{
						if (motorcycles[stop.motorcycle].Path().size() == 0)
							continue;
						auto h = (stop.locationInMotorcyclePath == 0 ? motorcycles[stop.motorcycle].Path().front()
							: mesh->opposite_halfedge_handle(motorcycles[stop.motorcycle].Path().back()));

						if (h == prevH && !motorcycles[stop.motorcycle].isDeactivated)
							prevActive = true;
						if (h == nextH && !motorcycles[stop.motorcycle].isDeactivated)
							nextActive = true;
					}
					return prevActive && nextActive;
				}
				return true;
			};

			if (!canDeactivateMotorcycleAtSingularity(firstH))
				continue;
			if (!canDeactivateMotorcycleAtSingularity(mesh->opposite_halfedge_handle(lastH)))
				continue;

			//If we reached this point, we are allowed to deactivate this motorcycle (and all motorcycles in its straight continuations)

			for (auto candidateIdx : arcsInStraightContinuation)
			{
				motorcycles[candidateIdx].isDeactivated = true;
				thisAttemptDeactivatedMotorcycles.push_back(candidateIdx);
				thisAttemptDeactivatedLength += motorcycles[candidateIdx].NonBoundaryEdges();
			}
			//Deactivating this motorcycle might allow more motorcycles to be deactivated; schedule them for the check
			for (auto& stop : visitedVertices.AccessAtToVertex(mesh->opposite_halfedge_handle(firstH)))
				candidatesToRemove.push_back(stop.motorcycle);
			for (auto& stop : visitedVertices.AccessAtToVertex(lastH))
				candidatesToRemove.push_back(stop.motorcycle);
		}

		if (thisAttemptDeactivatedLength > bestAttemptDeactivatedLength)
		{
			bestAttemptDeactivatedMotorcycles = std::move(thisAttemptDeactivatedMotorcycles);
			bestAttemptDeactivatedLength = thisAttemptDeactivatedLength;
		}
	}

	for (size_t i = 0; i < motorcycles.size(); ++i)
		motorcycles[i].isDeactivated = false;
	for (auto deactivated : bestAttemptDeactivatedMotorcycles)
		motorcycles[deactivated].isDeactivated = true;
}


template <typename PatchOutput>
void MotorcycleGraph::BuildPatchesFromArcs(const std::vector<size_t>& halfarcs, const std::vector<int>& halfarcColors, int firstArcOnFirstSide, int maxColor, const std::vector<HEMesh::HalfedgeHandle>& additionalBoundaries, PatchOutput&& patchOutput) const
{
	//Find the connected components of this patch
	std::vector<bool> arcHandled(halfarcs.size(), false);

	std::map<HEMesh::HalfedgeHandle, size_t> boundaryEdges;
	for (int i = 0; i < halfarcs.size(); ++i)
	{
		auto& arc = halfArcs[halfarcs[i]];
		for (auto segment : arc)
			boundaryEdges[MotorcycleHalfedge(segment)] = i;
	}
	for (auto h : additionalBoundaries)
		boundaryEdges[h] = (size_t)-1;

	for (int iArc = 0; iArc < halfarcs.size(); ++iArc)
	{
		if (arcHandled[iArc])
			continue;
		auto& arc = halfArcs[halfarcs[iArc]];

		std::vector<bool> isArcInThisConnectedComponent(halfarcs.size(), false);
		std::set<HEMesh::FaceHandle> facesInThisConnectedComponent;
		std::queue<HEMesh::FaceHandle> bfsQueue;
		std::vector<HEMesh::HalfedgeHandle> openBoundaryEdges;

		//Initialize BFS from this arc
		for (auto segment : arc)
		{
			auto h = MotorcycleHalfedge(segment);
			auto f = mesh->face_handle(h);
			if (!f.is_valid())
				continue;

			facesInThisConnectedComponent.insert(f);
			bfsQueue.push(f);
		}

		isArcInThisConnectedComponent[iArc] = true;
		arcHandled[iArc] = true;

		//Perform BFS to find the faces, arcs, and open boundaries of this connected component.
		while (!bfsQueue.empty())
		{
			auto f = bfsQueue.front();
			bfsQueue.pop();
			for (auto h : mesh->fh_range(f))
			{
				auto boundaryIt = boundaryEdges.find(h);
				if (boundaryIt != boundaryEdges.end())
				{
					if (boundaryIt->second == (size_t)-1)
					{
						openBoundaryEdges.push_back(h);
					}
					else
					{
						arcHandled[boundaryIt->second] = true;
						isArcInThisConnectedComponent[boundaryIt->second] = true;
					}
					continue;
				}

				auto opph = mesh->opposite_halfedge_handle(h);
				if (mesh->is_boundary(opph))
				{
					openBoundaryEdges.push_back(h);
					continue;
				}
				auto fn = mesh->face_handle(opph);
				if (facesInThisConnectedComponent.find(fn) == facesInThisConnectedComponent.end())
				{
					facesInThisConnectedComponent.insert(fn);
					bfsQueue.push(fn);
				}
			}
		}

		if (facesInThisConnectedComponent.empty())
			continue;

		//Reconstruct the patch sides
		std::vector<std::vector<size_t>> patchSides(maxColor);
		int i = firstArcOnFirstSide;
		do
		{
			if (isArcInThisConnectedComponent[i])
			{
				patchSides[halfarcColors[i]].push_back(halfarcs[i]);
			}

			++i;
			if (i >= halfarcs.size())
				i = 0;
		} while (i != firstArcOnFirstSide);

		//Build the patch
		std::forward<PatchOutput>(patchOutput)(std::move(patchSides), std::move(facesInThisConnectedComponent), std::move(openBoundaryEdges));
	}
}

void MotorcycleGraph::SplitNonRectangularPatches()
{
	nse::util::TimedBlock b("Splitting non-rectangular patches ..");
	std::cout << "Halfarcs before split: " << halfArcs.size() << std::endl;
	int patchCount = patches.size();
	for (int iPatch = 0; iPatch < patchCount; ++iPatch)
	{
		auto& patch = patches[iPatch];
		//holds for each of the two dimensions the paths that split along the dimension (index to motorcyclePaths)
		std::vector<size_t> splitPaths[2];		
		std::vector<float> splitSegmentsMultipliers[2];

		if (patch.PatchSides().size() != 4)
			continue;

		bool patchHasOpenBoundary = !patch.OpenBoundaryEdges().empty();
		if (patchHasOpenBoundary)
			continue;

		std::set<HEMesh::VertexHandle> verticesOnBoundary;
		std::map<HEMesh::VertexHandle, int> verticesOnBoundaryEdge[4];
		
		//Record the vertices on the boundary
		for (int i = 0; i < 4; ++i)
		{
			bool firstSegment = true;
			int nextVertexIdx = 0;
			for (int iArc = 0; iArc < patch.PatchSides()[i].size(); ++iArc)
			{
				auto arcIdx = patch.PatchSides()[i][iArc];
				auto& arc = halfArcs[arcIdx];

				auto h = MotorcycleHalfedge(*arc.begin());
				if (mesh->is_boundary(h) || mesh->is_boundary(mesh->opposite_halfedge_handle(h)))
					patchHasOpenBoundary = true;

				for (auto segment : arc)
				{
					auto h = MotorcycleHalfedge(segment);
					auto v = mesh->from_vertex_handle(h);
					verticesOnBoundary.insert(v);
					if (!firstSegment)
						verticesOnBoundaryEdge[i][v] = nextVertexIdx++;
					firstSegment = false;
				}
			}
		}

		if (patchHasOpenBoundary)
			continue;

		//Spawn motorcycles from every side into the patch
		std::vector<HEMesh::HalfedgeHandle> motorcycleStarts[4];
		std::vector<HEMesh::HalfedgeHandle> edgePaths[4];
		for (int side = 0; side < 4; ++side)
		{
			for (int iArc = 0; iArc < patch.PatchSides()[side].size(); ++iArc)
			{
				auto arcIdx = patch.PatchSides()[side][iArc];
				auto& arc = halfArcs[arcIdx];
				auto lastSegment = arc.LastPathSegment();
				for (auto segment : arc)
				{
					auto h = MotorcycleHalfedge(segment);
					edgePaths[side].push_back(h);
					if (iArc == patch.PatchSides()[side].size() - 1 && segment == lastSegment)
						continue;					
					CirculateForwardUntil<false>(h, *mesh, [&](HEMesh::HalfedgeHandle hInside)
					{
						auto v = mesh->to_vertex_handle(hInside);
						if (verticesOnBoundary.find(v) != verticesOnBoundary.end())
							return true;
						motorcycleStarts[side].push_back(hInside);
						return false;
					});
				}
			}
		}

		int patchSideForDirection[2];		
		std::vector<std::vector<HEMesh::HalfedgeHandle>> motorcyclePaths[2];
		for(int iSide = 0; iSide < 2; ++iSide)
		{
			int side = iSide;
			if (motorcycleStarts[side + 2].size() > motorcycleStarts[side].size())
				side = side + 2;
			patchSideForDirection[iSide] = side;

			auto oppSide = (side + 2) % 4;

			motorcyclePaths[iSide].resize(motorcycleStarts[side].size() + 2);
			motorcyclePaths[iSide].front() = edgePaths[(side + 3) % 4];
			motorcyclePaths[iSide].back() = edgePaths[(side + 1) % 4];

			const float minContinuationScore = 0.8f;
			//Run the motorcycles until they leave the patch
			for (int i = 0; i < motorcycleStarts[side].size(); ++i)
			{
				auto h = motorcycleStarts[side][i];
				motorcyclePaths[iSide][i + 1].push_back(h);
				while (verticesOnBoundary.find(mesh->to_vertex_handle(h)) == verticesOnBoundary.end())
				{
					PathContinuation continuation;
					auto canContinue = FindBestContinuation(*mesh, h, continuation);
					if(canContinue)
						canContinue = continuation.edge.is_valid();
					h = continuation.edge;
					if (!canContinue || continuation.score < minContinuationScore || !mesh->property(isOriginalEdgeProp, mesh->edge_handle(h)))
					{
						motorcyclePaths[iSide][i + 1].clear();
						break;
					}
					motorcyclePaths[iSide][i + 1].push_back(h);
				}				
				if (h.is_valid() && verticesOnBoundaryEdge[oppSide].find(mesh->to_vertex_handle(h)) == verticesOnBoundaryEdge[oppSide].end())
					motorcyclePaths[iSide][i + 1].clear();
			}
			//Remove all motorcycles that could not reach the other side
			motorcyclePaths[iSide].erase(std::remove_if(motorcyclePaths[iSide].begin(), motorcyclePaths[iSide].end(), [](const std::vector<HEMesh::HalfedgeHandle>& path) {return path.empty(); }), motorcyclePaths[iSide].end());
			
			//Use at least 10 motorcycles to achieve some kind of significance
			if (motorcyclePaths[iSide].size() < 10)
				continue;

			{
				//ensure that the motorcycles are in increasing order of the opposite side
				std::vector<int> targetIdx(motorcyclePaths[iSide].size());
				for (int i = 1; i < motorcyclePaths[iSide].size() - 1; ++i)
					targetIdx[i] = edgePaths[oppSide].size() - verticesOnBoundaryEdge[oppSide].at(mesh->to_vertex_handle(motorcyclePaths[iSide][i].back()));
				targetIdx[0] = 0;
				targetIdx.back() = edgePaths[oppSide].size();
				//find the longest increasing subsequence
				//algorithm from Wikipedia https://en.wikipedia.org/w/index.php?title=Longest_increasing_subsequence#Efficient_algorithms
				std::vector<int> P(motorcyclePaths[iSide].size());
				std::vector<int> M(motorcyclePaths[iSide].size());
				int L = 0;
				for (int i = 0; i < motorcyclePaths[iSide].size(); ++i)
				{
					int lo = 1;
					int hi = L;
					while (lo <= hi)
					{
						auto mid = std::ceil((lo + hi) / 2);
						if (targetIdx[M[mid - 1]] < targetIdx[i])
							lo = mid + 1;
						else
							hi = mid - 1;
					}
					int newL = lo;
					if(newL > 1)
						P[i] = M[newL - 2];
					M[newL - 1] = i;
					if (newL > L)
						L = newL;
				}
				std::vector<bool> useMotorcycle(motorcyclePaths[iSide].size(), false);
				int k = M[L - 1];
				for (int i = L - 1; i >= 0; --i)
				{
					useMotorcycle[k] = true;
					k = P[k];
				}
				for (int i = 0; i < motorcyclePaths[iSide].size(); ++i)
					if (!useMotorcycle[i])
						motorcyclePaths[iSide][i].clear();

				//ensure that the motorcycles start at different places on the current patch side
				HEMesh::VertexHandle lastStart;
				for (int i = 0; i < motorcyclePaths[iSide].size(); ++i)
				{
					if (motorcyclePaths[iSide][i].empty())
						continue;
					auto currentStart = mesh->from_vertex_handle(motorcyclePaths[iSide][i].front());
					if (lastStart.is_valid())
					{
						if (currentStart == lastStart)
							motorcyclePaths[iSide][i].clear();
					}
					lastStart = currentStart;
				}
				motorcyclePaths[iSide].erase(std::remove_if(motorcyclePaths[iSide].begin(), motorcyclePaths[iSide].end(), [](const std::vector<HEMesh::HalfedgeHandle>& path) {return path.empty(); }), motorcyclePaths[iSide].end());
			}

			//Measure the paths for all motorcycles
			std::vector<float> pathLengths(motorcyclePaths[iSide].size());			
			for (int i = 0; i < motorcyclePaths[iSide].size(); ++i)
			{
				switch (lengthMeasure)
				{
				case Topological:
					pathLengths[i] = motorcyclePaths[iSide][i].size();
					break;
				case Geometrical:
					pathLengths[i] = 0;
					for (auto h : motorcyclePaths[iSide][i])
						pathLengths[i] += mesh->calc_edge_length(h);
					break;
				}
				//std::cout << pathLengths[i] << std::endl;
			}			

			//Laplacian smoothing on path lengths
			for (int iteration = 0; iteration < 3; ++iteration)
			{
				std::vector<float> smoothedPathLengths(pathLengths.size());
				smoothedPathLengths[0] = pathLengths[0];
				smoothedPathLengths.back() = pathLengths.back();
				for (int i = 1; i < pathLengths.size() - 1; ++i)
				{
					smoothedPathLengths[i] = 0.5f * (pathLengths[i - 1] + pathLengths[i + 1]);
				}
				pathLengths = std::move(smoothedPathLengths);
			}		
			

			PatchSplittingProblem prob(pathLengths);

			const float immediateRoundThreshold = 0.01f;

			while (prob.HasUnfixedMultipliers())
			{
				prob.Optimize();

				//Start by rounding all multipliers that are integer
				for (int i = 0; i < prob.NumberOfMultipliers(); ++i)
				{					
					auto m = prob.GetMultiplier(i);
					auto rounded = std::ldexp(1.0f, (int)std::round(std::log2f(m)));
					if (std::abs(rounded - m) <= immediateRoundThreshold)
						prob.FixMultiplier(i, rounded);
				}

				if (!prob.HasUnfixedMultipliers())
					break;

				//Find the multiplier with the smallest increase in the energy
				float minEnergy = std::numeric_limits<float>::infinity();
				int minMultiplier;
				float minMultiplierFixedTo;
				for (int i = 0; i < prob.NumberOfMultipliers(); ++i)
				{
					if (prob.IsMultiplierFixed(i))
						continue;

					float m = prob.GetMultiplier(i);

					float roundedFloor = std::max(0.25f, std::min(std::ldexp(1.0f, (int)std::floor(std::log2f(m))), 4.0f));
					prob.FixMultiplier(i, roundedFloor);
					float energy = prob.Optimize();
					//std::cout << "Energy of rounding m" << i << " from " << m << " down to " << roundedFloor << ": " << energy << std::endl;
					if (energy < minEnergy)
					{
						minEnergy = energy;
						minMultiplierFixedTo = roundedFloor;
						minMultiplier = i;
					}

					float roundedCeiling = std::max(0.25f, std::min(std::ldexp(1.0f, (int)std::ceil(std::log2f(m))), 4.0f));
					prob.FixMultiplier(i, roundedCeiling);
					energy = prob.Optimize();
					//std::cout << "Energy of rounding m" << i << " from " << m << " up to " << roundedCeiling << ": " << energy << std::endl;
					if (energy < minEnergy)
					{
						minEnergy = energy;
						minMultiplierFixedTo = roundedCeiling;
						minMultiplier = i;
					}
					prob.UnfixMultiplier(i);
				}

				prob.FixMultiplier(minMultiplier, minMultiplierFixedTo);				
			}
			prob.Optimize();			

			//Record the subpatches with equal multipliers
			struct Subpatch
			{
				int startSegmentInclusive;
				int endSegmentExclusive;
				float multiplier;
			};
			std::vector<Subpatch> subpatches;
			subpatches.emplace_back(); 
			subpatches.back().startSegmentInclusive = -1; //the first segment has no multiplier (implicitly -1)
			subpatches.back().multiplier = 1;

			for (int i = 0; i < prob.NumberOfMultipliers(); ++i)
			{
				float m = prob.GetMultiplier(i);
				if (m != subpatches.back().multiplier)
				{
					subpatches.back().endSegmentExclusive = i;
					subpatches.emplace_back();
					subpatches.back().startSegmentInclusive = i;
					subpatches.back().multiplier = m;
				}
			}
			subpatches.back().endSegmentExclusive = prob.NumberOfMultipliers();

			//Eliminate small subpatches
			const int minSubpatchSize = 4;
			for (int i = 0; i < subpatches.size(); ++i)
			{
				auto& subpatch = subpatches[i];
				int subpatchSize = subpatch.endSegmentExclusive - subpatch.startSegmentInclusive;

				auto leftNeighborSubPatchIdx = i - 1;
				while (leftNeighborSubPatchIdx >= 0 && subpatches[leftNeighborSubPatchIdx].startSegmentInclusive == subpatches[leftNeighborSubPatchIdx].endSegmentExclusive)
					--leftNeighborSubPatchIdx;

				if (subpatchSize >= minSubpatchSize)
				{
					//merge patches with equal multipliers
					if (leftNeighborSubPatchIdx >= 0 && subpatches[leftNeighborSubPatchIdx].multiplier == subpatch.multiplier)
					{
						subpatches[leftNeighborSubPatchIdx].endSegmentExclusive = subpatch.endSegmentExclusive;
						subpatch.startSegmentInclusive = subpatch.endSegmentExclusive;
					}
					continue;
				}									

				if (leftNeighborSubPatchIdx < 0)
				{
					subpatches[i + 1].startSegmentInclusive = subpatch.startSegmentInclusive;						
					subpatch.endSegmentExclusive = subpatch.startSegmentInclusive;
					float adaptMultipliers = 1.0f / subpatches[i + 1].multiplier;
					for (int j = 0; j < subpatches.size(); ++j)
					{
						subpatches[j].multiplier *= adaptMultipliers;
						for (int k = std::max(0, subpatches[j].startSegmentInclusive); k < subpatches[j].endSegmentExclusive; ++k)
							prob.FixMultiplier(k, subpatches[j].multiplier);
					}
				}
				else if (i == subpatches.size() - 1)
				{
					subpatches[subpatches.size() - 2].endSegmentExclusive = subpatch.endSegmentExclusive;
					for (int j = subpatch.startSegmentInclusive; j != subpatch.endSegmentExclusive; ++j)
						prob.FixMultiplier(j, subpatches[subpatches.size() - 2].multiplier);
					subpatch.startSegmentInclusive = subpatch.endSegmentExclusive;
				}
				else
				{
					//Eliminate one segment at a time
					while(subpatch.startSegmentInclusive < subpatch.endSegmentExclusive)
					{							
						//try to push left edge to neighbor
						auto& leftNeighborPatch = subpatches[leftNeighborSubPatchIdx];
						prob.FixMultiplier(subpatch.startSegmentInclusive, leftNeighborPatch.multiplier);
						float energyLeft = prob.Optimize();
						prob.FixMultiplier(subpatch.startSegmentInclusive, subpatch.multiplier); //revert

						//try to push right edge to neighbor
						auto& rightNeighborPatch = subpatches[i + 1];
						prob.FixMultiplier(subpatch.endSegmentExclusive - 1, rightNeighborPatch.multiplier);
						float energyRight = prob.Optimize();

						if (energyLeft < energyRight)
						{
							prob.FixMultiplier(subpatch.endSegmentExclusive - 1, subpatch.multiplier); //revert
							prob.FixMultiplier(subpatch.startSegmentInclusive, leftNeighborPatch.multiplier);
							++leftNeighborPatch.endSegmentExclusive;
							++subpatch.startSegmentInclusive;
						}
						else
						{
							--rightNeighborPatch.startSegmentInclusive;
							--subpatch.endSegmentExclusive;
						}
					}
				}
			}
			//remove empty subpatches
			subpatches.erase(std::remove_if(subpatches.begin(), subpatches.end(), [](const Subpatch& subpatch) { return subpatch.startSegmentInclusive == subpatch.endSegmentExclusive; }), subpatches.end());

			splitSegmentsMultipliers[iSide].push_back(1);
			for (int i = 1; i < subpatches.size(); ++i)
			{
				splitPaths[iSide].push_back(subpatches[i].startSegmentInclusive);
				splitSegmentsMultipliers[iSide].push_back(subpatches[i].multiplier);
			}			
		} //for each of the two sides

		//Now realize the splits
		std::map<HEMesh::FaceHandle, size_t> faceToPatchIdx;
		if (!splitPaths[0].empty() || !splitPaths[1].empty())
		{
			std::vector<std::pair<size_t, size_t>> arcMerges;

			size_t motorcyclesStartSide0 = motorcycles.size();

			std::deque<size_t> halfarcsPerPatchSide[4];

			//record the vertices where motorcycles of side 0 pass through
			ManifoldnessAwareVertexMap<std::pair<int, int>> vertexVisitedByMotorcycleOfSide0(*mesh); //pair: (motorcycle, position)
			ManifoldnessAwareVertexMap<char> vertexVisitedByMotorcycleOfSide1(*mesh); //only on the patch edges, value if ignored
			for (int i = 0; i < splitPaths[0].size(); ++i)
			{
				auto& path = motorcyclePaths[0][splitPaths[0][i]];
				auto h = mesh->opposite_halfedge_handle(path.front());
				vertexVisitedByMotorcycleOfSide0.AccessOrCreateAtToVertex(h) = std::make_pair(i, 0);
				motorcycles.emplace_back();
				motorcycles.back().multiplier = splitSegmentsMultipliers[0][i + 1] / splitSegmentsMultipliers[0][i];
				for (int iSegment = 0; iSegment < path.size(); ++iSegment)
				{
					auto segment = path[iSegment];
#if __DEBUG
					std::pair<int, int>* entry;
					if (vertexVisitedByMotorcycleOfSide0.TryAccessAtToVertex(segment, entry))
						throw std::runtime_error("The cutting paths of the patch intersect.");
#endif
					vertexVisitedByMotorcycleOfSide0.AccessOrCreateAtToVertex(segment) = std::make_pair(motorcycles.size() - 1, iSegment + 1);
					motorcycles.back().AddToPath(segment, *mesh);
				}
			}

			size_t motorcyclesStartSide1 = motorcycles.size();

			//find the intersection point between motorcycles from side 0 and motorcycles of side 1

			for (int i = 0; i < splitPaths[1].size(); ++i)
			{
				auto& path = motorcyclePaths[1][splitPaths[1][i]];
				motorcycles.emplace_back();
				motorcycles.back().multiplier = splitSegmentsMultipliers[1][i + 1] / splitSegmentsMultipliers[1][i];
				vertexVisitedByMotorcycleOfSide1.AccessOrCreateAtToVertex(mesh->opposite_halfedge_handle(path.front()));
				vertexVisitedByMotorcycleOfSide1.AccessOrCreateAtToVertex(path.back());
				for (int iSegment = 0; iSegment < path.size(); ++iSegment)
				{
					auto h = path[iSegment];
					motorcycles.back().AddToPath(h, *mesh);
					const std::pair<int, int>* entry;
					if (vertexVisitedByMotorcycleOfSide0.TryAccessAtToVertex(h, entry))
					{
						motorcycles.back().crashesOnPath.insert(motorcycles.back().Path().size());						
						motorcycles[entry->first].crashesOnPath.insert(entry->second);
					}
				}

				auto mIdx = motorcycles.size() - 1;

				auto oldArcCount = halfArcs.size();
				ExtractArcsForMotorcycle(motorcycles.size() - 1, arcMerges, false);
				auto forwardEdge = (patchSideForDirection[1] == 1 ? 2 : 0);
				auto backwardEdge = (patchSideForDirection[1] == 1 ? 0 : 2);
				for (int i = oldArcCount; i < halfArcs.size(); i += 2)
				{
					halfarcsPerPatchSide[forwardEdge].push_back(i);
					halfarcsPerPatchSide[backwardEdge].push_front(i + 1);					
				}
			}

			//Add the side0 motorcycles
			for (int i = 0; i < splitPaths[0].size(); ++i)
			{
				auto oldArcCount = halfArcs.size();
				ExtractArcsForMotorcycle(motorcyclesStartSide0 + i, arcMerges, false);				
				auto forwardEdge = (patchSideForDirection[0] == 0 ? 1 : 3);
				auto backwardEdge = (patchSideForDirection[0] == 0 ? 3 : 1);
				for (int i = oldArcCount; i < halfArcs.size(); i += 2)
				{
					halfarcsPerPatchSide[forwardEdge].push_back(i);
					halfarcsPerPatchSide[backwardEdge].push_front(i + 1);
				}
			}

			//split the halfarcs of the patch edges if necessary
			for (int iSide = 0; iSide < 4; ++iSide)
			{
				for (int iArc = 0; iArc < patch.PatchSides()[iSide].size(); ++iArc)
				{
					auto arcIdx = patch.PatchSides()[iSide][iArc];
					halfarcsPerPatchSide[iSide].push_back(arcIdx);
					auto it = halfArcs[arcIdx].begin();
					while (it != halfArcs[patch.PatchSides()[iSide][iArc]].end())
					{
						auto h = MotorcycleHalfedge(*it);
						bool targetHasMotorcycle = false;
						if (iSide % 2 == 0)
						{
							const std::pair<int, int>* unused;
							targetHasMotorcycle = vertexVisitedByMotorcycleOfSide0.TryAccessAtToVertex(h, unused);
						}
						else
						{
							const char* unused;
							targetHasMotorcycle = vertexVisitedByMotorcycleOfSide1.TryAccessAtToVertex(h, unused);
						}

						if (targetHasMotorcycle)
						{
							SplitArcAfterSegment(patch.PatchSides()[iSide][iArc], it);
							break;
						}
						++it;
					}
				}
			}

			//filll the vector of all halfarcs
			std::vector<size_t> arcs;
			arcs.reserve(halfarcsPerPatchSide[0].size() + halfarcsPerPatchSide[1].size() + halfarcsPerPatchSide[2].size() + halfarcsPerPatchSide[3].size());
			std::vector<int> arcColors;
			arcColors.reserve(halfarcsPerPatchSide[0].size() + halfarcsPerPatchSide[1].size() + halfarcsPerPatchSide[2].size() + halfarcsPerPatchSide[3].size());

			for (int iSide = 0; iSide < 4; ++iSide)
			{
				for (auto arc : halfarcsPerPatchSide[iSide])
				{
					arcs.push_back(arc);
					arcColors.push_back(iSide);
				}
			}

			bool usedThisPatch = false;
			std::vector<HEMesh::HalfedgeHandle> empty;
			BuildPatchesFromArcs(arcs, arcColors, 0, 4, empty,
				[&](std::vector<std::vector<size_t>>&& patchSides, std::set<OpenMesh::FaceHandle>&& faces, std::vector<HEMesh::HalfedgeHandle>&& openBoundary)
			{			
				int patchIdx;
				if (!usedThisPatch)
				{
					patchIdx = iPatch;
					usedThisPatch = true;
				}
				else
				{
					patches.emplace_back(*this);
					patchIdx = patches.size() - 1;
				}
				for (auto f : faces)
					faceToPatchIdx[f] = patchIdx;
				patches[patchIdx].PrepareBuild(std::move(patchSides), std::move(faces), std::move(openBoundary));
				patches[patchIdx].Build(*mesh);
				AssignFaceIndexToHalfarcs(patchIdx);

			});						
		} //if we need to split
		else
		{
			for (auto f : patch.Faces())
				faceToPatchIdx[f] = iPatch;
		}

		//Record interior patch lengths
		for (int side = 0; side < 2; ++side)
		{
			for (auto& path : motorcyclePaths[side])
			{
				std::map<size_t, float> patchToLength;
				for (auto h : path)
				{
					float edgeLength;
					switch (lengthMeasure)
					{
					case Topological:
						edgeLength = 1;
						break;
					case Geometrical:
						edgeLength = mesh->calc_edge_length(h);
						break;
					}

					auto f = mesh->face_handle(h);
					size_t fPatch = -1;
					if (f.is_valid())
					{
						auto patchIt = faceToPatchIdx.find(f);
						if (patchIt != faceToPatchIdx.end())
						{
							patchToLength[patchIt->second] += edgeLength;
							fPatch = patchIt->second;
						}
					}

					auto oppf = mesh->opposite_face_handle(h);
					if (oppf.is_valid())
					{
						auto patchIt = faceToPatchIdx.find(oppf);
						if (patchIt != faceToPatchIdx.end() && patchIt->second != fPatch)
							patchToLength[patchIt->second] += edgeLength;
					}
				}

				for (auto& entry : patchToLength)
					patches[entry.first].AddInteriorPatchLength(side, entry.second);
			}			
		}

	} //for every patch
	std::cout << "Halfarcs after split: " << halfArcs.size() << std::endl;

	PatchesChanged();
}

void MotorcycleGraph::SplitArcAfterSegment(size_t iArc, HalfArc::PathSegmentIterator it)
{
	auto& arc = halfArcs[iArc];
	auto nextIt = it; ++nextIt;
	if (nextIt == arc.end())
		return; //there is nothing to split

	auto oppositeLocation = *it;
	oppositeLocation.goForward = !oppositeLocation.goForward;
	auto oppositeArcIdx = OppositeHalfarc(iArc);
	auto& oppositeArc = halfArcs[oppositeArcIdx];
	auto oppositeLocationIt = std::find(oppositeArc.begin(), oppositeArc.end(), oppositeLocation);

	//Create two new halfarcs and copy the properties
	HalfArc newArc(nextIt, arc.end()); 
	newArc.face = arc.face;
	newArc.multiplier = arc.multiplier;
	HalfArc newOppositeArc(oppositeArc.begin(), oppositeLocationIt);
	newOppositeArc.face = oppositeArc.face;
	newOppositeArc.multiplier = oppositeArc.multiplier;
	arc.SetEndExclusive(nextIt);
	oppositeArc.SetStartInclusive(oppositeLocationIt);

	if (arc.face != (size_t)-1)
		patches[arc.face].InsertArcAfter(halfArcs.size(), iArc);
	if (oppositeArc.face != (size_t)-1)
		patches[oppositeArc.face].InsertArcBefore(halfArcs.size() + 1, oppositeArcIdx);

	halfArcs.push_back(std::move(newArc));
	halfArcs.push_back(std::move(newOppositeArc));
}

void MotorcycleGraph::LoadExternal(std::vector<TexturePatch>&& patches, std::vector<HalfArc>&& halfarcs)
{
	this->patches = patches;
	this->halfArcs = halfarcs;

	for (auto& patch : this->patches)
		patch.Build(*mesh);

	PatchesChanged();
}

void MotorcycleGraph::ExtractPatches(ExtractionStatistics& statistics)
{
	nse::util::TimedBlock b("Extracting patches from motorcycles ..");

	patches.clear();

	if (!activeMotorcycles.empty())
	{
		std::cout << "The motorcycle graph still has active motorcycles. Cannot extract patches." << std::endl;
		return;
	}	
	
	ExtractArcs();

	std::vector<bool> halfArcUsed(halfArcs.size(), false);	

	bool verbose = false;
	
	struct PatchExtractionInfo
	{
		std::vector<size_t> halfArcsInPatch;
		std::vector<int> halfArcInPatchColors;
		int orientation;
		std::vector<HEMesh::HalfedgeHandle> additionalBoundaries;
	};
	std::vector<PatchExtractionInfo> extractionInfos;

	struct ColorZeroPatchInfo
	{
		int extractionInfoId;
		std::vector<size_t> skippedMotorcycles;

		ColorZeroPatchInfo(int extractionInfoId, std::vector<size_t>&& skippedMotorcycles)
			: extractionInfoId(extractionInfoId), skippedMotorcycles(skippedMotorcycles)
		{ }
	};
	std::vector<ColorZeroPatchInfo> colorZeroPatches;
	
	statistics.totalFaces = mesh->n_faces();

	//Extract patches by following the arcs around until they reach their starting point
	for(int iArc = 0; iArc < halfArcs.size(); ++iArc)
	{
		auto& arc = halfArcs[iArc];

		auto currentArc = iArc;		
				
		if (halfArcUsed[iArc])
			continue;

		auto firstHalfedge = MotorcycleHalfedge(arc.segments.front().location);
		if (mesh->is_boundary(firstHalfedge))
			continue;

		if (verbose)
		{
			auto startVertex = (arc.segments.front().location.goForward ?
				mesh->from_vertex_handle(motorcycles[arc.segments.front().location.motorcycle].Path()[arc.segments.front().location.pathSegment])
				: mesh->to_vertex_handle(motorcycles[arc.segments.front().location.motorcycle].Path()[arc.segments.front().location.pathSegment]));
			std::cout << "Start at (" << arc.segments.front().location.motorcycle << ", " << arc.segments.front().location.pathSegment
				<< ", " << arc.segments.front().location.goForward << "): vertex " << startVertex << " at " << mesh->point(startVertex)
				<< " going in direction " << mesh->calc_edge_vector(MotorcycleHalfedge(arc.segments.front().location)) << std::endl;
		}

		//We are creating a new patch
		extractionInfos.emplace_back();

		extractionInfos.back().orientation = 0;				

		std::vector<size_t> skippedMotorcycles;

		while (!halfArcUsed[currentArc])
		{
			extractionInfos.back().halfArcsInPatch.push_back(currentArc);
			extractionInfos.back().halfArcInPatchColors.push_back(extractionInfos.back().orientation);
			halfArcUsed[currentArc] = true;				

			//find the next arc
			auto lastLocation = halfArcs[currentArc].LastPathSegment();
			auto lastH = MotorcycleHalfedge(lastLocation);

			const VisitedEntry* visited;
			if (!visitedVertices.TryAccessAtToVertex(lastH, visited))
			{
				std::cout << "Cannot find visited vertex at " << mesh->point(mesh->to_vertex_handle(lastH)) << std::endl;
				return;
			}

			int colorChange;			
			auto newLocation = NextEdgeInPatch(lastLocation, colorChange, skippedMotorcycles, verbose);
			auto newH = MotorcycleHalfedge(newLocation);
			try
			{
				currentArc = firstHalfedgeToHalfArc.at(newH);
			}
			catch (...)
			{
				std::cout << "Halfedge does not correspond to a halfarc: ";
				PrintFullHalfedge(std::cout, newH, *mesh);
				std::cout << std::endl;
				return;
			}
			
			extractionInfos.back().orientation += colorChange;

			if (verbose)
			{
				std::cout << "Going over to motorcycle " << newLocation.motorcycle << " at segment " << newLocation.pathSegment << ", forward = " << newLocation.goForward << ", new color = " << extractionInfos.back().orientation << std::endl;
			}
		}

		if (verbose)
		{
			auto endLocation = halfArcs[currentArc].LastPathSegment();
			auto endVertex = (endLocation.goForward
				? mesh->from_vertex_handle(motorcycles[endLocation.motorcycle].Path()[endLocation.pathSegment]) :
				mesh->to_vertex_handle(motorcycles[endLocation.motorcycle].Path()[endLocation.pathSegment]));
			std::cout << "End after " << extractionInfos.back().orientation << " corners in total." << std::endl;
		}
		auto endHalfedge = MotorcycleHalfedge(halfArcs[currentArc].segments.front().location);
		
		if (endHalfedge != firstHalfedge || extractionInfos.back().orientation == 0)
		{			
			if (extractionInfos.back().orientation == 0)
			{
				std::cout << "Color-0 patch." << std::endl;
				colorZeroPatches.emplace_back(extractionInfos.size() - 1, std::move(skippedMotorcycles));
			}
			else
			{				
				std::cout << "Open patch." << std::endl;
				++statistics.openPatches;

				HEMesh::HalfedgeHandle lastH;
				for (auto arcIdx : extractionInfos.back().halfArcsInPatch)
					for (auto p : halfArcs[arcIdx])
					{
						lastH = MotorcycleHalfedge(p);
						PrintHalfedge(std::cout, lastH, *mesh);
					}
				std::cout << "(" << mesh->to_vertex_handle(lastH).idx() << ")" << std::endl;
				extractionInfos.pop_back();
			}					
		}
		else
		{
			//valid patch
			statistics.maxPatchDegree = std::max(statistics.maxPatchDegree, extractionInfos.back().orientation);
		}
	}

	if (!colorZeroPatches.empty())
	{		
		{
			nse::util::TimedBlock b("Reactivating motorcycles to avoid cylindrical patches ..");
			//There are cylindrical patches; we need to re-activate some motorcycles
			std::map<size_t, size_t> motorcycleToColorZeroPatchId;
			for (int i = 0; i < colorZeroPatches.size(); ++i)
			{
				auto& patch = colorZeroPatches[i];
				for (auto arcIdx : extractionInfos[patch.extractionInfoId].halfArcsInPatch)
					for (auto& location : halfArcs[arcIdx].segments)
						motorcycleToColorZeroPatchId[location.location.motorcycle] = i;
			}

			std::vector<bool> colorZeroPatchHandled(colorZeroPatches.size(), false);
			CirculateUntilMotorcycleCondition stopCondition(*this);
			for (int i = 0; i < colorZeroPatches.size(); ++i)
			{
				if (colorZeroPatchHandled[i])
					continue;
				colorZeroPatchHandled[i] = true;

				auto& patch = colorZeroPatches[i];

				std::queue<size_t> reactivateMotorcycleQueue;

				//greedily reactivate the first motorcycle
				reactivateMotorcycleQueue.push(patch.skippedMotorcycles.front());

				while (!reactivateMotorcycleQueue.empty())
				{
					auto motorcycleIdx = reactivateMotorcycleQueue.front();
					reactivateMotorcycleQueue.pop();
					auto& motorcycle = motorcycles[motorcycleIdx];
					if (!motorcycle.isDeactivated)
						continue;
					motorcycle.isDeactivated = false;

					if (motorcycle.straightContinuationStart != (size_t)-1)
						reactivateMotorcycleQueue.push(motorcycle.straightContinuationStart);
					if (motorcycle.straightContinuationEnd != (size_t)-1)
						reactivateMotorcycleQueue.push(motorcycle.straightContinuationEnd);

					const size_t* unused;

					//Check motorcycle target
					auto h = motorcycle.Path().back();
					if (motorcycle.straightContinuationEnd == (size_t)-1 && !motorcycle.collideWithPaddedBoundary && !vertexToMetaSingularityIndex.TryAccessAtToVertex(h, unused))
					{
						//There is no straight continuation and the target collision is not a singularity 
						//(which can not happen in the current formulation, but I'll leave it here anyway)					
						//Hence, the motorcycle must collide with another, reactivate this one, too.
						stopCondition.motorcyclesAtQueryTargetVertex = &visitedVertices.AccessAtToVertex(h);
						CirculateForwardUntil<true>(h, *mesh, stopCondition);
						reactivateMotorcycleQueue.push(stopCondition.result.motorcycle);
						
						auto otherPatchIt = motorcycleToColorZeroPatchId.find(stopCondition.result.motorcycle);
						if (otherPatchIt != motorcycleToColorZeroPatchId.end())
							colorZeroPatchHandled[otherPatchIt->second] = true;
					}
				}
			}
		}

		//try again
		ExtractPatches(statistics);
		return;
	}

	//do a BFS from each arc to find bridges between patches
	{
		nse::util::TimedBlock b("Finding bridges ..");
		struct PerFaceInfo
		{
			size_t belongsToExtractionUnit;
			PerFaceInfo() { }
			PerFaceInfo(size_t belongsToExtractionUnit)
				: belongsToExtractionUnit(belongsToExtractionUnit)
			{ }
		};
		std::map<HEMesh::FaceHandle, PerFaceInfo> perFaceInfos;
		std::set<HEMesh::HalfedgeHandle> boundaries;

		//initialize the BFS
		std::queue<HEMesh::FaceHandle> bfsQueue;
		for (int iExtractionUnit = 0; iExtractionUnit < extractionInfos.size(); ++iExtractionUnit)
		{
			auto& extractionUnit = extractionInfos[iExtractionUnit];
			for (size_t arcIdx : extractionUnit.halfArcsInPatch)
			{
				for (auto p : halfArcs[arcIdx])
				{
					auto h = MotorcycleHalfedge(p);
					boundaries.insert(h);
					auto f = mesh->face_handle(h);
					if (f.is_valid())
					{
						perFaceInfos[f].belongsToExtractionUnit = iExtractionUnit;
						bfsQueue.push(f);
					}
				}
			}
		}

		//do the BFS
		while (!bfsQueue.empty())
		{
			auto f = bfsQueue.front();
			bfsQueue.pop();

			auto& fInfo = perFaceInfos.at(f);

			for (auto h : mesh->fh_range(f))
			{
				if (boundaries.find(h) != boundaries.end())
					continue;
				auto oppF = mesh->opposite_face_handle(h);
				if (!oppF.is_valid())
					continue;
				auto inserted = perFaceInfos.emplace(oppF, fInfo);
				if (inserted.second)
					bfsQueue.push(oppF);
				else
				{
					//there was already an entry for this face
					if (inserted.first->second.belongsToExtractionUnit != fInfo.belongsToExtractionUnit)
					{
						//two different patches meet here, we are on a bridge
						extractionInfos[fInfo.belongsToExtractionUnit].additionalBoundaries.push_back(h);
						extractionInfos[inserted.first->second.belongsToExtractionUnit].additionalBoundaries.push_back(mesh->opposite_halfedge_handle(h));
					}
				}
			}
		}
	}

	//extract the actual patches

	statistics.patchCountPerNumberOfCorners.resize(statistics.maxPatchDegree + 1);
	statistics.patchSizesPerNumberOfCorners.resize(statistics.maxPatchDegree + 1);

	{
		nse::util::TimedBlock b("Building patches ..");
#pragma omp parallel for
		for (int i = 0; i < extractionInfos.size(); ++i)
		{
			auto& info = extractionInfos[i];

			int firstArcOnFirstSide = -1;
			for (int i = 0; i < info.halfArcsInPatch.size(); ++i)
			{
				auto& arc = halfArcs[info.halfArcsInPatch[i]];
				if (info.halfArcInPatchColors[i] == info.orientation && firstArcOnFirstSide == -1)
					firstArcOnFirstSide = i;
				info.halfArcInPatchColors[i] = info.halfArcInPatchColors[i] % info.orientation;
			}
			if (firstArcOnFirstSide == -1)
				firstArcOnFirstSide = 0;

			BuildPatchesFromArcs(info.halfArcsInPatch, info.halfArcInPatchColors, firstArcOnFirstSide, info.orientation, info.additionalBoundaries,
				[&](std::vector<std::vector<size_t>>&& patchSides, std::set<OpenMesh::FaceHandle>&& faces, std::vector<HEMesh::HalfedgeHandle>&& openBoundary)
			{
				int patchIdx;
#pragma omp critical
				{
					statistics.patchSizes.push_back(faces.size());
					patches.emplace_back(*this);
					patchIdx = patches.size() - 1;

					patches[patchIdx].PrepareBuild(std::move(patchSides), std::move(faces), std::move(openBoundary));
				}
				nse::data::atomicAdd(&statistics.patchCountPerNumberOfCorners[info.orientation], 1);
				nse::data::atomicAdd(&statistics.patchSizesPerNumberOfCorners[info.orientation], patches[patchIdx].Faces().size());
			});
		}
	}

	{
		nse::util::TimedBlock b("Finalizing patches ..");
#pragma omp parallel for
		for (int i = 0; i < patches.size(); ++i)
		{
			patches[i].Build(*mesh);
			AssignFaceIndexToHalfarcs(i);
		}
	}	

	PatchesChanged();
}

void MotorcycleGraph::OccupyEdge(HEMesh::EdgeHandle e, size_t motorcycle)
{
	auto h = mesh->halfedge_handle(e, 0);
	occupiedEdges[e] = motorcycle;
}

// -----------  CirculateUntilMotorcycleCondition  -----------


bool MotorcycleGraph::CirculateUntilMotorcycleCondition::operator()(HEMesh::HalfedgeHandle h)
{
	if (motorcyclesAtQueryTargetVertex == nullptr)
		return false;
	result = graph.EdgeToMotorcycleLocation(h, *motorcyclesAtQueryTargetVertex);
	if (result.motorcycle != (size_t)-1)
		return true;
	return false;
}

// -----------  StoppingMotorcycles  -----------

MotorcycleGraph::StoppingMotorcycles::StoppingMotorcycles(const HEMesh& mesh, HEMesh::HalfedgeHandle stopBeforeEdge, OpenMesh::EPropHandleT<bool> isOriginalEdgeProp)
	: stopBeforeEdge(stopBeforeEdge)
{
	if (!stopBeforeEdge.is_valid())
		return;

	useNonOriginalEdge = false;
	if (!mesh.is_boundary(stopBeforeEdge))
		newSingularityH[0] = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(stopBeforeEdge));
	else
		newSingularityH[0] = HEMesh::HalfedgeHandle(-1);

	if (!mesh.is_boundary(mesh.opposite_halfedge_handle(stopBeforeEdge)))
		newSingularityH[1] = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(stopBeforeEdge));
	else
		newSingularityH[1] = HEMesh::HalfedgeHandle(-1);

	for (int i = 0; i < 2; ++i)
	{
		if (!newSingularityH[i].is_valid())
			continue;
		if (!mesh.property(isOriginalEdgeProp, mesh.edge_handle(newSingularityH[i])))
			useNonOriginalEdge = true;
	}
}

void MotorcycleGraph::StoppingMotorcycles::TryToRoute(MotorcycleGraph& graph, bool tryAlternative)
{
	hasValidRoute = true;
	for (int i = 0; i < 2; ++i)
	{
		if (!newSingularityH[i].is_valid())
			continue;
		auto enteringIt = graph.patchSet->enteringEdgeToPatch.find(newSingularityH[i]);
		if (enteringIt != graph.patchSet->enteringEdgeToPatch.end())
		{
			//the new motorcycle is in a patch again
			newSingularityHRoutes[i] = graph.RouteThroughPatch(enteringIt->first, graph.patchSet->patches[enteringIt->second.patchIdx], *enteringIt);
			if (newSingularityHRoutes[i].type == NoRoute)
			{
				//try with the original direction
				if (tryAlternative)
					newSingularityHRoutes[i] = graph.RouteThroughPatch(stopBeforeEdge, graph.patchSet->patches[enteringIt->second.patchIdx], *enteringIt);
				if (!tryAlternative || newSingularityHRoutes[i].type == NoRoute)
				{
					failedStoppingMotorcycle = i;
					hasValidRoute = false;
				}
			}
		}
		else
			newSingularityHRoutes[i].type = NoRoute;
	}
}

void MotorcycleGraph::StoppingMotorcycles::Realize(MotorcycleGraph& graph, bool highPriority)
{
	if (newSingularityH[0].is_valid() && newSingularityH[1].is_valid())
	{
		auto idx = graph.AddMotorcyclePair({ { newSingularityH[0], newSingularityH[1] } }, highPriority);
		for (int i = 0; i < 2; ++i)
			graph.RealizeRoute(idx + i, newSingularityHRoutes[i]);
	}
	else
	{
		for (int i = 0; i < 2; ++i)
			if (newSingularityH[i].is_valid())
			{
				auto idx = graph.AddMotorcycle(newSingularityH[i], highPriority);
				graph.motorcycles[idx].startInPaddedBoundary = true;
				graph.RealizeRoute(idx, newSingularityHRoutes[i]);
			}
	}
}