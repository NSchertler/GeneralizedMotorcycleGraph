#pragma once

#include "common.h"

//Represents a singular vertex
struct SingularityInfo
{
	SingularityInfo(HEMesh::VertexHandle v, int valenceDefect)
		: vertexHandle(v), valenceDefect(valenceDefect), correspondingRegion((size_t)-1)
	{ }

	SingularityInfo(HEMesh::VertexHandle v, HEMesh::HalfedgeHandle contextOutgoingBoundary, int valenceDefect)
		: vertexHandle(v), contextOutgoingBoundaryEdge(contextOutgoingBoundary), valenceDefect(valenceDefect), correspondingRegion((size_t)-1)
	{ }

	//The vertex
	HEMesh::VertexHandle vertexHandle;

	//Only used for non-manifold vertices; the first outgoing edge in the respective surface patch.
	//This is a boundary halfedge
	HEMesh::HalfedgeHandle contextOutgoingBoundaryEdge;

	//The valence defect of the singularity (valence - 4)
	int valenceDefect;

	//index of the fenced region that this singularity belongs to or (size_t)-1
	size_t correspondingRegion;
};

//Represents a vertex that acts as a singularity (either an original singularity or an arbitrary vertex inside a fenced region)
class MetaSingularity
{
public:	
	//Priority of motorcycles that are generated from this meta singularity
	bool highPriority;

	//Returns the path of the motorcycle with the given orientation
	const std::vector<HEMesh::HalfedgeHandle>& GetEmanatingMotorcycle(int orientation) const;

	//Returns the path of the motorcycle with the given orientation. If the list of emanating
	//motorcycles is too short, it is extended.
	std::vector<HEMesh::HalfedgeHandle>& GetEmanatingMotorcycle(int orientation);

	//Sets if a given emanating motorcycle is active.
	void SetEmanatingMotorcycleActive(int orientation, bool active);

	//Returns if a given emanating motorcycle is active.
	bool GetEmanatingMotorcycleActive(int orientation) const;

	//Returns the orientation of the emanating motorcycle whose first halfedge is h
	//or (size_t)-1 if there is no such motorcycle.
	size_t FindOutgoingMotorcycle(HEMesh::HalfedgeHandle h) const;

	//Returns the number of emanating motorcycles.
	size_t Degree() const;

	//Adds a new emanating motorcycle at the end of this meta singularity.
	void AddMotorcycle(HEMesh::HalfedgeHandle direction, bool active);

	//Adds a new virtual motorcycle with no path at the end of this meta singularity. This
	//is only used to ensure correct ordering.
	void AddEmptyMotorcycle();

	//Adds a new emanating motorcycle at the beginning of this meta singularity.
	void AddMotorcycleAtFront(HEMesh::HalfedgeHandle direction, bool active);

private:
	//List of motorcycles that emanate from this meta singularity, one for each orientation.
	//Can be empty for a single edge orientation if the singuarity is on the boundary.
	//Sorted in counter-clockwise order. The second pair item determines if the motorcycle
	//is active.
	std::vector<std::pair<std::vector<HEMesh::HalfedgeHandle>, bool>> emanatingMotorcycles;
};
