#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>

#include "common.h"
#include "FencedRegion.h"
#include "PatchSet.h"
#include "MotorcycleGraph.h"

#include "Statistics.h"

#include <nsessentials/math/BoundingBox.h>
#include <nsessentials/util/UnionFind.h>

const Eigen::Vector4f valenceColors[] =
{
	Eigen::Vector4f(0.843137255f,	0.188235294f,	0.152941176f, 1.0f),
	Eigen::Vector4f(0.988235294f,	0.552941176f,	0.349019608f, 1.0f),
	Eigen::Vector4f(0.996078431f,	0.878431373f,	0.564705882f, 1.0f),
	Eigen::Vector4f(1.0f,			1.0f,			1.0f,		  1.0f),
	Eigen::Vector4f(0.878431373f,	0.952941176f,	0.97254902f , 1.0f),
	Eigen::Vector4f(0.568627451f,	0.749019608f,	0.858823529f, 1.0f),
	Eigen::Vector4f(0.270588235f,	0.458823529f,	0.705882353f, 1.0f),
};

//Returns the color for a given valence defect. A valence defect of 0 will produce a 
//neutral color.
static const Eigen::Vector4f GetValenceColor(int valenceDefect)
{
	size_t colorIndex = std::max<size_t>(0, std::min<size_t>(6, valenceDefect + 3));
	return valenceColors[colorIndex];
}

const int segmentColorCount = 9;
const int segmentColors[segmentColorCount][3] =
{
	{ 166, 206, 227 },
	{ 31,120,180 },
	{ 251,154,153 },
	/*{ 227,26,28 },*/
	{ 253,191,111 },
	{ 255,127,0 },
	{ 202,178,214 },
	{ 106,61,154 },
	{ 255,255,153 },
	{ 177,89,40 }
};

//Returns an arbitrary color for a given segment index. Neighboring segment
//indices will produce different colors.
static const Eigen::Vector4f GetSegmentColor(int segment)
{
	return Eigen::Vector4f(segmentColors[segment % segmentColorCount][0] / 255.0f, segmentColors[segment % segmentColorCount][1] / 255.0f, segmentColors[segment % segmentColorCount][2] / 255.0f, 1.0f);
}

//Represents the result of singularity classification
enum ClassificationResult
{
	//A non-degenerate regular fenced region could be found.
	RectilinearPatch,

	//A non-degenerate irregular fenced region could be found.
	NonRectilinearPatch,

	//No non-degenerate fenced region could be found.
	NoPatch
};


//Deprecated
enum MultiplierStrategy
{
	AllInvalid,
	IterativeRounding,
};

//Represents the strategy to calculate arc lengths during parametrization
enum ArclengthStrategy
{
	//ArclengthStrategySimple
	Simple,

	//ArclengthStrategyGurobiSeparate
	GurobiSeparate,

	//ArclengthStrategyGlobalGurobi
	GurobiGlobal,
};

//Holds all relevant data for calculating a motorcycle graph and parametrizing the patches.
class Data
{
public:
	Data();	

	//Resets the data object
	void Clear();

	//Loads a mesh from a file and prepares it for processing.
	//filename - the file to load
	//mergeTriangulatedQuads - try to merge neighboring triangles in the mesh if their union becomes near rectangular
	//mergeAngleTreshold - angle threshold for triangle merging (deviation from 90°) in radians
	void LoadMesh(const std::string & filename, bool mergeTriangulatedQuads, float mergeAngleThreshold);

	//Loads a mesh for the sole purpose of visualization. No processing can be performed on this mesh.
	//filename - the file to load
	//scaleTexCoords - determines if the read texture coordinates should be scaled up to match pixel coordinates
	void LoadMeshForVisualization(const std::string& filename, bool scaleTexCoords);

	//Saves the geometry of the current mesh to a PLY file.
	void SaveMeshPLY(const std::string& file) const;

	//Saves the geometry and parametrization of the current mesh to an OBJ file.
	//file - target file name
	//unsubdivide - removes the Catmull Clark subdivision where possible (if the original mesh was not a pure quad mesh)
	void SaveMeshOBJ(const std::string& file, bool unsubdivide = true) const;

	//Finds irregular vertices in the mesh. Also calls FindMetaSingularities() at the end.
	void FindSingularities();

	//Finds meta singularities based on the current set of irregular vertices and fenced regions. A meta
	//singularity is either an irregular vertex with no fenced region or an arbitrary vertex within a 
	//fenced region.
	void FindMetaSingularities();

	//Tries to grow minimal fenced region around every singularity. 
	void ClassifyAllSingularities();

	//Tries to merge neighboring fenced region in order to decrease the number of spawned motorcycles.
	void TryToMergePatches();

	//Calculates the motorcycle graph for the current set of meta singularities.
	void CalculateMotorcycleGraph(bool deactivateUnnecessaryMotorcycles = true);

	//Extracts rectangular patches from the current motorcycle graph.
	MotorcycleGraph::ExtractionStatistics ExtractPatches(bool splitNonRectangularPatches = true);
	
	//Calculates a parametrization for the current set of rectangular patches.
	//targetLengthPerEdge - the average parametric target length of edges; used to determine a global scaling 
	//                      factor between geometric and parametric domain.
	//parametrizationErrorThreshold - determines the upper threshold of the parametrization energy above which 
	//                                parametrization is re-attempted after relaxing the system. Details depend 
	//                                on the chosen parametrization strategy
	//discreteOptimizationTimeLimit - the time limit for discrete optimization steps if there are any.
	//arclengthStrategy - the strategy to use to determine parametric arc lengths.
	void CalculateParametrization(float targetLengthPerEdge, float parametrizationErrorThreshold, const std::chrono::steady_clock::duration& discreteOptimizationTimeLimit, ArclengthStrategy arclengthStrategy);
	
	//Packs the current parametrization into a rectangular domain for a given MIP level factor (positive power of 2).
	//Texture coordinates will be scaled to the [0, 1] range and the members textureWidth and textureHeight will be set.
	void PackTexture(int mipFactor);
	
	//Tries to grow a minimal fenced region around a single singularity.
	//v - the singularity to examine
	//patch - the resulting fenced region
	ClassificationResult ClassifySingularity(const SingularityInfo& v, FencedRegion& patch, bool verbose = true);
	
	

	//Prints information about the number and sizes of singularities and fenced regions to std::cout.
	void PrintSingularityStatistics() const;

	//Performs quantitative evaluation on the current parametrization.
	//logAreaRatios - Output variable. Statistics about logarithmic area ratios over all quads.
	//mips - Output variable. Statistics about MIPS energy over all triangles.
	//printStatistics - set to true in order to immediately print the gathered statistics.
	void EvaluateParametrization(Statistics& logAreaRatios, Statistics& mips, bool printStatistics);

	//TODO: fix
	//Performs remeshing of the current mesh based on the current parametrization (must be unpacked).
	//Saves the result in "remeshed.obj".
	void Remesh();

	//TODO: fix
	//Renders surface information represented in the mesh colors framework to a texture laid out
	//with the current packed parametrization.
	void RenderMeshColorsToTexture(const std::string& meshColorsFile);

	//Performs tangential smoothing of the current mesh geometry using OpenMesh.
	void TangentialSmooth();	

	//Saves the current fenced regions to a file.
	void SaveFencedRegions() const;		
	
	//Establishes colors for all texture patches and tries to ensure that neighboring patches
	//have different colors.
	void FindPatchColors();

	//Returns the object that stores texture coordinates for a given patch. Creates the storage if it does not exist.
	TextureCoordinatesStorage& AccessTexCoords(size_t patchIdx);

	//Returns the object that stores texture coordinates for a given patch. Throws an exception if it does not exist.
	const TextureCoordinatesStorage& AccessTexCoords(size_t patchIdx) const;

	//Generates texture coordinates storage objects for all currently existing patches.
	void GenerateTexCoordAccess();


	// --------  Constant accessors to data  --------

	//Returns the current motorcycle graph. Only available after calling CalculateMotorcycleGraph()
	const MotorcycleGraph* Motorcycles() const { return motorcycleGraph.get(); }

	//Returns the indices of all broken halfarc constraints (i.e. visible seams). Only available after
	//calling CalculateParametrization().
	const std::vector<size_t>& BrokenHalfarcs() const { return brokenHalfarcs; }

	//Returns the list of faces of the current mesh.
	const FaceList& Faces() const { return F; }

	//Returns the list of vertices of the current mesh.
	const Matrix3Xf& Vertices() const { return V; }

	//Returns the halfedge data structure of the current mesh.
	const HEMesh& Mesh() const { return mesh; }

	//Returns the bounding box of the current mesh.
	const nse::math::BoundingBox<float, 3>& MeshBoundingBox() const { return meshBoundingBox; }

	//Returns the current set of fenced regions. Only available after calling ClassifyAllSingularities().
	const PatchSet& FencedRegions() const { return fencedRegions; }

	//Returns if a given edge is an original edge (i.e. not created through subdivision).
	bool IsEdgeOriginal(HEMesh::EdgeHandle e) const { return mesh.property(isOriginalEdgeProp, e); }

	//Returns the current set of singularities.
	const std::vector<SingularityInfo>& Singularities() const { return singularities; }

	//Returns the current set of meta singularities (i.e. spawning points of motorcycles).
	const std::vector<MetaSingularity>& MetaSingularities() const { return metaSingularities; }

	//Returns the average edge length of the current mesh
	float AverageEdgeLength() const { return averageEdgeLength; }

	//Returns a structure that maps irregular vertices to their singularity information.
	const ManifoldnessAwareVertexMap<size_t>& VertexToSingularity() const { return vertexToSingularity; }

	//Returns a structure that maps vertices to their meta singularity information (if they are meta singularities)
	const ManifoldnessAwareVertexMap<size_t>& VertexToMetaSingularity() const { return vertexToMetaSingularity; }

	//Returns the color for a given patch. Only available after calling FindPatchColors().
	Eigen::Vector4f GetPatchColor(size_t patchIndex) const;
private:	

	//Represents an edge in the fenced region graph used for merging.
	struct GraphEdge
	{
		//Index of the source GraphNode of this edge.
		size_t from;
		
		//Index of the target GraphNode of this edge.
		size_t to;

		//Index of the last BFS node originating at the from GraphNode
		size_t fromBFSNode;
		
		//Index of the last BFS node originating at the to GraphNode
		size_t toBFSNode;

		//Length of the edge in number of faces
		size_t distance;

		//Temporary variable that determines if a merge operation uses this edge.
		bool enabled;
	};

	//Represents a node in the fenced region graph used for merging.
	struct GraphNode
	{
		//The index of the corresponding fenced region.
		size_t patchIdx;

		//List of all incident edges.
		//key = target patch id, value = idx of LoopGraphEdge
		std::map<size_t, size_t> outgoingEdges; 

		//Temporary variables that determine if the node has been processed.
		bool handled, handledForEvaluatingMergeQuality;

		GraphNode(size_t patchIdx)
			: patchIdx(patchIdx), handled(false)
		{ }
	};	

	//Mesh edge property that determines if an edge is original (i.e. not introduced through subdivision)
	OpenMesh::EPropHandleT<bool> isOriginalEdgeProp;

	//If Catmull Clark subdivision is performed after loading the mesh, represents information about
	//how the subdivision has been performed in order to undo it at the end.
	struct SubdivisionInfo
	{
		//The face that was generated after subdivision
		HEMesh::FaceHandle subdividedFace;

		//The corresponding corner of the original mesh that produced this face. This vertex is
		//also present in the subdivided mesh.
		HEMesh::VertexHandle originalCornerVertex;

		SubdivisionInfo() { }
		SubdivisionInfo(HEMesh::FaceHandle subdividedFace, HEMesh::VertexHandle originalCornerVertex)
			: subdividedFace(subdividedFace), originalCornerVertex(originalCornerVertex)
		{ }
	};
	
	//Loads a mesh into the structure, generating all auxiliary data structures
	void LoadPolygons(const Matrix3Xf& V, const FaceList& F, const std::vector<std::vector<unsigned int>>& subdivisionInfo);

	//For each fenced region, gathers the edges over which motorcycles can enter it and stores them in 
	//fencedRegions.enteringEdgeToPatch.
	void RecordPatchEnteringEdges();

	//Generates a new empty motorcycle graph and subscribes to the necessary events.
	void GenerateNewMotorcycleGraph();

	//Evaluates the cost for merging a subgraph of the fenced region graph.
	//Returns if merging is possible (i.e. if the union becomes regular).
	//nodeSubset - Indices into nodes that represent the nodes in the subgraph.
	//edges - edges of the fenced region graph
	//singularities - Output variable. Number of meta singularities resulting from this merge.
	//totalFacesInPatches - Output variable. Total number of faces covered by fenced region in the subgraph of this merge.
	bool EvaluateMergeQuality(const std::set<size_t>& nodeSubset, std::vector<GraphNode>& nodes, const std::vector<GraphEdge>& edges, size_t& singularities, size_t& totalFacesInPatches);

	//Maximum size of fenced region in number of faces.
	const int maxPatchSize = 100;// 400;

	//Determines if the texCoords storage objects need to be re-allocated.
	bool texCoordsDirty = false;

	//Objects for storing texture coordinates for every texture patch.
	std::vector<TextureCoordinatesStorage> texCoords;

	FaceList F;
	Matrix3Xf V;

	HEMesh mesh;
	nse::math::BoundingBox<float, 3> meshBoundingBox;

	//For each original face, a list of the subdivided faces.
	std::vector<std::vector<SubdivisionInfo>> subdivisionInfo;

	std::vector<SingularityInfo> singularities;
	ManifoldnessAwareVertexMap<size_t> vertexToSingularity;

	std::vector<MetaSingularity> metaSingularities;
	ManifoldnessAwareVertexMap<size_t> vertexToMetaSingularity;
	bool metaSingularitiesDirty = true;

	PatchSet fencedRegions;

	int textureWidth = 0, textureHeight = 0;

	std::unique_ptr<::MotorcycleGraph> motorcycleGraph;

	//Holds pairs (i -> (j, k)) that represent that fenced region i has been merged from
	//fenced region j and k.
	std::map<size_t, std::array<size_t, 2>> regionWasMergedFrom;

	float averageEdgeLength;
	std::vector<size_t> brokenHalfarcs;


	const LengthMeasure lengthMeasure = Geometrical;

	std::vector<int> patchColors;

	//Union find to extract the connected components of the mesh
	nse::util::UnionFind connectedComponentsUF;

	//vertex ids of the representatives	for every connected component
	std::vector<unsigned int> connectedComponents; 
};