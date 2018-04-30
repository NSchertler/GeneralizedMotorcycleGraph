#pragma once

#include "common.h"
#include "TexturePatch.h"
#include "MotorcycleGraph.h"

#include <nsessentials/math/LeastSquaresSystem.h>

#include <Eigen/Core>
#include <cmath>

#include <vector>
#include <set>

//We virtually triangulate quads as follows:
//          h2
//     v3  <--  v2
//  h3  |    /   | h1
//     v0  -->  v1
//          h0
//The halfedge h0 (v0->v1) is the halfedge referenced by the face.	

//Establishes unknown indices for every entry in texCoords. If separateSystemsForXY is set to true, the algorithm
//uses separate indices for the x-components and y-components. If an entry in texCoords is nan, it is considered unknown
//and an new contiguous index (starting at 0) is assigned. Otherwise, the index is set to -1.
extern Eigen::Vector2i SetupUnknownIndices(const TextureCoordinatesStorage& texCoords, std::vector<Eigen::Vector2i>& outIndices, bool separateSystemsForXY);

//Adds the cotangent matrix for the given set of triangles to the provided systems using the fixed texture coordinates from texCoords or the unknown indices
//from vertexIdxToUnkownIdx.
//system - pointers to the least squares systems for the x/y components. May point to the same system.
//outConstantEnergyTerm - Output variable. The constant energy term of the resulting energy that is not stored in the system.
//factor - weight for the energy terms
extern void AddCotanMatrix(const std::vector<std::array<HEMesh::VertexHandle, 3>>& triangles, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor = 1.0f);

//Adds the cotangent matrix for the given set of faces to the provided systems using the fixed texture coordinates from texCoords or the unknown indices
//from vertexIdxToUnkownIdx.
//system - pointers to the least squares systems for the x/y components. May point to the same system.
//outConstantEnergyTerm - Output variable. The constant energy term of the resulting energy that is not stored in the system.
//factor - weight for the energy terms
extern void AddCotanMatrix(const std::set<HEMesh::FaceHandle>& faces, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor = 1.0f);

//Adds the area matrix for the given set of boundary edges to the provided systems using the fixed texture coordinates from texCoords or the unknown indices
//from vertexIdxToUnkownIdx.
//system - pointers to the least squares systems for the x/y components. May point to the same system.
//outConstantEnergyTerm - Output variable. The constant energy term of the resulting energy that is not stored in the system.
//factor - weight for the energy terms
extern void AddAreaMatrix(std::set<HEMesh::HalfedgeHandle> boundary, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor = 1.0f);


//Function pointers that set part of the x/y coordinates to an axis aligned rectangle for each of the
//four rectangle sides. E.g. setupBoundaryConstraint[0] sets sets the bottom side, i.e. y=0.
static void(*setupBoundaryConstraint[4])(Eigen::Vector2f& uv, float sizeU, float sizeV) =
{
	[](Eigen::Vector2f& uv, float sizeU, float sizeV) { uv.y() = 0; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV) { uv.x() = sizeU; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV) { uv.y() = sizeV; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV) { uv.x() = 0; },
};

//Function pointers that set part of the x/y coordinates to a position on a side of an axis aligned rectangle.
//E.g. setupLineConstraint[0] sets sets the bottom side, i.e. x=position.
static void(*setupLineConstraint[4])(Eigen::Vector2f& uv, float sizeU, float sizeV, float position) =
{
	[](Eigen::Vector2f& uv, float sizeU, float sizeV, float position) { uv.x() = position; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV, float position) { uv.y() = position; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV, float position) { uv.x() = sizeU - position; },
	[](Eigen::Vector2f& uv, float sizeU, float sizeV, float position) { uv.y() = sizeV - position; },
};

//Sets up boundary constraints for the given patch.
//Use[Source/Interior/Target]LineConstraintForArcCallback: bool(int side, size_t iArc)
//texCoords - the texture coordinates for the given patch. This function will set the according coordinates of vertices on the boundary.
//patchSize - the parametric patch size
//parametricHalfarcLengths - the calculated list of parametric halfarc lengths
//mesh - the underlying mesh
//graph - the underlying motorcycle graph
//useSourceLineConstraintForArc - callback to determine if a line constraint should be used for the first vertex of a given arc.
//useTargetLineConstraintForArc - callback to determine if a line constraint should be used for the last vertex of a given arc.
//useInteriorLineConstraintForArc - callback to determine if line constraints should be used for a given arc, i.e. if the texels should be distributed uniformly along that arc.
template <typename UseSourceLineConstraintForArcCallback, typename UseTargetLineConstraintForArcCallback, typename UseInteriorLineConstraintForArcCallback>
void SetupBoundaryConstraints(TextureCoordinatesStorage& texCoords, float patchSize[2], const std::vector<float>& parametricHalfarcLengths, const HEMesh& mesh, const MotorcycleGraph& graph, UseSourceLineConstraintForArcCallback&& useSourceLineConstraintForArc, UseTargetLineConstraintForArcCallback&& useTargetLineConstraintForArc, UseInteriorLineConstraintForArcCallback&& useInteriorLineConstraintForArc)
{
	auto& patch = texCoords.Patch();
	//Iterate the four sides of the rectangular patch
	for (int side = 0; side < 4; ++side)
	{
		//sum of arc length along this patch side
		float summedArcLengths = 0;

		//Iterate the arcs that make up that patch side
		for (int iArc = 0; iArc < patch.PatchSides()[side].size(); ++iArc)
		{
			auto arcIdx = patch.PatchSides()[side][iArc];
			auto& arc = graph.Halfarcs()[arcIdx];

			//Set up the source line constraint if applicable
			if (std::forward<UseSourceLineConstraintForArcCallback>(useSourceLineConstraintForArc)(side, iArc))
			{
				auto h = graph.MotorcycleHalfedge(arc.segments.front().location);
				if (!mesh.is_boundary(mesh.opposite_halfedge_handle(h)))
					setupLineConstraint[side](texCoords.TexCoordAtFromVertex(h, mesh), patchSize[0], patchSize[1], summedArcLengths);
			}

			auto arcLength = arc.Length();
			bool firstSegment = true;
			//Iterate the vertices on the arc
			for (auto segment : arc)
			{
				auto h = graph.MotorcycleHalfedge(segment);
				if (mesh.is_boundary(mesh.opposite_halfedge_handle(h)))
					continue;

				//Set up the boundary constraint for the first vertex on the patch side
				if (iArc == 0 && firstSegment)
				{
					setupBoundaryConstraint[side](texCoords.TexCoordAtFromVertex(h, mesh), patchSize[0], patchSize[1]);
					firstSegment = false;
				}

				//Set up the boundary constraint for all other vertices on the patch side
				auto& toUV = texCoords.TexCoordAtToVertex(h, mesh);
				setupBoundaryConstraint[side](toUV, patchSize[0], patchSize[1]);

				summedArcLengths += parametricHalfarcLengths[arcIdx] / arcLength;

				//Set up the interior line constraints if applicable
				if (std::forward<UseInteriorLineConstraintForArcCallback>(useInteriorLineConstraintForArc)(side, iArc))
					setupLineConstraint[side](toUV, patchSize[0], patchSize[1], summedArcLengths);
			}

			//Set up the target line constraint if applicable
			if (std::forward<UseTargetLineConstraintForArcCallback>(useTargetLineConstraintForArc)(side, iArc))
			{
				auto h = graph.MotorcycleHalfedge(arc.LastPathSegment());
				if (!mesh.is_boundary(mesh.opposite_halfedge_handle(h)))
					setupLineConstraint[side](texCoords.TexCoordAtToVertex(h, mesh), patchSize[0], patchSize[1], summedArcLengths);
			}
		}
	}
}

//Calculates LSCM parametrization for the given patch.
extern float CalculateLSCM(TextureCoordinatesStorage& texCoords, const float patchSize[2], const HEMesh& mesh, const MotorcycleGraph& graph);

//Calculates scaffold map parametrization for the given patch, initialized with LSCM.
extern void CalculateScaffoldMap(TextureCoordinatesStorage& texCoords, const float patchSize[2], const HEMesh& mesh, const MotorcycleGraph& graph, float factor);