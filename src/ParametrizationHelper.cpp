#include "ParametrizationHelper.h"
#include <nsessentials/util/IndentationLog.h>

#include <igl/scaf.h>

Eigen::Matrix<int, 2, 1> SetupUnknownIndices(const TextureCoordinatesStorage& texCoords, std::vector<Eigen::Matrix<int, 2, 1>>& outIndices, bool separateSystemsForXY)
{
	Eigen::Matrix<int, 2, 1> _nextUnknown;
	_nextUnknown.setZero();
	int* nextUnknown[2] = { &_nextUnknown.x(), (separateSystemsForXY ? &_nextUnknown.y() : &_nextUnknown.x()) };
	for (int v = 0; v < texCoords.TexCoordCount(); ++v)
	{
		for (int i = 0; i < 2; ++i)
		{
			if (std::isnan(texCoords.TexCoord(v)[i]))
			{
				//Entry is unknown
				outIndices[v][i] = *nextUnknown[i];

				//Generate a new index
				++(*nextUnknown[i]);
			}
			else
			{
				//Entry is known
				outIndices[v][i] = -1;
			}
		}
	}
	return _nextUnknown;
}

void AddCoefficientOrFixed(int row, int column, float coefficient, float fixedValueOfCol, float fixedValueOfRow, nse::math::LeastSquaresSystem<1>& system, float& outConstantEnergyTerm)
{
	if(row != -1 && column != -1)
		system.lhs.coeffRef(row, column) += coefficient;
	else if (row == -1 && column == -1)
		outConstantEnergyTerm += coefficient * fixedValueOfRow * fixedValueOfCol;
	else if (column == -1)
		system.rhs.row(row)(0) -= coefficient * fixedValueOfCol;		
}

//Calculates the cotangent in the triangle v0 -- d0 --> v1 -- d1 --> v2 at v1.
float CalcCotangent(const OpenMesh::Vec3f& d0, const OpenMesh::Vec3f& d1)
{
	auto d2 = -(d0 + d1);
	float triangleArea = 0.5f * OpenMesh::cross(d0, -d1).norm();
	float cotangent = (d0.sqrnorm() + d1.sqrnorm() - d2.sqrnorm()) / (4.0f * triangleArea);

	return cotangent;
}

//calculates the cotangent of the angle opposite to the halfedge in the same triangle
float CalcHalfedgeCotangent(HEMesh::HalfedgeHandle h, const HEMesh& mesh)
{
	auto f = mesh.face_handle(h);
	if (!f.is_valid())
		return 0.0f;
	int halfedgeNumberInQuad = 0;
	auto testH = mesh.halfedge_handle(f);
	while (testH != h)
	{
		++halfedgeNumberInQuad;
		testH = mesh.next_halfedge_handle(testH);
	}
	OpenMesh::Vec3f d0, d1;
	switch (halfedgeNumberInQuad)
	{
	case 0:
	case 2:
	{
		auto next = mesh.next_halfedge_handle(h);
		d0 = mesh.calc_edge_vector(next);
		d1 = mesh.point(mesh.from_vertex_handle(h)) - mesh.point(mesh.to_vertex_handle(next));
		break;
	}
	case 1:
	case 3:
	{
		auto prev = mesh.prev_halfedge_handle(h);
		d0 = mesh.point(mesh.from_vertex_handle(prev)) - mesh.point(mesh.to_vertex_handle(h));
		d1 = mesh.calc_edge_vector(prev);
		break;
	}
	}
	return CalcCotangent(d0, d1);
}

void AddCotanWeight(size_t sourceIdx, size_t targetIdx, float weight, nse::math::LeastSquaresSystem<1>* system[2], 
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm)
{
	auto& sourceUnknownIdx = vertexIdxToUnknownIdx[sourceIdx];
	auto& targetUnknownIdx = vertexIdxToUnknownIdx[targetIdx];

	for (int i = 0; i < 2; ++i)
	{
		AddCoefficientOrFixed(sourceUnknownIdx[i], sourceUnknownIdx[i], weight, texCoords.TexCoord(sourceIdx)[i], texCoords.TexCoord(sourceIdx)[i], *system[i], outConstantEnergyTerm);
		AddCoefficientOrFixed(targetUnknownIdx[i], targetUnknownIdx[i], weight, texCoords.TexCoord(targetIdx)[i], texCoords.TexCoord(targetIdx)[i], *system[i], outConstantEnergyTerm);
		AddCoefficientOrFixed(sourceUnknownIdx[i], targetUnknownIdx[i], -weight, texCoords.TexCoord(targetIdx)[i], texCoords.TexCoord(sourceIdx)[i], *system[i], outConstantEnergyTerm);
		AddCoefficientOrFixed(targetUnknownIdx[i], sourceUnknownIdx[i], -weight, texCoords.TexCoord(sourceIdx)[i], texCoords.TexCoord(targetIdx)[i], *system[i], outConstantEnergyTerm);
	}
}

void AddCotanMatrix(const std::vector<std::array<HEMesh::VertexHandle, 3>>& triangles, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor)
{
	for (auto& tri : triangles)
	{
		for (int i = 0; i < 3; ++i)
		{
			int prev = (i + 2) % 3;
			int next = (i + 1) % 3;

			auto d0 = mesh.point(tri[i]) - mesh.point(tri[prev]);
			auto d1 = mesh.point(tri[next]) - mesh.point(tri[i]);

			float weight = CalcCotangent(d0, d1);
			AddCotanWeight(texCoords.Patch().IdOfInnerVertex(tri[next], mesh), texCoords.Patch().IdOfInnerVertex(tri[prev], mesh), factor * weight, system, texCoords, vertexIdxToUnknownIdx, outConstantEnergyTerm);
		}
	}
}

void AddCotanMatrix(const std::set<HEMesh::FaceHandle>& faces, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor)
{
	auto& patch = texCoords.Patch();

	for (auto f : faces)
	{
		for (auto h : mesh.fh_range(f))
		{
			auto weight = CalcHalfedgeCotangent(h, mesh);

			AddCotanWeight(patch.IdOfFromVertex(h, mesh), patch.IdOfToVertex(h, mesh), factor * weight, system, texCoords, vertexIdxToUnknownIdx, outConstantEnergyTerm);
		}

		//Add the quad diagonals
		auto h = mesh.halfedge_handle(f);
		auto target = patch.IdOfFromVertex(h, mesh);

		auto d0 = mesh.calc_edge_vector(h);
		h = mesh.next_halfedge_handle(h);
		auto source = patch.IdOfToVertex(h, mesh);

		auto d1 = mesh.calc_edge_vector(h);
		float weight = CalcCotangent(d0, d1);

		AddCotanWeight(source, target, factor * weight, system, texCoords, vertexIdxToUnknownIdx, outConstantEnergyTerm);


		h = mesh.next_halfedge_handle(h);
		target = patch.IdOfFromVertex(h, mesh);

		d0 = mesh.calc_edge_vector(h);
		h = mesh.next_halfedge_handle(h);
		source = patch.IdOfToVertex(h, mesh);

		d1 = mesh.calc_edge_vector(h);
		weight = CalcCotangent(d0, d1);

		AddCotanWeight(source, target, factor * weight, system, texCoords, vertexIdxToUnknownIdx, outConstantEnergyTerm);
	}
}

void AddAreaTermForEdge(HEMesh::HalfedgeHandle h, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh, 
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor)
{
	size_t vertexIndex[2] = { texCoords.Patch().IdOfFromVertex(h, mesh), texCoords.Patch().IdOfToVertex(h, mesh) };
	const Eigen::Vector2i* unknownIndex[2] = { &vertexIdxToUnknownIdx[vertexIndex[0]], &vertexIdxToUnknownIdx[vertexIndex[1]] };

	struct Entry
	{
		int firstVertex, firstCoordinate;
		int secondVertex, secondCoordinate;
		float coefficient;
	};

	const Entry entries[4] =
	{
		{ 0, 1, 1, 0, -0.5 },
		{ 1, 0, 0, 1, -0.5 },
		{ 0, 0, 1, 1, 0.5 },
		{ 1, 1, 0, 0, 0.5 },
	};

	for (auto& entry : entries)
	{
		AddCoefficientOrFixed((*unknownIndex[entry.firstVertex])[entry.firstCoordinate], (*unknownIndex[entry.secondVertex])[entry.secondCoordinate],
			factor * entry.coefficient, texCoords.TexCoord(vertexIndex[entry.secondVertex])[entry.secondCoordinate],
			texCoords.TexCoord(vertexIndex[entry.firstVertex])[entry.firstCoordinate], *system[entry.firstCoordinate], outConstantEnergyTerm);
	}
}

void AddAreaMatrix(std::set<HEMesh::HalfedgeHandle> boundary, nse::math::LeastSquaresSystem<1>* system[2], const HEMesh& mesh,
	const TextureCoordinatesStorage& texCoords, const std::vector<Eigen::Vector2i>& vertexIdxToUnknownIdx, float& outConstantEnergyTerm, float factor)
{
	for (auto e : boundary)
		AddAreaTermForEdge(e, system, mesh, texCoords, vertexIdxToUnknownIdx, outConstantEnergyTerm, factor);
}

float CalculateLSCM(TextureCoordinatesStorage& texCoords, const float patchSize[2], const HEMesh& mesh, const MotorcycleGraph& graph)
{	
	std::vector<Eigen::Matrix<int, 2, 1>> unknownIndexPerVertex(texCoords.TexCoordCount());

	std::set<HEMesh::HalfedgeHandle> patchBoundary;
	for (auto h : texCoords.Patch().OpenBoundaryEdges())
		patchBoundary.insert(h);

	for (int side = 0; side < 4; ++side)
	{
		for (int iArc = 0; iArc < texCoords.Patch().PatchSides()[side].size(); ++iArc)
		{
			auto arcIdx = texCoords.Patch().PatchSides()[side][iArc];
			auto& arc = graph.Halfarcs()[arcIdx];

			for (auto segment : arc)
			{
				auto h = graph.MotorcycleHalfedge(segment);
				patchBoundary.insert(h);
			}
		}
	}	

	//assign indices to all unknown variables
	auto numberOfUnknowns = SetupUnknownIndices(texCoords, unknownIndexPerVertex, false).x();
	if (numberOfUnknowns == 0)
		return 0;

	/*for (int i = 0; i < 2; ++i)
	{
	std::cout << "Starting dimension " << i << " at unknown index " << nextUnknown << std::endl;
	for (auto& v : patch.verticesToIndex)
	{
	if (std::isnan(patch.textureCoordinates[v.second][i]))
	{
	unknownIndexPerVertex[v.second][i] = nextUnknown;
	++nextUnknown;
	}
	else
	unknownIndexPerVertex[v.second][i] = -1;
	}
	}
	std::cout << nextUnknown << " unknowns, " << patch.verticesToIndex.size() * 2 - nextUnknown << " fixed variables" << std::endl;*/

	/*int s = 0;
	for (auto& side : patch.PatchSides())
	{
	std::cout << "Side " << (s++) << std::endl;
	for (auto iArc : side)
	{
	auto& arc = graph.Halfarcs()[iArc];
	std::cout << "Arc " << iArc << std::endl;
	std::cout << "Motorcycle " << arc.start.motorcycle << ", segment " << arc.start.pathSegment << ", forward = " << arc.start.goForward << std::endl;
	if (arc.length == 0)
	continue;
	auto v = optData.mesh->from_vertex_handle(graph.MotorcycleHalfedge(arc.start));
	auto index = patch.verticesToIndex[v];
	std::cout << patch.textureCoordinates[index].transpose() << " unknown index=" << unknownIndexPerVertex[index].transpose() << std::endl;
	for (int p = 0; p < arc.length; ++p)
	{
	v = optData.mesh->to_vertex_handle(graph.MotorcycleHalfedge(arc.GetSegment(p)));
	index = patch.verticesToIndex[v];
	std::cout << patch.textureCoordinates[index].transpose() << " unknown index=" << unknownIndexPerVertex[index].transpose() << std::endl;
	}
	}
	}*/

	nse::math::LeastSquaresSystem<1> system(numberOfUnknowns);
	nse::math::LeastSquaresSystem<1>* systemArray[2] = { &system, &system };
	system.lhs.reserve(Eigen::VectorXi::Constant(numberOfUnknowns, 9));
	float constantEnergyTerm = 0;

	//add the Laplacian to the system
	AddCotanMatrix(texCoords.Patch().Faces(), systemArray, mesh, texCoords, unknownIndexPerVertex, constantEnergyTerm);
	AddCotanMatrix(texCoords.Patch().FilledHoles(), systemArray, mesh, texCoords, unknownIndexPerVertex, constantEnergyTerm);

	//Add the area term to the system
	AddAreaMatrix(patchBoundary, systemArray, mesh, texCoords, unknownIndexPerVertex, constantEnergyTerm, -2.0f);

	Eigen::SimplicialLDLT<nse::math::LeastSquaresSystem<1>::MatrixType> solver;
	solver.compute(system.lhs);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "Matrix could not be decomposed: " << solver.info() << std::endl;
		return -1;
	}
	Eigen::Matrix<float, -1, 1> solution = solver.solve(system.rhs);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "System could not be solved: " << solver.info() << std::endl;
		return -1;
	}

	auto error = ((solution.transpose() * system.lhs * solution)(0) - (2 * solution.transpose() * system.rhs)(0) + constantEnergyTerm) / texCoords.TexCoordCount();
	//std::cout << "error for patch with vertex " 
	//	<< (mesh.from_vertex_handle(graph.MotorcycleHalfedge(*graph.Halfarcs()[patch.PatchSides().front().front()].begin()))) << " and " 
	//	<< (mesh.to_vertex_handle(graph.MotorcycleHalfedge(*graph.Halfarcs()[patch.PatchSides().front().front()].begin())))
	//	<< ": " << error << std::endl;

	for (size_t v = 0; v < texCoords.TexCoordCount(); ++v)
	{
		auto& unknownIndex = unknownIndexPerVertex[v];
		for (int i = 0; i < 2; ++i)
		{
			if (unknownIndex[i] != -1)
				texCoords.TexCoord(v)[i] = solution.row(unknownIndex[i])(0);
		}
	}

	return error;

	/*int s = 0;
	for (auto& side : patch.PatchSides())
	{
	std::cout << "Side " << (s++) << std::endl;
	for (auto iArc : side)
	{
	auto& arc = graph.Halfarcs()[iArc];
	std::cout << "Arc " << iArc << std::endl;
	std::cout << "Motorcycle " << arc.start.motorcycle << ", segment " << arc.start.pathSegment << ", forward = " << arc.start.goForward << std::endl;
	if (arc.length == 0)
	continue;
	auto v = optData.mesh->from_vertex_handle(graph.MotorcycleHalfedge(arc.start));
	std::cout << patch.textureCoordinates[patch.verticesToIndex[v]].transpose() << std::endl;
	for (int p = 0; p < arc.length; ++p)
	{
	v = optData.mesh->to_vertex_handle(graph.MotorcycleHalfedge(arc.GetSegment(p)));
	std::cout << patch.textureCoordinates[patch.verticesToIndex[v]].transpose() << std::endl;
	}
	}
	}*/
}

int InContextOf(int value, int context, int totalSize)
{
	if (value < context)
		return value + totalSize;
	else
		return value;
}

std::pair<int, int> InContextOf(const std::pair<int, int>& value, int context, int totalSize)
{
	auto ret = value;
	if (ret.second != -1 && ret.first < context)
	{
		ret.first += totalSize;
		ret.second += totalSize;
	}
	return ret;
}

#include <igl/writeOBJ.h>

void CalculateScaffoldMap(TextureCoordinatesStorage& texCoords, const float patchSize[2], const HEMesh& mesh, const MotorcycleGraph& graph, float factor)
{
	/*static int patchId = 0;
	std::cout << "Open patch " << (patchId++) << std::endl;
*/
	std::vector<int> fixedX, fixedY;
	for (int i = 0; i < texCoords.TexCoordCount(); ++i)
	{
		auto& tex = texCoords.TexCoord(i);
		if (!std::isnan(tex.x()))
			fixedX.push_back(i);
		if (!std::isnan(tex.y()))
			fixedY.push_back(i);
	}

	Eigen::VectorXi fixedXVec(fixedX.size()), fixedYVec(fixedY.size());
	memcpy(fixedXVec.data(), fixedX.data(), sizeof(int) * fixedX.size());
	memcpy(fixedYVec.data(), fixedY.data(), sizeof(int) * fixedY.size());

	//std::cout.imbue(std::locale());
	//std::cout << "Eigen::VectorXi fixedX(" << fixedX.size() << ");" << std::endl;
	//std::cout << "fixedX << ";
	//for (int i = 0; i < fixedX.size(); ++i)
	//{
	//	if (i != 0)
	//		std::cout << ", ";
	//	std::cout << fixedX[i];
	//}
	//std::cout << ";" << std::endl;
	//std::cout << "Eigen::VectorXi fixedY(" << fixedY.size() << ");" << std::endl;
	//std::cout << "fixedY << ";
	//for (int i = 0; i < fixedY.size(); ++i)
	//{
	//	if (i != 0)
	//		std::cout << ", ";
	//	std::cout << fixedY[i];
	//}
	//std::cout << ";" << std::endl;

	if (fixedY.size() == 0 || fixedX.size() == 0)
	{
		texCoords.ResetTextureCoordinates();
		return; //TODO: fix
	}

	std::map<HEMesh::VertexHandle, HEMesh::HalfedgeHandle> vertexToOutgoingFreeBoundary;
	{
		for (auto h : texCoords.Patch().OpenBoundaryEdges())
		{
			auto v = mesh.from_vertex_handle(h);
			vertexToOutgoingFreeBoundary[v] = h;
		}
	}

	std::vector<HEMesh::HalfedgeHandle> boundaryLoop;
	std::vector<int> cornerPositions(4, -1); //indices on the boundary loop
	std::vector<std::pair<int, int>> sideIntervals(4, std::make_pair(std::numeric_limits<int>::max(), -1)); //positions on the boundary loop

	{
		for (int side = 0; side < 4; ++side)
		{
			for (auto arcIdx : texCoords.Patch().PatchSides()[side])
			{
				auto firstP = *graph.Halfarcs()[arcIdx].begin();
				auto firstH = graph.MotorcycleHalfedge(firstP);
				auto& firstUV = texCoords.TexCoordAtFromVertex(firstH, mesh);
				if (!std::isnan(firstUV.x()) && !std::isnan(firstUV.y()))
					cornerPositions[side] = boundaryLoop.size();
				for (auto p : graph.Halfarcs()[arcIdx])
				{
					auto h = graph.MotorcycleHalfedge(p);
					sideIntervals[side].first = std::min(sideIntervals[side].first, (int)boundaryLoop.size());
					sideIntervals[side].second = std::max(sideIntervals[side].second, (int)boundaryLoop.size() + 1);
					boundaryLoop.push_back(h);
					std::map<HEMesh::VertexHandle, HEMesh::HalfedgeHandle>::iterator it;
					while ((it = vertexToOutgoingFreeBoundary.find(mesh.to_vertex_handle(h))) != vertexToOutgoingFreeBoundary.end())
					{
						h = it->second;
						boundaryLoop.push_back(h);
					}
				}

				auto lastP = graph.Halfarcs()[arcIdx].LastPathSegment();
				auto lastH = graph.MotorcycleHalfedge(lastP);
				auto& lastUV = texCoords.TexCoordAtToVertex(lastH, mesh);
				if (!std::isnan(lastUV.x()) && !std::isnan(lastUV.y()))
					cornerPositions[(side + 1) % 4] = sideIntervals[side].second;
			}
		}
	}

	//std::cout << "Boundary edges: " << boundaryLoop.size() << std::endl;
	//for (int i = 0; i < 4; ++i)
	//{
	//	std::cout << "Corner " << i << ": " << cornerPositions[i] << std::endl;
	//	std::cout << "Side Interval " << i << ": " << "[" << sideIntervals[i].first << ", " << sideIntervals[i].second << "]" << std::endl;
	//}

	{
		//find first specified corner
		int corner = 0;
		while (cornerPositions[corner] == -1)
			++corner;

		int unspecifiedCorners = 0;
		for (int i = 0; i < 4; ++i)
			if (cornerPositions[i] == -1)
				++unspecifiedCorners;

		while (unspecifiedCorners > 0)
		{
			//find next specified corner
			int nextKnownCorner = corner + 1;
			while (cornerPositions[nextKnownCorner % 4] == -1)
				++nextKnownCorner;

			int currentCornerPosition = cornerPositions[corner];
			int nextKnownCornerPosition = InContextOf(cornerPositions[nextKnownCorner % 4], currentCornerPosition, boundaryLoop.size());

			int edgesForDistribution = nextKnownCornerPosition - currentCornerPosition;			

			std::vector<float> segmentsRelativeSizePrefixSum(nextKnownCorner - corner);
			float totalRelativeSize = 0;
			for (int i = 0; i < segmentsRelativeSizePrefixSum.size(); ++i)
			{
				float relativeSize = patchSize[(corner + i) % 2];
				totalRelativeSize += relativeSize;
				segmentsRelativeSizePrefixSum[i] = relativeSize;
				if (i > 0)
					segmentsRelativeSizePrefixSum[i] += segmentsRelativeSizePrefixSum[i - 1];
			}

			int unknownCornerFrom = corner + 1;
			while (unknownCornerFrom < nextKnownCorner)
			{
				int unknownCornerFromPosition = InContextOf(cornerPositions[(unknownCornerFrom - 1 + 4) % 4], currentCornerPosition, boundaryLoop.size());				
				auto firstInterval = InContextOf(sideIntervals[(unknownCornerFrom - 1 + 4) % 4], currentCornerPosition, boundaryLoop.size());				
				
				int lowerPositionBoundIncl = std::max(unknownCornerFromPosition + 1, firstInterval.second);
				int upperBoundCornerIncl = unknownCornerFrom;
				int upperPositionBoundIncl;
				//find the upper bound for the segment
				while (true)
				{
					auto nextInterval = InContextOf(sideIntervals[upperBoundCornerIncl % 4], currentCornerPosition, boundaryLoop.size());
					if (nextInterval.second != -1)
					{
						upperPositionBoundIncl = nextInterval.first;
						break;
					}
					if (upperBoundCornerIncl + 1 == nextKnownCorner)
					{
						upperPositionBoundIncl = InContextOf(cornerPositions[nextKnownCorner % 4] - 1, currentCornerPosition, boundaryLoop.size());
						break;
					}
					++upperBoundCornerIncl;
				}				

				int cornersToDistribute = upperBoundCornerIncl + 1 - unknownCornerFrom;							

				for (int i = 0; i < cornersToDistribute; ++i)
				{
					int currentCorner = unknownCornerFrom + i;
					auto cornerPos = currentCornerPosition + (int)std::round(segmentsRelativeSizePrefixSum[currentCorner - corner - 1] * edgesForDistribution / (double)(totalRelativeSize));
					if (cornerPos < lowerPositionBoundIncl)
						cornerPos = lowerPositionBoundIncl;
					if (cornerPos > upperPositionBoundIncl - (cornersToDistribute - i - 1)) //reserve space for later corners
						cornerPos = upperPositionBoundIncl - (cornersToDistribute - i - 1);
					lowerPositionBoundIncl = cornerPos + 1;
					
					if (currentCorner >= 4)
						cornerPositions[currentCorner - 4] = cornerPos - boundaryLoop.size();
					else
						cornerPositions[currentCorner] = cornerPos;					
					--unspecifiedCorners;
				}
				unknownCornerFrom = upperBoundCornerIncl + 1;
			}
			corner = nextKnownCorner;
		}
	}

	/*std::cout << "After adaptation: " << std::endl;
	for (int i = 0; i < 4; ++i)
	{
		std::cout << "Corner " << i << ": " << cornerPositions[i] << std::endl;
		std::cout << "Side Interval " << i << ": " << "[" << sideIntervals[i].first << ", " << sideIntervals[i].second << "]" << std::endl;
	}*/

	//Fix boundary
	for (int side = 0; side < 4; ++side)
	{
		int cornerFrom = cornerPositions[side];
		if (cornerFrom < 0)
			cornerFrom += boundaryLoop.size();
		int cornerTo = cornerPositions[(side + 1) % 4];
		if (cornerTo < cornerFrom)
			cornerTo += boundaryLoop.size();
		for (int s = cornerFrom; s < cornerTo; ++s)
		{
			auto h = boundaryLoop[s % boundaryLoop.size()];
			if (s == cornerFrom)
			{
				auto& tex = texCoords.TexCoordAtFromVertex(h, mesh);
				//auto idx = patch.IdOfFromVertex(h, mesh);
				//std::cout << "Updating constraint from " << tex.transpose() << " to ";
				setupBoundaryConstraint[side](tex, patchSize[0], patchSize[1]);
				setupLineConstraint[side](tex, patchSize[0], patchSize[1], 0);
				//std::cout << h.idx() << "  " << idx << "   " << tex.transpose() << std::endl;
			}
			{
				auto& tex = texCoords.TexCoordAtToVertex(h, mesh);
				//auto idx = patch.IdOfToVertex(h, mesh);
				//std::cout << "Updating constraint from " << tex.transpose() << " to ";
				setupBoundaryConstraint[side](tex, patchSize[0], patchSize[1]);
				setupLineConstraint[side](tex, patchSize[0], patchSize[1], patchSize[side % 2] * (s + 1 - cornerFrom) / (float)(cornerTo - cornerFrom));
				//std::cout << h.idx() << "  " << idx << "   " << tex.transpose() << std::endl;
			}
		}
	}

	CalculateLSCM(texCoords, patchSize, mesh, graph);

	Eigen::MatrixXd V(texCoords.TexCoordCount(), 3);
	Eigen::MatrixXd UV_init(texCoords.TexCoordCount(), 2);
	Eigen::MatrixXi F(2 * texCoords.Patch().Faces().size(), 3);
	for (int i = 0; i < texCoords.TexCoordCount(); ++i)
	{
		V.row(i) = factor * ToEigenVector(mesh.point(texCoords.Patch().Vertex(i))).cast<double>();
		UV_init.row(i) = texCoords.TexCoord(i).cast<double>();
	}
	int i = 0;
	for (auto f : texCoords.Patch().Faces())
	{
		auto h = mesh.halfedge_handle(f);
		OpenMesh::HalfedgeHandle hs[4];
		for (int j = 0; j < 4; ++j)
		{
			hs[j] = h;
			h = mesh.next_halfedge_handle(h);
		}
		F.row(2 * i + 0) << texCoords.Patch().IdOfFromVertex(hs[0], mesh), texCoords.Patch().IdOfFromVertex(hs[1], mesh), texCoords.Patch().IdOfFromVertex(hs[2], mesh);
		F.row(2 * i + 1) << texCoords.Patch().IdOfFromVertex(hs[0], mesh), texCoords.Patch().IdOfFromVertex(hs[2], mesh), texCoords.Patch().IdOfFromVertex(hs[3], mesh);
		++i;
	}

	//Eigen::MatrixXd CN;
	//Eigen::MatrixXi FTC, FN;
	//igl::readOBJ("D:\\Users\\nicos\\CMakeBuilds\\ad48334b-103a-c139-a430-96b64dbecdcb\\build\\x64-Release\\test_cup.obj", V, UV_init, CN, F, FTC, FN);

	//Eigen::MatrixXd CN;
	//Eigen::MatrixXi FN;
	//igl::writeOBJ("patch.obj", V, F, CN, FN, UV_init, F);

	igl::SCAFData d;
	d.scaf_energy = d.slim_energy = igl::SLIMData::LOG_ARAP;
	igl::scaf_precompute(V, F, UV_init, d, 0);
	igl::scaf_solve(d, 50, fixedXVec, fixedYVec);
	
	for (int i = 0; i < texCoords.TexCoordCount(); ++i)
		texCoords.TexCoord(i) = d.w_uv.row(i).cast<float>();
}