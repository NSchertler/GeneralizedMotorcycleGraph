#include "common.h"
#include <nsessentials/util/TimedBlock.h>

#include <map>
#include <set>
#include <Eigen/Geometry>
#include <iostream>

int FaceValence(OpenMesh::VertexHandle v, const HEMesh& mesh)
{
	int valence = 0;
	for (auto f : mesh.vf_range(v))
		++valence;
	return valence;
}

bool IsToVertexSingularity(HEMesh::HalfedgeHandle h, const HEMesh& mesh)
{
	auto v = mesh.to_vertex_handle(h);

	if (mesh.is_manifold(v))
		return IsSingularityManifoldVertex(v, mesh);
	else
	{
		CirculateBackwardUntil<true>(h, mesh, [](HEMesh::HalfedgeHandle) { return false; });
		return IsSingularityNonManifoldVertex(h, mesh);
	}
}

bool IsSingularityManifoldVertex(OpenMesh::VertexHandle v, const HEMesh& mesh, int& valenceDefect)
{
	int faceValence = FaceValence(v, mesh);
	if (mesh.is_boundary(v))
	{
		valenceDefect = faceValence - 2;
		//all incident faces must be quads
		//for(auto f : mesh.vf_range(v))
		//	if(mesh.valence(f) != 4)
		//		return true;
		if (valenceDefect != 0)
			return true;
	}
	else
	{		
		valenceDefect = faceValence - 4;
		if (faceValence != 4)
			return true;
	}
	
	//check geometry
	//const float angleVarianceThreshold = 0.2f;
	//float n = 0;
	//float mean = 0;
	//float summedVariance = 0;
	//for (auto h : mesh.vih_range(v))
	//{
	//	if (mesh.is_boundary(h))
	//		continue;
	//	auto angle = mesh.calc_sector_angle(h);
	//	++n;
	//	auto newMean = mean + 1.0f / n * (angle - mean);
	//	summedVariance += (angle - mean) * (angle - newMean);
	//	mean = newMean;
	//}
	//if (summedVariance / n > angleVarianceThreshold)
	//	return true;

	return false;
}

bool IsSingularityManifoldVertex(OpenMesh::VertexHandle v, const HEMesh& mesh)
{
	int defect;
	return IsSingularityManifoldVertex(v, mesh, defect);
}

bool IsSingularityNonManifoldVertex(HEMesh::HalfedgeHandle contextOutgoingBoundary, const HEMesh& mesh, int& valenceDefect)
{
	auto h = mesh.opposite_halfedge_handle(contextOutgoingBoundary);
	int faces = 0;
	while (!mesh.is_boundary(h))
	{
		++faces;
		h = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(h));
	}
	valenceDefect = faces - 2;
	return valenceDefect != 0;
}

bool IsSingularityNonManifoldVertex(HEMesh::HalfedgeHandle contextOutgoingBoundary, const HEMesh& mesh)
{
	int defect;
	return IsSingularityNonManifoldVertex(contextOutgoingBoundary, mesh, defect);
}

size_t GetEdgeId(unsigned int v1, unsigned int v2, std::map<std::pair<unsigned int, unsigned int>, size_t>& vertexPairToEdge, size_t& nextEdgeId)
{
	auto lo = (v1 < v2 ? v1 : v2);
	auto hi = (v1 < v2 ? v2 : v1);
	auto inserted = vertexPairToEdge.insert(std::make_pair(std::make_pair(lo, hi), nextEdgeId));
	if (inserted.second)
	{
		++nextEdgeId;
		return nextEdgeId - 1;
	}
	else
		return inserted.first->second;
};

void CatmullClarkSubdivide(FaceList& F, Matrix3Xf& V, std::vector<std::vector<unsigned int>>& subdivisionInfo)
{
	nse::util::TimedBlock b("Catmull-Clark subdivision ..");

	size_t nextEdgeId = 0;
	std::map<std::pair<unsigned int, unsigned int>, size_t> vertexPairToEdge;	

	subdivisionInfo.resize(F.size());

	//establish the edges
	size_t totalFaceValence = 0;
	for (auto& f : F)
	{
		totalFaceValence += f.size();
		for (int i = 0; i < f.size() - 1; ++i)
			GetEdgeId(f[i], f[i + 1], vertexPairToEdge, nextEdgeId);
		GetEdgeId(f.back(), f[0], vertexPairToEdge, nextEdgeId);
	}

	size_t oldVertexCount = V.cols();
	V.conservativeResize(Eigen::NoChange, oldVertexCount + vertexPairToEdge.size() + F.size());
	FaceList newF(totalFaceValence);

	//Calculate vertices for edges
	for (auto& edgeEntry : vertexPairToEdge)
	{
		size_t index = edgeEntry.second;
		unsigned int v1 = edgeEntry.first.first;
		unsigned int v2 = edgeEntry.first.second;
		V.col(oldVertexCount + index) = 0.5f * (V.col(v1) + V.col(v2));
	}

	size_t nextFaceId = 0;
	//Calculate vertices for faces; calculate new faces
	for (int iFace = 0; iFace < F.size(); ++iFace)
	{
		auto& f = F.at(iFace);
		Eigen::Vector3f centroid; centroid.setZero();
		for (int i = 0; i < f.size(); ++i)
		{
			int next = (i + 1) % f.size();
			int prev = (i - 1 + f.size()) % f.size();
			centroid += V.col(f[i]);
			subdivisionInfo[iFace].push_back(nextFaceId);
			newF[nextFaceId].resize(4);
			newF[nextFaceId][0] = GetEdgeId(f[prev], f[i], vertexPairToEdge, nextEdgeId) + oldVertexCount;
			newF[nextFaceId][1] = f[i];
			newF[nextFaceId][2] = GetEdgeId(f[i], f[next], vertexPairToEdge, nextEdgeId) + oldVertexCount;
			newF[nextFaceId][3] = oldVertexCount + vertexPairToEdge.size() + iFace;
			++nextFaceId;
		}
		centroid *= 1.0f / f.size();
		V.col(oldVertexCount + vertexPairToEdge.size() + iFace) = centroid;
	}

	F = std::move(newF);
}

void MergeTriangulatedQuads(FaceList& F, const Matrix3Xf& V, float cosAngleThreshold)
{
	nse::util::TimedBlock b("Merging triangulated quads ..");

	size_t nextEdgeId = 0;
	std::map<std::pair<unsigned int, unsigned int>, size_t> vertexPairToEdge;
	std::map<size_t, std::vector<size_t>> edgeTriangleIncidences;

	//establish edge-face incidences
	for (int iFace = 0; iFace < F.size(); ++iFace)
	{
		auto& f = F.at(iFace);
		for (int i = 0; i < f.size(); ++i)
		{
			auto next = (i + 1) % f.size();
			auto edgeId = GetEdgeId(f[i], f[next], vertexPairToEdge, nextEdgeId);
			if(f.size() == 3)
				edgeTriangleIncidences[edgeId].push_back(iFace);
		}
	}

	for (auto& entry : edgeTriangleIncidences)
	{
		if (entry.second.size() == 2) //TODO: support non-manifold edges
		{
			auto& f1 = F[entry.second[0]];
			auto& f2 = F[entry.second[1]];

			if (f1.size() != 3 || f2.size() != 3)
				continue; //one of the faces has already been merged

			unsigned int f1DiagonalStartVertex, f2DiagonalStartVertex;
			for (int i = 0; i < 3; ++i)
			{
				auto next = (i + 1) % 3;
				if (GetEdgeId(f1[i], f1[next], vertexPairToEdge, nextEdgeId) == entry.first)
					f1DiagonalStartVertex = i;
				if (GetEdgeId(f2[i], f2[next], vertexPairToEdge, nextEdgeId) == entry.first)
					f2DiagonalStartVertex = i;
			}
			unsigned int quadV[] = {
				f1[(f1DiagonalStartVertex + 2) % 3],
				f1[f1DiagonalStartVertex],
				f2[(f2DiagonalStartVertex + 2) % 3],
				f2[f2DiagonalStartVertex] };

			if (f1[f1DiagonalStartVertex] == f2[f2DiagonalStartVertex])
				continue; //The faces are not oriented consistently

			//check the angles
			bool isQuad = true;
			for (int i = 0; i < 4; ++i)
			{
				auto next = (i + 1) % 4;
				auto prev = (i + 3) % 4;
				auto dir1 = (V.col(quadV[prev]) - V.col(quadV[i])).normalized();
				auto dir2 = (V.col(quadV[next]) - V.col(quadV[i])).normalized();
				float dot = std::abs(dir1.dot(dir2));
				if (1 - dot * dot < cosAngleThreshold * cosAngleThreshold)
					isQuad = false;
			}

			if (isQuad)
			{
				f1.resize(4);
				for (int i = 0; i < 4; ++i)
					f1[i] = quadV[i];
				f2.resize(0);
			}
		}
	}

	//delete empty faces
	size_t deleted = 0;
	for (int i = F.size() - 1; i >= 0; --i)
	{
		if (F[i].size() == 0)
		{
			F.erase(F.begin() + i);
			++deleted;
		}
	}	

	std::cout << "Deleted " << deleted << " triangles." << std::endl;
}

void ProjectToTangentSpace(OpenMesh::Vec3f& vec, const OpenMesh::Vec3f& normal)
{
	if (!std::isnan(normal[0]))
		vec -= OpenMesh::dot(vec, normal) * normal;
}

float ContinuationScoreToBoundary(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, const OpenMesh::Vec3f vertexNormal)
{
	auto currentVertex = mesh.to_vertex_handle(edge);

	//find the two boundary edges
	HEMesh::HalfedgeHandle incomingBoundary, outgoingBoundary;
	auto edgeVector = mesh.calc_edge_vector(edge).normalized();

	CirculateForwardUntil<false>(edge, mesh, [&](HEMesh::HalfedgeHandle h)
	{		
		if (mesh.is_boundary(h))
			outgoingBoundary = h;
		if (mesh.is_boundary(mesh.opposite_halfedge_handle(h)))
			incomingBoundary = h;

		return false;
	});

	if (incomingBoundary.is_valid())
	{
		//check continuation to boundary
		auto dIncomingOM = mesh.calc_edge_vector(incomingBoundary); ProjectToTangentSpace(dIncomingOM, vertexNormal);
		auto dOutgoingOM = mesh.calc_edge_vector(outgoingBoundary); ProjectToTangentSpace(dOutgoingOM, vertexNormal);
		dIncomingOM.normalize();
		dOutgoingOM.normalize();

		Eigen::Vector3f dIncoming(dIncomingOM[0], dIncomingOM[1], dIncomingOM[2]);
		Eigen::Vector3f dOutgoing(dOutgoingOM[0], dOutgoingOM[1], dOutgoingOM[2]);
		Eigen::Vector3f normal(vertexNormal[0], vertexNormal[1], vertexNormal[2]);

		//calculate the signed angle between the incoming and outgoing boundary edges
		float angle = std::atan2(dOutgoing.cross(dIncoming).dot(normal), dIncoming.dot(dOutgoing));
		Eigen::AngleAxis<float> rotation(angle / 2, normal);

		//calculate the virtual continuation as the angle bisector
		auto directionToBoundary = rotation * dOutgoing;

		//calculate the dot product
		return edgeVector[0] * directionToBoundary.x() + edgeVector[1] * directionToBoundary.y() + edgeVector[2] * directionToBoundary.z();
	}
	else
		return -std::numeric_limits<float>::infinity();
}

void FindPossibleContinuations(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, std::vector<PathContinuation>& outContinuations, const float minCosDeviationAngle, bool allowTurnsOnRegularVertices)
{
	auto currentVertex = mesh.to_vertex_handle(edge);
	outContinuations.clear();

	OpenMesh::Vec3f vertexNormal;
	CalcVertexNormalAngleWeights(mesh, currentVertex, vertexNormal);	

	auto dir1 = mesh.calc_edge_vector(edge);
	ProjectToTangentSpace(dir1, vertexNormal); dir1.normalize();

	//outSide is true if h1 comes before h2 in an ordering of the edges around vertexNormal (having angle less than 180°)
	auto calcDot = [&](HEMesh::HalfedgeHandle h2)
	{
		auto dir2 = mesh.calc_edge_vector(h2);
		ProjectToTangentSpace(dir2, vertexNormal);
		dir2.normalize();

		return OpenMesh::dot(dir1, dir2);
	};

	bool isSingularity = false;
	if (mesh.is_manifold(currentVertex))
		isSingularity = IsSingularityManifoldVertex(currentVertex, mesh);
	else
	{
		auto contextEdge = edge;
		CirculateBackwardUntil<true>(contextEdge, mesh, [](HEMesh::HalfedgeHandle) {return false; });
		isSingularity = IsSingularityNonManifoldVertex(contextEdge, mesh);
	}

	bool madeCycle = false;
	if (allowTurnsOnRegularVertices || isSingularity)
	{
		//irregular vertex or we want to allow turns everywhere
		//find the edges that have the allowed deviation from being straight

		auto h = edge;
		CirculateForwardUntil<false>(h, mesh, [&](HEMesh::HalfedgeHandle h)
		{
			float dot = calcDot(h);

			if (dot >= minCosDeviationAngle)
				outContinuations.emplace_back(h, dot);

			return false;
		});	

		auto toBoundaryScore = ContinuationScoreToBoundary(mesh, edge, vertexNormal);
		if (!std::isinf(toBoundaryScore) && toBoundaryScore > minCosDeviationAngle)
		{
			outContinuations.emplace_back(HEMesh::HalfedgeHandle(), toBoundaryScore);
		}		
	}
	else
	{
		//continuation over a regular vertex

		HEMesh::HalfedgeHandle h;
		if (RegularContinuation(mesh, edge, h))
		{
			auto dot = calcDot(h);

			outContinuations.emplace_back(h, dot);
		}
		else
		{
			auto toBoundaryScore = ContinuationScoreToBoundary(mesh, edge, vertexNormal);
			outContinuations.emplace_back(HEMesh::HalfedgeHandle(-1), toBoundaryScore);
		}
	}
}

bool FindBestContinuation(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, PathContinuation& outContinuation, const float minCosDeviationAngle, bool allowTurnsOnRegularVertices)
{
	//find the continuation with the best score

	std::vector<PathContinuation> continuations;
	FindPossibleContinuations(mesh, edge, continuations, minCosDeviationAngle, allowTurnsOnRegularVertices);
	int bestContinuation = -1;
	float bestScore = -std::numeric_limits<float>::infinity();
	for (int i = 0; i < continuations.size(); ++i)
		if (continuations[i].score > bestScore)
		{
			bestContinuation = i;
			bestScore = continuations[i].score;
		}
	if (bestContinuation == -1)
		return false;
	outContinuation = continuations[bestContinuation];
	return true;
}

bool RegularContinuation(const HEMesh& mesh, HEMesh::HalfedgeHandle edge, HEMesh::HalfedgeHandle& outContinuation)
{
	auto currentVertex = mesh.to_vertex_handle(edge);

	if (!mesh.is_boundary(currentVertex))
	{
		//regular non-boundary
		outContinuation = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(edge)));
		return true;
	}
	else
	{
		//regular boundary
		if (mesh.is_boundary(edge))
		{
			auto h = edge;
			CirculateBackwardUntil<true>(h, mesh, [](HEMesh::HalfedgeHandle) {return false; });
			outContinuation = h;
			return true;
		}
		else
		{
			auto nextOpposite = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(edge));
			if (!mesh.is_boundary(nextOpposite))
			{
				outContinuation = mesh.next_halfedge_handle(nextOpposite);
				return true;
			}
			else
				return false;
		}
	}
}

void PrintHalfedge(std::ostream& stream, HEMesh::HalfedgeHandle h, const HEMesh& mesh)
{
	stream << "(" << mesh.from_vertex_handle(h) << ") -- " << h << " --> ";
}

void PrintFullHalfedge(std::ostream& stream, HEMesh::HalfedgeHandle h, const HEMesh& mesh)
{
	stream << "(" << mesh.from_vertex_handle(h) << ") -- " << h << " --> (" << mesh.to_vertex_handle(h) << ")";
}

int TurnsBetweenEdges(HEMesh::HalfedgeHandle h1, HEMesh::HalfedgeHandle h2, const HEMesh& mesh, bool useMargin)
{
	if (useMargin && mesh.is_boundary(mesh.opposite_halfedge_handle(h1)) && mesh.is_boundary(mesh.opposite_halfedge_handle(h2)))
		return 0; //add regular margin around mesh

	//Circulate h1 until it is equal to h2 and count the hops
	int hops = 1;
	h1 = mesh.next_halfedge_handle(h1);
	while (h1 != h2)
	{
		h1 = mesh.opposite_halfedge_handle(h1);
		if (mesh.is_boundary(h1))
			++hops;
		h1 = mesh.next_halfedge_handle(h1);
		++hops;
	}
	return 2 - hops;
}

void CalcVertexNormalAngleWeights(const HEMesh& mesh, OpenMesh::VertexHandle _vh, OpenMesh::Vec3f& _n)
{
	//adaptation of OpenMesh::calc_vertex_normal_correct

	_n.vectorize(0.0);
	HEMesh::ConstVertexIHalfedgeIter cvih_it = mesh.cvih_iter(_vh);
	if (!cvih_it.is_valid())
	{//don't crash on isolated vertices
		return;
	}
	OpenMesh::Vec3f in_he_vec;
	mesh.calc_edge_vector(*cvih_it, in_he_vec);
	for (; cvih_it.is_valid(); ++cvih_it)
	{//calculates the sector normal defined by cvih_it and adds it to _n
		if (mesh.is_boundary(*cvih_it))
		{
			continue;
		}
		OpenMesh::HalfedgeHandle out_heh(mesh.next_halfedge_handle(*cvih_it));
		OpenMesh::Vec3f out_he_vec;
		mesh.calc_edge_vector(out_heh, out_he_vec);
		auto angle = acos(OpenMesh::sane_aarg(OpenMesh::dot(in_he_vec, out_he_vec)));
		_n += angle * cross(in_he_vec, out_he_vec).normalized();
		in_he_vec = out_he_vec;
		in_he_vec *= -1;//change the orientation
	}
	_n.normalize();
}