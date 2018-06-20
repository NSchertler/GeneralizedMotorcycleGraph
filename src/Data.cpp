#include "Data.h"

#include "meshio.h"
#include "MotorcycleGraph.h"
#include "MotorcycleOptionsInPatch.h"
#include "FencedRegionRepresentativeCalculator.h"

#include "Parametrization.h"
#include "TargetLengthStrategySimple.h"
#include "TargetLengthStrategyPatchAverage.h"
#include "MultiplierStrategyIterativeRounding.h"
#include "MultiplierStrategyAllInvalid.h"
#include "ArclengthStrategyGurobi.h"
#include "ArclengthStrategyGurobiSeparate.h"
#include "ArclengthStrategySimple.h"
#include "ParametrizationStrategyLSCMForInteriorScaffoldForBoundary.h"

#include <nanogui/opengl.h>
#include <nsessentials/gui/GLBuffer.h>
#include <nsessentials/gui/GLVertexArray.h>
#include "gui/ShaderPool.h"

#include <nsessentials/data/Serialization.h>
#include <nsessentials/util/TimedBlock.h>

#include <nsessentials/util/UnionFind.h>

#include "stb_image_write.h"

#include <RectangleBinPack/MaxRectsBinPack.h>

#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

Data::Data()
	: vertexToSingularity(mesh), vertexToMetaSingularity(mesh)
{ 
	mesh.add_property(isOriginalEdgeProp);
}

void Data::Clear()
{
	texCoords.clear();
	V.resize(3, 0);
	F.clear();
	mesh.clear();
	meshBoundingBox.reset();
	subdivisionInfo.clear();
	singularities.clear();
	vertexToSingularity.clear();
	metaSingularities.clear();
	vertexToMetaSingularity.clear();
	metaSingularitiesDirty = true;
	fencedRegions.clear();
	textureWidth = textureHeight = 0;	
	motorcycleGraph = nullptr;	
	regionWasMergedFrom.clear();
	brokenHalfarcs.clear();
	patchColors.clear();
	connectedComponentsUF.Clear();
	connectedComponents.clear();	
}

void Data::LoadPolygons(const Matrix3Xf& V, const FaceList& F, const std::vector<std::vector<unsigned int>>& subdivisionInfo)
{
	connectedComponentsUF.AddItems(V.cols());

	//Generating half-edge data structure
	{
		nse::util::TimedBlock b("Generating half-edge data structure ..");
		mesh.clear();
		mesh.reserve(V.cols(), 2 * V.cols(), F.size());
		for (int i = 0; i < V.cols(); ++i)
		{
			mesh.add_vertex(HEMesh::Point(V.coeff(0, i), V.coeff(1, i), V.coeff(2, i)));
		}

		meshBoundingBox.reset();
		std::vector<HEMesh::VertexHandle> face;
		std::set<size_t> faceVertices;
		std::vector<HEMesh::FaceHandle> faceHandles(F.size());
		for (int i = 0; i < F.size(); ++i)
		{
			auto& f = F[i];
			face.clear();
			faceVertices.clear();
			for (auto v : f)
			{
				face.push_back(mesh.vertex_handle(v));
				faceVertices.insert(v);

				meshBoundingBox.expand(V.col(v));

				connectedComponentsUF.Merge(f[0], v);
			}
			if (faceVertices.size() != f.size())
			{
				std::cout << "Face has duplicated vertices: ";
				for (auto v : f)
					std::cout << v << ", ";
				std::cout << std::endl;
			}
			faceHandles[i] = mesh.add_face(face);
		}

		this->subdivisionInfo.resize(subdivisionInfo.size());
		for (int i = 0; i < subdivisionInfo.size(); ++i)
		{
			for (auto j : subdivisionInfo[i])
			{
				auto f = faceHandles[j];
				this->subdivisionInfo[i].emplace_back(f, mesh.vertex_handle(F[j][1]));
			}
		}
	}

	//Find connected components
	for (unsigned int i = 0; i < V.cols(); ++i)
	{
		if (connectedComponentsUF.GetRepresentative(i) == i)
			connectedComponents.push_back(i);
	}

	//calculate average edge length
	averageEdgeLength = 0;
	for (auto e : mesh.edges())
		averageEdgeLength += mesh.calc_edge_length(e);
	averageEdgeLength /= mesh.n_edges();

	//Find singularities
	FindSingularities();
}

void Data::LoadMesh(const std::string& filename, bool mergeTriangulatedQuads, float mergeAngleThreshold)
{
	Clear();

	::LoadMesh(filename, V, F);

	if (mergeTriangulatedQuads)
		MergeTriangulatedQuads(F, V, std::cos(mergeAngleThreshold));

	bool hasNonQuads = false;
	for (auto& f : F)
		if (f.size() != 4)
		{
			hasNonQuads = true;
			break;
		}

	size_t originalVertices = V.cols();

	std::vector<std::vector<unsigned int>> subdivisionInfo;
	if (hasNonQuads)
		CatmullClarkSubdivide(F, V, subdivisionInfo);

	LoadPolygons(V, F, subdivisionInfo);

	//Set isOriginal property
	for (auto& e : mesh.edges())
	{
		auto h = mesh.halfedge_handle(e, 0);
		mesh.property(isOriginalEdgeProp, e) = mesh.from_vertex_handle(h).idx() < originalVertices || mesh.to_vertex_handle(h).idx() < originalVertices;
	}
}

void Data::LoadMeshForVisualization(const std::string& filename, bool scaleTexCoords)
{
	Clear();

	Matrix2Xf T;
	FaceList FT;
	::load_obj(filename, V, F, T, FT);

	std::vector<std::vector<unsigned int>> subdivisionInfo;
	LoadPolygons(V, F, subdivisionInfo);

	for (auto& e : mesh.edges())
		mesh.property(isOriginalEdgeProp, e) = true;

	GenerateNewMotorcycleGraph();

	//Find patches
	if (!FT.empty())
	{
		//Find connected components
		nse::util::UnionFind uf;
		uf.AddItems(T.cols());
		for (auto& ft : FT)
		{
			for (int i = 1; i < ft.size(); ++i)
				uf.Merge(ft[0], ft[i]);
		}

		//set patch indices
		std::map<unsigned int, size_t> repToPatchIdx;
		size_t nextPatchIdx = 0;
		for (int i = 0; i < uf.size(); ++i)
			if (uf.GetRepresentative(i) == i)
				repToPatchIdx[i] = nextPatchIdx++;

		std::vector<MotorcycleGraph::HalfArc> halfarcs;

		std::vector<std::set<HEMesh::HalfedgeHandle>> edgesInPatch(nextPatchIdx);

		std::vector<std::set<HEMesh::FaceHandle>> facesPerPatch(nextPatchIdx);
		for (int i = 0; i < FT.size(); ++i)
		{
			auto patchIdx = repToPatchIdx.at(uf.GetRepresentative(FT[i][0]));
			auto f = HEMesh::FaceHandle(i);
			facesPerPatch[patchIdx].insert(f);
			for (auto h : mesh.fh_range(f))
				edgesInPatch[patchIdx].insert(h);
		}

		//find boundary
		std::map<HEMesh::HalfedgeHandle, size_t> edgeToHalfarc;
		std::vector<std::vector<std::vector<size_t>>> sidesPerPatch(nextPatchIdx);
		for (int patchIdx = 0; patchIdx < nextPatchIdx; ++patchIdx)
		{
			sidesPerPatch[patchIdx].emplace_back();
			for (auto h : edgesInPatch[patchIdx])
			{
				auto opph = mesh.opposite_halfedge_handle(h);
				if (edgesInPatch[patchIdx].find(opph) == edgesInPatch[patchIdx].end())
				{
					//boundary
					size_t arcIdx;
					auto arcIt = edgeToHalfarc.find(h);
					if (arcIt == edgeToHalfarc.end())
					{
						edgeToHalfarc[h] = halfarcs.size();
						halfarcs.emplace_back();
						edgeToHalfarc[opph] = halfarcs.size();
						halfarcs.emplace_back();
						arcIdx = halfarcs.size() - 2;
					}
					else
						arcIdx = arcIt->second;

					halfarcs[arcIdx].face = patchIdx;
					sidesPerPatch[patchIdx].back().push_back(arcIdx);
				}
			}
		}

		std::vector<TexturePatch> patches;
		patches.reserve(nextPatchIdx);
		for (int i = 0; i < nextPatchIdx; ++i)
		{
			patches.emplace_back(*motorcycleGraph);
			std::vector<HEMesh::HalfedgeHandle> empty;
			patches.back().PrepareBuild(std::move(sidesPerPatch[i]), std::move(facesPerPatch[i]), std::move(empty));
		}

		motorcycleGraph->LoadExternal(std::move(patches), std::move(halfarcs));

		textureWidth = 1024;
		textureHeight = 1024;

		Eigen::Vector2f scale(1, 1);
		if (scaleTexCoords)
		{
			scale.x() = textureWidth / averageEdgeLength;
			scale.y() = textureHeight / averageEdgeLength;
		}
		for (int i = 0; i < FT.size(); ++i)
		{
			auto& ft = FT[i];
			auto f = HEMesh::FaceHandle(i);
			auto patchIdx = repToPatchIdx.at(uf.GetRepresentative(ft[0]));
			auto& patch = patches[patchIdx];

			auto h = mesh.halfedge_handle(f);
			auto& texCoords = AccessTexCoords(patchIdx);

			//Align the face halfedge to the face structure
			while (mesh.to_vertex_handle(h).idx() != F[i][0])
				h = mesh.next_halfedge_handle(h);

			for (int j = 0; j < ft.size(); ++j)
			{
				try
				{
					texCoords.TexCoordAtToVertex(h, mesh) = scale.cwiseProduct(T.col(ft[j]));
					h = mesh.next_halfedge_handle(h);
				}
				catch (...)
				{
				}
			}
		}
	}
}

void Data::TangentialSmooth()
{
	OpenMesh::Smoother::JacobiLaplaceSmootherT<HEMesh> smoother(mesh);
	smoother.initialize(OpenMesh::Smoother::JacobiLaplaceSmootherT<HEMesh>::Tangential, OpenMesh::Smoother::JacobiLaplaceSmootherT<HEMesh>::C1);
	smoother.smooth(100);
#pragma omp parallel for
	for (int i = 0; i < V.cols(); ++i)
		V.col(i) = ToEigenVector(mesh.point(mesh.vertex_handle(i)));
}

void Data::SaveMeshPLY(const std::string& file) const
{
	//find out which vertices are referenced
	size_t nextVertex = 0;
	std::map<HEMesh::VertexHandle, size_t> vertexToIndex;
	std::vector<HEMesh::VertexHandle> referencedVertices;
	for (auto f : mesh.faces())
	{
		for (auto v : mesh.fv_range(f))
		{
			auto inserted = vertexToIndex.insert(std::make_pair(v, nextVertex));
			if (inserted.second)
			{
				++nextVertex;
				referencedVertices.push_back(v);
			}
		}
	}

	std::ofstream ply(file, std::ios::binary);

	ply << "ply" << std::endl;
	ply << "format binary_little_endian 1.0" << std::endl;
	ply << "element vertex " << referencedVertices.size() << std::endl;
	ply << "property float32 x" << std::endl;
	ply << "property float32 y" << std::endl;
	ply << "property float32 z" << std::endl;
	ply << "element face " << mesh.n_faces() << std::endl;
	ply << "property list uint8 int32 vertex_index" << std::endl;
	ply << "end_header" << std::endl;


	for (auto v : referencedVertices)
	{
		auto& pos = mesh.point(v);
		ply.write(reinterpret_cast<const char*>(pos.data()), sizeof(OpenMesh::Vec3f));
	}

	std::vector<int32_t> faceIndices;
	for (auto f : mesh.faces())
	{
		uint8_t count = mesh.valence(f);
		faceIndices.resize(0);
		for (auto v : mesh.fv_range(f))
			faceIndices.push_back(vertexToIndex.at(v));
		ply.write(reinterpret_cast<const char*>(&count), 1);
		ply.write(reinterpret_cast<const char*>(faceIndices.data()), count * sizeof(int32_t));
	}
	ply.close();
}

void Data::SaveMeshOBJ(const std::string& file, bool unsubdivide) const
{
	//unsubdivide
	std::vector<bool> keepEdge(mesh.n_edges(), !unsubdivide);

	//Keep all original edges
	for (auto& e : mesh.edges())
	{
		if (mesh.property(isOriginalEdgeProp, e))
			keepEdge[e.idx()] = true;
	}

	//Keep all edges traversed by a motorcycle
	std::queue<HEMesh::EdgeHandle> unsubdivideQueue;
	for (auto& patch : motorcycleGraph->Patches())
	{
		for (auto& side : patch.PatchSides())
		{
			for (auto arcIdx : side)
			{
				auto& arc = motorcycleGraph->Halfarcs()[arcIdx];
				for (auto p : arc)
				{
					auto e = mesh.edge_handle(motorcycleGraph->MotorcycleHalfedge(p));
					unsubdivideQueue.push(e);
				}
			}
		}
	}

	while (!unsubdivideQueue.empty())
	{
		auto e = unsubdivideQueue.front();
		unsubdivideQueue.pop();
		if (keepEdge[e.idx()])
			continue;

		keepEdge[e.idx()] = true;
		continue; //do not propagate, allot T-junctions
		for (int i = 0; i < 2; ++i)
		{
			auto startH = mesh.halfedge_handle(e, i);
			if (!mesh.is_boundary(mesh.to_vertex_handle(startH)))
			{
				auto valence = mesh.valence(mesh.to_vertex_handle(startH));
				if (valence == 4)
					unsubdivideQueue.push(mesh.edge_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(startH)))));
				else
				{
					PathContinuation continuation;
					auto canContinue = FindBestContinuation(mesh, startH, continuation, -1.0f);
					assert(canContinue);
					unsubdivideQueue.push(mesh.edge_handle(continuation.edge));
				}
			}
		}
	}

	//Record a correspondence between mesh faces and texture patch indices
	std::vector<size_t> faceBelongsToPatch(mesh.n_faces(), (size_t)-1);
	for (int iPatch = 0; iPatch < motorcycleGraph->Patches().size(); ++iPatch)
	{
		auto& patch = motorcycleGraph->Patches()[iPatch];
		for (auto& face : patch.Faces())
		{
			faceBelongsToPatch[face.idx()] = iPatch;
		}
	}

	//If some of the subdivided faces do not exist (due to removal of non-manifold configurations),
	//keep all interior edges
	for (auto& subdivInfo : subdivisionInfo)
	{
		auto hasKilledFaces = false;
		for (auto& f : subdivInfo)
		{
			if (!f.subdividedFace.is_valid())
				hasKilledFaces = true;
		}
		if (hasKilledFaces)
		{
			for (auto& f : subdivInfo)
			{
				if (f.subdividedFace.is_valid())
				{
					for (auto e : mesh.fe_range(f.subdividedFace))
						keepEdge[e.idx()] = true;
				}
			}
		}
	}

	//pairs of patch idx and vertices (the from-vertices of the halfedges)
	std::vector<std::pair<size_t, std::vector<HEMesh::HalfedgeHandle>>> outputFaces; 

	//extract all faces from an original (unsubdivided) face represented by its boundary
	auto extractFaces = [&](const std::vector<HEMesh::HalfedgeHandle>& originalBoundary)
	{
		std::set<HEMesh::HalfedgeHandle> originalBoundaryHandled;
		//trace face outlines over all edges that are kept, starting from each of the original boundaries
		for (int i = 0; i < originalBoundary.size(); ++i)
		{
			auto startH = originalBoundary[i];
			if (originalBoundaryHandled.find(startH) != originalBoundaryHandled.end())
				continue;
			originalBoundaryHandled.insert(startH);
			auto h = startH;

			auto f = mesh.face_handle(h);
			auto patchIdx = faceBelongsToPatch[f.idx()];
			if (patchIdx == (size_t)-1)
				continue;
			auto& patch = motorcycleGraph->Patches()[patchIdx];
			//add a new empty output face
			outputFaces.emplace_back(patchIdx, std::vector<HEMesh::HalfedgeHandle>());
			do
			{
				//add the current halfedge to the current output face
				outputFaces.back().second.push_back(h);

				//advance the current halfedge
				h = mesh.next_halfedge_handle(h);
				
				//if we do not want to keep the current halfedge, advance it even more
				bool needNext = false;
				while (!keepEdge[mesh.edge_handle(h).idx()])
				{
					needNext = true;
					h = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(h));
					originalBoundaryHandled.insert(h);
				}
				if (needNext)
					h = mesh.next_halfedge_handle(h);
			} while (h != startH);
		}
	};

	for (auto &subdiv : subdivisionInfo)
	{
		//find the original boundary
		std::vector<HEMesh::HalfedgeHandle> originalBoundary;
		for (auto& f : subdiv)
		{
			if (!f.subdividedFace.is_valid())
				continue;
			for (auto h : mesh.fh_range(f.subdividedFace))
			{
				if (mesh.property(isOriginalEdgeProp, mesh.edge_handle(h)) && mesh.property(isOriginalEdgeProp, mesh.edge_handle(mesh.prev_halfedge_handle(h))))
					originalBoundary.push_back(h);
			}
		}

		extractFaces(originalBoundary);
	}

	//If the mesh has not been subdivided, start face extraction from every face
	if (subdivisionInfo.empty())
	{
		std::vector<HEMesh::HalfedgeHandle> originalBoundary(1);
		for (auto f : mesh.faces())
		{
			originalBoundary.front() = mesh.halfedge_handle(f);
			extractFaces(originalBoundary);
		}
	}

	std::map<HEMesh::VertexHandle, size_t> vertexToIndex;
	size_t nextVertex = 1;

	std::map<std::pair<size_t, HEMesh::HalfedgeHandle>, size_t> patchVertexPairToTextureIndex;
	size_t nextTextureIndex = 1;

	std::ofstream obj(file);
	obj << "OBJ" << std::endl;

	std::vector<Eigen::Vector4f> vboPositions;
	std::vector<Eigen::Vector4f> vboColors;
	std::vector<uint32_t> faceIndices;
	std::vector<uint32_t> wireframeIndices;

	for (auto& face : outputFaces)
	{
		auto& patch = motorcycleGraph->Patches()[face.first];
		auto& texCoordStorage = AccessTexCoords(face.first);
		auto& halfedges = face.second;
		std::vector<std::pair<size_t, size_t>> objPairs;
		for (auto& h : halfedges)
		{
			auto v = mesh.from_vertex_handle(h);
			size_t vertexIndex;
			auto insertedVertex = vertexToIndex.emplace(v, nextVertex);
			if (insertedVertex.second)
			{
				vertexIndex = nextVertex;
				++nextVertex;
				auto p = mesh.point(v);
				HEMesh::HalfedgeHandle removedIncidentEdge;
				for (auto h : mesh.voh_range(v))
				{
					auto e = mesh.edge_handle(h);
					if (!keepEdge[e.idx()])
						removedIncidentEdge = h;
				}
				//Project t-junction onto surface of the kept mesh (to avoid cracks and overlaps)
				if (removedIncidentEdge.is_valid())
				{
					p = 0.5f * (mesh.point(mesh.from_vertex_handle(mesh.prev_halfedge_handle(removedIncidentEdge)))
						+ mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(removedIncidentEdge)))));
				}
				obj << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			}
			else
				vertexIndex = insertedVertex.first->second;

			size_t textureIndex;
			auto insertedTexture = patchVertexPairToTextureIndex.emplace(std::make_pair(face.first, h), nextTextureIndex);
			if (insertedTexture.second)
			{
				textureIndex = nextTextureIndex;
				++nextTextureIndex;
				auto& t = texCoordStorage.TexCoordAtFromVertex(h, mesh);
				obj << "vt " << t.x() << " " << t.y() << std::endl;
				vboPositions.push_back(Eigen::Vector4f(t.x() * 2 - 1 + 1.0f / textureWidth, t.y() * 2 - 1 + 1.0f / textureHeight, 0, 1.0f));
				vboColors.push_back(GetSegmentColor(patchColors[face.first]));
			}
			else
				textureIndex = insertedTexture.first->second;

			objPairs.emplace_back(vertexIndex, textureIndex);
		}
		obj << "f";
		for (auto& pair : objPairs)
		{
			obj << " " << pair.first << "/" << pair.second;
		}
		for (int i = 0; i < objPairs.size(); ++i)
		{
			wireframeIndices.push_back(objPairs[i].second - 1);
			wireframeIndices.push_back(objPairs[(i + 1) % objPairs.size()].second - 1);
			if (i != 0)
			{
				faceIndices.push_back(objPairs.front().second - 1);
				faceIndices.push_back(objPairs[i].second - 1);
				faceIndices.push_back(objPairs[(i + 1) % objPairs.size()].second - 1);
			}
		}
			
		obj << std::endl;
	}

	obj.close();

	if(textureWidth != 0)
	{
		nse::util::TimedBlock b("Rendering texture ..");
	
		auto& shader = ShaderPool::Instance()->SimpleShader;
		shader.bind();
		nse::gui::GLVertexArray wireframeVAO;
		nse::gui::GLBuffer posBuffer(nse::gui::VertexBuffer), colorBuffer(nse::gui::VertexBuffer);
		nse::gui::GLBuffer wireframeIndexBuffer(nse::gui::IndexBuffer);
		wireframeVAO.generate();
		wireframeVAO.bind();

		Eigen::Matrix4f id;
		id.setIdentity();
		shader.setUniform("mvp", id);
		shader.setUniform("visualizeTexCoords", 0);		
		
		posBuffer.uploadData(vboPositions).bindToAttribute("position");
		colorBuffer.uploadData(vboColors).bindToAttribute("color");			
		
		GLuint textureFBO, textureColor;
		glGenFramebuffers(1, &textureFBO);
		glBindFramebuffer(GL_FRAMEBUFFER, textureFBO);

		glGenRenderbuffers(1, &textureColor);
		glBindRenderbuffer(GL_RENDERBUFFER, textureColor);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, textureWidth, textureHeight);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, textureColor);		

		glClearColor(1, 1, 1, 1);
		glClear(GL_COLOR_BUFFER_BIT);

		glViewport(0, 0, textureWidth, textureHeight);
		
		wireframeIndexBuffer.uploadData(faceIndices.size(), 1, sizeof(unsigned int), GL_UNSIGNED_INT, true, reinterpret_cast<uint8_t*>(faceIndices.data()));
		shader.setUniform("useUniformColor", 0);
		glDrawElements(GL_TRIANGLES, faceIndices.size(), GL_UNSIGNED_INT, 0);

		wireframeIndexBuffer.uploadData(wireframeIndices.size(), 1, sizeof(unsigned int), GL_UNSIGNED_INT, true, reinterpret_cast<uint8_t*>(wireframeIndices.data()));
		shader.setUniform("useUniformColor", 1);
		shader.setUniform("color", Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f));
		glDrawElements(GL_LINES, wireframeIndices.size(), GL_UNSIGNED_INT, 0);

		std::vector<unsigned char> pixels(textureWidth * textureHeight * 4);
		glReadPixels(0, 0, textureWidth, textureHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
		
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		glDeleteFramebuffers(1, &textureFBO);
		glDeleteRenderbuffers(1, &textureColor);

		stbi_write_png((file + ".png").c_str(), textureWidth, textureHeight, 4, pixels.data() + (textureHeight - 1) * textureWidth * 4 * sizeof(unsigned char), -textureWidth * 4 * sizeof(unsigned char));
		
	}
}

void Data::FindPatchColors()
{
	patchColors.clear();
	{
		nse::util::TimedBlock b("Finding patch colors ..");
		//Perform greedy graph coloring
		int orientation = 0;
		for (auto& patch : motorcycleGraph->Patches())
		{
			bool colorValid = false;
			int startColor = orientation;
			do
			{
				++orientation;
				if (orientation >= segmentColorCount)
					orientation = 0;
				colorValid = true;
				for (auto& side : patch.PatchSides())
				{
					for (auto arcIdx : side)
					{
						auto oppositeArc = motorcycleGraph->OppositeHalfarc(arcIdx);
						auto f = motorcycleGraph->Halfarcs()[oppositeArc].face;
						if (f >= patchColors.size())
							continue;
						if (patchColors[f] == orientation)
						{
							colorValid = false;
							break;
						}
					}
					if (!colorValid)
						break;
				}
			} while (!colorValid && orientation != startColor);
			patchColors.push_back(orientation);
		}
	}
}


//Triangle rasterization code from http://www.sunshine2k.de/coding/java/TriangleRasterization/TriangleRasterization.html
void FillBottomFlatTriangle(const Eigen::Vector2f& v1, const Eigen::Vector2f& v2, const Eigen::Vector2f& v3, std::vector<Eigen::Vector2i>& result, float epsilon)
{
	double invslope1 = (double)(v2.x() - v1.x()) / (v2.y() - v1.y());
	double invslope2 = (double)(v3.x() - v1.x()) / (v3.y() - v1.y());

	if (v2.x() > v3.x())
		std::swap(invslope1, invslope2);

	double scanlineYTo = std::max(v2.y(), v3.y()) + epsilon;
	int scanlineYFrom = (int)std::ceil(v1.y() - epsilon);

	double curx1 = v1.x() + invslope1 * (scanlineYFrom - v1.y());
	double curx2 = v1.x() + invslope2 * (scanlineYFrom - v1.y());

	for (int scanlineY = scanlineYFrom; scanlineY <= scanlineYTo; scanlineY++)
	{
		for (int x = (int)std::ceil(curx1 - epsilon); x <= curx2 + epsilon; ++x)
			result.push_back(Eigen::Vector2i(x, scanlineY));
		curx1 += invslope1;
		curx2 += invslope2;
	}
}

void FillTopFlatTriangle(const Eigen::Vector2f& v1, const Eigen::Vector2f& v2, const Eigen::Vector2f& v3, std::vector<Eigen::Vector2i>& result, float epsilon)
{
	double invslope1 = (double)(v3.x() - v1.x()) / (v3.y() - v1.y());
	double invslope2 = (double)(v3.x() - v2.x()) / (v3.y() - v2.y());

	if (v1.x() > v2.x())
		std::swap(invslope1, invslope2);

	double scanLineYTo = std::min(v1.y(), v2.y()) - epsilon;
	int scanLineYFrom = std::floor(v3.y() + epsilon);

	double curx1 = v3.x() + invslope1 * (scanLineYFrom - v3.y());
	double curx2 = v3.x() + invslope2 * (scanLineYFrom - v3.y());

	for (int scanlineY = scanLineYFrom; scanlineY >= scanLineYTo; scanlineY--)
	{
		for (int x = (int)std::ceil(curx1 - epsilon); x <= curx2 + epsilon; ++x)
			result.push_back(Eigen::Vector2i(x, scanlineY));
		curx1 -= invslope1;
		curx2 -= invslope2;
	}
}


void LatticePointsInTriangle(const Eigen::Vector2f& v1, const Eigen::Vector2f& v2, const Eigen::Vector2f& v3, std::vector<Eigen::Vector2i>& result, float epsilon)
{
	//sort vertices
	std::vector<const Eigen::Vector2f*> vSorted = { &v1, &v2, &v3 };
	std::sort(vSorted.begin(), vSorted.end(), [](const Eigen::Vector2f* v1, const Eigen::Vector2f* v2) {return v1->y() < v2->y(); });

	if (std::abs(vSorted[1]->y() - vSorted[2]->y()) < epsilon)
		FillBottomFlatTriangle(*vSorted[0], *vSorted[1], *vSorted[2], result, epsilon);
	else if (std::abs(vSorted[0]->y() - vSorted[1]->y()) < epsilon)
		FillTopFlatTriangle(*vSorted[0], *vSorted[1], *vSorted[2], result, epsilon);
	else
	{
		Eigen::Vector2f v4((vSorted[0]->x() + ((vSorted[1]->y() - vSorted[0]->y()) / (vSorted[2]->y() - vSorted[0]->y())) * (vSorted[2]->x() - vSorted[0]->x())), vSorted[1]->y());
		FillBottomFlatTriangle(*vSorted[0], *vSorted[1], v4, result, epsilon);
		FillTopFlatTriangle(*vSorted[1], v4, *vSorted[2], result, epsilon);
	}
}

template <typename T>
T Bilinear(const T& v0, const T& v1, const T& v2, const T& v3, const Eigen::Vector2f& param)
{
	return (1 - param.x()) * ((1 - param.y()) * v0 + param.y() * v2) + param.x() * ((1 - param.y()) * v1 + param.y() * v3);
}

float Cross2(const Eigen::Vector2f& v1, const Eigen::Vector2f& v2)
{
	return v1.x() * v2.y() - v1.y() * v2.x();
}

//Calculates the parameters (s, t) such that (1 - s) * [(1 - t) * v0 + t * v2] + s * [(1 - t) * v1 + t * v3] == searchPoint
//v0 = bottom left, v1 = bottom right, v2 = top left, v3 = top right
Eigen::Vector2f InverseBilinear(const Eigen::Vector2f& v0, const Eigen::Vector2f& v1, const Eigen::Vector2f& v2, const Eigen::Vector2f& v3, const Eigen::Vector2f& searchPoint)
{
	//maths from https://stackoverflow.com/a/812077/1210053
	float A = Cross2(v0 - searchPoint, v0 - v2);
	float B = 0.5f * (Cross2(v0 - searchPoint, v1 - v3) + Cross2(v1 - searchPoint, v0 - v2));
	float C = Cross2(v1 - searchPoint, v1 - v3);

	float s, t;

	float denom = A - 2 * B + C;
	if (std::abs(denom) < 0.00001f)
		s = A / (A - C);
	else
	{
		float root = sqrt(B * B - A * C);
		float s1 = ((A - B) + root) / denom;
		float s2 = ((A - B) - root) / denom;
		float s1DistanceToQuad, s2DistanceToQuad;
		if (s1 < 0)
			s1DistanceToQuad = -s1;
		else if (s1 > 1)
			s1DistanceToQuad = s1 - 1;
		else
			s1DistanceToQuad = 0;

		if (s2 < 0)
			s2DistanceToQuad = -s2;
		else if (s2 > 1)
			s2DistanceToQuad = s2 - 1;
		else
			s2DistanceToQuad = 0;

		if (s1DistanceToQuad < s2DistanceToQuad)
			s = s1;
		else
			s = s2;
	}

	float tDenom1 = (1 - s) * (v0.x() - v2.x()) + s * (v1.x() - v3.x());
	float tDenom2 = (1 - s) * (v0.y() - v2.y()) + s * (v1.y() - v3.y());
	if (std::abs(tDenom1) > std::abs(tDenom2))
		t = ((1 - s) * (v0.x() - searchPoint.x()) + s * (v1.x() - searchPoint.x())) / tDenom1;
	else
		t = ((1 - s) * (v0.y() - searchPoint.y()) + s * (v1.y() - searchPoint.y())) / tDenom2;
	return Eigen::Vector2f(s, t);
}


void Data::Remesh()
{
	//TODO: repair
	/*float epsilon = 0.01f;

	struct Vector2iComparer
	{
		bool operator()(const Eigen::Vector2i& first, const Eigen::Vector2i& second) const
		{
			if (first.x() != second.x())
				return first.x() < second.x();
			return first.y() < second.y();
		}
	};

	std::ofstream obj("remeshed.obj");
	obj << "OBJ" << std::endl;
	size_t nextVertexIndex = 1;

	for (auto& patch : motorcycleGraph->patches)
	{

		std::map<Eigen::Vector2i, size_t, Vector2iComparer> uvToVertexIndex;
		std::vector<Eigen::Vector2i> latticePoints;

		for (auto f : patch.faces)
		{
			latticePoints.clear();
			HEMesh::VertexHandle v[4];
			size_t vertexIndexInPatch[4];
			int i = 0;
			for (auto _v : mesh.fv_range(f))
			{
				v[i] = _v;
				vertexIndexInPatch[i] = patch.verticesToIndex.at(_v);
				++i;
			}
			LatticePointsInTriangle(patch.textureCoordinates[vertexIndexInPatch[0]], patch.textureCoordinates[vertexIndexInPatch[1]], patch.textureCoordinates[vertexIndexInPatch[2]], latticePoints, epsilon);
			LatticePointsInTriangle(patch.textureCoordinates[vertexIndexInPatch[0]], patch.textureCoordinates[vertexIndexInPatch[2]], patch.textureCoordinates[vertexIndexInPatch[3]], latticePoints, epsilon);
			for (auto& texCoord : latticePoints)
			{
				if (uvToVertexIndex.find(texCoord) != uvToVertexIndex.end())
					continue;
				auto param = InverseBilinear(patch.textureCoordinates[vertexIndexInPatch[0]], patch.textureCoordinates[vertexIndexInPatch[1]],
					patch.textureCoordinates[vertexIndexInPatch[3]], patch.textureCoordinates[vertexIndexInPatch[2]], texCoord.cast<float>());
				auto p = Bilinear(ToEigenVector(mesh.point(v[0])), ToEigenVector(mesh.point(v[1])), ToEigenVector(mesh.point(v[3])), ToEigenVector(mesh.point(v[2])), param);
				obj << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
				uvToVertexIndex[texCoord] = nextVertexIndex++;
			}
		}
		for (auto & entry : uvToVertexIndex)
		{
			auto rightIt = uvToVertexIndex.find(entry.first + Eigen::Vector2i(1, 0));
			if (rightIt == uvToVertexIndex.end())
				continue;
			auto topIt = uvToVertexIndex.find(entry.first + Eigen::Vector2i(0, 1));
			if (topIt == uvToVertexIndex.end())
				continue;
			auto diagonalIt = uvToVertexIndex.find(entry.first + Eigen::Vector2i(1, 1));
			if (diagonalIt == uvToVertexIndex.end())
				continue;
			obj << "f " << entry.second << " " << rightIt->second << " " << diagonalIt->second << " " << topIt->second << std::endl;
		}
	}
	obj.close();*/
}

void Data::RenderMeshColorsToTexture(const std::string& colorsFilename)
{
	//TODO: repair
	//std::vector<size_t> quadSubdivs;
	//std::vector<size_t> triSubdivs;

	//for (int i = 0; i < subdivisionInfo.size(); ++i)
	//{
	//	auto& originalFace = subdivisionInfo[i];

	//	if (originalFace.size() == 4)
	//		quadSubdivs.push_back(i);
	//	if (originalFace.size() == 3)
	//		triSubdivs.push_back(i);
	//}

	//GLuint textureFBO, textureColor;
	//glGenFramebuffers(1, &textureFBO);
	//glBindFramebuffer(GL_FRAMEBUFFER, textureFBO);
	//glGenRenderbuffers(1, &textureColor);
	//glBindRenderbuffer(GL_RENDERBUFFER, textureColor);

	//glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, textureWidth, textureHeight);
	//glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, textureColor);

	//nse::gui::GLVertexArray vao;
	//vao.generate();
	//vao.bind();

	//nse::gui::GLBuffer quadBuffer(nse::gui::ShaderStorageBuffer);
	//nse::gui::GLBuffer triBuffer(nse::gui::ShaderStorageBuffer);
	//nse::gui::GLBuffer colorBuffer(nse::gui::ShaderStorageBuffer);

	//struct QuadData
	//{
	//	Eigen::Vector2f cornerTexCoords[4];
	//	Eigen::Vector2f texCoords[9];
	//	uint32_t cPtr;
	//	uint32_t padding;

	//	QuadData()
	//	{
	//		for (int i = 0; i < 4; ++i)
	//			cornerTexCoords[i].setConstant(std::numeric_limits<float>::quiet_NaN());
	//		for (int i = 0; i < 4; ++i)
	//			texCoords[i].setConstant(std::numeric_limits<float>::quiet_NaN());
	//	}
	//};
	//struct TriData
	//{
	//	Eigen::Vector2f cornerTexCoords[3];
	//	Eigen::Vector2f edgeTexCoords[3];
	//	uint32_t cPtr;
	//	uint32_t padding;
	//};

	//std::vector<QuadData> quads;
	//std::vector<TriData> tris;
	//std::vector<Eigen::Vector4f> colors;

	//std::ifstream colorFile(colorsFilename, std::ios::binary);
	//uint32_t R;
	//colorFile.read(reinterpret_cast<char*>(&R), sizeof(uint32_t));
	//uint32_t quadCount, triCount;
	//colorFile.read(reinterpret_cast<char*>(&quadCount), sizeof(uint32_t));
	//colorFile.read(reinterpret_cast<char*>(&triCount), sizeof(uint32_t));
	//if (quadCount != quadSubdivs.size())
	//{
	//	std::cout << "There are " << quadCount << " quads in the color file but " << quadSubdivs.size() << " in the current mesh." << std::endl;
	//	return;
	//}
	//if (triCount != triSubdivs.size())
	//{
	//	std::cout << "There are " << triCount << " triangles in the color file but " << triSubdivs.size() << " in the current mesh." << std::endl;
	//	return;
	//}
	//auto texelsPerQuad = (R + 1) * (R + 1);
	//auto texelsPerTri = 3 * R*(2 + R) / 4 + 1; //assuming even R
	//colors.resize(quadCount * texelsPerQuad + triCount * texelsPerTri);
	//colorFile.read(reinterpret_cast<char*>(colors.data()), colors.size() * sizeof(Eigen::Vector4f));
	//colorFile.close();

	//std::vector<size_t> faceBelongsToPatch(mesh.n_faces(), (size_t)-1);
	//for (int iPatch = 0; iPatch < motorcycleGraph->patches.size(); ++iPatch)
	//{
	//	auto& patch = motorcycleGraph->patches[iPatch];
	//	for (auto& face : patch.faces)
	//	{
	//		faceBelongsToPatch[face.idx()] = iPatch;
	//	}
	//}

	//size_t nextCPtr = 0;
	//for (auto subdivInfo : quadSubdivs)
	//{
	//	auto& subdiv = subdivisionInfo[subdivInfo];

	//	std::vector<size_t> involvedPatches;
	//	for (auto& f : subdiv)
	//	{
	//		if (!f.subdividedFace.is_valid())
	//			continue;
	//		auto patch = faceBelongsToPatch.at(f.subdividedFace.idx());
	//		if (patch != (size_t)-1 && std::find(involvedPatches.begin(), involvedPatches.end(), patch) == involvedPatches.end())
	//			involvedPatches.push_back(patch);
	//	}

	//	auto cPtr = nextCPtr;
	//	nextCPtr += texelsPerQuad;

	//	for (auto patchIdx : involvedPatches)
	//	{
	//		auto& patch = motorcycleGraph->patches[patchIdx];

	//		quads.emplace_back();
	//		quads.back().cPtr = cPtr;

	//		Eigen::Vector2f center;
	//		Eigen::Vector2f edgeTexCoords[4];
	//		for (int i = 0; i < 4; ++i)
	//		{
	//			auto& f = subdiv[i];
	//			auto vIt = patch.verticesToIndex.find(f.originalCornerVertex);
	//			if (vIt != patch.verticesToIndex.end())
	//				quads.back().cornerTexCoords[i] = patch.textureCoordinates[vIt->second];

	//			if (f.subdividedFace.is_valid())
	//			{
	//				auto h = mesh.halfedge_handle(f.subdividedFace);
	//				while (mesh.to_vertex_handle(h) != f.originalCornerVertex)
	//					h = mesh.next_halfedge_handle(h);

	//				auto prevVIt = patch.verticesToIndex.find(mesh.from_vertex_handle(h));
	//				if (prevVIt != patch.verticesToIndex.end())
	//					edgeTexCoords[(i + 3) % 4] = patch.textureCoordinates[prevVIt->second];

	//				h = mesh.next_halfedge_handle(h);
	//				auto nextVIt = patch.verticesToIndex.find(mesh.to_vertex_handle(h));
	//				if (nextVIt != patch.verticesToIndex.end())
	//					edgeTexCoords[i] = patch.textureCoordinates[nextVIt->second];

	//				h = mesh.next_halfedge_handle(h);
	//				nextVIt = patch.verticesToIndex.find(mesh.to_vertex_handle(h));
	//				if (nextVIt != patch.verticesToIndex.end())
	//					center = patch.textureCoordinates[nextVIt->second];
	//			}
	//		}
	//		quads.back().texCoords[0] = quads.back().cornerTexCoords[0];
	//		quads.back().texCoords[1] = edgeTexCoords[0];
	//		quads.back().texCoords[2] = quads.back().cornerTexCoords[1];

	//		quads.back().texCoords[3] = edgeTexCoords[3];
	//		quads.back().texCoords[4] = center;
	//		quads.back().texCoords[5] = edgeTexCoords[1];

	//		quads.back().texCoords[6] = quads.back().cornerTexCoords[3];
	//		quads.back().texCoords[7] = edgeTexCoords[2];
	//		quads.back().texCoords[8] = quads.back().cornerTexCoords[2];
	//	}
	//}

	//for (auto subdivInfo : triSubdivs)
	//{
	//	auto& subdiv = subdivisionInfo[subdivInfo];

	//	std::vector<size_t> involvedPatches;
	//	for (auto& f : subdiv)
	//	{
	//		if (!f.subdividedFace.is_valid())
	//			continue;
	//		auto patch = faceBelongsToPatch.at(f.subdividedFace.idx());
	//		if (patch != (size_t)-1 && std::find(involvedPatches.begin(), involvedPatches.end(), patch) == involvedPatches.end())
	//			involvedPatches.push_back(patch);
	//	}

	//	auto cPtr = nextCPtr;
	//	nextCPtr += texelsPerTri;

	//	for (auto patchIdx : involvedPatches)
	//	{
	//		auto& patch = motorcycleGraph->patches[patchIdx];

	//		tris.emplace_back();
	//		tris.back().cPtr = cPtr;

	//		for (int i = 0; i < 3; ++i)
	//		{
	//			auto& f = subdiv[i];
	//			Eigen::Vector2f texCoord; texCoord.setConstant(std::numeric_limits<float>::quiet_NaN());
	//			auto vIt = patch.verticesToIndex.find(f.originalCornerVertex);
	//			if (vIt != patch.verticesToIndex.end())
	//				texCoord = patch.textureCoordinates[vIt->second];
	//			tris.back().cornerTexCoords[i] = texCoord;
	//		}
	//	}
	//}

	//colorBuffer.uploadData(sizeof(Eigen::Vector4f) * colors.size(), colors.data());
	//quadBuffer.uploadData(sizeof(QuadData) * quads.size(), quads.data());
	//triBuffer.uploadData(sizeof(TriData) * tris.size(), tris.data());
	//std::cout << sizeof(TriData) << std::endl;

	//glBindFramebuffer(GL_FRAMEBUFFER, textureFBO);
	//glClearColor(1, 1, 1, 0);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	//glViewport(0, 0, textureWidth, textureHeight);

	////render
	//for (int k = 0; k < 2; ++k)
	//{		
	//	auto& shader1 = ShaderPool::Instance()->MeshColorsQuadShader;
	//	shader1.bind();
	//	shader1.setUniform("overflow", k == 0 ? Eigen::Vector2f((float)texturePatchMargin / textureWidth, (float)texturePatchMargin / textureHeight) : Eigen::Vector2f(0.0f, 0.0f));
	//	shader1.setUniform("R", R);
	//	shader1.setUniform("useUniformColor", 0);
	//	shader1.setUniform("visualizeTexCoords", 0);

	//	quadBuffer.bindBufferBase(0);
	//	colorBuffer.bindBufferBase(1);
	//	glPatchParameteri(GL_PATCH_VERTICES, 4);
	//	glDrawArrays(GL_PATCHES, 0, 4 * quads.size());

	//	auto& shader2 = ShaderPool::Instance()->MeshColorsTriShader;
	//	shader2.bind();
	//	shader2.setUniform("R", R);
	//	shader2.setUniform("useUniformColor", 0);
	//	shader2.setUniform("visualizeTexCoords", 0);

	//	triBuffer.bindBufferBase(0);
	//	colorBuffer.bindBufferBase(1);
	//	glPatchParameteri(GL_PATCH_VERTICES, 3);
	//	glDrawArrays(GL_PATCHES, 0, 3 * tris.size());		
	//}
	////end render

	//std::vector<unsigned char> pixels(textureWidth * textureHeight * 4);
	//glReadPixels(0, 0, textureWidth, textureHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
	//
	//glBindFramebuffer(GL_FRAMEBUFFER, 0);
	//glDeleteFramebuffers(1, &textureFBO);
	//glDeleteRenderbuffers(1, &textureColor);

	//stbi_write_png("test.png", textureWidth, textureHeight, 4, pixels.data() + (textureHeight - 1) * textureWidth * 4 * sizeof(unsigned char), -textureWidth * 4 * sizeof(unsigned char));
}

void Data::FindSingularities()
{
	nse::util::TimedBlock b("Finding singularities ..");
	singularities.clear();
	vertexToSingularity.clear();
	for (auto& v : mesh.vertices())
	{
		if (!mesh.halfedge_handle(v).is_valid())
			continue; //unconnected vertex
		size_t edge_valence = 0;

		if (mesh.is_manifold(v))
		{
			int valenceDefect;
			if (IsSingularityManifoldVertex(v, mesh, valenceDefect))
			{
				singularities.push_back(SingularityInfo(v, valenceDefect));
				vertexToSingularity.AccessOrCreateAtManifoldVertex(v) = singularities.size() - 1;
			}
		}
		else
		{
			for (auto h : mesh.voh_range(v))
			{
				if (mesh.is_boundary(h))
				{
					//we have a new surface patch
					int valenceDefect;
					if (IsSingularityNonManifoldVertex(h, mesh, valenceDefect))
					{
						singularities.push_back(SingularityInfo(v, h, valenceDefect));
						vertexToSingularity.AccessOrCreateAtContextEdge(h) = singularities.size() - 1;
					}
				}
			}
		}
	}

	std::cout << singularities.size() << " singularities." << std::endl;

	FindMetaSingularities();

	PrintSingularityStatistics();
}

void Data::FindMetaSingularities()
{
	if (!metaSingularitiesDirty)
		return;
	metaSingularitiesDirty = false;

	nse::util::TimedBlock b("Setting up meta singularities ..");
	metaSingularities.clear();
	vertexToMetaSingularity.clear();

	//for every fenced region, find a representative and emanating motorcycles to make up the meta singularity
	for (auto& fencedRegion : fencedRegions.patches)
	{
		if (!fencedRegion.IsActive())
			continue;
		if (fencedRegion.IsRectilinear())
			continue;
		if (fencedRegion.Loops().front().Degree() < 2)
		{
			//Cannot trace motorcycles for degree < 2; use the original singularities
			fencedRegion.SetActive(false);
			continue;
		}

		//Try to find the emanating motorcycles
		DPFencedRegionRepresentativeCalculator calculator(mesh, fencedRegion, isOriginalEdgeProp);
		auto cost = calculator.CalculateRepresentative();

		if (std::isinf(cost))
		{
			//No representative could be found; use the original singularities
			fencedRegion.SetActive(false);
			continue;
		}
		else
		{
			metaSingularities.emplace_back();
			metaSingularities.back().highPriority = true;
			HEMesh::HalfedgeHandle outgoingEdge;
			for (int i = 0; i < fencedRegion.Loops()[0].Degree(); ++i)
			{
				bool active = calculator.FillWithSolution(i, metaSingularities.back().GetEmanatingMotorcycle(i));
				metaSingularities.back().SetEmanatingMotorcycleActive(i, active);

				if (metaSingularities.back().GetEmanatingMotorcycle(i).size() > 0)
					outgoingEdge = metaSingularities.back().GetEmanatingMotorcycle(i).front();
			}

			vertexToMetaSingularity.AccessOrCreateAtToVertex(mesh.opposite_halfedge_handle(outgoingEdge)) = metaSingularities.size() - 1;
		}
	}

	for (auto& singularity : singularities)
	{
		if (singularity.correspondingRegion != (size_t)-1 && fencedRegions.patches[singularity.correspondingRegion].IsActive())
			continue;

		auto& v = singularity.vertexHandle;
		metaSingularities.emplace_back();
		metaSingularities.back().highPriority = false;
		//find the context edge
		auto context = singularity.contextOutgoingBoundaryEdge.is_valid() ? singularity.contextOutgoingBoundaryEdge : mesh.halfedge_handle(v);
		
		auto h = context;
		//add the emanating motorcycles
		if (mesh.is_boundary(h))
			metaSingularities.back().AddMotorcycle(h, true);

		//emanate motorcycles along each of the incident edges
		h = mesh.opposite_halfedge_handle(h);
		auto circulationResult = CirculateForwardUntil<true>(h, mesh, [&](HEMesh::HalfedgeHandle incidentH)
		{
			metaSingularities.back().AddMotorcycleAtFront(incidentH, true);
			return false;
		});
		if (circulationResult == StoppedByBoundary)
			metaSingularities.back().AddEmptyMotorcycle(); //empty edge for motorcycle to padded boundary

		vertexToMetaSingularity.AccessOrCreateAtToVertex(mesh.opposite_halfedge_handle(context)) = metaSingularities.size() - 1;
	}


	//Check if we have at least one singularity per connected component. Create an artificial one otherwise.
	std::map<unsigned int, bool> ccHasSingularity;
	for (auto cc : connectedComponents)
		ccHasSingularity[cc] = false;

	for (auto sing : vertexToMetaSingularity)
	{
		auto v = mesh.from_vertex_handle(sing.first).idx();
		auto cc = connectedComponentsUF.GetRepresentative(v);
		ccHasSingularity[cc] = true;
	}

	for (unsigned int cc = 0; cc < connectedComponents.size(); ++cc)
	{
		if (ccHasSingularity[connectedComponents[cc]])
			continue;

		//this connected component does not have a singularity, create an artificial one
		//Use an arbitrary interior regular vertex of this connected component
		std::cout << "Adding an artificial singularity to connected component." << std::endl;
		
		for (int i = 0; i < mesh.n_halfedges(); ++i)
		{
			auto h = mesh.halfedge_handle(i);
			auto v = mesh.to_vertex_handle(h);

			if (connectedComponentsUF.GetRepresentative(v.idx()) != connectedComponents[cc])
				continue;

			if (!IsToVertexSingularity(h, mesh) && !mesh.is_boundary(v))
			{
				metaSingularities.emplace_back();
				metaSingularities.back().highPriority = true;
				int i = 0;
				for (auto oh : mesh.voh_range(v))
				{
					metaSingularities.back().GetEmanatingMotorcycle(i).push_back(oh);
					metaSingularities.back().SetEmanatingMotorcycleActive(i, true);
					++i;
				}

				vertexToMetaSingularity.AccessOrCreateAtToVertex(h) = metaSingularities.size() - 1;
				ccHasSingularity[cc] = true;
				break;
			}
		}
	}

	RecordPatchEnteringEdges();
}

void Data::ClassifyAllSingularities()
{
	nse::util::TimedBlock b("Classifying singularities ..");

	fencedRegions.clear();

#pragma omp parallel for
	for (int i = 0; i < singularities.size(); ++i)
	{
		singularities[i].correspondingRegion = (size_t)-1;
	}

	//#pragma omp parallel for reduction(+ : totalSingularitiesWithoutPatch)
	for (int i = 0; i < singularities.size(); ++i)
	{
		if (singularities[i].correspondingRegion != (size_t)-1)
			continue; //singularity is already covered by a patch
		FencedRegion patch(mesh, vertexToSingularity, false);
		auto classification = ClassifySingularity(singularities[i], patch, false);
		if (classification != NoPatch)
		{
			if (singularities[i].correspondingRegion == (size_t)-1)
//#pragma omp critical
			{
				if (singularities[i].correspondingRegion == (size_t)-1)
				{
					//There is no other patch that overlaps this singularity.
					size_t patchIndex = fencedRegions.patches.size();
					for (auto i : patch.CoveredSingularities())
					{
						auto& coveredSingularity = singularities[i];
						if (coveredSingularity.correspondingRegion != (size_t)-1)
							fencedRegions.patches[coveredSingularity.correspondingRegion].SetActive(false);
						coveredSingularity.correspondingRegion = patchIndex;
					}
					patch.SetActive(true);

					fencedRegions.patches.push_back(std::move(patch));
				}
			}
		}
	}

	//micro optimizations to patches
	std::set<HEMesh::FaceHandle> coveredFaces;
	for (auto& patch : fencedRegions.patches)
	{
		if (!patch.IsActive())
			continue;
		coveredFaces.insert(patch.IncludedFaces().begin(), patch.IncludedFaces().end());
	}

	//grow around concave corners if possible	
	for (int patchIdx = 0; patchIdx < fencedRegions.patches.size(); ++patchIdx)
	{
		auto& patch = fencedRegions.patches[patchIdx];

		if (!patch.IsActive())
			continue;
		auto& loop = patch.Loops().front();
		auto concave = loop.ConcaveCornersOutgoing(); //copy on purpose since we are modifying the loops
		for (auto concaveOutgoing : concave)
		{
			if (!concaveOutgoing.is_valid())
				continue;
			//std::cout << concaveOutgoing.idx() << " / " << mesh.n_halfedges() << std::endl;
			auto f = mesh.opposite_face_handle(concaveOutgoing);
			if (!f.is_valid())
				continue;
			bool isOk = true;
			//check if we are not adding a singularity
			for (auto h : mesh.fh_range(f))
			{
				if (IsToVertexSingularity(h, mesh))
					isOk = false;
			}
			isOk = isOk && !loop.AddFaceGeneratesHandle(f, mesh);
			if (!isOk)
				continue;
			if (f.is_valid() && coveredFaces.find(f) == coveredFaces.end())
			{
				patch.AddFace(f);
				coveredFaces.insert(f);
				patch.EstablishLoops();
			}
		}
	}	

	metaSingularitiesDirty = true;	

	PrintSingularityStatistics();
}

void Data::SaveFencedRegions() const
{
	auto f = fopen("fencedRegions.bin", "wb");
	nse::data::saveToFile(fencedRegions.patches, f);
	fclose(f);
}

void Data::RecordPatchEnteringEdges()
{
	fencedRegions.enteringEdgeToPatch.clear();
	for (int patchIdx = 0; patchIdx < fencedRegions.patches.size(); ++patchIdx)
	{
		auto& patch = fencedRegions.patches[patchIdx];
		if (!patch.IsActive())
			continue;
		//record the entering edges
		for (auto& l : patch.Loops())
		{
			for (int i = 0; i < l.Edges().size(); ++i)
			{
				auto prevLoopEdge = (i + l.Edges().size() - 1) % l.Edges().size();
				auto nextLoopEdge = (i + 1) % l.Edges().size();

				auto perpendicular1 = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(l.Edges()[i].edge));
				auto perpendicular2 = mesh.next_halfedge_handle(l.Edges()[i].edge);

				if (mesh.opposite_halfedge_handle(l.Edges()[prevLoopEdge].edge) != perpendicular1 || mesh.is_boundary(perpendicular1))
					fencedRegions.enteringEdgeToPatch[perpendicular1] = EnteringEdge(patchIdx, l.Edges()[i].orientation);
				if (l.Edges()[nextLoopEdge].orientation != l.Edges()[i].orientation && (l.Edges()[nextLoopEdge].edge != perpendicular2 || mesh.is_boundary(mesh.opposite_halfedge_handle(perpendicular2))))
					fencedRegions.enteringEdgeToPatch[perpendicular2] = EnteringEdge(patchIdx, l.Edges()[i].orientation);
			}
		}
	}
}

void Data::PrintSingularityStatistics() const
{
	size_t totalSingularitiesWithoutPatch = 0;
	for (auto& s : singularities)
	{
		if (s.correspondingRegion == (size_t)-1)
			++totalSingularitiesWithoutPatch;
	}

	size_t singularPatches = 0;
	size_t activeNonEmptyPatches = 0;
	size_t maxPatchSize = 0;
	for (auto& p : fencedRegions.patches)
	{
		if (!p.IsActive())
			continue;
		if (p.IncludedFaces().size() == 0)
			continue;
		++activeNonEmptyPatches;
		if (!p.IsRectilinear())
			++singularPatches;
		if (p.IncludedFaces().size() > maxPatchSize)
			maxPatchSize = p.IncludedFaces().size();
	}

	std::cout << singularities.size() << " singularities in total." << std::endl;
	std::cout << activeNonEmptyPatches << " fenced regions (including " << singularPatches << " singular regions)." << std::endl;
	std::cout << totalSingularitiesWithoutPatch << " singularities degenerate fenced region." << std::endl;
	std::cout << "Largest fenced region: " << maxPatchSize << " faces" << std::endl;
}

//returns true for real singularity and false for fake singularities
ClassificationResult Data::ClassifySingularity(const SingularityInfo& v, FencedRegion& patch, bool verbose)
{
	const int maxHoleLength = 40; //number of edges

	if (verbose)
		std::cout << "Classifying singularity " << v.vertexHandle.idx() << std::endl;

	int singularitySteps = 0;
	auto addFacesForSingularity = [&](const SingularityInfo& sing)
	{
		//add all faces incident to the given singularity
		if (sing.contextOutgoingBoundaryEdge.is_valid())
		{
			//singularity is non-manifold
			patch.FillMeshHole(sing.contextOutgoingBoundaryEdge, maxHoleLength);
			auto h = mesh.opposite_halfedge_handle(sing.contextOutgoingBoundaryEdge);
			CirculateForwardUntil<true>(h, mesh, [&](HEMesh::HalfedgeHandle h)
			{
				patch.AddFace(mesh.face_handle(h));
				return false;
			});
		}
		else
		{
			for (auto h : mesh.voh_range(sing.vertexHandle))
			{
				if (mesh.is_boundary(h))
					patch.FillMeshHole(h, maxHoleLength);
				else
					patch.AddFace(mesh.face_handle(h));
			}
		}
	};
	addFacesForSingularity(v);

	//grow the patch until there is no singularity on the boundary
	while (patch.HasSingularitiesOnBoundary())
	{
		if (patch.IncludedFaces().size() >= maxPatchSize)
		{
			if (verbose)
				std::cout << "Patch becomes too large." << std::endl;
			return NoPatch;
		}		

		auto& singularity = singularities[patch.NextSingularity()];
		addFacesForSingularity(singularity);
	}

	patch.EstablishLoops();

	if (patch.Loops().size() != 1)
	{
		if (verbose)
			std::cout << "Fenced region does not have a unique loop at " << mesh.point(v.vertexHandle) << "." << std::endl;
		return NoPatch;
	}

	auto& loop = patch.Loops()[0];
	if (loop.Degree() < 2)
	{
		if (verbose)
			std::cout << "Fenced region valence is smaller than 2." << std::endl;
		return NoPatch;
	}

	if (!loop.IsSingular())
		return RectilinearPatch;

	return NonRectilinearPatch;
}

bool Data::EvaluateMergeQuality(const std::set<size_t>& nodeSubset, std::vector<GraphNode>& nodes, const std::vector<GraphEdge>& edges, size_t& singularities, size_t& totalFacesInPatches)
{
	for (auto n : nodeSubset)
		nodes[n].handledForEvaluatingMergeQuality = false;

	singularities = 0;
	totalFacesInPatches = 0;

	//Find all connected components of the subgraph
	for (auto n : nodeSubset)
	{
		if (nodes[n].handledForEvaluatingMergeQuality)
			continue;

		//this is a new connected component

		int summedTurnCountDefect = 0;
		size_t nodesInConnectedComponent = 0;
		size_t sizeOfConnectedComponent = 0;
		std::vector<size_t> nodesOnFront;
		nodesOnFront.push_back(n);
		size_t countNegativeTurns = 0;
		size_t countPositiveTurns = 0;

		//perform BFS to explore entire connected component
		while (!nodesOnFront.empty())
		{
			auto currentNode = nodesOnFront.back();
			nodesOnFront.pop_back();
			if (nodes[currentNode].handledForEvaluatingMergeQuality)
				continue;
			nodes[currentNode].handledForEvaluatingMergeQuality = true;

			auto turnCountDefect = fencedRegions.patches[nodes[currentNode].patchIdx].Loops()[0].Degree() - 4;
			if (turnCountDefect > 0)
				++countPositiveTurns;
			else
				++countNegativeTurns;
			summedTurnCountDefect += turnCountDefect;
			++nodesInConnectedComponent;
			sizeOfConnectedComponent += fencedRegions.patches[nodes[currentNode].patchIdx].IncludedFaces().size();

			for (auto edge : nodes[currentNode].outgoingEdges)
			{
				if (edges[edge.second].enabled)
				{
					nodesOnFront.push_back(edge.first);
					if (currentNode < edge.first)
						sizeOfConnectedComponent += edges[edge.second].distance;
				}
			}
		}

		if (summedTurnCountDefect != 0 && nodesInConnectedComponent > 1)
			return false; //the resulting connected component is not regular and merges at least 2 fenced regions
		if (sizeOfConnectedComponent > maxPatchSize)
			return false; //the resulting connected component is too large
		if (countNegativeTurns != 1 && countPositiveTurns != 1)
			return false; //Both possible signs of the degree have more than one fenced region

		//if we leave the region unmerged ...
		if (summedTurnCountDefect != 0)
		{
			// ... record the number of singularities.
			auto& patch = fencedRegions.patches[nodes[n].patchIdx];
			if (patch.IsActive())
				++singularities;
			else
				singularities += patch.CoveredSingularities().size();
		}
		totalFacesInPatches += sizeOfConnectedComponent;
	}
	return true;
}

void Data::TryToMergePatches()
{
	FindMetaSingularities();

	nse::util::TimedBlock b("Trying to merge fenced regions ..");

	const size_t maxGrowRadius = 2;

	//Represents a face on the mesh that is traversed during BFS in order to build
	//the fenced region graph
	struct BFSNode
	{
		//the corresponding face
		HEMesh::FaceHandle face;

		//the distance from the originating fenced region in number of faces
		size_t distance;

		//index of the originating fenced region graph node
		size_t originatingSingularPatchNode;

		//index of the predecessor BFSNode
		size_t previousNode;

		BFSNode(HEMesh::FaceHandle face, size_t patch, size_t previousNode, size_t distance)
			: face(face), distance(distance), originatingSingularPatchNode(patch), previousNode(previousNode)
		{ }
	};

	//Edges of the fenced region graph
	std::vector<GraphEdge> singularPatchGraphEdges;
	//Nodes of the fenced region graph
	std::vector<GraphNode> patchNodes;

	//Perform multidirectional BFS in order to build the fenced region graph

	std::vector<BFSNode> bfsNodes;
	std::deque<size_t> bfsQueue; //indices in bfsNodes
	std::map<HEMesh::FaceHandle, size_t> faceToNodeIndex; //maps a mesh face to its BFSNode index

	int nextSingularPatchIdx = 0;
	//initialize BFS with singular loops
	for (auto iPatch = 0; iPatch < fencedRegions.patches.size(); ++iPatch)
	{
		bool isSingular = !fencedRegions.patches[iPatch].IsRectilinear();
		
		auto& p = fencedRegions.patches[iPatch];		
		if (isSingular)
			patchNodes.emplace_back(iPatch);
		for (auto f : p.IncludedFaces())
		{
			faceToNodeIndex[f] = (isSingular ? bfsNodes.size() : (size_t)-1);
			if (isSingular)
			{
				bfsQueue.push_back(bfsNodes.size());
				bfsNodes.emplace_back(f, nextSingularPatchIdx, -1, 0);
			}
		}
		if (isSingular)
			++nextSingularPatchIdx;
	}
	//make the neighbor faces of singularities unwalkable
	for (auto& sing : singularities)
	{
		if (sing.correspondingRegion != (size_t)-1)
			continue;
		auto hOut = (sing.contextOutgoingBoundaryEdge.is_valid() ? sing.contextOutgoingBoundaryEdge : mesh.halfedge_handle(sing.vertexHandle));
		auto hIn = mesh.opposite_halfedge_handle(hOut);
		CirculateForwardUntil<true>(hIn, mesh, [&](HEMesh::HalfedgeHandle h)
		{ 
			faceToNodeIndex[mesh.face_handle(h)] = (size_t)-1;
			return false;
		});
	}

	//Perform BFS to build fenced region graph
	while (!bfsQueue.empty())
	{
		//Take the next BFS node
		auto nodeIdx = bfsQueue.front();
		bfsQueue.pop_front();

		//only walk up to a given radius
		if (bfsNodes[nodeIdx].distance + 1 > maxGrowRadius)
			break;

		auto thisNodeTurncountDefect = fencedRegions.patches[patchNodes[bfsNodes[nodeIdx].originatingSingularPatchNode].patchIdx].Loops()[0].Degree() - 4;

		//Iterate all neighboring faces of the current face
		for (auto f : mesh.ff_range(bfsNodes[nodeIdx].face))
		{
			//create a new BFS node if the neighboring face has not yet been visited
			auto inserted = faceToNodeIndex.insert(std::make_pair(f, bfsNodes.size()));
			if (inserted.second)
			{
				//the neighboring node has not yet been visited
				bfsQueue.push_back(bfsNodes.size());
				size_t origin = bfsNodes[nodeIdx].originatingSingularPatchNode;
				size_t distance = bfsNodes[nodeIdx].distance;
				bfsNodes.emplace_back(f, origin, nodeIdx, distance + 1);
			}
			else
			{
				//the neighboring node has already been visited
				auto neighborNodeIdx = inserted.first->second;
				if (neighborNodeIdx != (size_t)-1) //if the face is walkable
				{
					auto& neighborNode = bfsNodes[neighborNodeIdx];
					auto neighborNodeTurncountDefect = fencedRegions.patches[patchNodes[neighborNode.originatingSingularPatchNode].patchIdx].Loops()[0].Degree() - 4;

					//If the originating region of the meeting nodes are different, this represents a new edge in the fenced region graph.
					//Consider only edges where the regions have different degree signs for performance reasons.
					if (neighborNode.originatingSingularPatchNode != bfsNodes[nodeIdx].originatingSingularPatchNode && ((thisNodeTurncountDefect > 0) ^ (neighborNodeTurncountDefect > 0)))
					{
						//Try to add the edge to the fenced region graph
						auto distanceInserted = patchNodes[bfsNodes[nodeIdx].originatingSingularPatchNode].outgoingEdges.insert(std::make_pair(neighborNode.originatingSingularPatchNode, singularPatchGraphEdges.size()));
						if (distanceInserted.second)
						{
							//this is a new fenced region graph edge
							patchNodes[neighborNode.originatingSingularPatchNode].outgoingEdges[bfsNodes[nodeIdx].originatingSingularPatchNode] = singularPatchGraphEdges.size();
							singularPatchGraphEdges.emplace_back();
							singularPatchGraphEdges.back().distance = bfsNodes[nodeIdx].distance + neighborNode.distance;
							singularPatchGraphEdges.back().from = neighborNode.originatingSingularPatchNode;
							singularPatchGraphEdges.back().fromBFSNode = neighborNodeIdx;
							singularPatchGraphEdges.back().to = bfsNodes[nodeIdx].originatingSingularPatchNode;
							singularPatchGraphEdges.back().toBFSNode = nodeIdx;
						}
					}
				}
			}
		}
	}

	//{
	//	nse::util::TimedBlock b("Generating graph visualization ..");
	//	//Generate patch graph visualization
	//	const std::string neatoPath = "\"E:\\Program Files (x86)\\Graphviz2.38\\bin\\neato\"";
	//	std::ofstream dot("layout.dot");
	//	dot << "graph G { node[style=filled];";
	//	for (int i = 0; i < loops.size(); ++i)
	//	{
	//		auto degree = patches.patches[loops[i].patchIdx].Loops()[loops[i].localLoopIdx].degree;
	//		auto color = (GetValenceColor(degree - 4) * 255.0f).cast<int>();
	//		dot << "n" << std::dec << i << " [label=\"" << i << "(" << degree << "; " << patches.patches[loops[i].patchIdx].IncludedFaces() << ")\";color=\"#" << std::hex << color.x() << color.y() << color.z() << "\"];";
	//	}
	//	dot << std::dec;
	//	for (int i = 0; i < loops.size(); ++i)
	//	{
	//		for (auto nl : loops[i].outgoingEdges)
	//			if (i < nl.first)
	//				dot << "n" << i << " -- n" << nl.first << " [label=\"" << loopGraphEdges[nl.second].distance << "\"];";
	//	}
	//	dot << "}";
	//	dot.close();
	//	std::string call = "\"" + neatoPath + " -Tpng layout.dot -O\"";
	//	system(call.c_str());
	//}

	//Merge patches to reduce number of singularities

	nse::util::UnionFind uf;
	uf.AddItems(fencedRegions.patches.size());
	auto oldRegionCount = fencedRegions.patches.size();

	//TODO: this can be done more efficiently because we know that there will be a turn count sign that will exist exactly once
	//Iterate all connected components
	for (int iLoop = 0; iLoop < patchNodes.size(); ++iLoop) //loop over all connected components
	{
		auto& l = patchNodes[iLoop];
		if (l.handled)
			continue;

		//find the connected component around l
		std::set<size_t> ccLoops;
		std::vector<size_t> ccEdges;

		std::deque<size_t> bfsQueue;
		bfsQueue.push_back(iLoop);
		ccLoops.insert(iLoop);
		while (!bfsQueue.empty())
		{
			auto loop = bfsQueue.front();
			bfsQueue.pop_front();
			for (auto edge : patchNodes[loop].outgoingEdges)
			{
				if (loop == singularPatchGraphEdges[edge.second].from)
					ccEdges.push_back(edge.second);
				if (ccLoops.find(edge.first) == ccLoops.end())
				{
					bfsQueue.push_back(edge.first);
					ccLoops.insert(edge.first);

					patchNodes[edge.first].handled = true;
				}
			}
		}

		//brute force all merge options	of this connected component

		uint32_t bestOption;
		size_t bestOptionSingularPatches = (size_t)-1;
		size_t bestOptionSize = (size_t)-1;
		uint32_t maxCheckOptions = 1 << 20;
		uint32_t checkOptions = 1 << ccEdges.size();
		if (ccEdges.size() > 32 || checkOptions > maxCheckOptions)
		{
			std::cout << "Too many edges in subgraph for merging. Simple brute forcing won't work. Restricting to a subset of options (" << (maxCheckOptions * 100.0 / std::ldexp(1.0, ccEdges.size())) << " %)." << std::endl;
			checkOptions = maxCheckOptions;
			//sort the edges by length
			std::sort(ccEdges.begin(), ccEdges.end(), [&](size_t e1, size_t e2) { return singularPatchGraphEdges[e1].distance < singularPatchGraphEdges[e2].distance; });
		}


		//The bits of option encode what edges to activate
		for (uint32_t option = 0; option < checkOptions; ++option)
		{
			for (int e = 0; e < ccEdges.size(); ++e)
				singularPatchGraphEdges[ccEdges[e]].enabled = (option & (1 << e)) != 0;

			size_t singularPatches, size;
			bool validOption = EvaluateMergeQuality(ccLoops, patchNodes, singularPatchGraphEdges, singularPatches, size);

			if (!validOption)
				continue;
			if (singularPatches > bestOptionSingularPatches)
				continue;
			if (singularPatches == bestOptionSingularPatches && size > bestOptionSize)
				continue;

			bestOption = option;
			bestOptionSingularPatches = singularPatches;
			bestOptionSize = size;
		}

		if (bestOption != 0)
		{
			//Merge all regions over active edges
			for (int e = 0; e < ccEdges.size(); ++e)
			{
				if (bestOption & (1 << e))
				{
					auto& edge = singularPatchGraphEdges[ccEdges[e]];
					auto currentNode = edge.fromBFSNode;

					auto rep1 = uf.GetRepresentative(patchNodes[edge.from].patchIdx);
					auto rep2 = uf.GetRepresentative(patchNodes[edge.to].patchIdx);

					if (rep1 == rep2)
						continue; //this should not happen

					fencedRegions.patches.emplace_back(mesh, vertexToSingularity, false);
					auto newRegionIdx = fencedRegions.patches.size() - 1;
					uf.AddItem();

					auto& mergedFrom = regionWasMergedFrom[newRegionIdx];
					mergedFrom[0] = rep1;
					mergedFrom[1] = rep2;

					uf.MergeWithPredefinedRoot(newRegionIdx, rep1);
					fencedRegions.patches.back().MergeWith(fencedRegions.patches[rep1]);
					uf.MergeWithPredefinedRoot(newRegionIdx, rep2);
					fencedRegions.patches.back().MergeWith(fencedRegions.patches[rep2]);

					fencedRegions.patches.back().SetActive(true);
					fencedRegions.patches[rep2].SetActive(false);
					fencedRegions.patches[rep1].SetActive(false);

					while (currentNode != (size_t)-1)
					{
						fencedRegions.patches.back().AddFace(bfsNodes[currentNode].face);
						currentNode = bfsNodes[currentNode].previousNode;
					}
					currentNode = edge.toBFSNode;
					while (currentNode != (size_t)-1)
					{
						fencedRegions.patches.back().AddFace(bfsNodes[currentNode].face);
						currentNode = bfsNodes[currentNode].previousNode;
					}
				}
			}
		}
	}

	//Update the covered singularities for merged regions
	for (auto patchIdx = oldRegionCount; patchIdx < fencedRegions.patches.size(); ++patchIdx)
	{
		auto& patch = fencedRegions.patches[patchIdx];
		if (patch.IncludedFaces().size() == 0 || !patch.IsActive())
			continue;
		for (auto iSing : patch.CoveredSingularities())
		{
			auto& sing = singularities[iSing];
			sing.correspondingRegion = patchIdx;
		}

		patch.EstablishLoops();
		if (patch.Loops().size() != 1)
		{
			std::cout << "Fenced region has " << patch.Loops().size() << " loops after merge. This should not happen" << std::endl;
			patch.SetActive(false); 
		}
		else if (patch.Loops()[0].Degree() != 4 || !patch.IsRectilinear())
		{
			std::cout << "Fenced region has turn count " << patch.Loops()[0].Degree() << " after merge at " << mesh.point(mesh.from_vertex_handle(patch.Loops()[0].Edges()[0].edge)) << ". This should not happen" << std::endl;
			patch.SetActive(false);
		}
	}
	
	metaSingularitiesDirty = true;	

	PrintSingularityStatistics();	
}

void Data::CalculateMotorcycleGraph(bool deactivateUnnecessaryMotorcycles)
{
	nse::util::TimedBlock b("Calculating motorcycle graph ..");

	while (true)
	{
		FindMetaSingularities();

		GenerateNewMotorcycleGraph();		

		{
			nse::util::TimedBlock b("Initializing motorcycles ..");
			for (auto& s : metaSingularities)
			{
				for (int d = 0; d < s.Degree(); ++d)
				{
					auto& motorcycle = s.GetEmanatingMotorcycle(d);
					motorcycleGraph->AddMotorcycleWithHistory(motorcycle, s.highPriority, !s.GetEmanatingMotorcycleActive(d));
				}
			}
		}

		motorcycleGraph->AdvanceToEnd();
		if (motorcycleGraph->Status() == MotorcycleGraph::Finished)
			break;

		if (motorcycleGraph->Status() == MotorcycleGraph::InternalFailure)
		{
			std::cout << "Internal failure while tracing motorcycles." << std::endl;
			break;
		}

		if (motorcycleGraph->Status() == MotorcycleGraph::CannotTraceThroughRegion)
		{
			std::cout << "Failed to trace through fenced region. Removing the region and retrying .." << std::endl;
			fencedRegions.patches[motorcycleGraph->ProblematicFencedRegion()].SetActive(false);
			std::deque<size_t> mergedFrom;
			//If the problematic fenced region is the result of a merge, reestablish the original regions first
			auto mergedFromIt = regionWasMergedFrom.find(motorcycleGraph->ProblematicFencedRegion());
			if (mergedFromIt != regionWasMergedFrom.end())
			{
				mergedFrom.push_back(mergedFromIt->second[0]);
				mergedFrom.push_back(mergedFromIt->second[1]);
			}
			while (!mergedFrom.empty())
			{
				auto patchIdx = mergedFrom.front();
				mergedFrom.pop_front();
				mergedFromIt = regionWasMergedFrom.find(patchIdx);
				if (mergedFromIt != regionWasMergedFrom.end())
				{
					mergedFrom.push_back(mergedFromIt->second[0]);
					mergedFrom.push_back(mergedFromIt->second[1]);
				}
				else
				{
					fencedRegions.patches[patchIdx].SetActive(true);
					for (auto iSing : fencedRegions.patches[patchIdx].CoveredSingularities())
					{
						auto& sing = singularities[iSing];
						sing.correspondingRegion = patchIdx;
					}
				}
			}

			metaSingularitiesDirty = true;
		}
	}

	if (motorcycleGraph->Status() == MotorcycleGraph::Finished && deactivateUnnecessaryMotorcycles)
		motorcycleGraph->DeactivateUnnecessaryMotorcycles();
}

MotorcycleGraph::ExtractionStatistics Data::ExtractPatches(bool splitNonRectangularPatches)
{
	MotorcycleGraph::ExtractionStatistics stats;
	motorcycleGraph->ExtractPatches(stats);	

	if (splitNonRectangularPatches)
		motorcycleGraph->SplitNonRectangularPatches();

	return stats;
}

void Data::CalculateParametrization(float targetLengthPerEdge, float parametrizationErrorThreshold, const std::chrono::steady_clock::duration& discreteOptimizationTimeLimit, ArclengthStrategy arclengthStrategy)
{	
	float parametrizationFactor = targetLengthPerEdge / averageEdgeLength;
	TargetLengthStrategyPatchAverage targetLengthStrategy(lengthMeasure, mesh, *motorcycleGraph, parametrizationFactor);
	//TargetLengthStrategySimple targetLengthStrategy;
	std::unique_ptr<IArclengthStrategy> arclengthStrategyObj;

	switch (arclengthStrategy)
	{
	case Simple:
		arclengthStrategyObj = std::make_unique<ArclengthStrategySimple>();
		break;
	case GurobiSeparate:
		arclengthStrategyObj = std::make_unique<ArclengthStrategyGurobiSeparate>();
		break;
	case GurobiGlobal:
		arclengthStrategyObj = std::make_unique<ArclengthStrategyGurobi>(discreteOptimizationTimeLimit);
		break;
	}
	ParametrizationStrategyLSCMForInteriorScaffoldForBoundary parametrizationStrategy(parametrizationFactor);

	GenerateTexCoordAccess();
	Parametrization p(mesh, *motorcycleGraph, motorcycleGraph->Patches(), texCoords, parametrizationErrorThreshold, &targetLengthStrategy, arclengthStrategyObj.get(), &parametrizationStrategy);

	p.run();
	brokenHalfarcs = p.BrokenHalfarcs();	
}

float TriangleMIPSEnergy(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2, const Eigen::Vector3f& p3, 
	const Eigen::Vector2f& uv1, const Eigen::Vector2f& uv2, const Eigen::Vector2f& uv3)
{
	Eigen::Vector3f d1 = p2 - p1; //=s-axis of 2D coordinate system
	Eigen::Vector3f d2 = p3 - p1;

	Eigen::Vector3f normal = d1.cross(d2).normalized();

	Eigen::Matrix2f jacobian;
	jacobian.col(0) = (uv2 - uv1);

	//t-axis of 2D coordinate system
	Eigen::Vector3f tAxis = normal.cross(d1);

	float d1LengthSq = d1.squaredNorm();
	float s3 = d1.dot(d2) / d1LengthSq;
	float t3 = tAxis.dot(d2) / d1LengthSq;

	jacobian.col(1) = (1.0f / t3) * (uv3 - uv1 - s3 * jacobian.col(0));

	float mips = (jacobian.transpose() * jacobian).trace() / (jacobian.determinant());
	return mips;
}

void Data::EvaluateParametrization(Statistics& logAreaRatios, Statistics& mips, bool printStatistics)
{
	size_t startIndex = 0;// logAreaRatios.size();
	logAreaRatios.Clear();
	mips.Clear();
		
	Eigen::Vector3f p[4];
	Eigen::Vector2f uv[4];
	for (auto i = 0; i < motorcycleGraph->Patches().size(); ++i)
	{
		auto& patch = motorcycleGraph->Patches()[i];
		auto& texCoordStorage = AccessTexCoords(i);
		for (auto f : patch.Faces())
		{
			//Calculate area and MIPS energy for each face
			OpenMesh::Vec3f area3D = OpenMesh::Vec3f(0, 0, 0);
			float area2D = 0;
			int i = 0;
			for (auto h : mesh.fh_range(f))
			{
				auto& p1 = mesh.point(mesh.from_vertex_handle(h));				
				auto& p2 = mesh.point(mesh.to_vertex_handle(h));
				area3D += OpenMesh::cross(p1, p2);
				p[i] = ToEigenVector(p1);

				auto& t1 = texCoordStorage.TexCoordAtFromVertex(h, mesh);
				auto& t2 = texCoordStorage.TexCoordAtToVertex(h, mesh);
				area2D += t1.x() * t2.y() - t1.y() * t2.x();
				uv[i] = t1;
				++i;
			}
			if (area2D < 1.0e-10)
				continue;
			if (area3D.norm() < 1.0e-10)
				continue;
			float logRatio = std::log(area2D / area3D.norm());
			logAreaRatios.AddDatum(logRatio);						

			auto mipsEnergy = TriangleMIPSEnergy(p[0], p[1], p[2], uv[0], uv[1], uv[2]);
			mips.AddDatum(mipsEnergy);
			mipsEnergy = TriangleMIPSEnergy(p[0], p[2], p[3], uv[0], uv[2], uv[3]);
			mips.AddDatum(mipsEnergy);
		}
	}
	size_t newRatiosCount = logAreaRatios.Size() - startIndex;
	if (newRatiosCount > 0)
	{		
		float medLogRatio = logAreaRatios.Median();
		
		logAreaRatios.SetGlobalShift(-medLogRatio);		

		if (printStatistics)
		{
			std::cout << "Parametrization statistics:" << std::endl;
			std::cout << "Minimum log area ratio: " << logAreaRatios.Min() << " (x" << std::exp(logAreaRatios.Min()) << ")" << std::endl;
			std::cout << "Maximum log area ratio: " << logAreaRatios.Max() << " (x" << std::exp(logAreaRatios.Max()) << ")" << std::endl;
			std::cout << "Average log area ratio: " << logAreaRatios.Average() << " (x" << std::exp(logAreaRatios.Average()) << ")" << std::endl;
			for (int i = 1; i < 20; ++i)
			{
				float percent = i * 0.05;
				auto percentile = logAreaRatios.Percentile(percent);
				std::cout << (100 * percent) << "-tile: " << percentile << " (x" << std::exp(percentile) << ")" << std::endl;
			}
			std::cout << "Standard deviation of log ratio: " << std::sqrt(logAreaRatios.Variance()) << " (x" << std::exp(std::sqrt(logAreaRatios.Variance())) << ")" << std::endl;

			std::cout << "10-tile MIPS: " << mips.Percentile(0.1f) << std::endl;
			std::cout << "Median MIPS: " << mips.Median() << std::endl;
			std::cout << "90-tile MIPS: " << mips.Percentile(0.9f) << std::endl;
		}
	}

	std::cout << "Total edges: " << mesh.n_edges() << std::endl;
	std::cout << "Broken arcs: " << brokenHalfarcs.size() << " out of " << motorcycleGraph->Halfarcs().size() / 2 << std::endl;
	int totalArcLength = 0, totalBrokenArclength = 0;
	for (auto& arc : motorcycleGraph->Halfarcs())
		for (auto p : arc)
			++totalArcLength;
	for(auto idx : brokenHalfarcs)
		for (auto p : motorcycleGraph->Halfarcs()[idx])
			++totalBrokenArclength;
	std::cout << "Broken arc length: " << totalBrokenArclength << " out of " << totalArcLength / 2 << std::endl;
}

void Data::PackTexture(int mipFactor)
{
	nse::util::TimedBlock b("Packing texture patches ..");

	//We place corners of texture patches at texel centers. If an application
	//needs them at texel corners, this packing procedure has to be adapted accordingly.

	std::vector<rbp::RectSize> patchSizes(motorcycleGraph->Patches().size()); //size on the targeted MIP level
	std::vector<std::pair<Eigen::Vector2f, Eigen::Vector2f>> minMax(motorcycleGraph->Patches().size());
	double area = 0;
	for (int i = 0; i < motorcycleGraph->Patches().size(); ++i)
	{
		auto& patch = motorcycleGraph->Patches()[i];
		auto& texCoordStorage = AccessTexCoords(i);

		texCoordStorage.GetMinMaxTextureCoordinates(minMax[i].first, minMax[i].second);
		auto& size = minMax[i].second - minMax[i].first;
		patchSizes[i].width = (int)std::ceil(size.x()) + 1; //two half pixels margin
		patchSizes[i].height = (int)std::ceil(size.y()) + 1; //two half pixels margin
		area += (double)patchSizes[i].width * patchSizes[i].height;
	}
	
	//Find the closest power of 2 texture width that would approximately fit the required number of texels
	double textureArea = area;
	textureWidth = (int)std::ldexp(1.0, std::round(std::log2(sqrt(textureArea))));
	
	//Use a texture width of at least 32 x 32
	if (textureWidth < 32)
		textureWidth = 32;
	textureHeight = std::ceil(textureArea / textureWidth);
	if (textureHeight < 32)
		textureHeight = 32;

	std::vector<rbp::Rect> placedRectangles;
	bool success;
	do
	{
		//Try packing
		std::cout << "Trying texture of size " << textureWidth << " x " << textureHeight << std::endl;
		rbp::MaxRectsBinPack bin(textureWidth, textureHeight);
		success = bin.Insert(patchSizes, placedRectangles, rbp::MaxRectsBinPack::RectBestShortSideFit);
		if (!success)
			//Increase attempted texture height by 1%
			textureHeight = (int)std::ceil(1.01 * textureHeight);
		else
			std::cout << "Occupancy: " << bin.Occupancy() << std::endl;

	} while (!success);

	//Modify the texture coordinates according to the packed rectangles
#pragma omp parallel for
	for (int i = 0; i < motorcycleGraph->Patches().size(); ++i)
	{
		bool rotate = patchSizes[i].width != placedRectangles[i].width;
		Eigen::Affine2f rotation = Eigen::Translation2f(placedRectangles[i].width, 0) * Eigen::Rotation2Df((float)M_PI / 2);
		auto& patch = motorcycleGraph->Patches()[i];
		auto& texCoordStorage = AccessTexCoords(i);

		for (int j = 0; j < texCoordStorage.TexCoordCount(); ++j)
		{
			auto& uv = texCoordStorage.TexCoord(j);
			uv -= minMax[i].first;
			uv += Eigen::Vector2f(0.5f, 0.5f);
			if (rotate)
				uv = rotation * uv;
			uv += Eigen::Vector2f(placedRectangles[i].x, placedRectangles[i].y);
			uv.x() /= textureWidth;
			uv.y() /= textureHeight;
		}
	}

	textureWidth *= mipFactor;
	textureHeight *= mipFactor;
}

void Data::GenerateNewMotorcycleGraph()
{
	motorcycleGraph = std::make_unique<::MotorcycleGraph>(mesh, fencedRegions, metaSingularities, vertexToMetaSingularity, isOriginalEdgeProp, lengthMeasure);
	texCoords.clear();
	texCoordsDirty = false;
	motorcycleGraph->PatchesChanged.Subscribe([this]()
	{
		texCoordsDirty = true;
	});
}

void Data::GenerateTexCoordAccess()
{
	if (texCoordsDirty)
	{
		texCoords.clear();
		for (auto& p : motorcycleGraph->Patches())
			texCoords.emplace_back(p);
		texCoordsDirty = false;
	}
}

Eigen::Vector4f Data::GetPatchColor(size_t patchIndex) const
{
	return GetSegmentColor(patchColors[patchIndex]);
}

TextureCoordinatesStorage& Data::AccessTexCoords(size_t patchIdx)
{
	GenerateTexCoordAccess();

	return texCoords[patchIdx];
}

const TextureCoordinatesStorage & Data::AccessTexCoords(size_t patchIdx) const
{
	if (texCoordsDirty)
		throw std::runtime_error("The texture coordinate storages are dirty. Call GenerateTexCoordAccess().");

	return texCoords[patchIdx];
}
