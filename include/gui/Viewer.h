#pragma once

#include <nsessentials/gui/AbstractViewer.h>

#include <nsessentials/gui/GLBuffer.h>
#include <nsessentials/gui/GLVertexArray.h>
#include <nsessentials/math/BoundingBox.h>

#include <chrono>

#include "common.h"
#include "Data.h"

//Viewer for the sample application
class Viewer : public nse::gui::AbstractViewer
{
public:
	Viewer();

	void drawContents();

	bool resizeEvent(const Eigen::Vector2i&);	

protected:
	bool mouseMotionHook(const Eigen::Vector2i & p, const Eigen::Vector2i & rel, int button, int modifiers);
	bool mouseButtonHook(const Eigen::Vector2i & p, int button, bool down, int modifiers);

private:	

	void SetupGUI();
	
	//Uploads vertex buffers, index buffers, and arrays to GPU.
	void UploadWireframeToGPU();
	void UploadMeshFaces();
	void UploadSingularities();	
	template <typename Iterator>
	void UploadFencedRegions(Iterator itBegin, Iterator itEnd, bool includeDeactivated = false);
	void UploadMotorcycles(const MotorcycleGraph& graph);
	void UploadPatches();
	void UploadBrokenArcs(const MotorcycleGraph& graph, const std::vector<size_t>& brokenArcs);

	//Takes a screenshot of the current view and saves it to a PNG file
	void TakeScreenshot(const std::string& filename, int screenshotWidth, int screenshotHeight);	

	//Reads the index from the picking framebuffer at the provided position.
	int32_t GetPickingIndex(const Eigen::Vector2i p) const;	
	
	void render(const Eigen::Matrix4f& mv, const Eigen::Matrix4f& proj, bool background = true);
	//Draws 2D text at a given 3D point if the point is not occluded by other parts of the scene
	void DrawTextIfNotHidden(const Eigen::Vector4f& pos, const Eigen::Matrix4f mvp, const std::string& text, float tolerance);
	
	nanogui::CheckBox* chkMergeTriangles;
	nanogui::Slider* angleThresholdSld;
	nanogui::CheckBox* chkFocusOnModel;
	nanogui::CheckBox* chkScaleTex;

	nanogui::CheckBox* chkShowSingularities;
	nanogui::Slider* sldSingularitySize;
	nanogui::CheckBox* chkShowMotorcycleGraph;
	nanogui::CheckBox* chkShowMotorcycleIndices;
	nanogui::CheckBox* chkShowWireframe;
	nanogui::Slider* sldWireframeWidth;
	nanogui::CheckBox* chkShowFaces;
	nanogui::CheckBox* chkShowFencedRegions;
	nanogui::CheckBox* chkShowFencedRegionIndices;
	nanogui::CheckBox* chkShowPatches;
	nanogui::CheckBox* chkShowUV;
	nanogui::ComboBox* cmbPatchColoring;
	nanogui::Slider* sldTubeWidth;
	nanogui::CheckBox* chkShowVertexIndices;
	nanogui::CheckBox* chkShowBrokenArcs;
	nanogui::Slider* sldParameterizationGridSize;

	nanogui::CheckBox* chkClassifySingularities;
	nanogui::CheckBox* chkMergeFencedRegion;
	nanogui::CheckBox* chkCalculateMotorcycleGraph;
	nanogui::CheckBox* chkDeactivateUnnecessaryMotorcycles;
	nanogui::CheckBox* chkSplitNonRectangularPatches;
	nanogui::CheckBox* chkExtractPatches;

	nanogui::Slider* sldParametricEdgeLength;
	nanogui::CheckBox* chkInvisibleSeams;
	nanogui::CheckBox* chkPackTexture;
		
	bool isFileDialogOpen = false;	

	GLuint fbo, fboSinglePixel, colorTexture, pickingTexture, pickingTextureSinglePixel, depthBuffer;

	nse::gui::GLBuffer meshVertices, meshIndices;
	nse::gui::GLVertexArray meshVAO;
	unsigned int indexCount;	

	nse::gui::GLBuffer wireframeIndices;
	nse::gui::GLVertexArray wireframeVAO;
	unsigned int wireframeIndexCount;

	nse::gui::GLBuffer singularityPositions, singularityColors;
	nse::gui::GLVertexArray singularityVAO;
	unsigned int singularityCount;

	nse::gui::GLVertexArray emptyVAO;

	nse::gui::GLBuffer fencedRegionIndices, fencedRegionFaceIndices;
	nse::gui::GLVertexArray fencedRegionsVAO, fencedRegionsFaceVAO;
	unsigned int fencedRegionIndexCount[3], fencedRegionFaceIndexCount[3]; //valence less than 4, 4, greater than 4

	nse::gui::GLBuffer motorcycleIndices;
	nse::gui::GLVertexArray motorcycleVAO;
	std::vector<unsigned int> motorcycleIndexCount;

	nse::gui::GLBuffer brokenArcsIndices;
	nse::gui::GLVertexArray brokenArcsVAO;
	unsigned int brokenArcsIndexCount;

	nse::gui::GLBuffer patchPositionsVBO, patchColorsVBO, patchTexCoordsVBO;
	nse::gui::GLVertexArray patchVAO;
	unsigned int patchIndexCount;	
	
	size_t selectedRegion;

	Data data;

	void LoadMesh(const std::string& filename);
	void LoadMeshForVisualization(const std::string& filename);
	void ExtractQuadLayout();
	void ProcessDirectory();

	std::map<size_t, size_t> primitiveIdToSingularityId;

	float parametrizationErrorThreshold;
	std::chrono::steady_clock::duration discreteOptimizationTimeLimit;

	int mipLevel;
};