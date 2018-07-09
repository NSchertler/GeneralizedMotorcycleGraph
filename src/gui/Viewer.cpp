#include "gui/Viewer.h"
#include "gui/ShaderPool.h"

#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/button.h>
#include <nanogui/slider.h>
#include <nanogui/label.h>
#include <nanogui/textbox.h>
#include <nanogui/checkbox.h>
#include <nanogui/combobox.h>
#include <nanogui/messagedialog.h>

#include <nsessentials/data/FileHelper.h>

#include <nsessentials/util/TimedBlock.h>
#include <nsessentials/util/UnionFind.h>
#include <iostream>
#include <set>
#include <deque>
#include <queue>
#include <utility>
#include <fstream>

#include "stb_image_write.h"

#define USE_MULTISAMPLING

//#define LAPTOP_COMPAT

Viewer::Viewer()
	: AbstractViewer("Generalized Motorcycle Graph", 1280, 800, 
#ifdef USE_MULTISAMPLING
		4
#else
		1
#endif
	),
	meshVertices(nse::gui::VertexBuffer), meshIndices(nse::gui::IndexBuffer),
	singularityPositions(nse::gui::VertexBuffer), singularityColors(nse::gui::VertexBuffer), singularityCount(0),
	wireframeIndices(nse::gui::IndexBuffer),
	fencedRegionIndices(nse::gui::IndexBuffer),
	fencedRegionFaceIndices(nse::gui::IndexBuffer),
	motorcycleIndices(nse::gui::IndexBuffer), motorcycleIndexCount(0),
	patchPositionsVBO(nse::gui::VertexBuffer), patchColorsVBO(nse::gui::VertexBuffer), patchTexCoordsVBO(nse::gui::VertexBuffer), patchIndexCount(0),
	selectedRegion(-1),
	brokenArcsIndices(nse::gui::IndexBuffer)
{
	ShaderPool::Instance()->CompileAll();

	SetupGUI();
	
	glGenFramebuffers(1, &fbo);
	glGenFramebuffers(1, &fboSinglePixel);
	glGenTextures(1, &colorTexture);
	glGenTextures(1, &pickingTexture);
	glGenTextures(1, &pickingTextureSinglePixel);
	glGenRenderbuffers(1, &depthBuffer);	

	emptyVAO.generate();

	for (int i = 0; i < 3; ++i)
	{
		fencedRegionIndexCount[i] = 0;
		fencedRegionFaceIndexCount[i] = 0;
	}
}

bool Viewer::resizeEvent(const Eigen::Vector2i& size)
{
	AbstractViewer::resizeEvent(size);

	GLenum status;

#ifndef LAPTOP_COMPAT
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	auto textureTarget = (nSamples == 1 ? GL_TEXTURE_2D : GL_TEXTURE_2D_MULTISAMPLE);

	glBindTexture(textureTarget, colorTexture);
	if (nSamples == 1)
		glTexImage2D(textureTarget, 0, GL_RGBA8, width(), height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	else
		glTexImage2DMultisample(textureTarget, nSamples, GL_RGBA8, width(), height(), true);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, textureTarget, colorTexture, 0);

	glBindTexture(textureTarget, pickingTexture);
	if (nSamples == 1)
		glTexImage2D(textureTarget, 0, GL_R32I, width(), height(), 0, GL_RED_INTEGER, GL_INT, nullptr);
	else
		glTexImage2DMultisample(textureTarget, nSamples, GL_R32I, width(), height(), true);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, textureTarget, pickingTexture, 0);
	glBindTexture(textureTarget, 0);

	glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
	if (nSamples == 1)
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width(), height());
	else
		glRenderbufferStorageMultisample(GL_RENDERBUFFER, nSamples, GL_DEPTH24_STENCIL8, width(), height());
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE)
	{
		std::cout << "Framebuffer is incomplete: " << status << std::endl;
	}
#else
	fbo = 0;
#endif

	glBindFramebuffer(GL_FRAMEBUFFER, fboSinglePixel);

	glBindTexture(GL_TEXTURE_2D, pickingTextureSinglePixel);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32I, 1, 1, 0, GL_RED_INTEGER, GL_INT, nullptr);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pickingTextureSinglePixel, 0);
	
	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE)
	{
		std::cout << "Single-pixel framebuffer is incomplete: " << status << std::endl;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	return true;
}

nanogui::Slider* AddWidgetWidthSlider(nanogui::Widget* parent, const std::string& caption, const std::pair<float, float>& range, float defaultValue, nanogui::Widget*& outWidget)
{
	outWidget = new nanogui::Widget(parent);
	outWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 6));
	new nanogui::Label(outWidget, caption);

	auto slider = new nanogui::Slider(outWidget);
	slider->setFixedWidth(100);
	slider->setValue(defaultValue);
	slider->setRange(range);

	return slider;
}

nanogui::Slider* AddLabeledSlider(nanogui::Widget* parent, const std::string& caption, const std::pair<float, float>& range, float defaultValue)
{
	nanogui::Widget* widget;
	return AddWidgetWidthSlider(parent, caption, range, defaultValue, widget);
}

nanogui::Slider* AddLabeledSlider(nanogui::Widget* parent, const std::string& caption, const std::pair<float, float>& range, float defaultValue, nanogui::TextBox*& out_label)
{
	nanogui::Widget* widget;
	auto slider = AddWidgetWidthSlider(parent, caption, range, defaultValue, widget);	

	out_label = new nanogui::TextBox(widget);
	out_label->setFixedSize(Eigen::Vector2i(50, 25));

	return slider;
}

void Viewer::SetupGUI()
{
	auto ctx = nvgContext();

	auto mainWindow = new nanogui::Window(this, "Regular Mesh Texturing");
	mainWindow->setPosition(Eigen::Vector2i(15, 15));
	mainWindow->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 4, 4));

	auto loadingOptionsBtn = new nanogui::PopupButton(mainWindow, "Loading Options");
	loadingOptionsBtn->popup()->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 10, 4));

	nanogui::TextBox* angleThresholdLabel;
	angleThresholdSld = AddLabeledSlider(loadingOptionsBtn->popup(), "Angle Threshold:", std::make_pair(1, 45), 20, angleThresholdLabel);
	angleThresholdSld->setCallback([angleThresholdLabel](float value)
	{
		std::stringstream ss;
		ss.precision(2);
		ss << value << "\302\260"; //\302\260 = °
		angleThresholdLabel->setValue(ss.str());
	});
	angleThresholdSld->callback()(angleThresholdSld->value());
	chkMergeTriangles = new nanogui::CheckBox(loadingOptionsBtn->popup(), "Merge Triangulated Quads");

	chkFocusOnModel = new nanogui::CheckBox(loadingOptionsBtn->popup(), "Focus On Model");	chkFocusOnModel->setChecked(true);
	chkScaleTex = new nanogui::CheckBox(loadingOptionsBtn->popup(), "Scale Texture Coordinates");	chkScaleTex->setChecked(true);

	auto loadVisBtn = new nanogui::Button(loadingOptionsBtn->popup(), "Load Mesh for Vis");
	loadVisBtn->setCallback([&]()
	{
		isFileDialogOpen = true;
		std::string filename = nanogui::file_dialog({ { "obj", "OBJ" } }, false);
		if (filename != "")
			LoadMeshForVisualization(filename);
		isFileDialogOpen = false;
	});



	auto loadBtn = new nanogui::Button(mainWindow, "Load Mesh");
	loadBtn->setCallback([&]()
	{
		isFileDialogOpen = true;
		std::string filename = nanogui::file_dialog({ {"ply", "PLY"}, {"obj", "OBJ"} }, false);
		if (filename != "")
			LoadMesh(filename);
		isFileDialogOpen = false;
	});

	//auto processDirectoryBtn = new nanogui::Button(mainWindow, "Process Directory");
	//processDirectoryBtn->setCallback([this]() { ProcessDirectory(); });

	auto saveBtn = new nanogui::Button(mainWindow, "Save Untextured Mesh");
	saveBtn->setCallback([&]()
	{
		isFileDialogOpen = true;
		std::string filename = nanogui::file_dialog({ { "ply", "PLY" } }, true);		
		if (!filename.empty())
		{
			if (filename.substr(filename.length() - 4) != ".ply")
				filename.append(".ply");
			if (filename != "")
				data.SaveMeshPLY(filename);
		}
		isFileDialogOpen = false;
	});

	auto saveBtnObj = new nanogui::Button(mainWindow, "Save Textured Mesh");
	saveBtnObj->setCallback([&]()
	{
		isFileDialogOpen = true;
		std::string filename = nanogui::file_dialog({ { "obj", "OBJ" } }, true);
		if (!filename.empty())
		{
			if (filename.substr(filename.length() - 4) != ".obj")
				filename.append(".obj");

			if (filename != "")
			{
				data.SaveMeshOBJ(filename);
				glViewport(0, 0, width(), height());
			}
		}
		isFileDialogOpen = false;
	});

	auto displayOptionsBtn = new nanogui::PopupButton(mainWindow, "Display Options");
	displayOptionsBtn->popup()->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 10, 4));

	chkShowFaces = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Faces");	chkShowFaces->setChecked(true);
	chkShowWireframe = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Wireframe");	chkShowWireframe->setChecked(true);
	sldWireframeWidth = AddLabeledSlider(displayOptionsBtn->popup(), "Wireframe Thickness", std::make_pair(0.1f, 15.0f), 1.0f);
	sldWireframeWidth->setFixedWidth(250);
	chkShowSingularities = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Singularities");	chkShowSingularities->setChecked(true);
	sldSingularitySize = AddLabeledSlider(displayOptionsBtn->popup(), "Singularity Size", std::make_pair(0.1f, 10.0f), 1.0f);
	sldSingularitySize->setCallback([this](float) { UploadSingularities(); });
	sldSingularitySize->setFixedWidth(250);
	chkShowMotorcycleGraph = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Motorcycle Graph");	chkShowMotorcycleGraph->setChecked(true);
	chkShowMotorcycleIndices = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Motorcycle Indices");
	chkShowFencedRegions = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Fenced Regions");	chkShowFencedRegions->setChecked(true);
	chkShowFencedRegionIndices = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Fenced Region Indices");
	chkShowPatches = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Patches");	chkShowPatches->setChecked(true);
	chkShowUV = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show UV");	chkShowUV->setChecked(true);
	chkShowBrokenArcs = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Broken Arc Constraints");
	chkShowVertexIndices = new nanogui::CheckBox(displayOptionsBtn->popup(), "Show Vertex Indices");
	cmbPatchColoring = new nanogui::ComboBox(displayOptionsBtn->popup(), { "Color by Patch Index", "Color by Patch Degree" });
	cmbPatchColoring->setCallback([this](bool) { this->UploadPatches(); });
	sldTubeWidth = AddLabeledSlider(displayOptionsBtn->popup(), "Motorcycle Graph Thickness", std::make_pair(0.1f, 20.0f), 1.0f);
	sldTubeWidth->setFixedWidth(250);
	sldParameterizationGridSize = AddLabeledSlider(displayOptionsBtn->popup(), "Param Grid Size", std::make_pair(-1.0f, 8.0f), 2.0f);

	(new nanogui::Widget(mainWindow))->setHeight(10); //spacing

	auto quadLayoutOptionsBtn = new nanogui::PopupButton(mainWindow, "Quad Layout Options");
	quadLayoutOptionsBtn->popup()->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 10, 4));

	chkClassifySingularities = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Classify All Singularities");
	chkClassifySingularities->setChecked(true);

	chkMergeFencedRegion = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Merge Fenced Regions");
	chkMergeFencedRegion->setChecked(true);

	chkCalculateMotorcycleGraph = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Calculate Motorcycle Graph");
	chkCalculateMotorcycleGraph->setChecked(true);

	chkDeactivateUnnecessaryMotorcycles = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Deactivate Unnecessary Motorcycles");
	chkDeactivateUnnecessaryMotorcycles->setChecked(true);

	chkExtractPatches = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Extract Patches");
	chkExtractPatches->setChecked(true);

	chkSplitNonRectangularPatches = new nanogui::CheckBox(quadLayoutOptionsBtn->popup(), "Split Non-Rectangular Patches");
	chkSplitNonRectangularPatches->setChecked(true);

	auto quadLayoutBtn = new nanogui::Button(mainWindow, "Calculate Quad Layout");
	quadLayoutBtn->setCallback([this]() { ExtractQuadLayout(); });

	(new nanogui::Widget(mainWindow))->setHeight(10); //spacing

	auto parametrizationOptionsBtn = new nanogui::PopupButton(mainWindow, "Parametrization Options");
	parametrizationOptionsBtn->popup()->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 10, 4));

	nanogui::TextBox* txtParametricEdgeLength;
	sldParametricEdgeLength = AddLabeledSlider(parametrizationOptionsBtn->popup(), "Parametric Edge Length:", std::make_pair(-3, 3), std::log(10.0f), txtParametricEdgeLength);
	sldParametricEdgeLength->setCallback([txtParametricEdgeLength](float value)
	{
		auto length = std::exp(value);
		std::stringstream ss;
		ss.precision(2);
		ss << length;
		txtParametricEdgeLength->setValue(ss.str());
	});
	sldParametricEdgeLength->callback()(sldParametricEdgeLength->value());

	nanogui::TextBox* txtParametrizationErrorThreshold;
	auto sldParametrizationErrorThreshold = AddLabeledSlider(parametrizationOptionsBtn->popup(), "Parametrization Error Threshold:", std::make_pair(1, 100), 50, txtParametrizationErrorThreshold);
	sldParametrizationErrorThreshold->setCallback([this, txtParametrizationErrorThreshold](float value)
	{
		parametrizationErrorThreshold = value;
		if (value == 100)
			parametrizationErrorThreshold = std::numeric_limits<float>::infinity();
		
		std::stringstream ss;
		ss.precision(2);
		ss << parametrizationErrorThreshold;
		txtParametrizationErrorThreshold->setValue(ss.str());
	});
	sldParametrizationErrorThreshold->callback()(sldParametrizationErrorThreshold->value());

	nanogui::TextBox* txtOptimizationTimeLimit;
	auto sldOptimizationTimeLimit = AddLabeledSlider(parametrizationOptionsBtn->popup(), "Optimization Time Limit:", std::make_pair(1.0 / 60.0, 120.0), 5.0, txtOptimizationTimeLimit);
	sldOptimizationTimeLimit->setCallback([this, txtOptimizationTimeLimit](float value)
	{
		std::stringstream ss;
		ss << std::fixed;
		ss.precision(1);

		if (value == 120)
		{
			discreteOptimizationTimeLimit = std::chrono::steady_clock::duration::max();
			ss << "unlimited";
		}
		else
		{
			discreteOptimizationTimeLimit = std::chrono::seconds((int)std::round(value * 60));
			ss << value << " min";			
		}
		
		txtOptimizationTimeLimit->setValue(ss.str());
	});
	sldOptimizationTimeLimit->callback()(sldOptimizationTimeLimit->value());
	txtOptimizationTimeLimit->setFixedWidth(90);
	sldOptimizationTimeLimit->setFixedWidth(250);

	nanogui::TextBox* txtMipLevel;
	auto sldMipLevel = AddLabeledSlider(parametrizationOptionsBtn->popup(), "Texture Layout for MIP Level up to:", std::make_pair(0, 8), 0, txtMipLevel);
	sldMipLevel->setCallback([this, txtMipLevel](float value)
	{
		mipLevel = (int)std::round(value);		

		std::stringstream ss;
		ss << mipLevel;
		txtMipLevel->setValue(ss.str());
	});
	sldMipLevel->callback()(sldMipLevel->value());	

	chkInvisibleSeams = new nanogui::CheckBox(parametrizationOptionsBtn->popup(), "Invisible Seams");	
#ifndef WITH_GUROBI
	chkInvisibleSeams->setEnabled(false);
	chkInvisibleSeams->setCaption("Invisible Seams are not supported because Gurobi is not available.");
#endif // !WITH_GUROBI


	chkPackTexture = new nanogui::CheckBox(parametrizationOptionsBtn->popup(), "Pack Texture");

	auto parameterizationBtn = new nanogui::Button(mainWindow, "Calculate Parametrization");
	parameterizationBtn->setCallback([this]() 
	{
		if (data.Motorcycles() == nullptr || data.Motorcycles()->Patches().size() == 0)
		{
			new nanogui::MessageDialog(this, nanogui::MessageDialog::Type::Warning, "Parametrization",
				"A quad layout is needed to calculate the parametrization.");
			return;
		}

		int mipFactor = 1 << mipLevel;
		ArclengthStrategy arclengthStrategy = chkInvisibleSeams->checked() ? ArclengthStrategy::GurobiGlobal : ArclengthStrategy::Simple;
		data.CalculateParametrization(std::exp(sldParametricEdgeLength->value()) / mipFactor, parametrizationErrorThreshold, discreteOptimizationTimeLimit, arclengthStrategy);
		Statistics ratios, mips;
		data.EvaluateParametrization(ratios, mips, true);
		UploadPatches();
		UploadBrokenArcs(*data.Motorcycles(), data.BrokenHalfarcs());

		if (chkPackTexture->checked())
			data.PackTexture(mipFactor);
	});	

	/*(new nanogui::Widget(mainWindow))->setHeight(10); //spacing

	auto renderMeshColorsTextureBtn = new nanogui::Button(mainWindow, "Render Mesh Colors to Texture");
	renderMeshColorsTextureBtn->setCallback([this]() 
	{
		isFileDialogOpen = true;
		std::string colorsFilename = nanogui::file_dialog({ { "color", "Mesh Colors" } }, false);
		if (!colorsFilename.empty())
		{
			data.RenderMeshColorsToTexture(colorsFilename);
			glViewport(0, 0, width(), height());
		}
		isFileDialogOpen = false;
	});	

	auto remeshBtn = new nanogui::Button(mainWindow, "Remesh");
	remeshBtn->setCallback([this]() { data.Remesh(); });

	auto tangentialSmoothBtn = new nanogui::Button(mainWindow, "Tangential Smooth");
	tangentialSmoothBtn->setCallback([this]()
	{
		data.TangentialSmooth(); 
		meshVertices.uploadData(data.Vertices());
	});*/

	auto focusOnPointBtn = new nanogui::Button(displayOptionsBtn->popup(), "Focus on Point");
	focusOnPointBtn->setCallback([this]()
	{
		auto window = new nanogui::Window(this, "Focus on Point");		
		window->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 4, 4));

		auto widget = new nanogui::Widget(window);
		widget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Fill, 4, 4));

		auto& focusPoint = this->camera().GetFocusPoint();
		auto txtX = new nanogui::TextBox(widget, std::to_string(focusPoint.x())); txtX->setEditable(true);
		auto txtY = new nanogui::TextBox(widget, std::to_string(focusPoint.y())); txtY->setEditable(true);
		auto txtZ = new nanogui::TextBox(widget, std::to_string(focusPoint.z())); txtZ->setEditable(true);
		
		auto btnFocus = new nanogui::Button(window, "Focus");
		btnFocus->setCallback([this, txtX, txtY, txtZ, window]()
		{
			this->updateFocus(nullptr);
			try
			{
				this->camera().FocusOnPoint(Eigen::Vector3f(std::stof(txtX->value()), std::stof(txtY->value()), std::stof(txtZ->value())));
				this->removeChild(window);
			}
			catch (...)
			{
				std::cerr << "Cannot update camera focus point." << std::endl;
			}
		});

		performLayout(nvgContext());
		window->setPosition(Eigen::Vector2i(this->width() / 2 - window->width() / 2, this->height() / 2 - window->height() / 2) );

		performLayout(nvgContext());
	});

	auto focusOnVertexBtn = new nanogui::Button(displayOptionsBtn->popup(), "Focus on Vertex");
	focusOnVertexBtn->setCallback([this]()
	{
		auto window = new nanogui::Window(this, "Focus on Point");
		window->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Fill, 4, 4));		

		auto txt = new nanogui::TextBox(window, "0"); txt->setEditable(true);		

		auto btnFocus = new nanogui::Button(window, "Focus");
		btnFocus->setCallback([this, txt, window]()
		{
			this->updateFocus(nullptr);
			try
			{
				auto v = data.Mesh().vertex_handle(std::stoi(txt->value()));
				auto& p = data.Mesh().point(v);
				this->camera().FocusOnPoint(Eigen::Vector3f(p[0], p[1], p[2]));
				this->removeChild(window);
			}
			catch (...)
			{
				std::cerr << "Cannot update camera focus point." << std::endl;
			}
		});

		performLayout(nvgContext());
		window->setPosition(Eigen::Vector2i(this->width() / 2 - window->width() / 2, this->height() / 2 - window->height() / 2));

		performLayout(nvgContext());
	});

	//auto screenshotBtn = new nanogui::Button(mainWindow, "Take Screenshot");
	//screenshotBtn->setCallback([this]()
	//{
	//	const int screenshotWidth = 4096;
	//	const int screenshotHeight = screenshotWidth * height() / width();

	//	isFileDialogOpen = true;
	//	std::string filename = nanogui::file_dialog({ { "png", "PNG" } }, true);
	//	if (!filename.empty())
	//	{
	//		if (filename.substr(filename.length() - 4) != ".png")
	//			filename.append(".png");
	//		TakeScreenshot(filename, screenshotWidth, screenshotHeight);
	//	}
	//	isFileDialogOpen = false;
	//});

	//auto makeTurntableBtn = new nanogui::Button(mainWindow, "Render Turntable");
	//makeTurntableBtn->setCallback([this]()
	//{			
	//	camera().MakeHorizontal();

	//	int frames = 30 /*fps*/ * 5 /*seconds per turn*/;
	//	auto rot = Eigen::Quaternionf(Eigen::AngleAxisf(2 * M_PI / frames, Eigen::Vector3f::UnitY()));
	//	for (int i = 0; i < frames; ++i)
	//	{			
	//		TakeScreenshot(std::string("turntable_") + std::to_string(i) + ".png", 1920, 1080);
	//		camera().RotateAroundFocusPointGlobal(rot);
	//	}
	//});

	performLayout(ctx);
}

void Viewer::TakeScreenshot(const std::string& filename, int screenshotWidth, int screenshotHeight)
{
	GLuint textureFBO, textureColor, depthBuffer;
	glGenFramebuffers(1, &textureFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, textureFBO);

	glGenRenderbuffers(1, &textureColor);
	glBindRenderbuffer(GL_RENDERBUFFER, textureColor);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, screenshotWidth, screenshotHeight);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, textureColor);

	glGenRenderbuffers(1, &depthBuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH32F_STENCIL8, screenshotWidth, screenshotHeight);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);

	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glViewport(0, 0, screenshotWidth, screenshotHeight);
	Eigen::Matrix4f model, view, proj;
	_camera.ComputeCameraMatrices(view, proj, (float)screenshotWidth / screenshotHeight);

	render(view, proj, false);

	std::vector<unsigned char> pixels(screenshotWidth * screenshotHeight * 4);
	glReadPixels(0, 0, screenshotWidth, screenshotHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

	glViewport(0, 0, width(), height());
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glDeleteFramebuffers(1, &textureFBO);
	glDeleteRenderbuffers(1, &textureColor);
	glDeleteRenderbuffers(1, &depthBuffer);

	stbi_write_png(filename.c_str(), screenshotWidth, screenshotHeight, 4, pixels.data() + (screenshotHeight - 1) * screenshotWidth * 4 * sizeof(unsigned char), -screenshotWidth * 4 * sizeof(unsigned char));
}

void Viewer::LoadMesh(const std::string& filename)
{
	indexCount = 0;
	wireframeIndexCount = 0;
	singularityCount = 0;
	for (int i = 0; i < 3; ++i)
	{
		fencedRegionIndexCount[i] = 0;
		fencedRegionFaceIndexCount[i] = 0;
	}	
	motorcycleIndexCount.clear();
	patchIndexCount = 0;
	brokenArcsIndexCount = 0;	

	data.LoadMesh(filename, chkMergeTriangles->checked(), angleThresholdSld->value() * (float)M_PI / 180.0f);
		
	if (chkFocusOnModel->checked())
		camera().FocusOnBBox(data.MeshBoundingBox());

	{
		nse::util::TimedBlock b("Uploading data.Mesh() to GPU ..");
		ShaderPool::Instance()->FlatMeshShader.bind();
		meshVAO.generate();
		meshVAO.bind();

		meshVertices.uploadData(data.Vertices()).bindToAttribute("position");

		UploadMeshFaces();

		meshVAO.unbind();
	}

	UploadWireframeToGPU();	
	UploadSingularities();
}

void Viewer::LoadMeshForVisualization(const std::string& filename)
{
	indexCount = 0;
	wireframeIndexCount = 0;
	singularityCount = 0;
	for (int i = 0; i < 3; ++i)
	{
		fencedRegionIndexCount[i] = 0;
		fencedRegionFaceIndexCount[i] = 0;
	}
	motorcycleIndexCount.clear();
	patchIndexCount = 0;
	brokenArcsIndexCount = 0;

	data.LoadMeshForVisualization(filename, chkScaleTex->checked());

	if(chkFocusOnModel->checked())
		camera().FocusOnBBox(data.MeshBoundingBox());

	{
		nse::util::TimedBlock b("Uploading data.Mesh() to GPU ..");
		ShaderPool::Instance()->FlatMeshShader.bind();
		meshVAO.generate();
		meshVAO.bind();

		meshVertices.uploadData(data.Vertices()).bindToAttribute("position");

		UploadMeshFaces();

		meshVAO.unbind();
	}

	UploadWireframeToGPU();
	UploadSingularities();
	UploadPatches();
}

void Viewer::ExtractQuadLayout()
{
	nse::util::TimedBlock b("Calculating quad layout ..");

	if (chkClassifySingularities->checked())
		data.ClassifyAllSingularities();

	if (chkMergeFencedRegion->checked())
		data.TryToMergePatches();

	if (chkCalculateMotorcycleGraph->checked())
		data.CalculateMotorcycleGraph(chkDeactivateUnnecessaryMotorcycles->checked());

	UploadFencedRegions(data.FencedRegions().patches.begin(), data.FencedRegions().patches.end());
	UploadSingularities();
	if(data.Motorcycles() != nullptr)
		UploadMotorcycles(*data.Motorcycles());	

	if (data.Motorcycles() == nullptr || data.Motorcycles()->Status() != MotorcycleGraph::Finished)
		return;

	if (chkExtractPatches->checked())
	{
		auto stats = data.ExtractPatches(chkSplitNonRectangularPatches->checked());

		b.closeBlock();

		std::sort(stats.patchSizes.begin(), stats.patchSizes.end());
		std::cout << stats;

		UploadPatches();
		
		int nonOriginalEdges = 0;
		for (auto& patch : data.Motorcycles()->Patches())
		{
			for (auto& side : patch.PatchSides())
			{
				for (auto arcIdx : side)
				{
					auto& arc = data.Motorcycles()->Halfarcs()[arcIdx];
					for (auto p : arc)
					{
						auto e = data.Mesh().edge_handle(data.Motorcycles()->MotorcycleHalfedge(p));
						if (!data.IsEdgeOriginal(e))
							++nonOriginalEdges;
					}
				}
			}
		}
		std::cout << "Used " << nonOriginalEdges / 2 << " non-original edges." << std::endl; //every edge is counted twice
	}
}

void Viewer::ProcessDirectory()
{
	bool evaluateOriginalMCG = true;
	bool evaluateUnmergedRegions = false;
	bool evaluateSimpleParametrization = true;
	bool evaluateInvisibleSeamParametrization = true;
	bool renderSnippet = false;
	bool packTexture = true;
	bool saveResult = false;

	std::string filename = nanogui::file_dialog({
		{ "obj", "Wavefront OBJ" },
		{ "ply", "Stanford PLY" },		
	}, false);

	std::string directory = nse::data::parent_path(filename);
	std::vector<std::string> files;
	nse::data::files_in_dir(directory, files);
	
	//for (parametrizationErrorThreshold = 30; parametrizationErrorThreshold <= 90; parametrizationErrorThreshold += 10)
	{
		std::vector<std::string> modelsWithNonRectangles;
		std::vector<std::string> modelsWithFailedTracing;
		std::vector<std::string> modelsWithTimeLimitExceeded;
		Statistics logAreaRatios, mips;

		std::ofstream csv("stats_ParamError" + std::to_string((int)parametrizationErrorThreshold) + ".csv");
		csv << "file;vertices;faces;edges;singularities;";
		if (evaluateOriginalMCG)
			csv << "Patches With Simple MCG;";

		csv << "Classification Time;";

		if (evaluateUnmergedRegions)
			csv << "Meta Singularities; Patches With Classification;";

		csv << "Merge Time; Meta Singularities; Trace Time; Extract Time; Patches With Merge; Non-original edges;Original edges;";

		if (evaluateSimpleParametrization)
			csv << "Simple Parametrization Time; Ratio Percentiles 5-95;;;;;;;;;;;;;;;;;;;MIPS avg; MIPS 90-percentile;";
		if (evaluateInvisibleSeamParametrization)
			csv << "Invisible Seam Parametrization Time; Ratio Percentiles 5-95;;;;;;;;;;;;;;;;;;;MIPS avg; MIPS 90-percentile; Broken Arcs; BrokenArcLength; TotalArcs; TotalArcLength;";

		if (packTexture)
			csv << "Texture Packing Time; ";
		csv << std::endl;

		for (auto& entry : files)
		{
			if (!nse::data::is_directory(entry))
			{
				if ((nse::data::extension(entry) == ".obj" || nse::data::extension(entry) == ".ply") && !nse::data::str_ends_with(entry, ".textured.obj"))
					data.LoadMesh(entry, false, 0);
				else
					continue;

				csv << entry << ";" << data.Vertices().cols() << ";" << data.Faces().size() << ";" << data.Mesh().n_edges() << ";" << data.Singularities().size() << ";";

				nse::util::TimedBlock b("Processing file " + entry + " ..");

				MotorcycleGraph::ExtractionStatistics stats;
				try
				{
					if (evaluateOriginalMCG)
					{
						data.CalculateMotorcycleGraph();
						stats.Clear();
						stats = data.ExtractPatches(false);
						csv << data.Motorcycles()->Patches().size() << ";";
					}

					{
						nse::util::TimedBlock b("Classify ..");
						data.ClassifyAllSingularities();
						csv << b.time() << ";";
					}

					if (evaluateUnmergedRegions)
					{
						data.FindMetaSingularities();
						csv << data.MetaSingularities().size() << ";";
						data.CalculateMotorcycleGraph();
						stats = data.ExtractPatches(false);
						csv << data.Motorcycles()->Patches().size() << ";";
					}

					{
						nse::util::TimedBlock b("Merge ..");
						data.TryToMergePatches();
						csv << b.time() << ";";
					}


					{
						nse::util::TimedBlock b("Tracing ..");
						data.CalculateMotorcycleGraph(true);
						auto t = b.time();
						csv << data.MetaSingularities().size() << ";";
						csv << t << ";";
					}

					if (data.Motorcycles()->Status() != MotorcycleGraph::Finished)
					{
						modelsWithFailedTracing.push_back(entry);
						csv << std::endl;
						continue;
					}

					stats.Clear();
					{
						nse::util::TimedBlock b("Extracting ..");
						stats = data.ExtractPatches(true);
						csv << b.time() << ";";
					}
					csv << data.Motorcycles()->Patches().size() << ";";

					bool hasNonRectangles = false;
					for (int i = 0; i < stats.patchCountPerNumberOfCorners.size(); ++i)
					{
						if (i != 4 && stats.patchCountPerNumberOfCorners[i] > 0)
							hasNonRectangles = true;
					}
					if (hasNonRectangles)
						modelsWithNonRectangles.push_back(entry);

					int nonOriginalEdges = 0;
					int originalEdges = 0;
					for (auto& patch : data.Motorcycles()->Patches())
					{
						for (auto& side : patch.PatchSides())
						{
							for (auto arcIdx : side)
							{
								auto& arc = data.Motorcycles()->Halfarcs()[arcIdx];
								for (auto p : arc)
								{
									auto e = data.Mesh().edge_handle(data.Motorcycles()->MotorcycleHalfedge(p));
									if (!data.IsEdgeOriginal(e))
										++nonOriginalEdges;
									else
										++originalEdges;
								}
							}
						}
					}
					csv << nonOriginalEdges / 2 << ";" << originalEdges << ";";


					if (evaluateSimpleParametrization)
					{
						{
							nse::util::TimedBlock b("Simple Parametrization ..");
							data.CalculateParametrization(10.0f, parametrizationErrorThreshold, discreteOptimizationTimeLimit, Simple);
							csv << b.time() << ";";
						}
						logAreaRatios.Clear();
						data.EvaluateParametrization(logAreaRatios, mips, false);
						for (int i = 1; i < 20; ++i)
						{
							float percent = i * 0.05;
							auto percentile = logAreaRatios.Percentile(percent);
							csv << std::exp(percentile) << ";";
						}
						csv << mips.Average() << ";" << mips.Percentile(0.9f) << ";";
					}

					if (renderSnippet)
					{
						camera().FocusOnBBox(data.MeshBoundingBox());

						ShaderPool::Instance()->FlatMeshShader.bind();
						meshVAO.generate();
						meshVAO.bind();

						meshVertices.uploadData(data.Vertices()).bindToAttribute("position");

						UploadMeshFaces();

						meshVAO.unbind();

						UploadPatches();
						TakeScreenshot(entry + ".png", 512, 512 * height() / width());
					}

					if (evaluateInvisibleSeamParametrization)
					{
						{
							nse::util::TimedBlock b("Invisible Seam Parametrization ..");
							data.CalculateParametrization(10.0f, parametrizationErrorThreshold, discreteOptimizationTimeLimit, GurobiGlobal);
							csv << b.time() << ";";
						}
						logAreaRatios.Clear();
						data.EvaluateParametrization(logAreaRatios, mips, false);
						for (int i = 1; i < 20; ++i)
						{
							float percent = i * 0.05;
							float percentile = logAreaRatios.Percentile(percent);
							csv << std::exp(percentile) << ";";
						}
						csv << mips.Average() << ";" << mips.Percentile(0.9f) << ";";

						int totalArcLength = 0, totalBrokenArclength = 0;
						for (auto& arc : data.Motorcycles()->Halfarcs())
							for (auto p : arc)
								++totalArcLength;
						for (auto idx : data.BrokenHalfarcs())
							for (auto p : data.Motorcycles()->Halfarcs()[idx])
								++totalBrokenArclength;
						csv << data.BrokenHalfarcs().size() << ";" << totalBrokenArclength << ";" << data.Motorcycles()->Halfarcs().size() / 2 << ";" << totalArcLength / 2 << ";";
					}

					if (packTexture)
					{
						nse::util::TimedBlock b("Packing ..");
						data.PackTexture(1 << mipLevel);
						csv << b.time() << ";";
					}

					//try
					//{
					//	data.CalculateParametrization(10.0f, AllInvalid, GurobiGlobal);
					//}
					//catch (std::runtime_error& e)
					//{
					//	std::cout << e.what() << std::endl;
					//	modelsWithTimeLimitExceeded.push_back(entry);
					//	continue;
					//}

					//if(packTexture)
					//	data.PackTexture();

					if (saveResult)
						data.SaveMeshOBJ(entry + ".textured.obj");
				}
				catch (std::runtime_error& error)
				{
					std::cout << "Runtime Error: " << error.what() << std::endl;
				}
				catch (...)
				{
					std::cout << "Unknown error." << std::endl;
				}
				csv << std::endl;

				data.Clear();
			}
		}
		csv.close();
		std::cout << "Models with non-rectangles: " << modelsWithNonRectangles.size() << std::endl;
		for (auto& entry : modelsWithNonRectangles)
			std::cout << "  " << entry << std::endl;
		std::cout << "Models with failed tracing: " << modelsWithFailedTracing.size() << std::endl;
		for (auto& entry : modelsWithFailedTracing)
			std::cout << "  " << entry << std::endl;
		std::cout << "Models with time limit exceeded: " << modelsWithTimeLimitExceeded.size() << std::endl;
		for (auto& entry : modelsWithTimeLimitExceeded)
			std::cout << "  " << entry << std::endl;
	}
}

//EmitVertexFunctor: void(const HEMesh::HalfedgeHandle[3]) //the to-vertices of the halfedges are the triangle corners
template <typename EmitTriangleFunctor> 
void TriangulateMeshFace(HEMesh::FaceHandle f, const HEMesh& mesh, EmitTriangleFunctor&& emitTriangle)
{
	OpenMesh::HalfedgeHandle base;
	for (auto h : mesh.fh_range(f))
	{
		if (base.idx() == -1)
		{
			base = h;
			continue;
		}
		auto nextH = mesh.next_halfedge_handle(h);
		if (nextH == base)
			break;
		else
		{
			HEMesh::HalfedgeHandle triangle[3] = { base, h, nextH };
			std::forward<EmitTriangleFunctor>(emitTriangle)(triangle);
		}
	}
}

void Viewer::UploadMeshFaces()
{
	std::vector<unsigned int> indices;
	indices.reserve(4 * data.Faces().size());
	for (auto f : data.Mesh().faces()) //for each face
	{
		TriangulateMeshFace(f, data.Mesh(), [&](const HEMesh::HalfedgeHandle h[3])
		{
			indices.push_back(data.Mesh().to_vertex_handle(h[0]).idx());
			indices.push_back(data.Mesh().to_vertex_handle(h[1]).idx());
			indices.push_back(data.Mesh().to_vertex_handle(h[2]).idx());
		});		
	}
	indexCount = (unsigned int)indices.size();
	meshIndices.uploadData(indexCount, 1, sizeof(unsigned int), GL_UNSIGNED_INT, true, reinterpret_cast<uint8_t*>(indices.data()));
	meshIndices.bind();
}

void Viewer::UploadWireframeToGPU()
{
	nse::util::TimedBlock b("Uploading wireframe to GPU ..");
	std::vector<unsigned int> indices;
	for (auto& e : data.Mesh().edges())
	{
		if (!data.IsEdgeOriginal(e))
			continue;

		auto h0 = data.Mesh().halfedge_handle(e, 0);
		auto h1 = data.Mesh().halfedge_handle(e, 1);

		auto v0 = data.Mesh().to_vertex_handle(h0);
		auto v1 = data.Mesh().to_vertex_handle(h1);

		indices.push_back(v0.idx());
		indices.push_back(v1.idx());
	}
	wireframeIndexCount = (unsigned int)indices.size();

	ShaderPool::Instance()->SimpleShader.bind();
	wireframeVAO.generate();
	wireframeVAO.bind();
	meshVertices.bindToAttribute("position");
	wireframeIndices.uploadData(indices.size(), 1, sizeof(unsigned int), GL_UNSIGNED_INT, true, reinterpret_cast<uint8_t*>(indices.data()));
	wireframeVAO.unbind();
}

void Viewer::UploadSingularities()
{
	std::vector<Eigen::Vector4f> singPos, singCol;
	primitiveIdToSingularityId.clear();

	auto addVertex = [&](HEMesh::VertexHandle v, bool realSingularity, HEMesh::HalfedgeHandle contextOutgoingEdge, int valenceDefect)
	{
		float radius = 0;
		int edgeValence = 0;
		for (auto h : data.Mesh().voh_range(v))
		{
			radius += data.Mesh().calc_edge_length(h);
			++edgeValence;
		}
		radius /= edgeValence;
		radius *= (realSingularity ? 0.3f : 0.1f) * sldSingularitySize->value();
		auto p = data.Mesh().point(v);
		if (contextOutgoingEdge.is_valid())
			p = p + data.Mesh().calc_edge_vector(contextOutgoingEdge).normalized() * radius;
		singPos.emplace_back();
		for (int j = 0; j < 3; ++j)
			singPos.back()(j) = p[j];
		singPos.back()(3) = radius;
		singCol.push_back(realSingularity ? GetValenceColor(valenceDefect) : Eigen::Vector4f(0.1f, 0.8f, 0.1f, 1.0f));
	};
	
	for (int i = 0; i < data.Singularities().size(); ++i)
	{
		if (data.MetaSingularities().size() == 0)
		{
			addVertex(data.Singularities()[i].vertexHandle, true, data.Singularities()[i].contextOutgoingBoundaryEdge, data.Singularities()[i].valenceDefect);
			primitiveIdToSingularityId[singPos.size() - 1] = i;
		}
		else
		{
			bool realSingularity = false;
			const size_t* ptr;
			if (!realSingularity)
			{
				if (data.Singularities()[i].contextOutgoingBoundaryEdge.is_valid())
					realSingularity = data.VertexToMetaSingularity().TryAccessAtContextEdge(data.Singularities()[i].contextOutgoingBoundaryEdge, ptr);
				else
					realSingularity = data.VertexToMetaSingularity().TryAccessAtManifoldVertex(data.Singularities()[i].vertexHandle, ptr);
			}

			if (!realSingularity)
			{
				addVertex(data.Singularities()[i].vertexHandle, realSingularity, data.Singularities()[i].contextOutgoingBoundaryEdge, data.Singularities()[i].valenceDefect);
				primitiveIdToSingularityId[singPos.size() - 1] = i;
			}
		}
	}

	for (auto& s : data.MetaSingularities())
	{
		auto valence = s.Degree();
		for (int color = 0; color < valence; ++color)
		{
			auto& motorcycle = s.GetEmanatingMotorcycle(color);
			if (motorcycle.size() > 0)
			{
				auto e = motorcycle.front();
				HEMesh::HalfedgeHandle context;
				if (!data.Mesh().is_manifold(data.Mesh().from_vertex_handle(e)))
				{
					context = data.Mesh().opposite_halfedge_handle(e);
					CirculateBackwardUntil<true>(context, data.Mesh(), [&](HEMesh::HalfedgeHandle) {return false; });
				}
				addVertex(data.Mesh().from_vertex_handle(e), true, context, valence - 4);
				const size_t* singIdx;
				if(data.VertexToMetaSingularity().TryAccessAtToVertex(data.Mesh().opposite_halfedge_handle(e), singIdx))
					primitiveIdToSingularityId[singPos.size() - 1] = *singIdx;
				break;
			}
		}
	}

	ShaderPool::Instance()->SphereShader.bind();
	singularityVAO.generate();
	singularityVAO.bind();
	singularityPositions.uploadData(singPos).bindToAttribute("positionRadius");
	singularityColors.uploadData(singCol).bindToAttribute("color");
	singularityVAO.unbind();
	singularityCount = singPos.size();
}

int32_t Viewer::GetPickingIndex(const Eigen::Vector2i p) const
{	
	if (p.x() >= 0 && p.y() >= 0 && p.x() < width() && p.y() < height())
	{
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fboSinglePixel);
		glDrawBuffer(GL_COLOR_ATTACHMENT0);
		glReadBuffer(GL_COLOR_ATTACHMENT1);

		//copy the interesting pixel to a non-multisample frame buffer
		auto readX = p.x();
		auto readY = height() - 1 - p.y();
		glBlitFramebuffer(readX, readY, readX + 1, readY + 1, 0, 0, 1, 1, GL_COLOR_BUFFER_BIT, GL_NEAREST);

		//copy the pixel to CPU memory
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fboSinglePixel);
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		int32_t index;
		glReadPixels(0, 0, 1, 1, GL_RED_INTEGER, GL_INT, &index);

		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		return index;
	}
	else
		return -1;
}

bool Viewer::mouseMotionHook(const Eigen::Vector2i & p, const Eigen::Vector2i & rel, int button, int modifiers)
{
	return false;
}

bool Viewer::mouseButtonHook(const Eigen::Vector2i & p, int button, bool down, int modifiers)
{		
	if (isFileDialogOpen)
		return true;
	if (!down && button == 0 && modifiers == GLFW_MOD_CONTROL)
	{
		int pickedSingularity = -1;
		auto singIt = primitiveIdToSingularityId.find(GetPickingIndex(p));
		if (singIt != primitiveIdToSingularityId.end())
			pickedSingularity = singIt->second;
		if (pickedSingularity != -1)
		{
			FencedRegion calculatedPatch(data.Mesh(), data.VertexToSingularity());
			const FencedRegion* patchPtr;
			selectedRegion = (size_t)-1;
			if (data.Singularities()[pickedSingularity].correspondingRegion != (size_t)-1)
			{
				selectedRegion = data.Singularities()[pickedSingularity].correspondingRegion;
				patchPtr = &data.FencedRegions().patches[selectedRegion];
			}
			else
			{
				std::cout << "Calculating new patch" << std::endl;
				auto classificationResult = data.ClassifySingularity(data.Singularities()[pickedSingularity], calculatedPatch);
				patchPtr = &calculatedPatch;
				if (classificationResult == NoPatch)
				{
					std::cout << "Cannot grow patch." << std::endl;
					for (int i = 0; i < 3; ++i)
					{
						fencedRegionIndexCount[i] = 0;
						fencedRegionFaceIndexCount[i] = 0;
					}
				}
			}

			if (!patchPtr->IsRectilinear())
			{
				//MotorcycleGraph data.graph(data.Mesh(), fencedRegions, metaSingularities, vertexToMetaSingularity);
				//InitializeMotorcyclesForSingularPatch(*patchPtr, data.graph);
				//UploadMotorcycles(data.graph);
			}

			//clicked a singularity
			std::cout << "Corresponding patch of singularity: " << data.Singularities()[pickedSingularity].correspondingRegion << std::endl;
			for (auto& l : patchPtr->Loops())
				std::cout << "Loop with turn count of " << l.Degree() << std::endl;

			if(patchPtr->IsRectilinear())
			{
				std::cout << "Patch is rectilinear and covers " << patchPtr->CoveredSingularities().size() << " singularities." << std::endl;
			}
			else
			{
				std::cout << "Patch is not rectilinear and covers " << patchPtr->CoveredSingularities().size() << " singularities." << std::endl;
			}

			UploadFencedRegions(patchPtr, patchPtr + 1, true);

			return true;
		}
	}
	return false;
}

template <typename Iterator>
void Viewer::UploadFencedRegions(Iterator itBegin, Iterator itEnd, bool includeDeactivated)
{
	nse::util::TimedBlock b("Uploading fenced regions ..");

	for (int i = 0; i < 3; ++i)
	{
		fencedRegionIndexCount[i] = 0;
		fencedRegionFaceIndexCount[i] = 0;
	}	

	Eigen::Matrix<unsigned int, 2, Eigen::Dynamic> indices[3];
	std::vector<unsigned int> faceIndices[3];
	for (auto it = itBegin; it != itEnd; ++it)
	{
		const FencedRegion& patch = *it;
		if (!patch.IsActive() && !includeDeactivated)
			continue;
		int matrixIndex;
		if (patch.Loops()[0].Degree() < 4)
			matrixIndex = 0;
		else if (patch.Loops()[0].Degree() == 4)
			matrixIndex = 1;
		else
			matrixIndex = 2;

		auto patchIndexMatrix = patch.GenerateBoundaryIndices();
		indices[matrixIndex].conservativeResize(Eigen::NoChange, indices[matrixIndex].cols() + patchIndexMatrix.cols());
		indices[matrixIndex].rightCols(patchIndexMatrix.cols()) = patchIndexMatrix;

		for (auto f : patch.IncludedFaces())
		{
			TriangulateMeshFace(f, data.Mesh(), [&](HEMesh::HalfedgeHandle h[3])
			{
				for (int i = 0; i < 3; ++i)
					faceIndices[matrixIndex].push_back(data.Mesh().to_vertex_handle(h[i]).idx());
			});
		}
		//auto patchFaceIndexData = patch.GenerateFaceIndices();
		//faceIndices[matrixIndex].insert(faceIndices[matrixIndex].end(), patchFaceIndexData.begin(), patchFaceIndexData.end());
	}

	for (int i = 0; i < 3; ++i)
	{
		fencedRegionIndexCount[i] = indices[i].size();
		fencedRegionFaceIndexCount[i] = faceIndices[i].size();
	}	

	indices[0].conservativeResize(Eigen::NoChange, indices[0].cols() + indices[1].cols());
	indices[0].rightCols(indices[1].cols()) = indices[1];
	indices[0].conservativeResize(Eigen::NoChange, indices[0].cols() + indices[2].cols());
	indices[0].rightCols(indices[2].cols()) = indices[2];
	faceIndices[0].insert(faceIndices[0].end(), faceIndices[1].begin(), faceIndices[1].end());
	faceIndices[0].insert(faceIndices[0].end(), faceIndices[2].begin(), faceIndices[2].end());

	ShaderPool::Instance()->SimpleShader.bind();
	fencedRegionsVAO.generate();
	fencedRegionsVAO.bind();
	meshVertices.bind();
	meshVertices.bindToAttribute("position");
	fencedRegionIndices.uploadData(indices[0]);
	fencedRegionsVAO.unbind();

	fencedRegionsFaceVAO.generate();
	fencedRegionsFaceVAO.bind();
	meshVertices.bind();
	meshVertices.bindToAttribute("position");
	fencedRegionFaceIndices.uploadData(faceIndices[0].size(), 1, sizeof(unsigned int), GL_UNSIGNED_INT, true, reinterpret_cast<uint8_t*>(faceIndices[0].data()));
	fencedRegionsFaceVAO.unbind();
}

void Viewer::UploadPatches()
{
	nse::util::TimedBlock b("Uploading patches ..");

	data.FindPatchColors();

	std::vector<Eigen::Vector4f> positions;
	std::vector<Eigen::Vector4f> colors;
	std::vector<Eigen::Vector2f> texCoords;
	int patchIndex = 0;
	for (int i = 0; i < data.Motorcycles()->Patches().size(); ++i)
	{
		const TexturePatch& patch = data.Motorcycles()->Patches()[i];
		const TextureCoordinatesStorage& texCoordStorage = data.AccessTexCoords(i);

		Eigen::Vector4f color;
		switch (cmbPatchColoring->selectedIndex())
		{
		case 0:
			color = data.GetPatchColor(patchIndex++);
			color.w() = 0.7f;
			break;
		case 1:
			color = GetValenceColor(patch.Corners() - 4);
			break;
		}		
		for (HEMesh::FaceHandle f : patch.Faces())
		{
			TriangulateMeshFace(f, data.Mesh(), [&](const HEMesh::HalfedgeHandle h[3])
			{
				for (int i = 0; i < 3; ++i)
				{
					auto v = data.Mesh().to_vertex_handle(h[i]);
					auto pos = data.Mesh().point(v);
					positions.push_back(Eigen::Vector4f(pos[0], pos[1], pos[2], 1));
					colors.push_back(color);
					texCoords.push_back(texCoordStorage.TexCoordAtToVertex(h[i], data.Mesh()));
				}
			});
		}
		/*for (auto& tri : patch.filledHoles)
		{
			for (int i = 0; i < 3; ++i)
			{
				auto pos = data.Mesh().point(tri[i]);
				positions.push_back(Eigen::Vector4f(pos[0], pos[1], pos[2], 1));
				colors.push_back(color);
				texCoords.push_back(patch.TexCoordAtInnerVertex(tri[i], data.Mesh()));
			}
		}*/
	}
	ShaderPool::Instance()->SimpleShader.bind();
	patchVAO.generate();
	patchVAO.bind();
	patchPositionsVBO.uploadData(positions).bindToAttribute("position");
	patchColorsVBO.uploadData(colors).bindToAttribute("color");
	patchTexCoordsVBO.uploadData(texCoords).bindToAttribute("texCoords");
	patchVAO.unbind();
	patchIndexCount = positions.size();
}

void Viewer::UploadMotorcycles(const MotorcycleGraph& graph)
{
	nse::util::TimedBlock b("Uploading motorcycles ..");

	//Upload to GPU
	motorcycleIndexCount.clear();
	std::vector<unsigned int> indices;
	for (auto& cycle : graph.Motorcycles())
	{
		motorcycleIndexCount.push_back(2 * cycle.Path().size());
		for (auto h : cycle.Path())
		{
			indices.push_back(data.Mesh().from_vertex_handle(h).idx());
			indices.push_back(data.Mesh().to_vertex_handle(h).idx());
		}
	}

	ShaderPool::Instance()->SimpleShader.bind();
	motorcycleVAO.generate();
	motorcycleVAO.bind();
	meshVertices.bind();
	meshVertices.bindToAttribute("position");
	motorcycleIndices.uploadData(sizeof(unsigned int) * indices.size(), indices.data());
	motorcycleVAO.unbind();
}

void Viewer::UploadBrokenArcs(const MotorcycleGraph& graph, const std::vector<size_t>& brokenArcs)
{
	brokenArcsIndexCount = 0;
	
	std::vector<unsigned int> indices;
	for (auto iArc : brokenArcs)
	{
		auto& arc = graph.Halfarcs()[iArc];		
		for (auto segment : arc)
		{
			auto h = graph.MotorcycleHalfedge(segment);
			indices.push_back(data.Mesh().from_vertex_handle(h).idx());
			indices.push_back(data.Mesh().to_vertex_handle(h).idx());
		}
	}
	
	ShaderPool::Instance()->SimpleShader.bind();
	brokenArcsVAO.generate();
	brokenArcsVAO.bind();
	meshVertices.bind();
	meshVertices.bindToAttribute("position");
	brokenArcsIndices.uploadData(sizeof(unsigned int) * indices.size(), indices.data());
	brokenArcsVAO.unbind();
	brokenArcsIndexCount = indices.size();
}

void Viewer::drawContents()
{
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	Eigen::Matrix4f model, view, proj;
	_camera.ComputeCameraMatrices(view, proj);

	render(view, proj);
	
#ifndef LAPTOP_COMPAT
	glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glDrawBuffer(GL_BACK);
	glBlitFramebuffer(0, 0, width(), height(), 0, 0, width(), height(), GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT, GL_NEAREST);
	
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
#endif
}

void Viewer::DrawTextIfNotHidden(const Eigen::Vector4f& pos, const Eigen::Matrix4f mvp, const std::string& text, float tolerance)
{
	Eigen::Vector4f projectedPosition = mvp * pos;
	projectedPosition /= projectedPosition.w();
	projectedPosition.x() = (projectedPosition.x() + 1) / 2 * width();
	projectedPosition.y() = (1 - projectedPosition.y()) / 2 * height();

	if (projectedPosition.x() < 0 || projectedPosition.x() >= width()
		|| projectedPosition.y() < 0 || projectedPosition.y() >= height())
		return;

	Eigen::Vector3f reprojectedPosition;
	float depth = get3DPosition(Eigen::Vector2i((int)projectedPosition.x(), (int)projectedPosition.y()), reprojectedPosition);

	if ((pos.head<3>() - reprojectedPosition).norm() <= tolerance)
		nvgText(mNVGContext, projectedPosition.x(), projectedPosition.y(), text.c_str(), nullptr);
}

void Viewer::render(const Eigen::Matrix4f & mv, const Eigen::Matrix4f & proj, bool background)
{
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
	glBlendEquationSeparate(GL_FUNC_ADD, GL_MAX);

	Eigen::Matrix4f mvp = proj * mv;

#ifndef LAPTOP_COMPAT
	GLenum bufs[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT0 + 1 };
	if (background)
	{		
		glDrawBuffers(2, bufs);

		ShaderPool::Instance()->ClearShader.bind();
		glDepthMask(GL_FALSE);
		emptyVAO.bind();
		glDrawArrays(GL_TRIANGLES, 0, 3);
		emptyVAO.unbind();
		glDepthMask(GL_TRUE);

		glDrawBuffers(1, bufs);
	}
#endif

	Eigen::Vector3f cameraPosition = mv.topLeftCorner<3, 3>().transpose() * (-mv.topRightCorner<3, 1>());

	if (data.Vertices().cols() > 0)
	{
		if (chkShowFaces->checked())
		{
			auto& shader = ShaderPool::Instance()->FlatMeshShader;
			shader.bind();
			shader.setUniform("mv", mv);
			shader.setUniform("mvp", mvp);
			meshVAO.bind();
			glDrawElements(GL_TRIANGLES, indexCount, GL_UNSIGNED_INT, 0);
			meshVAO.unbind();
		}

		if (chkShowWireframe->checked())
		{
			/*auto& shader2 = ShaderPool::Instance()->SimpleShader;
			shader2.bind();
			shader2.setUniform("mvp", mvp);
			shader2.setUniform("color", Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f));
			wireframeVAO.bind();
			glDepthRange(0.0, 0.999995f);
			glDrawElements(GL_LINES, wireframeIndexCount, GL_UNSIGNED_INT, 0);
			glDepthRange(0.0, 1.0f);*/

			auto radius = 0.03f * data.AverageEdgeLength() * sldWireframeWidth->value();

			wireframeVAO.bind();
			Eigen::Vector4f color(0, 0, 0, 1);
			auto& shader = ShaderPool::Instance()->CylinderShader;
			shader.bind();
			shader.setUniform("mvp", mvp);
			shader.setUniform("radius", radius);
			shader.setUniform("cameraPosition", cameraPosition);
			shader.setUniform("color", color);
			shader.setUniform("shaded", 0);
			glDrawElements(GL_LINES, wireframeIndexCount, GL_UNSIGNED_INT, 0);

			auto& shader2 = ShaderPool::Instance()->SphereShader;
			shader2.bind();
			shader2.setUniform("mv", mv);
			shader2.setUniform("p", proj);
			shader2.setUniform("constantRadius", 1u);
			shader2.setUniform("radius", radius);
			shader2.setUniform("constantColor", 1u);
			shader2.setUniform("color", color);
			shader2.setUniform("shaded", 0);
			glDrawElements(GL_POINTS, wireframeIndexCount, GL_UNSIGNED_INT, 0);

			wireframeVAO.unbind();
		}
	}

	if (motorcycleIndexCount.size() > 0 && chkShowMotorcycleGraph->checked() && data.Motorcycles() != nullptr)
	{
		float globalRadius = 0.05f * data.AverageEdgeLength() * sldTubeWidth->value();

		unsigned int start = 0;
		for (int i = 0; i < motorcycleIndexCount.size(); ++i)
		{
			Eigen::Vector4f color = GetSegmentColor(i);
			float radius = globalRadius * (data.Motorcycles()->Motorcycles()[i].isDeactivated ? 0.3f : 1.0f);

			auto& shader = ShaderPool::Instance()->CylinderShader;
			shader.bind();
			shader.setUniform("mvp", mvp);
			shader.setUniform("radius", radius);			
			shader.setUniform("cameraPosition", cameraPosition);
			shader.setUniform("color", color);
			shader.setUniform("shaded", 1);

			motorcycleVAO.bind();
			glDrawElements(GL_LINES, motorcycleIndexCount[i], GL_UNSIGNED_INT, (void*)(start * 4));

			auto& shader2 = ShaderPool::Instance()->SphereShader;
			shader2.bind();
			shader2.setUniform("mv", mv);
			shader2.setUniform("p", proj);
			shader2.setUniform("constantRadius", 1u);
			shader2.setUniform("radius", radius);
			shader2.setUniform("constantColor", 1u);
			shader2.setUniform("color", color);
			shader2.setUniform("shaded", 1);
			glDrawElements(GL_POINTS, motorcycleIndexCount[i], GL_UNSIGNED_INT, (void*)(start * 4));
			motorcycleVAO.unbind();

			start += motorcycleIndexCount[i];
		}
	}

	if (brokenArcsIndexCount > 0 && chkShowBrokenArcs->checked())
	{
		Eigen::Vector4f color(1.0f, 0.0f, 0.0f, 1.0f);
		float radius = 0.10f * data.AverageEdgeLength() * sldTubeWidth->value();

		auto& shader = ShaderPool::Instance()->CylinderShader;
		shader.bind();
		shader.setUniform("mvp", mvp);
		shader.setUniform("radius", radius);
		shader.setUniform("cameraPosition", cameraPosition);
		shader.setUniform("color", color);
		shader.setUniform("shaded", 1);

		brokenArcsVAO.bind();
		glDrawElements(GL_LINES, brokenArcsIndexCount, GL_UNSIGNED_INT, 0);

		auto& shader2 = ShaderPool::Instance()->SphereShader;
		shader2.bind();
		shader2.setUniform("mv", mv);
		shader2.setUniform("p", proj);
		shader2.setUniform("constantRadius", 1u);
		shader2.setUniform("radius", radius);
		shader2.setUniform("constantColor", 1u);
		shader2.setUniform("color", color);
		shader2.setUniform("shaded", 1);
		glDrawElements(GL_POINTS, brokenArcsIndexCount, GL_UNSIGNED_INT, 0);
		brokenArcsVAO.unbind();
	}

	if (singularityCount > 0 && chkShowSingularities->checked())
	{		
#ifndef LAPTOP_COMPAT
		glDrawBuffers(2, bufs);
#endif

		auto& shader = ShaderPool::Instance()->SphereShader;
		shader.bind();
		shader.setUniform("mv", mv);
		shader.setUniform("p", proj);
		shader.setUniform("constantRadius", 0u);
		shader.setUniform("constantColor", 0u);
		shader.setUniform("shaded", 1);
		singularityVAO.bind();
		glDrawArrays(GL_POINTS, 0, (GLsizei)singularityCount);
		singularityVAO.unbind();

#ifndef LAPTOP_COMPAT
		glDrawBuffers(1, bufs);
#endif
	}

	if (patchIndexCount > 0 && chkShowPatches->checked())
	{
		auto& shader = ShaderPool::Instance()->SimpleShader;
		shader.bind();
		shader.setUniform("mvp", mvp);
		shader.setUniform("useUniformColor", 0);
		shader.setUniform("visualizeTexCoords", chkShowUV->checked() ? 1 : 0);
		shader.setUniform("texCoordScale", std::pow(2.0f, -sldParameterizationGridSize->value()));

		glDepthRange(0.0, 0.999999f);

		patchVAO.bind();
		glDrawArrays(GL_TRIANGLES, 0, patchIndexCount);
		patchVAO.unbind();

		glDepthRange(0.0, 1.0f);
	}

	if (fencedRegionIndexCount[0] + fencedRegionIndexCount[1] + fencedRegionIndexCount[2] > 0 && chkShowFencedRegions->checked())
	{				
		float radius = 0.02f * data.AverageEdgeLength() * sldTubeWidth->value();
		
		auto& shader = ShaderPool::Instance()->SimpleShader;
		shader.bind();
		shader.setUniform("mvp", mvp);
		shader.setUniform("useUniformColor", 1);
		shader.setUniform("visualizeTexCoords", 0);

		//(0.933, 0.839, 0.686
		//const Eigen::Vector4f colors[] = { Eigen::Vector4f(0.1f, 0.7f, 0.1f, 0.4f), Eigen::Vector4f(0.7f, 0.1f, 0.1f, 0.4f) };
		Eigen::Vector4f colors[] = { GetValenceColor(-2), Eigen::Vector4f(0.4f, 0.9f, 0.4f, 1.0f), GetValenceColor(2) };
		for (int i = 0; i < 3; ++i)
			colors[i].w() = 0.6f;

		glDepthRange(0.0, 0.999998f);

		fencedRegionsFaceVAO.bind();
		unsigned int offset = 0;
		for (int i = 0; i < 3; ++i)
		{
			if (fencedRegionFaceIndexCount[i] > 0)
			{
				shader.setUniform("color", colors[i]);
				glDrawElements(GL_TRIANGLES, fencedRegionFaceIndexCount[i], GL_UNSIGNED_INT, (void*)(offset * 4));
			}
			offset += fencedRegionFaceIndexCount[i];
		}		

		glDepthRange(0.0, 1.0f);

		fencedRegionsVAO.bind();
		
		auto& shader2 = ShaderPool::Instance()->CylinderShader;
		shader2.bind();
		shader2.setUniform("mvp", mvp);
		shader2.setUniform("radius", radius);
		shader2.setUniform("cameraPosition", cameraPosition);
		shader2.setUniform("shaded", 1);

		offset = 0;
		for (int i = 0; i < 3; ++i)
		{			
			if (fencedRegionIndexCount[i] > 0)
			{
				auto color = colors[i]; color.w() = 1;
				shader2.setUniform("color", color);
				glDrawElements(GL_LINES, fencedRegionIndexCount[i], GL_UNSIGNED_INT, (void*)(offset * 4));
			}
			offset += fencedRegionIndexCount[i];
		}
		
		auto& shader3 = ShaderPool::Instance()->SphereShader;
		shader3.bind();
		shader3.setUniform("mv", mv);
		shader3.setUniform("p", proj);
		shader3.setUniform("constantRadius", 1u);
		shader3.setUniform("radius", radius);
		shader3.setUniform("constantColor", 1u);
		shader3.setUniform("shaded", 1);
		
		offset = 0;
		for (int i = 0; i < 3; ++i)
		{
			if (fencedRegionIndexCount[i] > 0)
			{
				auto color = colors[i]; color.w() = 1;
				shader3.setUniform("color", color);
				glDrawElements(GL_POINTS, fencedRegionIndexCount[i], GL_UNSIGNED_INT, (void*)(offset * 4));
			}
			offset += fencedRegionIndexCount[i];
		}		

		fencedRegionsVAO.unbind();		
	}	

	bool drawRegionIndices = data.FencedRegions().patches.size() > 0 && chkShowFencedRegionIndices->checked();
	bool drawMotorcycleIndices = data.Motorcycles() != nullptr && chkShowMotorcycleIndices->checked();
	bool drawVertexIndices = data.Vertices().cols() > 0 && chkShowVertexIndices->checked();

	if (drawRegionIndices || drawMotorcycleIndices || drawVertexIndices)
	{
#ifndef LAPTOP_COMPAT
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		glDrawBuffer(GL_BACK);
		glBlitFramebuffer(0, 0, width(), height(), 0, 0, width(), height(), GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT, GL_NEAREST);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
#endif

		nvgBeginFrame(mNVGContext, mSize[0], mSize[1], mPixelRatio);
		nvgFontSize(mNVGContext, 14.0f);
		nvgFontFace(mNVGContext, "sans-bold");
		nvgTextAlign(mNVGContext, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
		nvgFillColor(mNVGContext, { 0.1f, 0.1f, 0.1f, 1.0f });
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

		float tolerance = 0.5 * data.AverageEdgeLength();

		if (drawVertexIndices)
		{
			for (int i = 0; i < data.Mesh().n_vertices(); ++i)
			{
				Eigen::Vector4f pos; pos.setZero();
				
				auto p = data.Mesh().point(data.Mesh().vertex_handle(i));
				pos.x() = p[0];
				pos.y() = p[1];
				pos.z() = p[2];
				pos.w() = 1;

				DrawTextIfNotHidden(pos, mvp, std::to_string(i), tolerance);
			}
		}


		if (drawRegionIndices)
		{
			for (int i = 0; i < data.FencedRegions().patches.size(); ++i)
			{
				Eigen::Vector4f pos; pos.setZero();
				if (data.FencedRegions().patches[i].Loops().size() == 0)
					continue;
				auto& e = data.FencedRegions().patches[i].Loops()[0].Edges()[0];

				auto p = data.Mesh().point(data.Mesh().from_vertex_handle(e.edge));
				pos.x() = p[0];
				pos.y() = p[1];
				pos.z() = p[2];
				pos.w() = 1;

				DrawTextIfNotHidden(pos, mvp, std::to_string(i), tolerance);
			}			
		}

		if (drawMotorcycleIndices)
		{
			auto& motorcycles = data.Motorcycles()->Motorcycles();
			for (int i = 0; i < motorcycles.size(); ++i)
			{
				if (motorcycles[i].Path().size() == 0)
					continue;
				auto h = motorcycles[i].Path()[motorcycles[i].Path().size() / 2];
				Eigen::Vector4f pos; pos.setZero();
				auto p = data.Mesh().point(data.Mesh().from_vertex_handle(h));
				if (motorcycles[i].Path().size() % 2 == 1)
					p = 0.6f * p +  0.4f * data.Mesh().point(data.Mesh().to_vertex_handle(h));
				pos.x() = p[0];
				pos.y() = p[1];
				pos.z() = p[2];
				pos.w() = 1;

				DrawTextIfNotHidden(pos, mvp, std::to_string(i), tolerance);
			}
		}

		nvgEndFrame(mNVGContext);
	}
}