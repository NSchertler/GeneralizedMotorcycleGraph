#pragma once

#include <nsessentials/gui/GLShader.h>

//Singleton. Holds all shaders during the runtime of the application.
class ShaderPool
{
private:
	static ShaderPool* _instance;
	ShaderPool();

public:
	static ShaderPool* Instance();
	void CompileAll();

	//Shader that renders a mesh using flat shading
	nse::gui::GLShader FlatMeshShader;

	//Shader that renders a sphere using raycasting
	nse::gui::GLShader SphereShader;

	//Shader that renders objects with a constant color without illumination
	nse::gui::GLShader SimpleShader;

	//Shader used to clear the render target
	nse::gui::GLShader ClearShader;

	//Shader that renders cylinders using raycasting
	nse::gui::GLShader CylinderShader;

	//Shader used to render mesh colors to texture
	nse::gui::GLShader MeshColorsQuadShader;

	//Shader used to render mesh colors to texture
	nse::gui::GLShader MeshColorsTriShader;
};