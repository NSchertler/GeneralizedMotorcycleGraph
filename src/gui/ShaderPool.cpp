#include "gui/ShaderPool.h"

#include "glsl.h"

ShaderPool* ShaderPool::_instance(nullptr);

ShaderPool::ShaderPool()
{
}

ShaderPool* ShaderPool::Instance()
{
	if (_instance == nullptr)
		_instance = new ShaderPool();
	return _instance;
}

void ShaderPool::CompileAll()
{	
	FlatMeshShader.init("FlatMeshShader", std::string((const char*)mesh_flat_vert, mesh_flat_vert_size), std::string((const char*)blinnphong_frag, blinnphong_frag_size), std::string((const char*)flat_shading_geom, flat_shading_geom_size));	
	SphereShader.init("Sphere Shader", std::string((const char*)sphere_vert, sphere_vert_size), std::string((const char*)sphere_frag, sphere_frag_size), std::string((const char*)sphere_geom, sphere_geom_size));
	SimpleShader.init("Simple Shader", std::string((const char*)simple_vert, simple_vert_size), std::string((const char*)constant_color_frag, constant_color_frag_size));
	ClearShader.init("Clear Shader", std::string((const char*)clear_vert, clear_vert_size), std::string((const char*)clear_frag, clear_frag_size));
	CylinderShader.init("Cylinder Shader", std::string((const char*)cylinder_vert, cylinder_vert_size), std::string((const char*)cylinder_frag, cylinder_frag_size), std::string((const char*)cylinder_geom, cylinder_geom_size));

	//MeshColorsQuadShader.initWithTessellation("MeshColorsQuadShader", std::string((const char*)void_vert, void_vert_size), std::string((const char*)mesh_colors_quad_tcs, mesh_colors_quad_tcs_size), std::string((const char*)mesh_colors_quad_tes, mesh_colors_quad_tes_size), std::string((const char*)constant_color_frag, constant_color_frag_size));
	//MeshColorsTriShader.initWithTessellation("MeshColorsTriShader", std::string((const char*)void_vert, void_vert_size), std::string((const char*)mesh_colors_tri_tcs, mesh_colors_tri_tcs_size), std::string((const char*)mesh_colors_tri_tes, mesh_colors_tri_tes_size), std::string((const char*)constant_color_frag, constant_color_frag_size));
}