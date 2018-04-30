cmake_minimum_required (VERSION 3.1)
project(RegularMeshTexturing)

FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

FUNCTION(JOIN VALUES GLUE OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
ENDFUNCTION()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

set(LIBS)

# Enable folders for projects in Visual Studio
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

#add_definitions(-DREDUCED_TIMINGS)

# Sanitize build environment for static build with C++11
if (MSVC)
  add_definitions (/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (/D "__TBB_NO_IMPLICIT_LINKAGE")

  add_definitions (-DNOMINMAX )
  add_definitions(/D_USE_MATH_DEFINES)

  # Parallel build on MSVC (all targets)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

  if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")

    # Disable Eigen vectorization for Windows 32 bit builds (issues with unaligned access segfaults)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DEIGEN_DONT_ALIGN")
  endif()

  # Static build
  set(CompilerFlags
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
  foreach(CompilerFlag ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
  endforeach()
elseif(APPLE)
  # Try to auto-detect a suitable SDK
  execute_process(COMMAND bash -c "xcodebuild -version -sdk | grep MacOSX | grep Path | head -n 1 | cut -f 2 -d ' '" OUTPUT_VARIABLE CMAKE_OSX_SYSROOT)
  string(REGEX REPLACE "(\r?\n)+$" "" CMAKE_OSX_SYSROOT "${CMAKE_OSX_SYSROOT}")
  string(REGEX REPLACE "^.*X([0-9.]*).sdk$" "\\1" CMAKE_OSX_DEPLOYMENT_TARGET "${CMAKE_OSX_SYSROOT}")
endif()

set (CMAKE_CXX_STANDARD 14)

# Build NanoGUI
set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_SHARED  OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_USE_GLAD ON CACHE BOOL " " FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui ext_build/nanogui)
set_property(TARGET nanogui nanogui-obj glfw glfw_objects PROPERTY FOLDER "dependencies")
include_directories(ext/nanogui/include
					${NANOGUI_EXTRA_INCS})
set(LIBS ${LIBS} nanogui ${NANOGUI_EXTRA_LIBS})
set(NSE_WITH_NANOGUI ON CACHE BOOL " " FORCE)
set(NSE_WITH_EIGEN ON CACHE BOOL " " FORCE)
add_definitions(${NANOGUI_EXTRA_DEFS})

#NS Essentials
SET(NSE_SSBO_SUPPORT ON CACHE BOOL " " FORCE)
SET(NSE_TESSELLATION_SUPPORT ON CACHE BOOL " " FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/nsessentials ext_build/nsessentials)
set_property(TARGET nsessentials PROPERTY FOLDER "dependencies")
include_directories(ext/nsessentials/include)
set(LIBS nsessentials ${LIBS})
add_definitions(${NSE_EXTRA_DEFS})

#OpenMesh
set(BUILD_APPS OFF CACHE BOOL " " FORCE)
set(OPENMESH_BUILD_PYTHON_BINDINGS OFF CACHE BOOL " " FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/OpenMesh ext_build/OpenMesh)
set(LIBS ${LIBS} OpenMeshCore OpenMeshTools)
include_directories(ext/OpenMesh/src)
set_property(TARGET OpenMeshCore OpenMeshTools PROPERTY FOLDER "dependencies")

find_package(OpenMP)
if(OPENMP_FOUND)
  add_definitions(-DOPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#Gurobi
find_package(Gurobi)
if(GUROBI_FOUND)	
	include_directories(${GUROBI_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${GUROBI_LIBRARY} optimized ${GUROBI_CXX_LIBRARY_RELEASE} debug ${GUROBI_CXX_LIBRARY_DEBUG})	
	add_definitions (/D "WITH_GUROBI")
endif()

#nlopt
#find_package(NLopt)
#if(NLOPT_INCLUDE_DIRS)
#	include_directories(${NLOPT_INCLUDE_DIRS})
#	set(LIBS ${LIBS} NLopt::nlopt)
#	add_definitions (/D "WITH_NLOPT")
#endif()

#RectPack
add_subdirectory(ext/RectangleBinPack ext_build/RectangleBinPack)
set(LIBS ${LIBS} RectangleBinPack)
include_directories(ext/RectangleBinPack/include)
set_property(TARGET RectangleBinPack PROPERTY FOLDER "dependencies")

include_directories(ext/libigl)
include_directories(ext/scaffold-map)
add_library(scaffold-map STATIC ext/scaffold-map/igl/scaf.cpp)
set(LIBS ${LIBS} scaffold-map)
set_property(TARGET scaffold-map PROPERTY FOLDER "dependencies")

#triangle
add_subdirectory(ext/libigl/external/triangle)
set(LIBS ${LIBS} triangle)
include_directories(ext/libigl/external/triangle)
set_property(TARGET triangle PROPERTY FOLDER "dependencies")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  # RPLY mesh I/O library
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/rply  
)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/ext/stb)

set(GLSL_FILES	mesh_flat.vert flat_shading.geom blinnphong.frag
				sphere.vert sphere.geom sphere.frag
				cylinder.vert cylinder.geom cylinder.frag
				simple.vert constant_color.frag
				clear.vert clear.frag
				void.vert
				mesh_colors_quad.tcs mesh_colors_quad.tes
				mesh_colors_tri.tcs mesh_colors_tri.tes)

PREPEND(GLSL_FILES "${CMAKE_CURRENT_SOURCE_DIR}/glsl/" ${GLSL_FILES})

source_group(glsl FILES ${GLSL_FILES})

JOIN("${GLSL_FILES}" "," glsl_string)
set(bin2c_cmdline
    -DOUTPUT_C=glsl.cpp
    -DOUTPUT_H=glsl.h
    "-DINPUT_FILES=${glsl_string}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/resources/bin2c.cmake")

add_custom_command(
  OUTPUT glsl.cpp glsl.h
  COMMAND ${CMAKE_COMMAND} ARGS ${bin2c_cmdline}
  DEPENDS ${GLSL_FILES}
  COMMENT "Running bin2c"
  PRE_BUILD VERBATIM)

include_directories( ${CMAKE_CURRENT_BINARY_DIR}  )

SET(THIRD_PARTY 	ext/rply/rply.c)
source_group(third_party FILES ${THIRD_PARTY})

add_executable(RegularMeshTexturing MACOSX_BUNDLE
	glsl.cpp
	src/main.cpp	
	src/common.cpp include/common.h
	src/gui/Viewer.cpp include/gui/Viewer.h
	src/meshio.cpp include/meshio.h
	src/gui/ShaderPool.cpp include/gui/ShaderPool.h
	src/MotorcycleGraph.cpp include/MotorcycleGraph.h
	src/FencedRegion.cpp include/FencedRegion.h
	src/TexturePatch.cpp include/TexturePatch.h
	src/FencedRegionRepresentativeCalculator.cpp include/FencedRegionRepresentativeCalculator.h
	include/MotorcycleOptionsInPatch.h
	src/Parametrization.cpp include/Parametrization.h
	src/ParametrizationData.cpp include/ParametrizationData.h
	include/TargetLengthStrategySimple.h
	src/MultiplierStrategyIterativeRounding.cpp include/MultiplierStrategyIterativeRounding.h
	src/ArclengthStrategyGurobi.cpp include/ArclengthStrategyGurobi.h
	src/ArclengthStrategyGurobiSeparate.cpp include/ArclengthStrategyGurobiSeparate.h
	src/ParametrizationStrategyLSCMForInteriorScaffoldForBoundary.cpp include/ParametrizationStrategyLSCMForInteriorScaffoldForBoundary.h
	src/ParametrizationHelper.cpp include/ParametrizationHelper.h
	src/PolygonTriangulation.cpp include/PolygonTriangulation.h
	src/Data.cpp include/Data.h
	src/Singularity.cpp include/Singularity.h
	src/Statistics.cpp include/Statistics.h
	include/ManifoldnessAwareVertexMap.h
	src/stbi_impl.cpp
	${GLSL_FILES}
	${THIRD_PARTY})

target_link_libraries(RegularMeshTexturing ${LIBS})