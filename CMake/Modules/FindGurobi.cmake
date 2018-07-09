#### Based on https://github.com/joschu/trajopt/blob/master/cmake/modules/FindGUROBI.cmake


# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

find_path(GUROBI_INCLUDE_DIRS 
          NAMES gurobi_c++.h
          PATHS ENV GUROBI_HOME
          PATH_SUFFIXES include)

set(RELEASE_SUFFIX)
set(DEBUG_SUFFIX)
if(MSVC)
	if(MSVC_VERSION EQUAL 1900)
		set(RELEASE_SUFFIX mt2015)
		set(DEBUG_SUFFIX md2015)
	elseif(MSVC_VERSION EQUAL 1910 OR MSVC_VERSION EQUAL 1911)
		set(RELEASE_SUFFIX mt2017)
		set(DEBUG_SUFFIX md2017)
	endif()
endif()		 
		 
find_library( GUROBI_LIBRARY 
              NAMES		    
				gurobi50 
				gurobi51
				gurobi52
				gurobi55
				gurobi56
				gurobi60
				gurobi65
				gurobi70
				gurobi75
              PATHS ENV GUROBI_HOME
			  PATH_SUFFIXES lib)

find_library( GUROBI_CXX_LIBRARY_RELEASE
              NAMES gurobi_c++${RELEASE_SUFFIX}
              PATHS ENV GUROBI_HOME
			  PATH_SUFFIXES lib)
			  
find_library( GUROBI_CXX_LIBRARY_DEBUG
              NAMES gurobi_c++${DEBUG_SUFFIX}
              PATHS ENV GUROBI_HOME
			  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GUROBI_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                  GUROBI_CXX_LIBRARY_RELEASE GUROBI_CXX_LIBRARY_DEBUG GUROBI_INCLUDE_DIRS)