/*
	@author Wenzel Jakob
	@author Nico Schertler
*/

#pragma once

#include "common.h"
#include <nsessentials/util/IndentationLog.h>
#include <nsessentials/data/FileHelper.h>

#include <fstream>

#ifndef _WIN32
#include <libgen.h>
#endif


extern void LoadMesh(const std::string &filename, Matrix3Xf& V, FaceList& F);

extern void load_obj(const std::string &filename, Matrix3Xf &V, FaceList& F, Matrix2Xf &T, FaceList& FT);