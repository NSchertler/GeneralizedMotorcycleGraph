cmake_minimum_required (VERSION 2.8)

project("RectangleBinPack")

include_directories(include/RectangleBinPack)

add_library(RectangleBinPack STATIC			
			src/GuillotineBinPack.cpp
			src/MaxRectsBinPack.cpp
			src/Rect.cpp
			src/ShelfBinPack.cpp
			src/ShelfNextFitBinPack.cpp
			src/SkylineBinPack.cpp)