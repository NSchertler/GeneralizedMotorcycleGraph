#pragma once

#include <vector>
#include <array>

#include <Eigen/Core>

//Triangulates a 3D polygon according to the DP algorithm presented in 
//  Filling gaps in the boundary of a polyhedron, Barequet and Sharir, Computer Aided Geometric Design 1995
extern std::vector<std::array<int, 3>> TriangulatePolygon(const std::vector<Eigen::Vector3f>& polygon);