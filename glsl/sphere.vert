#version 330

/*
	This file is part of the implementation for the technical paper

		Field-Aligned Online Surface Reconstruction
		Nico Schertler, Marco Tarini, Wenzel Jakob, Misha Kazhdan, Stefan Gumhold, Daniele Panozzo
		ACM TOG 36, 4, July 2017 (Proceedings of SIGGRAPH 2017)

	Use of this source code is granted via a BSD-style license, which can be found
	in License.txt in the repository root.

	@author Nico Schertler
*/

in vec4 positionRadius;
in vec4 color;

out vData
{
	vec4 pos;
	vec4 color;
	float radius;
} vertex;

uniform mat4 mv;
uniform bool constantRadius;
uniform float radius;

void main() 
{
	vec4 p = vec4(positionRadius.xyz, 1);
	vertex.pos = mv * p;
	vertex.radius = (constantRadius ? radius : positionRadius.w * length(mv[0].xyz));
	vertex.color = color;
}
