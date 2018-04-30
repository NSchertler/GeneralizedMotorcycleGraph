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

uniform mat4 mv;
uniform mat4 mvp;

in vec3 position;

out GEOM_IN
{
	vec4 pos;
	vec4 color;	
} vOut;


void main() 
{	
	vec4 p = vec4(position, 1);
	vOut.pos = mv * p;
	gl_Position = mvp * p;
	vOut.color = vec4(0.933, 0.839, 0.686, 1);
}