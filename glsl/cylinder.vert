#version 330

in vec3 position;
in vec4 color;

out vData
{
	vec3 pos;
	vec4 color;
} vertex;

uniform mat4 mv;

void main() 
{	
	vertex.pos = position;
	vertex.color = color;
}
