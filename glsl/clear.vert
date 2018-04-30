#version 330 core

const vec2 vertices[3] = vec2[]( vec2(-1.0, 1.0), vec2(-1.0, -3.0), vec2(3.0, 1.0));

out vec4 clipPos;

void main()
{
    clipPos = vec4(vertices[gl_VertexID], 0.0, 1.0);
	gl_Position = clipPos;
}