#version 330

in vec4 position;
in vec4 color;
in vec2 texCoords;

out vec4 vertexColor;
out vec2 vertexTexCoords;

uniform mat4 mvp;

void main()
{
	gl_Position = mvp * position;
	vertexColor = color;
	vertexTexCoords = texCoords;
}