#version 330 compatibility

layout(location=0) out vec4 color;
layout(location=1) out int id;

in vec4 clipPos;

void main()
{
	float t = (clipPos.y + 1) / 2;
	color = (1 - t) * vec4(0.1, 0.1f, 0.3f, 1.00) + t * vec4(0.5, 0.6, 0.8, 1.0);	
	id = -1;
}