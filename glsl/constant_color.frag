#version 330

in vec4 vertexColor;
in vec2 vertexTexCoords;

uniform vec4 color;
uniform bool useUniformColor;
uniform bool visualizeTexCoords;
uniform float texCoordScale;

out vec4 result;

void main()
{
	result = (useUniformColor ? color : vertexColor);
	if(visualizeTexCoords)
	{
		//result.xyz *= fract(texCoordScale * vertexTexCoords.y);
		int cellId = (int(texCoordScale * vertexTexCoords.x) + int(texCoordScale * vertexTexCoords.y));
		if(cellId % 2 == 0)
			;//result.xyz *= 0.5 + 0.5 * fract(texCoordScale * vertexTexCoords.x);
		else
			result.xyz *= 0.3;
			//result.xyz *= 0.3 + 0.2 * fract(texCoordScale * vertexTexCoords.y);
	}
}