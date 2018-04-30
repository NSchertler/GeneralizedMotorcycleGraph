#version 330

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

in vData
{
	vec3 pos;
	vec4 color;
} vertex[];

out fData
{
	vec2 cameraInGlyphSpace;
	vec3 posInGlyphSpace;
	float halfHeight;
	vec3 cylinderCenter;
	mat3 glyphSpaceToWorldSpace;
} outputV;

uniform float radius;

uniform vec3 cameraPosition;
uniform mat4 mvp;

void main(void)
{
	//cylinder axis
	vec3 x = vertex[1].pos - vertex[0].pos;
	outputV.halfHeight = length(x) * 0.5;
	x = normalize(x);
	outputV.cylinderCenter = 0.5 * (vertex[0].pos + vertex[1].pos);
	vec3 centerToCamera = cameraPosition - outputV.cylinderCenter;

	outputV.cameraInGlyphSpace.x =  dot(x, centerToCamera);
	vec3 cameraProjectedOnAxis = outputV.cameraInGlyphSpace.x * x;
	
	vec3 cameraRejected = centerToCamera - cameraProjectedOnAxis;
	outputV.cameraInGlyphSpace.y = length(cameraRejected);
	if(outputV.cameraInGlyphSpace.y <= radius)
		return; //camera is inside elongated cylinder -> outside is not visible

	//construct coordinate system
	vec3 y = normalize(cameraRejected);
	vec3 z = cross(x, y);
	outputV.glyphSpaceToWorldSpace = mat3(x, y, z);
	
	float quadLowerX = - outputV.halfHeight;
	float quadUpperX = outputV.halfHeight;
	
	if(outputV.cameraInGlyphSpace.x < quadLowerX)
		quadLowerX = -outputV.halfHeight + radius * (outputV.cameraInGlyphSpace.x + outputV.halfHeight) / outputV.cameraInGlyphSpace.y;;
	if(outputV.cameraInGlyphSpace.x > quadUpperX)
		quadUpperX = outputV.halfHeight + radius * (outputV.cameraInGlyphSpace.x - outputV.halfHeight) / outputV.cameraInGlyphSpace.y;

	outputV.posInGlyphSpace = vec3(quadLowerX, radius, -radius);
	gl_Position = mvp * vec4(outputV.glyphSpaceToWorldSpace * outputV.posInGlyphSpace + outputV.cylinderCenter, 1);
	EmitVertex();
	outputV.posInGlyphSpace = vec3(quadLowerX, radius, radius);
	gl_Position = mvp * vec4(outputV.glyphSpaceToWorldSpace * outputV.posInGlyphSpace + outputV.cylinderCenter, 1);
	EmitVertex();
	outputV.posInGlyphSpace = vec3(quadUpperX, radius, -radius);
	gl_Position = mvp * vec4(outputV.glyphSpaceToWorldSpace * outputV.posInGlyphSpace + outputV.cylinderCenter, 1);
	EmitVertex();
	outputV.posInGlyphSpace = vec3(quadUpperX, radius, radius);
	gl_Position = mvp * vec4(outputV.glyphSpaceToWorldSpace * outputV.posInGlyphSpace + outputV.cylinderCenter, 1);
	EmitVertex();
}