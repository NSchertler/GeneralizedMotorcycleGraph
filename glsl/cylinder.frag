#version 330
#extension GL_ARB_conservative_depth: enable

in fData
{
	vec2 cameraInGlyphSpace;
	vec3 posInGlyphSpace;
	float halfHeight;
	vec3 cylinderCenter;
	mat3 glyphSpaceToWorldSpace;
} frag;

uniform float radius;

uniform mat4 mvp;

out vec4 result;
layout (depth_greater) out float gl_FragDepth;

uniform vec4 color;
uniform bool shaded;

void main()
{
	vec3 camPos = vec3(frag.cameraInGlyphSpace, 0);
	vec3 dir = frag.posInGlyphSpace - camPos;
	//quadratic equation
	float a = dir.y * dir.y + dir.z * dir.z;
	float b = 2 * frag.cameraInGlyphSpace.y * dir.y;
	float c = frag.cameraInGlyphSpace.y * frag.cameraInGlyphSpace.y - radius * radius;

	b /= a;
	c /= a;

	float discriminant = 0.25 * b * b - c;
	if(discriminant < 0)
		discard;

	float t = -b / 2 - sqrt(discriminant);

	vec3 intersectionPos = camPos + t * dir;
	if(intersectionPos.x < -frag.halfHeight || intersectionPos.x > frag.halfHeight)
		discard;
	
		vec3 normal = intersectionPos;
		normal.x = 0;
		normal = normalize(normal);
		dir = normalize(dir);

		vec4 proj = mvp * vec4(frag.glyphSpaceToWorldSpace * intersectionPos + frag.cylinderCenter, 1);
		gl_FragDepth = ((gl_DepthRange.diff * proj.z / proj.w) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

	if(shaded)
	{
		float nDotL = dot(-dir, normal);
		vec3 col = color.rgb * nDotL + vec3(0.5, 0.5, 0.5) * pow(nDotL, 120);

		result = vec4(col, color.a);
	}
	else
		result = color;
}