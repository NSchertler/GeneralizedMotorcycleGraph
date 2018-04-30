#version 330

layout(location=0) out vec4 result;
layout(location=1) out int id;

in fData
{
	vec4 toPixel;
	vec4 cam;
	vec4 color;
	float radius;
} frag;

uniform mat4 p;
uniform bool constantColor;
uniform vec4 color;
uniform bool shaded;

void main(void)
{
	vec3 v = frag.toPixel.xyz - frag.cam.xyz;
	vec3 e = frag.cam.xyz;
	float ev = dot(e, v);
	float vv = dot(v, v);
	float ee = dot(e, e);
	float rr = frag.radius * frag.radius;

	float radicand = ev * ev - vv * (ee - rr);
	if(radicand < 0)
		discard;
	
	float rt = sqrt(radicand);
	

	float lambda = max(0, (-ev - rt) / vv);
	float lambda2 = (-ev + rt) / vv;
	if(lambda2 < lambda)
		discard;

	vec3 hit = lambda * v;
	vec3 normal = (frag.cam.xyz + hit) / frag.radius;
	
	vec4 proj = p * vec4(hit, 1);
	gl_FragDepth = ((gl_DepthRange.diff * proj.z / proj.w) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

	vec4 myColor = (constantColor ? color : frag.color);

	if(shaded)
	{
		vec3 vNormalized = -normalize(v);
		float nDotL = dot(vNormalized, normal);
		vec3 c = myColor.rgb * nDotL + vec3(0.5, 0.5, 0.5) * pow(nDotL, 120);	
		result = vec4(c, myColor.a);
	}
	else
		result = myColor;	

	id = gl_PrimitiveID;
}