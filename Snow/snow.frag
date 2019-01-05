#version 330 core

in vec3 pos;
in vec3 pColor;
uniform mat4 mView;
uniform mat4 projection;
uniform float pointRadius;

out vec4 fragColor;

float rand(vec2 co){
  return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void main() {
	//calculate normal
	vec3 normal;
	normal.xy = gl_PointCoord * 2.0 - 1.0;
	float r2 = dot(normal.xy, normal.xy);
	
	if (r2 > 1.0) {
		discard;
	}
	
	//fragColor = vec4(1);
	//fragColor =vec4(0.8667,0.6627,0.4,1);
	//fragColor =vec4(0.5,0.5,1.0,1);
	fragColor=vec4(pColor,1.0f);
}