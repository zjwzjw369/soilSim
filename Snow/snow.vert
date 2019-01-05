#version 330 core

layout (location = 0) in vec3 vertexPos;
layout (location = 1) in vec3 vertexColor;

uniform mat4 projection;
uniform mat4 mView;
uniform float pointRadius;
uniform float pointScale;

out vec3 pos;
out vec3 pColor;

void main() {
	vec4 viewPos = mView * vec4(vertexPos, 1.0);
    gl_Position = projection * viewPos;
	pos = viewPos.xyz;
	gl_PointSize = pointScale*0.5f * (pointRadius / gl_Position.w);
	pColor=vertexColor;
}