#version 330 core

layout (location = 0) in vec3 vertexPos;

uniform mat4 projection;
uniform mat4 mView;
uniform float pointRadius;
uniform float pointScale;

out vec3 pos;

void main() {
	vec4 viewPos = mView * vec4(vertexPos, 1.0);
    gl_Position = projection * viewPos;
	pos = viewPos.xyz;
	gl_PointSize = pointScale*0.5f * (pointRadius / gl_Position.w);
}