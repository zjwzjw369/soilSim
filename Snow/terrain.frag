#version 400 core
uniform sampler2D texture0;
in vec2 coord;
out vec4 color;

void main() {
    color = texture2D(texture0,coord);
	//color=vec4(1.0,1.0,1.0,1.0);
}