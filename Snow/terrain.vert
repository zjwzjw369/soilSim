#version 400 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;
uniform mat4 mView;
uniform mat4 projection;
uniform mat4 model;

out vec2 coord;

void main() {
    gl_Position = projection * mView * model *vec4(position, 1.0);

    coord = texCoord;
}  