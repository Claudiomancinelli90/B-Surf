#version 330
layout(location = 0) in vec3 vposition;
out vec3 position;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;

void main() {
    position = (frame * vec4(vposition, 1)).xyz;
    gl_Position = projection * view * vec4(position, 1);
}