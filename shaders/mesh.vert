#version 330
layout(location = 0) in vec3 vposition;
layout(location = 1) in vec3 vnormal;
layout(location = 2) in float vfield;
out vec3  position;
out vec3  normal;
out float field;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;

void main() {
    position    = (frame * vec4(vposition, 1)).xyz;
    normal      = (frame * vec4(vnormal, 0)).xyz;
    field       = vfield;
    gl_Position = projection * view * vec4(position, 1);
}