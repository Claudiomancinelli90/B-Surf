#version 330
layout(location = 0) in vec3 vposition;
out vec3 Position;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;

void main() {
  Position    = (frame * vec4(vposition, 1)).xyz;
  gl_Position = projection * view * vec4(Position, 1);
}