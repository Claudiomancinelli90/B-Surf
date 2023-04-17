#version 330
layout(location = 0) in vec3 vposition;
layout(location = 1) in vec3 vcolor;
out vec3 Position;
out vec3 Color;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;

void main() {
  Position    = (frame * vec4(vposition, 1)).xyz;
  Color       = vcolor;
  gl_Position = projection * view * vec4(Position, 1);
}