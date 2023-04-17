#version 330

layout(location = 0) in vec3 positions;
layout(location = 1) in vec3 centers;

uniform mat4  frame;
uniform mat4  view;
uniform mat4  projection;
uniform vec3 eye;
uniform float size;

out vec3 Position;
out vec3 Center;

// main function
void main() {
  // copy values
  Position = positions;
  float dist = length(eye - centers);
  Position *= atan(dist) * size;
  Center   = centers;

  // world projection
  // Position = (frame * vec4(Position, 1)).xyz;
  Position += centers;
  // clip
  gl_Position = projection * view * vec4(Position, 1);
}