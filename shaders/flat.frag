#version 330
in vec2 Pos;
in vec3 Color;

out vec4 result;

void main() { result = vec4(Color, 1); }
