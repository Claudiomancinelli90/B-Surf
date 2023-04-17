#version 330
out vec4 result;
in vec3  Position;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 color;

vec3 gamma(vec3 c) { return pow(clamp(c, 0, 1), vec3(1 / 2.2)); }

void main() { result = vec4(gamma(color), 1); }