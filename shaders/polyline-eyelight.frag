#version 330

out vec4 result;
in vec3  Position;
in vec3  Normal;

uniform mat4  frame;
uniform mat4  view;
uniform mat4  projection;
uniform vec3  color;
uniform vec3  eye;
uniform float gamma;
uniform float exposure;

void main() {
  vec3 radiance = vec3(0.0);
  vec3 N        = normalize(Normal);

  radiance = color;

  // final color correction
  radiance = pow(radiance * pow(2, exposure), vec3(1 / gamma));

  result = vec4(radiance, 1);
}