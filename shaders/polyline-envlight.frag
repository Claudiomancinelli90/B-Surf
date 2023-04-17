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

uniform samplerCube envlight_irradiance;
uniform samplerCube envlight_reflection;
uniform sampler2D   envlight_brdflut;

vec3 sample_prefiltered_refleciton(vec3 incoming, float roughness) {
  int   MAX_REFLECTION_LOD = 5;
  float lod                = sqrt(roughness) * MAX_REFLECTION_LOD;
  return textureLod(envlight_reflection, incoming, lod).rgb;
}

void main() {
  vec3 radiance = vec3(0.0);
  vec3 N        = normalize(Normal);

  vec3 V = normalize(eye - Position);
  vec3 L = normalize(vec3(1, 1, -1));
  vec3 H = normalize(V + L);

  // diffuse
  radiance += color * textureLod(envlight_irradiance, N, 0).rgb;

  // specular
  vec3  incoming   = normalize(reflect(-V, N));
  float roughness  = 0.5;
  vec3  specular   = vec3(0.1);
  vec3  reflection = sample_prefiltered_refleciton(incoming, roughness);
  vec2  env_brdf =
      texture(envlight_brdflut, vec2(max(dot(N, V), 0.0), roughness)).rg;
  radiance += reflection * (specular * env_brdf.x + env_brdf.y);

  radiance *= 3;

  // final color correction
  radiance = pow(radiance * pow(2, exposure), vec3(1 / gamma));

  result = vec4(radiance, 1);
}