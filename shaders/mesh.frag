#version 330
out vec4 result;
in vec3  position;
in vec3  normal;
in float field;

uniform mat4 frame;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 color;

vec3 gamma(vec3 c) { return pow(clamp(c, 0, 1), vec3(1 / 2.2)); }

void main() {
    vec3 c = color * (normal.y * 0.5 + 0.5);
    // c      = mix(vec3(1, 0, 0), vec3(0, 0, 1), field);
    c.xyz += 0.1 * sin(10000 * field);
    // result.xyz += cos(1000 * field);
    // result.xyz *= 0;
    result = vec4(gamma(c), 1);
}