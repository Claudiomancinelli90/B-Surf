#version 330
layout(location = 0) in vec2 vpos;
layout(location = 1) in vec3 vcolor;
out vec2      Pos;
out vec3      Color;
uniform mat2  frame;
uniform vec2  center;
uniform float ratio = 1;

void main() {
    Pos = inverse(frame) * (vpos + center);
    Pos.x /= ratio;
    Color       = vcolor;
    gl_Position = vec4(Pos.xy, 0, 1);
}