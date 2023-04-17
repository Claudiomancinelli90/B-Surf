#version 330
layout(location = 0) in vec2 vpos;
out vec2 Pos;

void main() {
    Pos         = vpos;
    gl_Position = vec4(vpos.xy, 0, 1);
}