#version 330

layout(location = 0) in vec3 positions;
layout(location = 1) in vec3 normals;
layout(location = 2) in vec3 instance_from;
layout(location = 3) in vec3 instance_to;

uniform mat4 frame;

uniform mat4 view;
uniform mat4 projection;

uniform vec3  eye;
uniform float size;

out vec3 Position;
out vec3 Normal;

// main function
void main() {
  // copy values
  Position = positions;
  Normal   = normals;
  float dist = length(eye - instance_from);
  Position.xy *= atan(dist) * size;

  // world projection
  Position = (frame * vec4(Position, 1)).xyz;
  vec3 dir = instance_to - instance_from;
  if (dir != vec3(0)) {
    vec3 up = abs(dir.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent   = normalize(cross(up, dir));
    vec3 bitangent = normalize(cross(dir, tangent));

    mat3 mat;
    mat[2]   = dir;
    mat[0]   = tangent;
    mat[1]   = bitangent;
    Position = mat * Position;
    Normal   = mat * Normal;
  }
  Position += instance_from;

  // clip
  gl_Position = projection * view * vec4(Position, 1);
}