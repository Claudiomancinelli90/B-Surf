#include "gpu.h"

// #include <graphics/common.h>
// #include <graphics/commonio.h>
#include <yocto/yocto_common.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <cassert>

//

#include <algorithm>
#include <cstdarg>
#include <deque>
#include <mutex>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "ext/glad/glad.h"


#ifdef _WIN32
#undef near
#undef far
#endif

namespace gpu {
using namespace yocto;

bool init_opengl() {
  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error("Cannot initialize OpenGL context.");
  glDepthFunc(GL_LEQUAL);
  return true;
}

void check_error() {
  // auto error = glGetError();
  // if (error != GL_NO_ERROR) printf("gl error: %d (%x)\n", error, error);
  assert(glGetError() == GL_NO_ERROR);
}

void clear_framebuffer(const vec4f& color, bool clear_depth) {
  assert(glGetError() == GL_NO_ERROR);
  glClearColor(color.x, color.y, color.z, color.w);
  if (clear_depth) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
  } else {
    glClear(GL_COLOR_BUFFER_BIT);
  }
}

void set_viewport(const vec4i& viewport) {
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);
}

void set_wireframe(bool enabled) {
  if (enabled)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void set_blending(bool enabled) {
  if (enabled) {
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
  } else {
    glDisable(GL_BLEND);
  }
}

void set_depth_test(DepthTest flag) {
  int glflag;
  switch (flag) {
    case DepthTest::never: {
      glflag = GL_NEVER;
      break;
    }
    case DepthTest::less: {
      glflag = GL_LESS;
      break;
    }
    case DepthTest::equal: {
      glflag = GL_EQUAL;
      break;
    }
    case DepthTest::lequal: {
      glflag = GL_LEQUAL;
      break;
    }
    case DepthTest::greater: {
      glflag = GL_GREATER;
      break;
    }
    case DepthTest::notequal: {
      glflag = GL_NOTEQUAL;
      break;
    }
    case DepthTest::gequal: {
      glflag = GL_GEQUAL;
      break;
    }
    case DepthTest::always: {
      glflag = GL_ALWAYS;
      break;
    }
  }
  glDepthFunc(glflag);
}

void set_point_size(int size) { glPointSize(size); }

void load_shader_code(Shader& shader) {
  auto error = string{};
  load_text(shader.vertex_filename, shader.vertex_code, error);
  load_text(shader.fragment_filename, shader.fragment_code, error);
}

Shader make_shader_from_file(const string& vertex_filename,
    const string& fragment_filename, bool abort_on_error) {
  Shader shader;
  shader.vertex_filename   = vertex_filename;
  shader.fragment_filename = fragment_filename;
  load_shader_code(shader);
  init_shader(shader, abort_on_error);
  return shader;
}

// bool init_shader(Shader& shader, bool abort_on_error) {
//   delete_shader(shader);
//   return init_shader(shader, shader.vertex_code.c_str(),
//       shader.fragment_code.c_str(), abort_on_error);
// }

bool init_shader(Shader& shader, bool abort_on_error) {
  check_error();
  const char* vertex   = shader.vertex_code.data();
  const char* fragment = shader.fragment_code.data();
  int         errflags;
  char        errbuf[10000];

  // create vertex
  check_error();
  auto vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_id, 1, &vertex, NULL);
  glCompileShader(vertex_id);
  glGetShaderiv(vertex_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(vertex_id, 10000, 0, errbuf);
    errbuf[6] = '\n';
    printf("\n*** VERTEX SHADER COMPILATION: %s\n %s\n",
        shader.vertex_filename.c_str(), errbuf);
    if (abort_on_error) {
      throw std::runtime_error("shader compilation failed\n");
    }
    return false;
    glDeleteProgram(vertex_id);
  }
  check_error();

  // create fragment
  check_error();
  auto fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_id, 1, &fragment, NULL);
  glCompileShader(fragment_id);
  glGetShaderiv(fragment_id, GL_COMPILE_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(fragment_id, 10000, 0, errbuf);
    errbuf[6] = '\n';
    printf("\n*** FRAGMENT SHADER COMPILATION: %s\n%s\n",
        shader.fragment_filename.c_str(), errbuf);
    if (abort_on_error) {
      throw std::runtime_error("shader compilation failed\n");
    }
    glDeleteProgram(fragment_id);
    return false;
  }
  check_error();

  // create shader
  check_error();
  auto shader_id = glCreateProgram();
  glAttachShader(shader_id, vertex_id);
  glAttachShader(shader_id, fragment_id);
  glLinkProgram(shader_id);
  glValidateProgram(shader_id);
  glGetProgramiv(shader_id, GL_LINK_STATUS, &errflags);
  if (!errflags) {
    glGetShaderInfoLog(fragment_id, 10000, 0, errbuf);
    //    errbuf[6] = '\n';
    printf("\n*** SHADER LINKING %s\n", errbuf);
    if (abort_on_error) {
      throw std::runtime_error("shader linking failed\n");
    }
    glDeleteProgram(shader_id);
    return false;
  }

  delete_shader(shader);
  shader.shader_id   = shader_id;
  shader.vertex_id   = vertex_id;
  shader.fragment_id = fragment_id;
  check_error();
  return true;
}

void delete_shader(Shader& shader) {
  glDeleteProgram(shader.shader_id);
  glDeleteShader(shader.vertex_id);
  glDeleteShader(shader.fragment_id);
  shader.shader_id   = 0;
  shader.vertex_id   = 0;
  shader.fragment_id = 0;
}

void init_texture(Texture& texture, const vec2i& size, bool as_float,
    bool as_srgb, bool linear, bool mipmap) {
  if (texture) delete_texture(texture);
  check_error();
  glGenTextures(1, &texture.id);
  texture.size     = size;
  texture.mipmap   = mipmap;
  texture.is_srgb  = as_srgb;
  texture.is_float = as_float;
  glBindTexture(GL_TEXTURE_2D, texture.id);
  if (as_float) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  } else if (as_srgb) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, size.x, size.y, 0, GL_RGBA,
        GL_UNSIGNED_BYTE, nullptr);
  } else {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, size.x, size.y, 0, GL_RGBA,
        GL_FLOAT, nullptr);
  }
  if (mipmap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        (linear) ? GL_LINEAR : GL_NEAREST);
  }
  check_error();
}

void init_texture(Texture& texture, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap) {
  init_texture(texture, img.imsize(), as_float, false, linear, mipmap);
  update_texture(texture, img, mipmap);
}

void init_texture(Texture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap) {
  init_texture(texture, img.imsize(), false, as_srgb, linear, mipmap);
  update_texture(texture, img, mipmap);
}

void update_texture(Texture& texture, const image<vec4f>& img, bool mipmap) {
  check_error();
  glBindTexture(GL_TEXTURE_2D, texture.id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.imsize().x, img.imsize().y,
      GL_RGBA, GL_FLOAT, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  check_error();
}

// void update_texture_region(Texture& texture, const image<vec4f>& img,
//                            const image_region& region, bool mipmap) {
//     check_error();
//     glBindTexture(GL_TEXTURE_2D, texture.id);
//     auto clipped = image<vec4f>{};
//     get_region(clipped, img, region);
//     glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y,
//                     region.size().x, region.size().y, GL_RGBA, GL_FLOAT,
//                     clipped.data());
//     if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
//     check_error();
// }

void update_texture(Texture& texture, const image<vec4b>& img, bool mipmap) {
  check_error();
  glBindTexture(GL_TEXTURE_2D, texture.id);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.imsize().x, img.imsize().y,
      GL_RGBA, GL_UNSIGNED_BYTE, img.data());
  if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  check_error();
}

// void update_texture_region(Texture& texture, const image<vec4b>& img,
//                            const image_region& region, bool mipmap) {
//     check_error();
//     glBindTexture(GL_TEXTURE_2D, texture.id);
//     auto clipped = image<vec4b>{};
//     get_region(clipped, img, region);
//     glTexSubImage2D(GL_TEXTURE_2D, 0, region.min.x, region.min.y,
//                     region.size().x, region.size().y, GL_RGBA,
//                     GL_UNSIGNED_BYTE, clipped.data());
//     if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
//     check_error();
// }

void delete_texture(Texture& texture) {
  if (!texture) return;
  glDeleteTextures(1, &texture.id);
  texture.id   = 0;
  texture.size = zero2i;
}

void init_arraybuffer(Arraybuffer& buffer, const void* data, bool dynamic) {
  check_error();
  glGenBuffers(1, &buffer.id);
  auto flag = buffer.is_index ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
  glBindBuffer(flag, buffer.id);
  glBufferData(flag, buffer.num * buffer.elem_size, data,
      (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  check_error();
}

void delete_arraybuffer(Arraybuffer& buffer) {
  if (!buffer) return;
  glDeleteBuffers(1, &buffer.id);
  buffer.id        = 0;
  buffer.elem_size = 0;
  buffer.num       = 0;
}

void bind_shader(const Shader& shader) {
  assert(shader.shader_id);
  glUseProgram(shader.shader_id);
}
void unbind_shader() { glUseProgram(0); }

int get_uniform_location(const Shader& shader, const char* name) {
  return glGetUniformLocation(shader.shader_id, name);
}

void set_uniform(int location, int value) {
  check_error();
  glUniform1i(location, value);
  check_error();
}

void set_uniform(int location, const vec2i& value) {
  check_error();
  glUniform2i(location, value.x, value.y);
  check_error();
}

void set_uniform(int location, const vec3i& value) {
  check_error();
  glUniform3i(location, value.x, value.y, value.z);
  check_error();
}

void set_uniform(int location, const vec4i& value) {
  check_error();
  glUniform4i(location, value.x, value.y, value.z, value.w);
  check_error();
}

void set_uniform(int location, float value) {
  check_error();
  glUniform1f(location, value);
  check_error();
}

void set_uniform(int location, const vec2f& value) {
  check_error();
  glUniform2f(location, value.x, value.y);
  check_error();
}

void set_uniform(int location, const vec3f& value) {
  check_error();
  glUniform3f(location, value.x, value.y, value.z);
  check_error();
}

void set_uniform(int location, const vec4f& value) {
  check_error();
  glUniform4f(location, value.x, value.y, value.z, value.w);
  check_error();
}

void set_uniform(int location, const mat2f& value) {
  check_error();
  glUniformMatrix2fv(location, 1, false, &value.x.x);
  check_error();
}

void set_uniform(int location, const mat4f& value) {
  check_error();
  glUniformMatrix4fv(location, 1, false, &value.x.x);
  check_error();
}

void set_uniform(int location, const frame3f& value) {
  check_error();
  glUniformMatrix4x3fv(location, 1, false, &value.x.x);
  check_error();
}

void set_uniform(int location, const float* values, int num_values) {
  check_error();
  glUniform1fv(location, num_values, values);
  check_error();
}
void set_uniform(int location, const vec3f* values, int num_values) {
  check_error();
  glUniform3fv(location, num_values, &values[0].x);
  check_error();
}

void set_uniform_texture(int location, const Texture& texture, int unit) {
  check_error();
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture.id);
  glUniform1i(location, unit);
  check_error();
}

void set_uniform_texture(
    const Shader& shader, const char* name, const Texture& texture, int unit) {
  set_uniform_texture(get_uniform_location(shader, name), texture, unit);
}

void set_uniform_texture(
    int location, int locatiom_on, const Texture& texture, int unit) {
  check_error();
  if (texture.id) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture.id);
    glUniform1i(location, unit);
    glUniform1i(locatiom_on, 1);
  } else {
    glUniform1i(locatiom_on, 0);
  }
  check_error();
}

void set_uniform_texture(const Shader& shader, const char* name,
    const char* name_on, const Texture& texture, int unit) {
  set_uniform_texture(get_uniform_location(shader, name),
      get_uniform_location(shader, name_on), texture, unit);
}

int get_vertex_attribute_location(const Shader& shader, const char* name) {
  return glGetAttribLocation(shader.shader_id, name);
}

void set_vertex_attribute(int location, float value) {
  glVertexAttrib1f(location, value);
}
void set_vertex_attribute(int location, const vec2f& value) {
  glVertexAttrib2f(location, value.x, value.y);
}
void set_vertex_attribute(int location, const vec3f& value) {
  glVertexAttrib3f(location, value.x, value.y, value.z);
}
void set_vertex_attribute(int location, const vec4f& value) {
  glVertexAttrib4f(location, value.x, value.y, value.z, value.w);
}

void set_vertex_attribute(int location, const Arraybuffer& buffer) {
  check_error();
  assert(buffer.id);
  glBindBuffer(GL_ARRAY_BUFFER, buffer.id);
  glEnableVertexAttribArray(location);
  glVertexAttribPointer(
      location, buffer.elem_size / 4, GL_FLOAT, false, 0, nullptr);
  check_error();
}

void set_vertexattrib(
    const Shader& shader, const char* name, const vec3f& value) {
  auto location = get_vertex_attribute_location(shader, name);
  set_vertex_attribute(location, value);
}

void draw_point_array(const Arraybuffer& buffer) {
  glDrawArrays(GL_POINTS, 0, buffer.num);
}

void draw_line_array(const Arraybuffer& buffer) {
  glDrawArrays(GL_LINES, 0, buffer.num);
}

void draw_triangle_array(const Arraybuffer& buffer) {
  glDrawArrays(GL_TRIANGLES, 0, buffer.num);
}

void draw_point_primitives(const Arraybuffer& buffer) {
  check_error();
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.id);
  glDrawElements(GL_POINTS, buffer.num, GL_UNSIGNED_INT, nullptr);
  check_error();
}

void draw_line_primitives(const Arraybuffer& buffer) {
  check_error();
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.id);
  glDrawElements(GL_LINES, buffer.num * 2, GL_UNSIGNED_INT, nullptr);
  check_error();
}

void draw_triangle_primitives(const Arraybuffer& buffer) {
  check_error();
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer.id);
  glDrawElements(GL_TRIANGLES, buffer.num * 3, GL_UNSIGNED_INT, nullptr);
  check_error();
}

void draw_line_strip(const Arraybuffer& buffer) {
  glDrawArrays(GL_LINE_STRIP, 0, buffer.num);
}

void draw_triangle_strip(const Arraybuffer& buffer) {
  glDrawArrays(GL_TRIANGLE_STRIP, 0, buffer.num);
}

void draw_image(const Texture& texture, int win_width, int win_height,
    const vec2f& image_center, float image_scale) {
  static Shader gl_prog = {};
  //  static Arraybuffer gl_texcoord  = {};
  //  static Arraybuffer gl_triangles = {};
  static Shape quad;

  // initialization
  if (!gl_prog) {
    gl_prog.vertex_code   = R"(
            #version 330
            in vec2 texcoord;
            out vec2 frag_texcoord;
            uniform vec2 window_size, image_size;
            uniform vec2 image_center;
            uniform float image_scale;
            void main() {
                vec2 pos = (texcoord - vec2(0.5,0.5)) * image_size * image_scale + image_center;
                gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
                frag_texcoord = texcoord;
            }
        )";
    gl_prog.fragment_code = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform sampler2D txt;
            void main() {
                // frag_color = texture(txt, frag_texcoord);
                frag_color = vec4(1,0,0,1);
            }
        )";
    init_shader(gl_prog, true);

    init_shape(quad);
    add_vertex_attribute(quad, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}});
  }

  // draw
  check_error();
  bind_shader(gl_prog);
  set_uniform_texture(gl_prog, "txt", texture, 0);
  set_uniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_uniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_uniform(gl_prog, "image_center", image_center);
  set_uniform(gl_prog, "image_scale", image_scale);
  //  set_vertexattrib(gl_prog, "texcoord", gl_texcoord, zero2f);
  draw_shape(quad);
  unbind_shader();
  check_error();
}

void draw_image_background(const Texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size) {
  static Shader      gl_prog      = {};
  static Arraybuffer gl_texcoord  = {};
  static Arraybuffer gl_triangles = {};

  // initialization
  if (!gl_prog) {
    gl_prog.vertex_code   = R"(
            #version 330
            in vec2 texcoord;
            out vec2 frag_texcoord;
            uniform vec2 window_size, image_size, border_size;
            uniform vec2 image_center;
            uniform float image_scale;
            void main() {
                vec2 pos = (texcoord - vec2(0.5,0.5)) * (image_size + border_size*2) * image_scale + image_center;
                gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
                frag_texcoord = texcoord;
            }
        )";
    gl_prog.fragment_code = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform vec2 image_size, border_size;
            uniform float image_scale;
            void main() {
                ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
                ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
                ivec2 tile = tilecoord / 16;
                if(imcoord.x <= 0 || imcoord.y <= 0 || 
                    imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
                else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
                else frag_color = vec4(0.3,0.3,0.3,1);
            }
        )";
    init_shader(gl_prog, true);
    init_arraybuffer(
        gl_texcoord, vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    init_arraybuffer(gl_triangles, vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
  }

  // draw
  bind_shader(gl_prog);
  set_uniform(
      gl_prog, "window_size", vec2f{(float)win_width, (float)win_height});
  set_uniform(gl_prog, "image_size",
      vec2f{(float)texture.size.x, (float)texture.size.y});
  set_uniform(
      gl_prog, "border_size", vec2f{(float)border_size, (float)border_size});
  set_uniform(gl_prog, "image_center", image_center);
  set_uniform(gl_prog, "image_scale", image_scale);
  set_vertex_attribute(gl_prog, "texcoord", gl_texcoord);
  draw_triangle_primitives(gl_triangles);
  unbind_shader();
}

void init_shape(Shape& shape) {
  check_error();
  glGenVertexArrays(1, &shape.id);
  glBindVertexArray(shape.id);
  check_error();
}

void delete_shape(Shape& shape) {
  for (auto& attribute : shape.vertex_attributes) {
    delete_arraybuffer(attribute);
  }
  delete_arraybuffer(shape.primitives);
  glDeleteVertexArrays(1, &shape.id);
  shape = {};
}

void bind_shape(const Shape& shape) {
  check_error();
  glBindVertexArray(shape.id);
  check_error();
}

Shape make_points_shape(const vector<vec3f>& positions) {
  auto shape = Shape{};
  init_shape(shape);
  add_vertex_attribute(shape, positions);
  shape.type = Shape::type::points;
  return shape;
}

Shape make_polyline_shape(
    const vector<vec3f>& positions, const vector<vec3f>& normals) {
  auto shape = Shape{};
  init_shape(shape);
  add_vertex_attribute(shape, positions);
  if (normals.size()) {
    add_vertex_attribute(shape, normals);
  }
  shape.type = Shape::type::lines;
  return shape;
}

Shape make_quad_shape() {
  auto shape = Shape{};
  init_shape(shape);
  add_vertex_attribute(
      shape, vector<vec2f>{{-1, -1}, {1, -1}, {-1, 1}, {1, 1}});
  shape.type = Shape::type::triangles;
  return shape;
}

Shape make_regular_polygon_shape(int num_sides) {
  auto shape = Shape{};
  init_shape(shape);
  auto positions   = vector<vec2f>(num_sides + 1);
  auto triangles   = vector<vec3i>(num_sides);
  positions.back() = {0, 0};
  positions[0]     = {1.0f / (2 * sinf(pif / num_sides)), 0};
  triangles[0]     = {0, 1, num_sides};
  for (int i = 1; i < num_sides; ++i) {
    auto& tr     = triangles[i];
    tr.x         = num_sides;
    tr.y         = i;
    tr.z         = (i + 1) % num_sides;
    float angle  = (2 * pif * i) / num_sides;
    positions[i] = {yocto::cos(angle), yocto::sin(angle)};
    positions[i] /= 2 * yocto::sin(pif / num_sides);
    // yocto::tan(pif * 2 - 2 * pif / (num_sides * 2));
  }
  // auto values   = vector<float>(positions.size(), 1.0f);
  // values.back() = 0;
  add_vertex_attribute(shape, positions);
  // add_vertex_attribute(shape, values);
  init_primitives(shape, triangles);
  return shape;
}

Shape make_mesh_shape(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals) {
  auto shape = Shape{};
  init_shape(shape);
  add_vertex_attribute(shape, positions);
  add_vertex_attribute(shape, normals);
  init_primitives(shape, triangles);
  return shape;
}

Shape make_vector_field_shape(
    const vector<vec3f>& vector_field, const vector<vec3f>& from, float scale) {
  assert(vector_field.size() == from.size());
  auto shape = Shape{};
  init_shape(shape);
  auto size      = vector_field.size();
  auto positions = vector<vec3f>(size * 2);

  for (int i = 0; i < size; i++) {
    auto to              = from[i] + scale * vector_field[i];
    positions[i * 2]     = from[i];
    positions[i * 2 + 1] = to;
  }
  add_vertex_attribute(shape, positions);

  auto elements = vector<vec2i>(size);
  for (int i = 0; i < elements.size(); i++) {
    elements[i] = {2 * i, 2 * i + 1};
  }
  init_primitives(shape, elements);
  return shape;
}

#define SURFACE_OFFSET 0.00002f

Shape make_vector_field_shape(const vector<vec3f>& vector_field,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    float scale) {
  assert(vector_field.size() == triangles.size());

  auto froms               = vector<vec3f>(vector_field.size());
  auto vector_field_lifted = vector_field;

  for (int i = 0; i < triangles.size(); i++) {
    auto x      = positions[triangles[i].x];
    auto y      = positions[triangles[i].y];
    auto z      = positions[triangles[i].z];
    auto normal = triangle_normal(x, y, z);
    auto center = (x + y + z) / 3;
    froms[i]    = center + SURFACE_OFFSET * normal;
    vector_field_lifted[i] += SURFACE_OFFSET * normal;
  }
  return make_vector_field_shape(vector_field, froms, scale);
}

void draw_shape(const Shape& shape) {
  // @SPEED: This is for extra-safety, but may have
  //         an impact with many draw calls
  if (!shape.id) return;
  if (shape.vertex_attributes.empty()) return;

  bind_shape(shape);
  auto& primitives = shape.primitives;

  if (primitives) {
    if (shape.type == Shape::type::points) {
      draw_point_primitives(primitives);
    } else if (shape.type == Shape::type::lines) {
      draw_line_primitives(primitives);
    } else if (shape.type == Shape::type::triangles) {
      draw_triangle_primitives(primitives);
    }
  } else {
    auto& positions = shape.vertex_attributes[0];
    if (shape.is_strip) {
      if (!positions) return;
      if (shape.type == Shape::type::points) draw_point_array(positions);
      if (shape.type == Shape::type::lines) draw_line_strip(positions);
      if (shape.type == Shape::type::triangles) draw_triangle_strip(positions);
    } else {
      if (shape.type == Shape::type::points) draw_point_array(positions);
      if (shape.type == Shape::type::lines) draw_line_array(positions);
      if (shape.type == Shape::type::triangles) draw_triangle_array(positions);
    }
    check_error();
  }
}

Camera make_lookat_camera(const vec3f& from, const vec3f& to, const vec3f& up) {
  auto camera  = Camera{};
  camera.frame = lookat_frame(from, to, {0, 1, 0});
  camera.focus = length(from - to);
  return camera;
}

mat4f make_view_matrix(const Camera& camera) {
  return frame_to_mat(inverse(camera.frame));
}

mat4f make_projection_matrix(const Camera& camera, const vec2i& viewport) {
  auto camera_aspect = (float)viewport.x / (float)viewport.y;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(camera.film.x / (camera_aspect * 2 * camera.lens)))
          : (2 * yocto::atan(camera.film.x / (2 * camera.lens)));
  return perspective_mat(camera_yfov, camera_aspect, camera.near, camera.far);
}

Rendertarget make_render_target(
    const vec2i& size, bool as_float, bool as_srgb, bool linear, bool mipmap) {
  auto target = Rendertarget{};
  glGenFramebuffers(1, &target.frame_buffer);
  glBindFramebuffer(GL_FRAMEBUFFER, target.frame_buffer);

  // create a color attachment texture
  auto& texture       = target.texture.id;
  target.texture.size = size;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, size.x, size.y, 0, GL_RGB,
      GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glFramebufferTexture2D(
      GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

  // create render buffer for depth and stencil
  glGenRenderbuffers(1, &target.render_buffer);
  glBindRenderbuffer(GL_RENDERBUFFER, target.render_buffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, size.x, size.y);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  // bind frame buffer and render buffer
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT,
      GL_RENDERBUFFER, target.render_buffer);
  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // done
  return target;
}

void bind_render_target(const Rendertarget& target) {
  glBindFramebuffer(GL_FRAMEBUFFER, target.frame_buffer);
}

void unbind_render_target() { glBindFramebuffer(GL_FRAMEBUFFER, 0); }

image<vec4f> capture_screenshot() {
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  int x      = viewport[0];
  int y      = viewport[1];
  int width  = viewport[2];
  int height = viewport[3];

  auto img = image<vec4f>({width, height});
  // char *data = (char*) malloc((size_t) (width * height * 3)); // 3 components
  // (R, G, B)

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(x, y, width, height, GL_RGBA, GL_FLOAT, img.data());

  return img;
}

#if 0
// init image shader
void init_image_shader(Image& glimage) {
  if (glimage.shader) return;
  auto vert =
      R"(
      #version 330
      in vec2 texcoord;
      out vec2 frag_texcoord;
      uniform vec2 window_size, image_size;
      uniform vec2 image_center;
      uniform float image_scale;
      void main() {
          vec2 pos = (texcoord - vec2(0.5,0.5)) * image_size * image_scale + image_center;
          gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
          frag_texcoord = texcoord;
      }
      )";
#if 0
  auto vert = R"(
          #version 330
          in vec2 texcoord;
          out vec2 frag_texcoord;
          uniform vec2 window_size, image_size, border_size;
          uniform vec2 image_center;
          uniform float image_scale;
          void main() {
              vec2 pos = (texcoord - vec2(0.5,0.5)) * (image_size + border_size*2) * image_scale + image_center;
              gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
              frag_texcoord = texcoord;
          }
      )";
#endif
  auto frag =
      R"(
      #version 330
      in vec2 frag_texcoord;
      out vec4 frag_color;
      uniform sampler2D txt;
      void main() {
          frag_color = texture(txt, frag_texcoord);
      }
      )";
#if 0
    auto frag = R"(
            #version 330
            in vec2 frag_texcoord;
            out vec4 frag_color;
            uniform vec2 image_size, border_size;
            uniform float image_scale;
            void main() {
                ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
                ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
                ivec2 tile = tilecoord / 16;
                if(imcoord.x <= 0 || imcoord.y <= 0 || 
                    imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
                else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
                else frag_color = vec4(0.3,0.3,0.3,1);
            }
        )";
#endif

  init_shader(glimage.shader, vert, frag);
  init_shape(glimage.shape);
  add_vertex_attribute(
      glimage.shape, vector<vec2f>{{1, 0}, {1, 1}, {0, 0}, {0, 1}});
  glimage.shape.type = Shape::type::triangles;
}

// update image data
void update_image(
    Image& glimage, const image<vec4f>& img, bool linear, bool mipmap) {
  init_image_shader(glimage);
  if (!glimage.texture) {
    init_texture(glimage.texture, img, false, linear, mipmap);
  } else if (glimage.texture.size != img.size() ||
             glimage.texture.linear != linear ||
             glimage.texture.mipmap != mipmap) {
    delete_texture(glimage.texture);
    init_texture(glimage.texture, img, false, linear, mipmap);
  } else {
    update_texture(glimage.texture, img, mipmap);
  }
}
void update_image(
    Image& glimage, const image<vec4b>& img, bool linear, bool mipmap) {
  init_image_shader(glimage);
  if (!glimage.texture) {
    init_texture(glimage.texture, img, false, linear, mipmap);
  } else if (glimage.texture.size != img.size() ||
             glimage.texture.linear != linear ||
             glimage.texture.mipmap != mipmap) {
    delete_texture(glimage.texture);
    init_texture(glimage.texture, img, false, linear, mipmap);
  } else {
    update_texture(glimage.texture, img, mipmap);
  }
}

void update_image_region(
    Image& glimage, const image<vec4f>& img, const image_region& region) {
  if (!glimage.texture) throw std::runtime_error("glimage is not initialized");
  update_texture_region(glimage.texture, img, region, glimage.texture.mipmap);
}
void update_image_region(
    Image& glimage, const image<vec4b>& img, const image_region& region) {
  if (!glimage.texture) throw std::runtime_error("glimage is not initialized");
  update_texture_region(glimage.texture, img, region, glimage.texture.mipmap);
}

// draw image
void draw_image(Image& glimage, const draw_image_params& params) {
  check_error();
  set_viewport(params.framebuffer);
  clear_framebuffer(params.background);
  bind_shader(glimage.shader);
  set_uniform_texture(glimage.shader, "txt", glimage.texture, 0);
  set_uniform(glimage.shader, "window_size",
      vec2f{(float)params.window.x, (float)params.window.y});
  set_uniform(glimage.shader, "image_size",
      vec2f{(float)glimage.texture.size.x, (float)glimage.texture.size.y});
  set_uniform(glimage.shader, "image_center", params.center);
  set_uniform(glimage.shader, "image_scale", params.scale);
  draw_shape(glimage.shape);
  unbind_shader();
  check_error();
}

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL FUNCTIONS
// -----------------------------------------------------------------------------

// Initialize an OpenGL scene
void make_scene(Scene& glscene) {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

  static const char* vertex =
      R"(
        #version 330

        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;                    // shape transform
        uniform mat4 shape_xform_invtranspose;       // shape transform
        uniform float shape_normal_offset;           // shape normal offset

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        out vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
        out vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
        out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
        out vec4 color;                 // [to fragment shader] vertex color
        out vec4 tangsp;                // [to fragment shader] vertex tangent space

        // main function
        void main() {
            // copy values
            pos = vert_pos;
            norm = vert_norm;
            tangsp = vert_tangsp;

            // normal offset
            if(shape_normal_offset != 0) {
                pos += shape_normal_offset * norm;
            }

            // world projection
            pos = (shape_xform * vec4(pos,1)).xyz;
            norm = (shape_xform_invtranspose * vec4(norm,0)).xyz;
            tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
        }
        )";

  static const char* fragment =
      R"(
        #version 330

        float pif = 3.14159265;

        uniform bool eyelight;         // eyelight shading
        uniform vec3 lamb;             // ambient light
        uniform int lnum;              // number of lights
        uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
        uniform vec3 lpos[16];         // light positions
        uniform vec3 lke[16];          // light intensities

        void evaluate_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(ltype[lid] == 0) {
                // compute point light color at pos
                cl = lke[lid] / pow(length(lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lpos[lid]-pos);
            }
            else if(ltype[lid] == 1) {
                // compute light color
                cl = lke[lid];
                // compute light direction
                wi = normalize(lpos[lid]);
            }
        }

        vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
            vec3 n, vec3 wi, vec3 wo) {
            if(etype == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(rs*rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(etype == 1) {
                return ((1+dot(wo,wi))/2) * kd/pif;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pif;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pif)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pif;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pif)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        uniform int elem_type;
        uniform bool elem_faceted;
        uniform vec4 highlight;   // highlighted color

        uniform int mat_type;          // material type
        uniform vec3 mat_ke;           // material ke
        uniform vec3 mat_kd;           // material kd
        uniform vec3 mat_ks;           // material ks
        uniform float mat_rs;          // material rs
        uniform float mat_op;          // material op

        uniform bool mat_ke_txt_on;    // material ke texture on
        uniform sampler2D mat_ke_txt;  // material ke texture
        uniform bool mat_kd_txt_on;    // material kd texture on
        uniform sampler2D mat_kd_txt;  // material kd texture
        uniform bool mat_ks_txt_on;    // material ks texture on
        uniform sampler2D mat_ks_txt;  // material ks texture
        uniform bool mat_rs_txt_on;    // material rs texture on
        uniform sampler2D mat_rs_txt;  // material rs texture
        uniform bool mat_op_txt_on;    // material op texture on
        uniform sampler2D mat_op_txt;  // material op texture

        uniform bool mat_norm_txt_on;    // material norm texture on
        uniform sampler2D mat_norm_txt;  // material norm texture

        uniform bool mat_double_sided;   // double sided rendering

        uniform mat4 shape_xform;              // shape transform
        uniform mat4 shape_xform_invtranspose; // shape transform

        bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                           out vec3 kd, out vec3 ks, out float rs, out float op) {
            if(mat_type == 0) {
                ke = mat_ke;
                kd = vec3(0,0,0);
                ks = vec3(0,0,0);
                op = 1;
                return false;
            }

            ke = color.xyz * mat_ke;
            kd = color.xyz * mat_kd;
            ks = color.xyz * mat_ks;
            rs = mat_rs;
            op = color.w * mat_op;

            vec4 ke_txt = (mat_ke_txt_on) ? texture(mat_ke_txt,texcoord) : vec4(1,1,1,1);
            vec4 kd_txt = (mat_kd_txt_on) ? texture(mat_kd_txt,texcoord) : vec4(1,1,1,1);
            vec4 ks_txt = (mat_ks_txt_on) ? texture(mat_ks_txt,texcoord) : vec4(1,1,1,1);
            vec4 rs_txt = (mat_rs_txt_on) ? texture(mat_rs_txt,texcoord) : vec4(1,1,1,1);
            vec4 op_txt = (mat_op_txt_on) ? texture(mat_op_txt,texcoord) : vec4(1,1,1,1);

            // get material color from textures and adjust values
            if(mat_type == 1) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt.y;
                rs = rs*rs;
                op *= op_txt.x * kd_txt.w;
            } else if(mat_type == 2) {
                ke *= ke_txt.xyz;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(mat_type == 3) {
                ke *= ke_txt.xyz;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                float gs = (1 - rs) * ks_txt.w;
                rs = 1 - gs;
                rs = rs*rs;
                op *= kd_txt.w;
            }

            return true;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!mat_norm_txt_on) return norm;
            vec3 tangu = normalize((shape_xform * vec4(normalize(tangsp.xyz),0)).xyz);
            vec3 tangv = normalize(cross(norm, tangu));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 texture = 2 * texture(mat_norm_txt,texcoord).xyz - 1;
            texture.y = -texture.y;
            return normalize( tangu * texture.x + tangv * texture.y + norm * texture.z );
        }

        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        uniform vec3 cam_pos;          // camera position
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        uniform float exposure; 
        uniform float gamma;

        out vec4 frag_color;      

        vec3 triangle_normal(vec3 pos) {
            vec3 fdx = dFdx(pos); 
            vec3 fdy = dFdy(pos); 
            return normalize((shape_xform * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
        }

        // main
        void main() {
            // view vector
            vec3 wo = normalize(cam_pos - pos);

            // prepare normals
            vec3 n;
            if(elem_faceted) {
                n = triangle_normal(pos);
            } else {
                n = normalize(norm);
            }

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(mat_double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
            bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

            // exit if needed
            if(brdf_op < 0.005) discard;

            // check const color
            if(elem_type == 0) {
                frag_color = vec4(brdf_ke,brdf_op);
                return;
            }

            // emission
            vec3 c = brdf_ke;

            // check early exit
            if(brdf_kd != vec3(0,0,0) || brdf_ks != vec3(0,0,0)) {
                // eyelight shading
                if(eyelight) {
                    vec3 wi = wo;
                    c += pif * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lamb * brdf_kd;
                    // foreach light
                    for(int lid = 0; lid < lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        evaluate_light(lid, pos, cl, wi);
                        c += cl * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                    }
                }
            }

            // final color correction
            c = pow(c * pow(2,exposure), vec3(1/gamma));

            // highlighting
            if(highlight.w > 0) {
                if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
                    c = highlight.xyz * highlight.w + c * (1-highlight.w);
            }

            // output final color by setting gl_FragColor
            frag_color = vec4(c,brdf_op);
        }
        )";
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

  // load shader
  init_shader(glscene.shader, vertex, fragment);
}

// Draw a shape
void draw_instance(
    Scene& state, const Instance& instance, const draw_scene_params& params) {
  if (instance.shape < 0 || instance.shape > state.shapes.size()) return;
  if (instance.material < 0 || instance.material > state.materials.size())
    return;

  auto& shape    = state.shapes[instance.shape];
  auto& material = state.materials[instance.material];

  if (!shape.id) return;

  set_uniform(state.shader, "shape_xform", mat4f(instance.frame));
  set_uniform(state.shader, "shape_xform_invtranspose",
      transpose(mat4f(inverse(instance.frame, params.non_rigid_frames))));
  set_uniform(state.shader, "shape_normal_offset", 0.0f);
  set_uniform(state.shader, "highlight",
      instance.highlighted ? vec4f{1, 1, 0, 1} : zero4f);

  auto mtype = 2;
  if (material.gltf_textures) mtype = 3;
  set_uniform(state.shader, "mat_type", mtype);
  set_uniform(state.shader, "mat_ke", material.emission);
  set_uniform(state.shader, "mat_kd", material.diffuse);
  set_uniform(state.shader, "mat_ks", vec3f{material.metallic});
  set_uniform(state.shader, "mat_rs", material.roughness);
  set_uniform(state.shader, "mat_op", material.opacity);
  //set_uniform(state.shader, "mat_double_sided", (int)params.double_sided);
  set_uniform(state.shader, "mat_double_sided", 1);
  if (material.emission_map >= 0) {
    set_uniform_texture(state.shader, "mat_ke_txt", "mat_ke_txt_on",
        state.textures.at(material.emission_map), 0);
  } else {
    set_uniform_texture(
        state.shader, "mat_ke_txt", "mat_ke_txt_on", Texture{}, 0);
  }
  if (material.diffuse_map >= 0) {
    set_uniform_texture(state.shader, "mat_kd_txt", "mat_kd_txt_on",
        state.textures.at(material.diffuse_map), 1);
  } else {
    set_uniform_texture(
        state.shader, "mat_kd_txt", "mat_kd_txt_on", Texture{}, 1);
  }
  if (material.metallic_map >= 0) {
    set_uniform_texture(state.shader, "mat_ks_txt", "mat_ks_txt_on",
        state.textures.at(material.metallic_map), 2);
  } else {
    set_uniform_texture(
        state.shader, "mat_ks_txt", "mat_ks_txt_on", Texture{}, 2);
  }
  if (material.roughness_map >= 0) {
    set_uniform_texture(state.shader, "mat_rs_txt", "mat_rs_txt_on",
        state.textures.at(material.roughness_map), 3);
  } else {
    set_uniform_texture(
        state.shader, "mat_rs_txt", "mat_rs_txt_on", Texture{}, 3);
  }
  if (material.normal_map >= 0) {
    set_uniform_texture(state.shader, "mat_norm_txt", "mat_norm_txt_on",
        state.textures.at(material.normal_map), 5);
  } else {
    set_uniform_texture(
        state.shader, "mat_norm_txt", "mat_norm_txt_on", Texture{}, 5);
  }

  draw_shape(shape);

#if 0
    if ((vbos.gl_edges && edges && !wireframe) || highlighted) {
        enable_culling(false);
        check_error();
        set_uniform(state.shader, "mtype"), 0);
        glUniform3f(glGetUniformLocation(state.shader, "ke"), 0, 0, 0);
        set_uniform(state.shader, "op"), material.op);
        set_uniform(state.shader, "shp_normal_offset"), 0.01f);
        check_error();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.gl_edges);
        glDrawElements(GL_LINES, vbos.triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        check_error();
    }
#endif
  if (params.edges) throw std::runtime_error("edges are momentarily disabled");
}

// Display a scene
void draw_scene(
    Scene& state, const vec4i& viewport, const draw_scene_params& params) {
  auto& glcamera      = state.cameras.at(params.camera);
  auto  camera_aspect = (float)viewport.z / (float)viewport.w;
  auto  camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(
                     glcamera.film / (camera_aspect * 2 * glcamera.lens)))
          : (2 * yocto::atan(glcamera.film / (2 * glcamera.lens)));
  auto camera_view = mat4f(inverse(glcamera.frame));
  auto camera_proj = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  clear_framebuffer(params.background);
  set_viewport(viewport);

  bind_shader(state.shader);
  set_uniform(state.shader, "cam_pos", glcamera.frame.o);
  set_uniform(state.shader, "cam_xform_inv", camera_view);
  set_uniform(state.shader, "cam_proj", camera_proj);
  set_uniform(state.shader, "eyelight", (int)params.eyelight);
  set_uniform(state.shader, "exposure", params.exposure);
  set_uniform(state.shader, "gamma", params.gamma);

  if (!params.eyelight) {
    set_uniform(state.shader, "lamb", zero3f);
    set_uniform(state.shader, "lnum", (int)state.lights.size());
    for (auto i = 0; i < state.lights.size(); i++) {
      auto is = std::to_string(i);
      set_uniform(
          state.shader, ("lpos[" + is + "]").c_str(), state.lights[i].position);
      set_uniform(
          state.shader, ("lke[" + is + "]").c_str(), state.lights[i].emission);
      set_uniform(state.shader, ("ltype[" + is + "]").c_str(),
          (int)state.lights[i].type);
    }
  }

  if (params.wireframe) set_wireframe(true);
  check_error();

  for (auto& instance : state.instances) {
    draw_instance(state, instance, params);
  }

  unbind_shader();
  if (params.wireframe) set_wireframe(false);
}
#endif
}  // namespace gpu
