#ifndef _REALTIME_GPU_
#define _REALTIME_GPU_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

// #include <graphics/image.h>
// #include <graphics/math.h>
#include <yocto/yocto_image.h>

#include <cassert>
#include <functional>   // std::function
#include <type_traits>  // std::is_floating_point

namespace gpu {
using namespace yocto;

// Commands to setup the opengl context and issue gpu operations.
bool init_opengl();
void check_error();
void clear_framebuffer(const vec4f& color, bool clear_depth = true);
void set_viewport(const vec4i& viewport);
void set_wireframe(bool enabled);
void set_blending(bool enabled);
void set_point_size(int size);

enum struct DepthTest {
  never,     // Never passes.
  less,      // Passes if the incoming depth is less than the stored depth.
  equal,     // Passes if the incoming depth is equal to the stored depth.
  lequal,    // Passes if the incoming depth is less than or equal to the stored
             // depth.
  greater,   // Passes if the incoming depth is greater than the stored depth.
  notequal,  // Passes if the incoming depth is not equal to the stored depth.
  gequal,    // Passes if the incoming depth is greater than or equal to the
             // stored depth.
  always     // Always passes.
};
void set_depth_test(DepthTest flag);

/* IMPORTANT: ************************************************************

  The following data structures help keep track of the resources stored in
  the gpu. We don't employ C++ RAII, so RESOURCES MUST BE DELETED MANUALLY
  by callind delete_<type>().

***************************************************************************/

// A Texture is an image stored on the gpu. We store the id assigned to the
// texture upon creation so that it can be bound and keep track of how it must
// be sampled inside the shaders using it.
struct Texture {
  uint  id       = 0;
  vec2i size     = {0, 0};
  bool  mipmap   = false;
  bool  linear   = false;
  bool  is_srgb  = false;
  bool  is_float = false;
        operator uint() const { return id; }
};

void init_texture(Texture& texture, const vec2i& size, bool as_float,
    bool as_srgb, bool linear, bool mipmap);

void init_texture(Texture& texture, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap);
void init_texture(Texture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap);

void update_texture(Texture& texture, const image<vec4f>& img, bool mipmap);
void update_texture(Texture& texture, const image<vec4b>& img, bool mipmap);

// void update_texture_region(Texture& texture, const image<vec4f>& img,
//                            const image_region& region, bool mipmap);
// void update_texture_region(Texture& texture, const image<vec4b>& img,
//                            const image_region& region, bool mipmap);

void delete_texture(Texture& texture);

// An Arraybuffer is an array stored on the gpu. Here we keep track of its id,
// the number of elements, the size of each element and whether it contains
// floating point values or integers (to represent indices).
struct Arraybuffer {
  uint id        = 0;
  int  num       = 0;
  int  elem_size = 0;
  bool is_index  = false;
       operator bool() const { return (bool)id; }
};
void init_arraybuffer(Arraybuffer& buffer, const void* data, bool dynamic);
void init_arraybuffer(
    Arraybuffer& buffer, const vector<float>& array, bool dynamic);

template <typename T>
void init_arraybuffer(
    Arraybuffer& buffer, const vector<T>& array, bool dynamic = false) {
  buffer           = Arraybuffer{};
  buffer.num       = size(array);
  buffer.elem_size = sizeof(T);
  if constexpr (std::is_floating_point<T>::value) {
    buffer.is_index = false;
  } else {
    buffer.is_index = !std::is_floating_point<typeof(array[0][0])>::value;
  }
  init_arraybuffer(buffer, array.data(), dynamic);
}
void delete_arraybuffer(Arraybuffer& buffer);

// A Shader is a program running on the gpu. We store the source code
// of the shader, the path of the file containing the shader source for
// hot-reloading and the ids assigned to the shader once compiled so that it can
// be bound later.
struct Shader {
  string vertex_code;
  string fragment_code;
  string vertex_filename;
  string fragment_filename;
  uint   shader_id   = 0;
  uint   vertex_id   = 0;
  uint   fragment_id = 0;
         operator uint() const { return shader_id; }
};

bool init_shader(Shader& shader, bool abort_on_error = false);
void delete_shader(Shader& shader);

void   load_shader_code(Shader& shader);
Shader make_shader_from_file(const string& vertex_filename,
    const string& fragment_filename, bool abort_on_error = true);

void bind_shader(const Shader& shader);
void unbind_shader();

// A Uniform is a parameter used by a shader that can be updated from the cpu. A
// uniform value is identified by the name to which is bound in the shader.
template <typename Type>
struct Uniform {
  const char* name;
  Type        value;
  Uniform(const char* n, const Type& v) : name(n), value(v) {}
};
int get_uniform_location(const Shader& shader, const char* name);

void set_uniform(int location, int value);
void set_uniform(int location, const vec2i& value);
void set_uniform(int location, const vec3i& value);
void set_uniform(int location, const vec4i& value);
void set_uniform(int location, float value);
void set_uniform(int location, const vec2f& value);
void set_uniform(int location, const vec3f& value);
void set_uniform(int location, const vec4f& value);
void set_uniform(int location, const mat2f& value);
void set_uniform(int location, const mat3f& value);
void set_uniform(int location, const mat4f& value);
void set_uniform(int location, const frame3f& value);

inline void set_uniform(const Shader& shader) {}

template <typename T>
inline void set_uniform(
    const Shader& shader, const char* name, const T& value) {
  set_uniform(get_uniform_location(shader, name), value);
}

void set_uniform(int location, const float* value, int num_values);
void set_uniform(int location, const vec3f* value, int num_values);

template <typename T>
inline void set_uniform(
    const Shader& shader, const char* name, const T* values, int num_values) {
  set_uniform(get_uniform_location(shader, name), values, num_values);
}

void set_uniform_texture(int location, const Texture& texture, int unit);
void set_uniform_texture(
    const Shader& shader, const char* name, const Texture& texture, int unit);
void set_uniform_texture(
    int location, int location_on, const Texture& texture, int unit);
void set_uniform_texture(const Shader& shader, const char* name,
    const char* name_on, const Texture& texture, int unit);

int  get_vertex_attribute_location(const Shader& shader, const char* name);
void set_vertex_attribute(int location, const Arraybuffer& buffer);
void set_vertex_attribute(int location, float value);
void set_vertex_attribute(int location, const vec2f& value);
void set_vertex_attribute(int location, const vec3f& value);
void set_vertex_attribute(int location, const vec4f& value);

template <typename Type>
void set_vertex_attribute(
    const Shader& shader, const char* name, const Type& buffer) {
  auto location = get_vertex_attribute_location(shader, name);
  set_vertex_attribute(location, buffer);
}

template <typename Type>
void set_uniform(const Shader& shader, const Uniform<Type>& u) {
  set_uniform(shader, u.name, u.value);
}

template <typename Type, typename... Args>
void set_uniform(const Shader& shader, const Uniform<Type>& u,
    const Uniform<Args>&... args) {
  set_uniform(shader, u);
  set_uniform(shader, args...);
}

// A Shape is an entity stored on the gpu that can be drawn on the screen. A
// Shape is always a discrete collection of primitives of the same type, where
// supported types are points, lines and triangles. It holds a list of vertex
// attributs, which are arrays in the gpu all the same size, storing
// information about the vertices of the primitves. We probably use the
// convention that the first vertex attribute is the position of the
// vertices. The "elements" array must contain indicse or be empty.
// If not empty, the elemente array tells the gpu how to group the vertices to
// create the primivies. If empty, the shape is drawn as a "strip".
struct Shape {
  vector<Arraybuffer> vertex_attributes = {};
  Arraybuffer         primitives        = {};
  uint                id                = 0;
  bool                is_strip          = true;

  enum struct type { points, lines, triangles };
  type type = type::triangles;
};

void init_shape(Shape& shape);
void delete_shape(Shape& shape);
void bind_shape(const Shape& shape);
void draw_shape(const Shape& shape);

void draw_point_array(const Arraybuffer& buffer);
void draw_line_array(const Arraybuffer& buffer);
void draw_triangle_array(const Arraybuffer& buffer);
void draw_point_primitives(const Arraybuffer& buffer);
void draw_line_primitives(const Arraybuffer& buffer);
void draw_triangle_primitives(const Arraybuffer& buffer);
void draw_line_strip(const Arraybuffer& buffer);
void draw_triangle_strip(const Arraybuffer& buffer);

template <typename Type>
void draw_points(const vector<Type>& points) {
  static auto array = Arraybuffer{};
  init_arraybuffer(array, points, false);
  draw_point_array(array);
}

template <typename T>
int set_vertex_attribute(Shape& shape, int index, const vector<T>& data) {
  assert(index < shape.vertex_attributes.size());
  // check all attributes have same size
  bind_shape(shape);
  // @Speed: update instead of delete
  delete_arraybuffer(shape.vertex_attributes[index]);
  init_arraybuffer(shape.vertex_attributes[index], data);
  set_vertex_attribute(index, shape.vertex_attributes[index]);
  return index;
}

template <typename T>
int set_vertex_attribute(Shape& shape, int index, const T& data) {
  assert(index < shape.vertex_attributes.size());
  // check all attributes have same size
  bind_shape(shape);
  // @Speed: update instead of delete
  //  delete_arraybuffer(shape.vertex_attributes[index]);
  //  init_arraybuffer(shape.vertex_attributes[index], data);
  // int elem_size = sizeof(T) / sizeof(float);
  set_vertex_attribute(index, data);

  return index;
}

template <typename T>
int add_vertex_attribute(Shape& shape, const vector<T>& data) {
  assert(shape.vertex_attributes.empty() ||
         shape.vertex_attributes[0].num == data.size());
  bind_shape(shape);
  int   index     = shape.vertex_attributes.size();
  auto& attribute = shape.vertex_attributes.emplace_back();
  init_arraybuffer(attribute, data);
  set_vertex_attribute(index, attribute);
  return index;
}

template <typename T>
int add_vertex_attribute(Shape& shape, const T& data) {
  bind_shape(shape);
  int index = shape.vertex_attributes.size();
  shape.vertex_attributes.push_back({});
  //    init_arraybuffer(attribute, data);
  set_vertex_attribute(index, data);
  return index;
}

template <typename T>
void init_primitives(Shape& shape, const vector<T>& data) {
  check_error();
  bind_shape(shape);
  init_arraybuffer(shape.primitives, data);
  int elem_size = sizeof(T) / sizeof(int);
  if (elem_size == 1) shape.type = Shape::type::points;
  if (elem_size == 2) shape.type = Shape::type::lines;
  if (elem_size == 3) shape.type = Shape::type::triangles;
  check_error();
}

Shape make_points_shape(const vector<vec3f>& positions);
Shape make_polyline_shape(
    const vector<vec3f>& position, const vector<vec3f>& normals = {});
Shape make_quad_shape();
Shape make_regular_polygon_shape(int num_sides);
Shape make_mesh_shape(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals);
Shape make_vector_field_shape(const vector<vec3f>& vector_field,
    const vector<vec3f>& from, float scale = 0.01);
Shape make_vector_field_shape(const vector<vec3f>& vector_field,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    float scale = 0.001);

template <typename... Args>
void draw_shape(
    const Shape& shape, const Shader& shader, const Uniform<Args>&... args) {
  bind_shader(shader);
  set_uniform(shader, args...);
  draw_shape(shape);
}

// A render target is a texture in the gpu onto which an image can be rendered
// with the provided draw calls.
struct Rendertarget {
  Texture texture;
  uint    frame_buffer;
  uint    render_buffer;
};

Rendertarget make_render_target(
    const vec2i& size, bool as_float, bool as_srgb, bool linear, bool mipmap);

void         bind_render_target(const Rendertarget& target);
void         unbind_render_target();
image<vec4f> capture_screenshot();

// A Camera holds the information required to compute the view-projection matrix
// needed for 3d rendering.
struct Camera {
  frame3f frame = identity3x4f;
  float   lens  = 0.050;
  vec2f   film  = {0.036, 0.024};
  float   near  = 0.001;
  float   far   = 10000;
  float   focus = flt_max;
};

Camera make_lookat_camera(
    const vec3f& from, const vec3f& to, const vec3f& up = {0, 1, 0});
mat4f make_view_matrix(const Camera& camera);
mat4f make_projection_matrix(const Camera& camera, const vec2i& viewport);

/*
struct Image {
  Texture texture = {};
  Shader  shader  = {};
  Shape   shape   = {};

  vec2i size() const { return texture.size; }
};
struct draw_image_params {
  vec2i window      = {512, 512};
  vec4i framebuffer = {0, 0, 512, 512};
  vec2f center      = {0, 0};
  float scale       = 1;
  bool  fit         = true;
  bool  checker     = true;
  float border_size = 2;
  vec4f background  = {0.15f, 0.15f, 0.15f, 1.0f};
};
void update_image(Image& glimage, const image<vec4f>& img, bool linear = false,
    bool mipmap = false);
void update_image(Image& glimage, const image<vec4b>& img, bool linear = false,
    bool mipmap = false);
void update_image_region(
    Image& glimage, const image<vec4f>& img, const image_region& region);
void update_image_region(
    Image& glimage, const image<vec4b>& img, const image_region& region);
void draw_image(Image& glimage, const draw_image_params& params);
*/

// TODO: We don't enforce any shading model here since this code is for just for
// drawing pixels. Maybe we can keep a simple shading model associated with
// provided shader functions.
// struct Material {
//   vec3f emission      = zero3f;
//   vec3f diffuse       = zero3f;
//   vec3f specular      = zero3f;
//   float metallic      = 0;
//   float roughness     = 0;
//   float opacity       = 1;
//   int   emission_map  = -1;
//   int   diffuse_map   = -1;
//   int   specular_map  = -1;
//   int   metallic_map  = -1;
//   int   roughness_map = -1;
//   int   normal_map    = -1;
//   bool  gltf_textures = false;
// };

// No assumption on how a scene is represented
/*
// Opengl instance group
struct Instance {
  frame3f frame       = identity3x4f;
  int     shape       = 0;
  int     material    = 0;
  bool    highlighted = false;
};

// Opengl light
struct Light {
  vec3f position = zero3f;
  vec3f emission = zero3f;
  int   type     = 0;
};

// Opengl scene
struct Scene {
  vector<Camera>   cameras   = {};
  vector<Instance> instances = {};
  vector<Shape>    shapes    = {};
  vector<Material> materials = {};
  vector<Texture>  textures  = {};
  vector<Light>    lights    = {};
  Shader           shader    = {};
};

// Draw options
struct draw_scene_params {
  int   camera           = 0;
  int   resolution       = 1280;
  bool  wireframe        = false;
  bool  edges            = false;
  float edge_offset      = 0.01f;
  bool  eyelight         = false;
  float exposure         = 0;
  float gamma            = 2.2f;
  vec3f ambient          = {0, 0, 0};
  bool  double_sided     = true;
  bool  non_rigid_frames = true;
  float near             = 0.01f;
  float far              = 10000.0f;
  vec4f background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void make_scene(Scene& scene);

// Draw an OpenGL scene
void draw_scene(
    Scene& state, const vec4i& viewport, const draw_scene_params& params);
*/

// template <int N>
// Render_target make_render_targets(array<N, Texture>&
// textures);

}  // namespace gpu

#endif
