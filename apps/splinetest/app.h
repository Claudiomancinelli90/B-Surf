#pragma once
#include <realtime/gpu.h>
#include <splinesurf/spline.h>
#include <splinesurf/splineio.h>
#include <yocto/yocto_bvh.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>
#include <yocto_gui/yocto_window.h>

#include <unordered_set>

#include "flipout/flipout.h"
#include "vector_heat/vector_heat.h"
using namespace yocto;

#define PROFILE 0

struct Gui_Line {
  geodesic_path path      = {};
  vector<vec3f> positions = {};
  float         length    = 0;
  ogl_shape*    shape     = new ogl_shape{};
};

void set_line_shape(Gui_Line& line, const bezier_mesh& mesh);

struct Gui_Curve {
  bezier_tree             tree;
  shade_shape*            shape         = nullptr;
  vector<vec3f>           positions     = {};
  vector<mesh_point>      samples       = {};
  vector<float>           sampling_rate = {};
  vector<vec2i>           offsets       = {};
  std::array<Gui_Line, 2> tangents      = {};
};

inline std::unordered_set<int> make_set(size_t n) {
  auto result = std::unordered_set<int>();
  result.reserve(n);
  for (int i = 0; i < n; i++) {
    result.insert(i);
  }
  return result;
}

struct Gui_Spline {
  vector<mesh_point> control_points = {};
  vector<bool>       is_smooth      = {};
  vector<Gui_Curve>  curves         = {};

  vector<vec2f>           local_coords            = {};
  mat2f                   frame                   = {};
  mesh_point              center                  = {};
  std::unordered_set<int> curves_to_update        = {};
  ogl_shape*              anchor_points_shape     = new ogl_shape{};
  int                     anchor_points_shape_new = 0;
  vec3f                   color                   = {1, 0, 0};
};

inline void clear_spline(Gui_Spline& spline) {
  for (auto& curve : spline.curves) {
    clear_shape(curve.shape);
    clear_shape(curve.tangents[0].shape);
    clear_shape(curve.tangents[1].shape);
    curve.shape             = nullptr;
    curve.tangents[0].shape = nullptr;
    curve.tangents[1].shape = nullptr;
  }
  clear_shape(spline.anchor_points_shape);
  spline.control_points.clear();
  spline.curves.clear();
}

void update_control_points_shape(Gui_Spline& spline, const bezier_mesh& mesh);
vector<vector<vec3f>> split_broken_curve(
    const vector<vec3f>& positions, const float& threshold);
void update_curve_shape(Gui_Spline& spline, const int curve_id,
    const bezier_mesh& mesh, const bezier_params& params);
void update_curve_shape(Gui_Spline& spline, int curve_id,
    const bezier_mesh& mesh, const bezier_params& params,
    const bool quadric = false, const bool use_vector_heat = false);

inline bezier_segment get_control_polygon(
    const Gui_Spline& spline, int curve_id) {
  auto polygon = bezier_segment{};
  for (int i = 0; i < 4; ++i) {
    polygon[i] = spline.control_points[curve_id * 3 + i];
  }
  return polygon;
}

enum struct editing_context : uint {
  is_doing_nothing,
  // is_rotating,
  is_editing_existing_curve,
  is_creating_new_curve,
  is_creating_first_curve,
};

static string context_names[4] = {
    "is_doing_nothing",
    "is_editing_existing_curve",
    "is_creating_new_curve",
    "is_creating_first_curve",
};

struct Gui_Input {
  int        selected_control_point      = -1;
  int        active_control_point        = -1;
  mesh_point new_point                   = {};
  mesh_point first_handle                = {};
  float      drag_start_time             = 0;
  mesh_point drag_start                  = {};
  mesh_point drag_current                = {};
  ogl_shape* temp_points                 = new ogl_shape{};
  bool       enforce_tangents            = true;
  vec2f      active_control_point_offset = {};
  int        selected_spline             = 0;
  float      parameter_t                 = 0.0;
  mesh_point eval_point                  = {};
  mesh_point eval_point_cheap            = {};
  bool       translating                 = false;
  bool       rotating                    = false;
  bool       scaling                     = false;
};

struct Added_Path {
  geodesic_path   path;
  vector<vec3f>   positions;
  vec3f           color;
  float           radius;
  shade_instance* instance;
};

struct Added_Points {
  vector<mesh_point> points;
  vec3f              color;
  float              radius;
  shade_instance*    instance;
};

struct App {
  // TEMPORARY DATA
  vector<Added_Path*>   added_paths         = {};
  vector<Added_Points*> added_points        = {};
  string                exported_scene_name = "scene.ply";
  int                   hue                 = 240;
  int                   sat                 = 100;
  int                   value               = 100;
  int                   level               = 0;
  int                   OLR_steps           = 1;
  bool                  export_edges        = false;
  bool                  weighted_average    = false;
  bool                  use_vector_heat     = false;
  bool                  show_construction   = false;
  bool                  show_wa_curve       = true;
  bool                  show_WA_samples     = false;
  bool                  show_flipout_curve  = false;
  bool                  show_olr            = false;
  bool                  quadric_curve       = false;
  bool                  show_gradient       = false;
  vector<mesh_point>    control_points      = {};
  vector<vec3f>         WA_samples          = {};
  Added_Points*         WA_samples_shape    = {};
  Added_Path*           WA_curve_shape      = {};
  vector<vec3f>         WA_curve            = {};

  int                type_of_average   = 0;
  int                control_polygon_h = 240;
  int                control_polygon_s = 100;
  int                control_polygon_v = 100;
  int                first_level_h     = 252;
  int                first_level_s     = 89;
  int                first_level_v     = 100;
  int                second_level_h    = 324;
  int                second_level_s    = 75;
  int                second_level_v    = 90;
  int                third_level_h     = 0;
  int                third_level_s     = 0;
  int                third_level_v     = 0;
  vector<ogl_shape*> temp_points       = {};
  int                temp_levels       = -1;
  ogl_shape*         eval_point_shape  = new ogl_shape{};
  Svg                svg               = {};
  mesh_point         xxx_point         = {};

  vector<mesh_point>    eval_points    = {};
  vector<mesh_point>    bad_points     = {};
  vector<float>         bad_t          = {};
  int                   curr_bad_t     = 0;
  vector<vector<vec3f>> karcher_grads  = {};
  vector<vector<float>> karcher_fields = {};
  vec2i                 window_size    = {};
  bool                  started        = false;
  int                   playback_tick  = 0;
  bool                  playback       = false;

  bezier_mesh mesh = {};

  float           scale_factor         = 1.f;
  bool            enable_jumps         = false;
  bool            cut_after_first_jump = false;
  editing_context context              = editing_context::is_doing_nothing;
  float           svg_size             = 0.3;

  struct {
    mat4f view            = identity4x4f;
    mat4f projection      = identity4x4f;
    mat4f projection_view = identity4x4f;
  } matrices;

  bool              recording_input = false;
  vector<gui_input> input_record    = {};

  struct Editing_State {
    vector<Gui_Spline> splines = {};
    Gui_Input          input   = {};
  };
  Editing_State state = {{Gui_Spline{}}, {}};

  int                   editing_history_count = 0;
  vector<Editing_State> editing_history       = {};
  int                   history_index         = 0;

  const vector<Gui_Spline>& splines() const { return state.splines; }
  vector<Gui_Spline>&       splines() { return state.splines; }
  const Gui_Input&          input() const { return state.input; }
  Gui_Input&                input() { return state.input; }

  const Gui_Spline& spline() const {
    assert(input().selected_spline != -1);
    return state.splines[input().selected_spline];
  }
  Gui_Spline& spline() {
    assert(input().selected_spline != -1);
    return state.splines[input().selected_spline];
  }

  void commit_state() {
    editing_history.resize(history_index + 1);
    editing_history[history_index] = state;
    history_index += 1;
    printf("%s: %d\n", "commit_state()", history_index);
  }
  void undo() {
    if (history_index > 0) {
      history_index -= 1;
      state = editing_history[history_index];
    }
  }
  void redo() {
    if (history_index < editing_history.size() - 1) {
      history_index += 1;
      state = editing_history[history_index];
    }
  }

  bezier_params bezier_params = {};
  string        filename      = "data/mesh.obj";
  string        testname      = "tests/loop-bunny.json";
  string        comparison    = "mesh";
  string        scene_name    = "data/mesh.obj_gamma";

  float         line_size          = 1.5;
  float         curve_size         = 1;
  float         vector_size        = 0.01;
  vector<vec3f> vector_field       = {};
  int           vector_size_factor = 1;
  mesh_point    last_placed_point  = {};

  float time_of_last_click = -1;

  float angle = 0;

  string error = "";

  bool show_edges           = false;
  bool show_points          = true;
  bool show_control_points  = false;
  bool envlight             = false;
  bool show_control_polygon = true;

  ogl_shape edges_shape = {};

  gui_widget*     widget          = new gui_widget{};
  shade_scene*    scene           = nullptr;
  shade_material* spline_material = nullptr;
  shade_material* mesh_material   = nullptr;
  shade_shape*    mesh_shape      = nullptr;
  shade_camera*   camera          = {};
  shade_params    shade_params{};
  float           camera_focus;
  shape_bvh       bvh = {};

  // Data stored on the gpu for rendering.
  std::unordered_map<string, ogl_shape>   shapes;
  std::unordered_map<string, ogl_program> shaders;
  std::unordered_map<string, gpu::Shape>  gpu_shapes;
  std::unordered_map<string, gpu::Shader> gpu_shaders;
};

inline void clear_all_splines(App& app) {
  for (auto& spline : app.splines()) {
    clear_spline(spline);
  }
}

#include <thread>
template <typename F>
inline void parallel_for(int size, F&& f) {
  auto num_threads = min(size, 16);
  auto threads     = vector<std::thread>(num_threads);
  auto batch_size  = (size + num_threads - 1) / num_threads;

  auto batch = [&](int k) {
    int from = k * batch_size;
    int to   = min(from + batch_size, size);
    for (int i = from; i < to; i++) f(i);
  };

  for (int k = 0; k < num_threads; k++) {
    threads[k] = std::thread(batch, k);
  }
  for (int k = 0; k < num_threads; k++) {
    threads[k].join();
  }
}

template <typename F>
inline void serial_for(int size, F&& f) {
  for (int i = 0; i < size; i++) {
    f(i);
  }
}

inline int get_selected_anchor_point(const App& app) {
  auto id = app.input().selected_control_point;
  if (id != -1 && id % 3 == 1) id -= 1;
  if (id != -1 && id % 3 == 2) id += 1;
  return id;
}
inline vec2i get_selected_handle_points(const App& app) {
  auto id     = app.input().selected_control_point;
  auto result = vec2i{-1, -1};
  if (id % 3 != 0) return result;

  if (id - 1 >= 0) result.x = id - 1;
  if (id + 1 < app.spline().control_points.size()) result.y = id + 1;
  return result;
}

void init_bvh(App& app);

shade_camera _make_framing_camera(const vector<vec3f>& positions);

void init_camera(App& app, const vec3f& from = vec3f{0, 0.5, 1.5},
    const vec3f& to = {0, 0, 0});

vector<vec3f> make_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);

void init_bvh(App& app);

ray3f camera_ray(const App& app, vec2f mouse);

vec2f screenspace_from_worldspace(App& app, const vec3f& position);

mesh_point intersect_mesh(const App& app, vec2f mouse);

void init_gpu(App& app, bool envlight);

void delete_app(App& app);

bool load_program(ogl_program* program, const string& vertex_filename,
    const string& fragment_filename);

void set_points_shape(ogl_shape* shape, const vector<vec3f>& positions);
void set_points_shape(ogl_shape* shape, const bezier_mesh& mesh,
    const vector<mesh_point>& points);
void set_mesh_shape(ogl_shape* shape, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals);
void set_polyline_shape(ogl_shape* shape, const vector<vec3f>& positions);

void update_all_splines(App& app);
void export_scene(App& app, const string& filename, const bool export_edges,
    const vector<vec3f>& spline_colors = {});
inline void update_camera_info(App& app, const gui_input& input) {
  auto& camera   = *app.camera;
  auto  viewport = input.framebuffer_viewport;
  camera.aspect  = (float)viewport.z / (float)viewport.w;

  auto camera_yfov =
      (camera.aspect >= 0)
          ? (2 * yocto::atan(camera.film / (camera.aspect * 2 * camera.lens)))
          : (2 * yocto::atan(camera.film / (2 * camera.lens)));

  app.matrices.view       = frame_to_mat(inverse(camera.frame));
  app.matrices.projection = perspective_mat(
      camera_yfov, camera.aspect, app.shade_params.near, app.shade_params.far);
  app.matrices.projection_view = app.matrices.projection * app.matrices.view;
}

void save_editing(const App& app, const string& filename);
bool load_editing(App& app, const string& filename);

void        export_scene(App& app, const string& filename);
void        update_glvector_field(App& app, const vector<vec3f>& vector_field,
           const float& scale, const string& name);
Added_Path* add_path_shape(
    App& app, const geodesic_path& path, float radius, const vec3f& color);
Added_Path*   add_path_shape(App& app, const vector<vec3f>& positions,
      float radius, const vec3f& color, const bool jumps = false);
Added_Points* add_points_shape(App& app, const vector<mesh_point>& points,
    float radius, const vec3f& color);
Added_Points* add_points_shape(
    App& app, const vector<vec3f>& points, float radius, const vec3f& color);
void update_path_shape(shade_shape* shape, const bezier_mesh& mesh,
    const geodesic_path& path, const float& radius,
    const float& threshold = flt_max);
void update_path_shape(shade_shape* shape, const bezier_mesh& mesh,
    const vector<vec3f>& positions, const float& radius,
    const float& threshold = flt_max);
void update_path_shape(Gui_Curve& curve, const bezier_mesh& mesh,
    const float& radius, const float& treshold,
    const bool cut_after_first = false);
void update_points_shape(
    shade_shape* shape, const vector<vec3f>& positions, float radius);
void update_points_shape(shade_shape* shape, const bezier_mesh& mesh,
    const vector<mesh_point>& points, float radius);

inline vector<mesh_point> karcher_average(const bezier_mesh& mesh,
    const vector<mesh_point>& points, const int number_of_subdivision,
    const bool use_vector_heat) {
  return bezier_karcher_final(mesh.solver, mesh.triangles, mesh.positions,
      mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles, mesh.total_angles,
      mesh.Grad, points, number_of_subdivision, *mesh.flipout, use_vector_heat);
}