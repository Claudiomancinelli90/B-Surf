#pragma once
#include <yocto/yocto_color.h>
#include <yocto/yocto_common.h>  // append...
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deque>
#include <iostream>

#include "karcher.h"
using namespace yocto;

// TODO(giacomo): is this necessary?
inline bool operator==(const mesh_point& a, const mesh_point& b) {
  return (a.face == b.face) && (a.uv == b.uv);
}

inline bool check_point(const mesh_point& point) {
  assert(point.face != -1);
  assert(point.uv.x >= 0);
  assert(point.uv.y >= 0);
  assert(point.uv.x <= 1);
  assert(point.uv.y <= 1);
  return true;
}

namespace flipout {
struct flipout_mesh;
}
struct bezier_mesh {
  vector<vec3i>          triangles   = {};
  vector<vec3i>          adjacencies = {};
  vector<vec3f>          positions   = {};
  vector<vec3f>          normals     = {};
  vector<vec2f>          texcoords   = {};
  flipout::flipout_mesh* flipout     = nullptr;
  // Additional data for experimental algorithms
  geodesic_solver                solver       = {};
  vector<vector<int>>            v2t          = {};
  vector<vector<float>>          angles       = {};
  vector<float>                  total_angles = {};
  dual_geodesic_solver           dual_solver  = {};
  Eigen::SparseMatrix<double, 1> Grad;
  float                          avg_edge_length = 0.f;
};

enum struct spline_algorithm {
  de_casteljau_uniform = 0,
  de_casteljau_adaptive,
  de_casteljau_classic,
  subdivision_uniform,
  subdivision_adaptive,
  karcher,
  flipout
};
const auto spline_algorithm_names = vector<string>{
    "de_casteljau_uniform", "de_casteljau_adaptive", "de_casteljau_classic",
    "subdivision_uniform", "subdivision_adaptive", "karcher", "flipout"

};

const auto levels_names = vector<string>{
    "control polygon", "first level", "second level", "third level", "curve"};

const auto average_names = vector<string>{
    "Karcher", "Euclidean", "Geodesic", "Straight line"};

struct bezier_params {
  spline_algorithm algorithm      = spline_algorithm::subdivision_uniform;
  int              subdivisions   = 8;
  float            precision      = 0.1;
  float            min_curve_size = 0.001;
  int              max_depth      = 10;
  bool             parallel       = false;
};

using bezier_segment = array<mesh_point, 4>;

using quadratic_bezier_segment = array<mesh_point, 3>;

inline bool check_segment(const bezier_segment& segment) {
  assert(check_point(segment[0]));
  assert(check_point(segment[1]));
  assert(check_point(segment[2]));
  assert(check_point(segment[3]));
  return true;
}

struct bezier_polygon {
  bezier_segment               segment;
  std::array<geodesic_path, 3> lines;
};

struct bezier_curve {
  vector<vector<bezier_polygon>> levels;
};

struct bezier_node {
  std::array<mesh_point, 4>    points      = {};
  std::array<geodesic_path, 3> lines       = {};
  int                          parent      = -1;
  int                          children[2] = {-1, -1};
  float                        t_start     = 0;
  float                        t_end       = 1;
};
struct spline_node {
  std::array<mesh_point, 4>    points  = {};
  std::array<geodesic_path, 3> lines   = {};
  vec2f                        t       = {};
  int                          depth   = 0;
  bool                         is_good = false;
};

struct bezier_tree {
  vector<bezier_node> nodes = {};
  int                 depth = 0;
};

inline void add_children(bezier_tree& tree, int parent) {
  auto id = (int)tree.nodes.size();
  tree.nodes.push_back({});
  tree.nodes.push_back({});

  tree.nodes[parent].children[0] = id;
  tree.nodes[parent].children[1] = id + 1;
  tree.nodes[id].parent          = parent;
  tree.nodes[id + 1].parent      = parent;
}

pair<bezier_segment, bezier_segment> subdivide_bezier_polygon(
    const bezier_mesh& mesh, const bezier_segment& input, float t = 0.5);

void subdivide_bezier(
    const bezier_mesh& mesh, bezier_curve& curve, bezier_params params);

void subdivide_bezier_tree(const bezier_mesh& mesh, bezier_tree& tree,
    bezier_params params, float t = 0.5);

inline vector<mesh_point> bezier_tree_points(
    const bezier_mesh& mesh, const bezier_tree& tree, int depth) {
  // TODO(giacomo): slow!
  auto from   = (1 << (depth - 1)) - 1;
  auto to     = (1 << depth) - 1;
  auto points = vector<mesh_point>((to - from) * 3 + 1);

  auto i = 0;
  for (auto node = from; node < to; node++) {
    points[i++] = tree.nodes[node].points[0];
    points[i++] = tree.nodes[node].points[1];
    points[i++] = tree.nodes[node].points[2];
  }
  // printf("[temp] from: %d\n", from);
  // printf("[temp] to: %d\n", to);
  // printf("[temp] size: %zd\n", points.size());
  // printf("[temp] guess: %d\n", (1 << tree.depth) * 3);
  // 0 | 1 2 | 3 4 5 6 | 7 8 9 10 11 12 13 14 | 15
  points[i] = tree.nodes.back().points.back();
  return points;
}

std::array<bezier_segment, 2> insert_point(
    const bezier_mesh& mesh, const bezier_segment& polygon, float t0);

std::array<quadratic_bezier_segment, 2> insert_point(
    const bezier_mesh& mesh, quadratic_bezier_segment& polygon, float t0);

std::array<bezier_segment, 2> insert_point_spline(const bezier_mesh& mesh,
    const bezier_segment& polygon, const float& t0,
    const bezier_params& params);

vector<mesh_point> bezier(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params);

vector<mesh_point> bezier_adaptive(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params);
vector<mesh_point> bezier_adaptive_simple_claudio(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params);
vector<mesh_point> bezier_uniform(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params);
vector<mesh_point> bezier_uniform(const bezier_mesh& mesh,
    const quadratic_bezier_segment&                  control_points,
    const bezier_params&                             params);
vector<mesh_point> bezier_adaptive_simple(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params);
vector<mesh_point> OLR(const bezier_mesh& mesh,
    const vector<mesh_point>& control_polygon, const int step);
vector<mesh_point> weighted_average(const bezier_mesh& mesh,
    const vector<mesh_point>& rectangle, const vector<vector<float>>& weights);

vector<mesh_point> de_casteljau_classic(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    vector<int>& badones, const bool jumps = false);

std::pair<vector<geodesic_path>, vector<vector<mesh_point>>>
dc_classic_construction(
    const bezier_mesh& mesh, const bezier_segment& polygon, const float& t);

std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
dc_construction(
    bezier_mesh& mesh, const bezier_segment& polygon, const float& t0);

std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
LR_algorithm(bezier_mesh& mesh, const bezier_segment& control_points);

std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
RDC_algorithm(bezier_mesh& mesh, const bezier_segment& control_points);

geodesic_path compute_geodesic_path(const bezier_mesh& mesh,
    const mesh_point& start, const mesh_point& end, int thread_id = 0);

geodesic_path my_compute_geodesic_path(const bezier_mesh& mesh,
    const mesh_point& start, const mesh_point& end, int thread_id = 0);

geodesic_path straightest_geodesic_biermann(const bezier_mesh& mesh,
    const mesh_point& start, const vec2f& dir, const float& path_len);

pair<vector<vec3f>, vector<vec3f>> import_curve(
    const bezier_mesh& mesh, const string& filename);

vector<bezier_segment> import_control_points(
    const bezier_mesh& mesh, const string& filename);

vector<vec3f> generate_polyline_from_positions(
    const bezier_mesh& mesh, const vector<vec3f>& pos, const int seed);

inline mesh_point make_point_from_vert(
    const bezier_mesh& mesh, const int vid, const int tid = -1) {
  return point_from_vert(mesh.triangles, mesh.v2t, vid);
}
inline mesh_point force_point_on_vert(
    const bezier_mesh& mesh, const mesh_point& p) {
  auto vid = vert_from_point(mesh.triangles, p);
  return make_point_from_vert(mesh, vid);
}
inline geodesic_path straightest_geodesic_biermann(
    const bezier_mesh& mesh, const mesh_point& start, const vec2f& coord) {
  auto len = length(coord);
  return straightest_geodesic_biermann(mesh, start, coord / len, len);
}

inline geodesic_path straightest_path(const bezier_mesh& mesh,
    const mesh_point& start, const vec2f& direction, float length) {
  return straightest_path(mesh.triangles, mesh.positions, mesh.adjacencies,
      start, direction, length);
}

inline vector<mesh_point> bezier_karcher_test(const bezier_mesh& mesh,
    const bezier_segment& control_points, const int number_of_subdivision,
    vector<int>& jumps) {
  return karcher_test(mesh.dual_solver, mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad,
      {control_points[0], control_points[1], control_points[2],
          control_points[3]},
      number_of_subdivision, jumps);
}
inline mesh_point weighted_average(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const float& t) {
  return weighted_average(mesh.solver, mesh.triangles, mesh.positions,
      mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles, mesh.total_angles,
      mesh.Grad, control_points, t);
}
inline geodesic_path straightest_path(
    const bezier_mesh& mesh, const mesh_point& start, const vec2f& coord) {
  auto len = length(coord);
  return straightest_path(mesh.triangles, mesh.positions, mesh.adjacencies,
      start, coord / len, len);
}

inline mesh_point eval_path_midpoint(
    const bezier_mesh& mesh, const geodesic_path& path) {
  return eval_path_point(
      path, mesh.triangles, mesh.positions, mesh.adjacencies, 0.5);
}

inline vec3f eval_position(const bezier_mesh& mesh, const mesh_point& point) {
  return eval_position(mesh.triangles, mesh.positions, point);
}

inline vec3f eval_normal(const bezier_mesh& mesh, const mesh_point& point) {
  return eval_normal(mesh.triangles, mesh.positions, point);
}

inline mesh_point eval_path_point(
    const bezier_mesh& mesh, const geodesic_path& path, float t) {
  return eval_path_point(
      path, mesh.triangles, mesh.positions, mesh.adjacencies, t);
}

inline vector<vec3f> path_positions(
    const bezier_mesh& mesh, const geodesic_path& path) {
  return path_positions(path, mesh.triangles, mesh.positions, mesh.adjacencies);
}

inline vector<vec3f> path_positions(
    const bezier_mesh& mesh, const mesh_path& path) {
  auto positions = vector<vec3f>(path.points.size());
  for (int i = 0; i < positions.size(); i++) {
    positions[i] = eval_position(mesh, path.points[i]);
  }
  return positions;
}

// TODO(giacomo): migrate everything to mesh_path
// inline vector<vec3f> path_positions_new(
//     const bezier_mesh& mesh, const geodesic_path& path) {
//   auto mp = convert_mesh_path(mesh.triangles, mesh.adjacencies, path);
//   return path_positions(mesh, mp);
// }

// TODO(giacomo): put in yocto_mesh (explode params)
inline vec2f tangent_path_direction(
    const bezier_mesh& mesh, const geodesic_path& path, bool start = true) {
  auto find = [](const vec3i& vec, int x) {
    for (int i = 0; i < size(vec); i++)
      if (vec[i] == x) return i;
    return -1;
  };

  auto direction = vec2f{};

  if (start) {
    auto start_tr = triangle_coordinates(
        mesh.triangles, mesh.positions, path.start);

    if (path.lerps.empty()) {
      direction = interpolate_triangle(
          start_tr[0], start_tr[1], start_tr[2], path.end.uv);
    } else {
      auto x    = path.lerps[0];
      auto k    = find(mesh.adjacencies[path.strip[0]], path.strip[1]);
      direction = lerp(start_tr[k], start_tr[(k + 1) % 3], x);
    }
  } else {
    auto end_tr = triangle_coordinates(
        mesh.triangles, mesh.positions, path.end);
    if (path.lerps.empty()) {
      direction = interpolate_triangle(
          end_tr[0], end_tr[1], end_tr[2], path.start.uv);
    } else {
      auto x = path.lerps.back();
      auto k = find(
          mesh.adjacencies[path.strip.rbegin()[0]], path.strip.rbegin()[1]);
      direction = lerp(end_tr[k], end_tr[(k + 1) % 3], 1 - x);
    }
  }
  return normalize(direction);
}

inline geodesic_path continue_path(
    const bezier_mesh& mesh, const geodesic_path& path, float length) {
  auto direction = tangent_path_direction(mesh, path);
  if (length < 0) direction = -direction;
  return straightest_path(mesh.triangles, mesh.positions, mesh.adjacencies,
      path.start, direction, yocto::abs(length));
};

// TODO: organize wrappers
// TODO: organize wrappers
// TODO: organize wrappers
// TODO: organize wrappers
// TODO: organize wrappers
// TODO: organize wrappers

mesh_point eval_bezier_point_cheap(
    const bezier_mesh& mesh, const bezier_segment& polygon, float t);

mesh_point eval_bezier_point(const bezier_mesh& mesh,
    const bezier_segment& segment, float t0, float t_start = 0,
    float t_end = 1);

mesh_point eval_spline_point(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    const float& t0);

inline vector<mesh_point> eval_spline_point(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    int num_points) {
  auto result = vector<mesh_point>(num_points);
  for (int i = 0; i < num_points; i++) {
    result[i] = eval_spline_point(mesh, polygon, params, float(i) / num_points);
  }
  return result;
}

inline vector<geodesic_path> make_polyline(
    const bezier_mesh& mesh, const vector<mesh_point>& points) {
  auto result = vector<geodesic_path>(points.size() - 1);
  for (int i = 0; i < points.size() - 1; i++) {
    auto a    = points[i];
    auto b    = points[i + 1];
    result[i] = compute_geodesic_path(mesh, a, b);
  }
  return result;
}

inline vector<vec3f> make_polyline_positions(
    const bezier_mesh& mesh, const vector<mesh_point>& points) {
  auto result = vector<vec3f>();  // TODO: reserve what?
  for (int i = 0; i < points.size() - 1; i++) {
    auto a    = points[i];
    auto b    = points[i + 1];
    auto path = compute_geodesic_path(mesh, a, b);
    append(result, path_positions(mesh, path));
  }
  return result;
}

inline float path_length(const bezier_mesh& mesh, const geodesic_path& path) {
  auto positions = path_positions(
      path, mesh.triangles, mesh.positions, mesh.adjacencies);

  auto result = 0.0f;
  for (int i = 1; i < positions.size(); ++i)
    result += length(positions[i] - positions[i - 1]);

  return result;
}

inline void print(const bezier_segment& segment) {
  for (int i = 0; i < 4; i++) {
    auto p = segment[i];
    printf("{\"face\": %d, \"uv\": [%f, %f]}}\n", p.face, p.uv.x, p.uv.y);
  }
}

vector<mesh_point> spline_subdivision_uniform(const bezier_mesh& mesh,
    const array<mesh_point, 4>& control_points, int num_subdivisions);

vector<mesh_point> spline_subdivision_uniform(const bezier_mesh& mesh,
    const quadratic_bezier_segment& control_points, int num_subdivisions);

vector<mesh_point> spline_subdivision_adaptive(const bezier_mesh& mesh,
    const array<mesh_point, 4>& polygon, const bezier_params& params);

bool point_are_close(const mesh_point& a, const mesh_point& b);

inline double avg_edge_length(
    const vector<vec3f>& positions, const vector<vec3i>& triangles) {
  double sum = 0;

  for (int i = 0; i < triangles.size(); ++i) {
    sum += length(positions[triangles[i].x] - positions[triangles[i].y]);
    sum += length(positions[triangles[i].z] - positions[triangles[i].y]);
    sum += length(positions[triangles[i].x] - positions[triangles[i].z]);
  }

  return sum / (positions.size() * 6);
}

/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

#if 0
#include "karcher.h"

inline mesh_point transport_vector(const bezier_mesh& mesh,
    const geodesic_path& path, const mesh_point& point) {
  return transport_vector(mesh.solver, mesh.angles, mesh.total_angles,
      mesh.triangles, mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t,
      path, point);
}

inline vector<mesh_point> bezier_karcher_full_gradient(
    const bezier_mesh& mesh, const vector<mesh_point>& control_points) {
  return bezier_karcher_full_gradient(mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points);
}
inline vector<mesh_point> bezier_karcher(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, splits& composite_bezier,
    const bool all_in_once = true, const vec2i& knot = {-1, -1}) {
  return bezier_karcher(mesh.dual_solver, mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points, number_of_subdivision,
      type_of_field, composite_bezier, all_in_once, knot);
}
inline vector<mesh_point> bezier_karcher(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field) {
  return bezier_karcher(mesh.dual_solver, mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points, number_of_subdivision,
      type_of_field);
}
inline vector<mesh_point> bezier_karcher_test(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<vector<vec3f>>& directions,
    vector<vector<mesh_point>>& steps, vector<pair<mesh_point, int>>& badones) {
  return bezier_karcher_test(mesh.dual_solver, mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points, number_of_subdivision,
      type_of_field, verts, tids, directions, steps, badones);
}
inline vector<mesh_point> bezier_karcher_almost_gradients(
    const bezier_mesh& mesh, const vector<mesh_point>& control_points,
    const int number_of_subdivision, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<pair<mesh_point, float>>& badones) {
  return bezier_karcher_almost_gradients(mesh.dual_solver, mesh.solver,
      mesh.triangles, mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t,
      mesh.angles, mesh.total_angles, control_points, number_of_subdivision,
      verts, tids, badones);
}

inline vector<mesh_point> bezier_karcher_tet(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const int number_of_subdivision) {
  return bezier_karcher_tet(mesh.solver, mesh.dual_solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points, number_of_subdivision);
}
inline vector<mesh_point> bezier_karcher_tet_subdivision(
    const bezier_mesh& mesh, const vector<mesh_point>& control_points,
    const int number_of_subdivision) {
  return bezier_karcher_tet_subdivision(mesh.solver, mesh.dual_solver,
      mesh.triangles, mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t,
      mesh.angles, mesh.total_angles, mesh.Grad, control_points,
      number_of_subdivision);
}
inline vector<mesh_point> bezier_karcher_bisection(const bezier_mesh& mesh,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<vector<mesh_point>>& final_control_polygons,
    vector<vector<mesh_point>>& bad_control_polygons,
    vector<int>&                type_of_badones) {
  return bezier_karcher_bisection(mesh.dual_solver, mesh.solver, mesh.triangles,
      mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t, mesh.angles,
      mesh.total_angles, mesh.Grad, control_points, number_of_subdivision,
      type_of_field, final_control_polygons, bad_control_polygons,
      type_of_badones);
}
inline vector<mesh_point> bezier_karcher_bisection_hybrid(
    const bezier_mesh& mesh, const vector<mesh_point>& control_points,
    const int number_of_subdivision, const int type_of_field,
    vector<vector<mesh_point>>& final_control_polygons,
    vector<vector<mesh_point>>& bad_control_polygons,
    vector<int>&                type_of_badones) {
  return bezier_karcher_bisection_hybrid(mesh.solver, mesh.dual_solver,
      mesh.triangles, mesh.positions, mesh.adjacencies, mesh.normals, mesh.v2t,
      mesh.angles, mesh.total_angles, mesh.Grad, control_points,
      number_of_subdivision, type_of_field, final_control_polygons,
      bad_control_polygons, type_of_badones);
}
#else

#endif
