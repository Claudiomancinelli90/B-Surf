#include "spline.h"

#include <logging.h>

#include <deque>
#include <mutex>

#include "strip.h"

static auto       arenas            = vector<strip_arena>{};
static auto       arenas_count      = 0;
static auto       free_arenas       = vector<int>{};
static auto       free_arenas_count = 0;
static std::mutex mutex;

int get_free_arena_id(size_t size) {
  mutex.lock();
  if (arenas.empty()) arenas.resize(64);
  if (free_arenas.empty()) free_arenas.resize(64);

  auto id = -1;
  if (free_arenas_count == 0) {
    id = (int)arenas.size();
    init_arena(arenas[arenas_count], size);
    arenas_count += 1;
  } else {
    id = free_arenas[free_arenas_count - 1];
    free_arenas_count -= 1;
  }
  mutex.unlock();
  return id;
}

void release_arena(int id) {
  mutex.lock();
  free_arenas[free_arenas_count++] = id;
  mutex.unlock();
}
bool point_are_close(const mesh_point& a, const mesh_point& b) {
  if (a.face == b.face) {
    if (yocto::abs(a.uv.x - b.uv.x) < 1e-10 &&
        yocto::abs(a.uv.y - b.uv.y) < 1e-10)
      return true;
  }
  return false;
}
// void clean_bary(mesh_point& sample) {
//   auto bary   = sample.uv;
//   auto coords = vector<pair<float, int>>{
//       {1 - bary.x - bary.y, 0}, {bary.x, 1}, {bary.y, 2}};
//   sort(coords.begin(), coords.end());
//   vec3f bary3d = get_bary(bary);
//   if (coords[0].first < 0 && coords[1].first > 0) {
//     bary3d[coords[0].second] = 0;
//     bary3d[coords[2].second] = 1 - bary3d[coords[1].second] -
//                                bary3d[coords[0].second];
//   } else if (coords[0].first < 0) {
//     bary3d[coords[0].second] = 0;
//     bary3d[coords[1].second] = 0;
//     bary3d[coords[2].second] = 1;
//   }
//   // } else if (coords[2].second > 1) {
//   //   bary3d[coords[0].second] = 1e-5;
//   //   bary3d[coords[1].second] = 1e-5;
//   //   bary3d[coords[2].second] = 0.99998;
//   // }

//   sample = {sample.face, {bary3d.y, bary3d.z}};
// }
geodesic_path my_compute_geodesic_path(const bezier_mesh& mesh,
    const mesh_point& start, const mesh_point& end, int thread_id) {
  // profile_function();

  vector<int> parents;
  auto        strip = get_strip(mesh.solver, mesh.triangles, mesh.positions,
             mesh.adjacencies, mesh.v2t, mesh.angles, end, start, parents);
  auto        path  = geodesic_path{};
  auto [cleaned, new_start, new_end] = cleaned_strip(
      mesh.triangles, mesh.positions, mesh.adjacencies, strip, start, end);

  if (new_start.face == new_end.face) {
    path.start = new_start;
    path.end   = new_end;
    path.strip = {new_start.face};
    return path;
  }

  path = shortest_path(mesh.triangles, mesh.positions, mesh.adjacencies,
      new_start, new_end, cleaned);
  return path;
}
geodesic_path compute_geodesic_path(const bezier_mesh& mesh,
    const mesh_point& start, const mesh_point& end, int thread_id) {
  // profile_function();
  check_point(start);
  check_point(end);
  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end   = end;
    path.strip = {start.face};
    return path;
  }
  auto strip = compute_strip_tlv(mesh, end.face, start.face);
  // vector<int> parents;
  // auto        strip = get_strip(mesh.solver, mesh.triangles, mesh.positions,
  //     mesh.adjacencies, mesh.v2t, mesh.angles, end, start, parents);
  // auto [cleaned, new_start, new_end] = cleaned_strip(
  //     mesh.triangles, mesh.positions, mesh.adjacencies, strip, start, end);

  // if (new_start.face == new_end.face) {
  //   path.start = new_start;
  //   path.end   = new_end;
  //   path.strip = {new_start.face};
  //   return path;
  // }

  path = shortest_path(
      mesh.triangles, mesh.positions, mesh.adjacencies, start, end, strip);
  return path;
}

inline mesh_point geodesic_midpoint(
    const bezier_mesh& mesh, const mesh_point& start, const mesh_point& end) {
  // profile_function();

  if (start.face == end.face) {
    return mesh_point{start.face, (start.uv + end.uv) * 0.5};
  }

  auto path     = compute_geodesic_path(mesh, start, end);
  auto midpoint = eval_path_midpoint(
      path, mesh.triangles, mesh.positions, mesh.adjacencies);
  //  assert(check_point(midpoint));
  return midpoint;
}

inline mesh_point geodesic_lerp(const bezier_mesh& mesh,
    const mesh_point& start, const mesh_point& end, float t) {
  // profile_function();

  if (start.face == end.face) {
    return mesh_point{start.face, lerp(start.uv, end.uv, t)};
  }

  auto path  = compute_geodesic_path(mesh, start, end);
  auto point = eval_path_point(
      path, mesh.triangles, mesh.positions, mesh.adjacencies, t);
  assert(check_point(point));
  return point;
}

inline mesh_point geodesic_lerp(const bezier_mesh& mesh, const mesh_point& a,
    const mesh_point& b, const mesh_point& c, float t0, float t1) {
  // den := (1-t0-t1) + t0 = 1 - t1;
  auto t  = t0 / (1 - t1);
  auto ab = geodesic_lerp(mesh, a, b, t);
  return geodesic_lerp(mesh, ab, c, t1);
}
pair<int, float> binary_search(const vector<float>& v, float t) {
  auto L = 0, R = (int)v.size(), entry = -1, m = 0;
  auto factor = 0.0;
  while (L < R) {
    m = (int)floor((L + R) / 2);
    if (t < v[m])
      R = m;
    else
      L = m + 1;
  }
  entry = R - 1;
  if (entry == v.size() - 1) {
    assert(t >= v[entry - 1] && t <= v[entry]);
    auto width = (v[entry] - v[entry - 1]);
    factor     = (t - v[entry - 1]) / width;
  } else if (t >= v[entry] && t <= v[entry + 1]) {
    auto width = (v[entry + 1] - v[entry]);
    factor     = (t - v[entry]) / width;
  } else
    assert(false);

  return {entry, factor};
}
vec3f path_pos(const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path, int entry) {
  auto u   = path.lerps[entry];
  auto eid = get_edge(triangles, positions, adjacencies, path.strip[entry],
      path.strip[entry + 1]);
  if (eid.x < 0) {
    assert(path.strip.size() == 1);
    return eval_position(triangles, positions, path.end);
  }
  auto p0 = positions[eid.x];
  auto p1 = positions[eid.y];
  return (1 - u) * p0 + u * p1;
}
mesh_point eval_geodesic_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, const vector<float>& path_parameter_t,
    const float& t) {
  // WHY DOUBLE IN THE INTERFACE?
  auto a     = mesh_point{};
  auto entry = -1;
  if (path.strip.size() == 1) {
    entry               = 0;
    auto p              = eval_position(triangles, positions, path.start);
    auto q              = eval_position(triangles, positions, path.end);
    auto start          = (p + q) / 2;
    auto tid            = path.end.face;
    auto [inside, bary] = point_in_triangle(triangles, positions, tid, start);

    assert(inside);

    a = {tid, bary};
  } else {
    auto i = binary_search(path_parameter_t, t);
    auto p = zero3f, q = zero3f;
    auto tid = -1;
    entry    = i.first;
    if (entry == 0) {
      tid = path.start.face;
      p   = eval_position(triangles, positions, path.start);
      q   = path_pos(triangles, positions, adjacencies, path, 0);
    } else if (entry == path_parameter_t.size() - 1) {
      return path.end;

    } else {
      p   = path_pos(triangles, positions, adjacencies, path, i.first - 1);
      q   = (entry == path.strip.size() - 1)
                ? eval_position(triangles, positions, path.end)
                : path_pos(triangles, positions, adjacencies, path, i.first);
      tid = path.strip[i.first];
    }
    auto start = q - p;
    start *= i.second;
    start += p;

    auto [inside, bary] = point_in_triangle(triangles, positions, tid, start);
    if (!inside) {
      tid                    = path.strip[i.first - 1];
      std::tie(inside, bary) = point_in_triangle(
          triangles, positions, tid, start);
      entry = i.first - 1;
    }
    assert(inside);
    a = {tid, bary};
  }

  return a;
}
using bezier_segment = array<mesh_point, 4>;

pair<bezier_segment, bezier_segment> subdivide_bezier_polygon(
    const bezier_mesh& mesh, const bezier_segment& input, float t) {
  auto Q0 = geodesic_lerp(mesh, input[0], input[1], t);
  auto Q1 = geodesic_lerp(mesh, input[1], input[2], t);
  auto Q2 = geodesic_lerp(mesh, input[2], input[3], t);
  auto R0 = geodesic_lerp(mesh, Q0, Q1, t);
  auto R1 = geodesic_lerp(mesh, Q1, Q2, t);
  auto S  = geodesic_lerp(mesh, R0, R1, t);
  return {{input[0], Q0, R0, S}, {S, R1, Q2, input[3]}};
}
pair<quadratic_bezier_segment, quadratic_bezier_segment>
subdivide_bezier_polygon(
    const bezier_mesh& mesh, const quadratic_bezier_segment& input, float t) {
  auto Q0 = geodesic_lerp(mesh, input[0], input[1], t);
  auto Q1 = geodesic_lerp(mesh, input[1], input[2], t);
  auto S  = geodesic_lerp(mesh, Q0, Q1, t);
  return {{input[0], Q0, S}, {S, Q1, input[2]}};
}

inline pair<bezier_polygon, bezier_polygon> subdivide_bezier(
    const bezier_polygon& input, const bezier_mesh& mesh) {
  auto& segment = input.segment;
  auto  Q0      = eval_path_midpoint(mesh, input.lines[0]);
  auto  Q1      = eval_path_midpoint(mesh, input.lines[1]);
  auto  Q2      = eval_path_midpoint(mesh, input.lines[2]);

  auto R0                = geodesic_midpoint(mesh, Q0, Q1);
  auto R1                = geodesic_midpoint(mesh, Q1, Q2);
  auto S                 = geodesic_midpoint(mesh, R0, R1);
  auto result            = pair<bezier_polygon, bezier_polygon>{};
  result.first.lines[0]  = compute_geodesic_path(mesh, segment[0], Q0);
  result.first.lines[1]  = compute_geodesic_path(mesh, Q0, R0);
  result.first.lines[2]  = compute_geodesic_path(mesh, R0, S);
  result.first.segment   = {segment[0], Q0, R0, S};
  result.second.lines[0] = compute_geodesic_path(mesh, S, R1);
  result.second.lines[1] = compute_geodesic_path(mesh, R1, Q2);
  result.second.lines[2] = compute_geodesic_path(mesh, Q2, segment[3]);
  result.second.segment  = {S, R1, Q2, segment[1]};
  return result;
  // return {{P0, Q0, R0, S}, {S, R1, Q2, P3}};
}

void subdivide_bezier_node(
    bezier_tree& tree, int node, const bezier_mesh& mesh, float t) {
  add_children(tree, node);
  auto& polygon = tree.nodes[node];
  auto  Q0      = eval_path_point(mesh, polygon.lines[0], t);
  auto  Q1      = eval_path_point(mesh, polygon.lines[1], t);
  auto  Q2      = eval_path_point(mesh, polygon.lines[2], t);
  assert(check_point(Q0));
  assert(check_point(Q1));
  assert(check_point(Q2));

  auto R0 = geodesic_lerp(mesh, Q0, Q1, t);
  auto R1 = geodesic_lerp(mesh, Q1, Q2, t);
  auto S  = geodesic_lerp(mesh, R0, R1, t);
  assert(check_point(R0));
  assert(check_point(R1));
  assert(check_point(S));

  auto& left  = tree.nodes[tree.nodes[node].children[0]];
  left.points = {polygon.points[0], Q0, R0, S};
  assert(check_segment(left.points));
  left.lines[0] = compute_geodesic_path(mesh, polygon.points[0], Q0);
  left.lines[1] = compute_geodesic_path(mesh, Q0, R0);
  left.lines[2] = compute_geodesic_path(mesh, S, R0);
  left.t_start  = polygon.t_start;
  left.t_end    = lerp(polygon.t_start, polygon.t_end, t);
  auto& right   = tree.nodes[tree.nodes[node].children[1]];
  right.points  = {S, R1, Q2, polygon.points[3]};
  assert(check_segment(right.points));
  right.lines[0] = compute_geodesic_path(mesh, S, R1);
  right.lines[1] = compute_geodesic_path(mesh, R1, Q2);
  right.lines[2] = compute_geodesic_path(mesh, polygon.points[3], Q2);
  right.t_start  = left.t_end;
  right.t_end    = polygon.t_end;
  // return {{P0, Q0, R0, S}, {S, R1, Q2, P3}};
}

inline bool is_right_child(const bezier_tree& tree, int node) {
  assert(tree.nodes[node].parent != -1);
  return tree.nodes[tree.nodes[node].parent].children[1] == node;
}

inline bool is_left_child(const bezier_tree& tree, int node) {
  assert(tree.nodes[node].parent != -1);
  return tree.nodes[tree.nodes[node].parent].children[0] == node;
}

bool is_control_polygon_unfoldable(
    const bezier_mesh& mesh, const bezier_segment& segment) {
  if (segment[0].face != segment[1].face) return false;
  if (segment[1].face != segment[2].face) return false;
  if (segment[2].face != segment[3].face) return false;
  return true;

  // auto faces = std::vector<int>(1);
  // faces[0]   = P0.face;
  // if (find_in_vector(faces, P1.face) == -1) faces.push_back(P1.face);
  // if (find_in_vector(faces, P2.face) == -1) faces.push_back(P2.face);
  // if (find_in_vector(faces, P3.face) == -1) faces.push_back(P3.face);
  // if (faces.size() == 1) return true;
  // if (faces.size() == 2) {
  //   if (find_in_vector(mesh.adjacencies[faces[0]], faces[1]) != -1) {
  //     return true;
  //   }
  // }
  // return false;

  // for (int i = 1; i < 4; ++i) {
  //   if (segment[i - 1].face == segment[i].face) continue;
  //   auto adj = mesh.adjacencies[segment[i - 1].face];
  //   if (find_in_vector(adj, segment[i].face) == -1) return false;
  // }
  // return true;
}

mesh_point eval_bezier_point_cheap(
    const bezier_mesh& mesh, const bezier_segment& polygon, float t) {
  auto Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t);
  auto Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t);
  auto Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t);
  auto R0 = geodesic_lerp(mesh, Q0, Q1, t);
  auto R1 = geodesic_lerp(mesh, Q1, Q2, t);
  return geodesic_lerp(mesh, R0, R1, t);
}

// static bezier_segment bezier_polygon_from_spline_polygon(
//     const bezier_mesh& mesh, const bezier_segment& polygon) {
//   auto p0 = geodesic_lerp(
//       mesh, polygon[0], polygon[1], polygon[2], 2 / 3.f, 1 / 6.f);
//   auto p1 = geodesic_lerp(mesh, polygon[1], polygon[2], 1 / 3.f);
//   auto p2 = geodesic_lerp(mesh, polygon[1], polygon[2], 2 / 3.f);
//   auto p3 = geodesic_lerp(
//       mesh, polygon[1], polygon[2], polygon[3], 2 / 3.f, 1 / 6.f);
//   return {p0, p1, p2, p3};
// }

mesh_point eval_bezier_point(const bezier_mesh& mesh,
    const bezier_segment& segment, float t0, float t_start, float t_end) {
  auto points = segment;

  // Recursively subdivide, shrinking the active control polygon around the
  // point of interest, until it is unfoldable.
  while (true) {
    auto Q0    = geodesic_lerp(mesh, points[0], points[1], 0.5);
    auto Q1    = geodesic_lerp(mesh, points[1], points[2], 0.5);
    auto Q2    = geodesic_lerp(mesh, points[2], points[3], 0.5);
    auto R0    = geodesic_lerp(mesh, Q0, Q1, 0.5);
    auto R1    = geodesic_lerp(mesh, Q1, Q2, 0.5);
    auto S     = geodesic_lerp(mesh, R0, R1, 0.5);
    auto mid_t = (t_start + t_end) / 2;
    if (t0 < mid_t) {
      points[1] = Q0;
      points[2] = R0;
      points[3] = S;
      t_end     = mid_t;
    } else {
      points[0] = S;
      points[1] = R1;
      points[2] = Q2;
      t_start   = mid_t;
    }

    if (is_control_polygon_unfoldable(mesh, points)) break;
  }

  // Exact evaulation.
  float t = (t0 - t_start) / (t_end - t_start);
  return eval_bezier_point_cheap(mesh, points, t);
}

bezier_segment from_spline_to_bezier(const bezier_mesh& mesh,
    const vector<float>& interval, const bezier_segment& polygon) {
  vector<float>  deltas(5);
  bezier_segment bezier_polygon;
  for (auto i = 1; i < 6; ++i) {
    deltas[i] = interval[i] - interval[i - 1];
  }
  auto a = (pow(deltas[2], 2)) /
           ((deltas[0] + deltas[1] + deltas[2]) * (deltas[1] + deltas[2]));
  auto b = (pow(deltas[2], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  auto c = (pow(deltas[2], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[2] + deltas[3]));
  auto d = (pow(deltas[2], 2)) /
           ((deltas[2] + deltas[3] + deltas[4]) * (deltas[3] + deltas[4]));
  auto e = (deltas[1] * deltas[2]) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  auto f = (pow(deltas[1], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  if (a != 0 && 1 - a - f != 0 && f != 0)
    bezier_polygon[0] = geodesic_lerp(
        mesh, polygon[0], polygon[1], polygon[2], 1 - a - f, f);
  else if (a != 0 && 1 - a != 0)
    bezier_polygon[0] = geodesic_lerp(mesh, polygon[0], polygon[1], 1 - a);
  else if (a != 0 && f != 0)
    bezier_polygon[0] = geodesic_lerp(mesh, polygon[0], polygon[2], f);
  else if (a != 0)
    bezier_polygon[0] = polygon[0];
  else if (f != 0 && 1 - f != 0)
    bezier_polygon[0] = geodesic_lerp(mesh, polygon[1], polygon[2], 1 - f);
  else if (f != 0)
    bezier_polygon[0] = polygon[2];
  else
    bezier_polygon[0] = polygon[1];

  if (1 - e - f != 0 && e + f != 0)
    bezier_polygon[1] = geodesic_lerp(mesh, polygon[1], polygon[2], e + f);
  else if (1 - e - f != 0)
    bezier_polygon[1] = polygon[1];
  else
    bezier_polygon[1] = polygon[2];

  if (1 - 2 * e - b - f != 0 && b + 2 * e + f != 0)
    bezier_polygon[2] = geodesic_lerp(
        mesh, polygon[1], polygon[2], b + 2 * e + f);
  else if (1 - 2 * e - b - f != 0)
    bezier_polygon[2] = polygon[1];
  else
    bezier_polygon[2] = polygon[2];

  if (1 - 3 * e - 2 * b + c - f != 0 && 2 * b - c - d + 3 * e + f != 0 &&
      d != 0)
    bezier_polygon[3] = geodesic_lerp(
        mesh, polygon[1], polygon[2], polygon[3], 2 * c - c - d + 3 * e + f, d);
  else if (1 - 3 * e - 2 * b + c - f != 0 && 2 * b - c + 3 * e + f != 0)
    bezier_polygon[3] = geodesic_lerp(
        mesh, polygon[1], polygon[2], 2 * b - c + 3 * e + f);
  else if (1 - 3 * e - 2 * b + c - f != 0 && d != 0)
    bezier_polygon[3] = geodesic_lerp(mesh, polygon[1], polygon[2], d);
  else if (1 - 3 * e - 2 * b + c - f != 0)
    bezier_polygon[3] = polygon[1];
  else if (2 * b - c - d + 3 * e + f != 0 && d != 0)
    bezier_polygon[3] = geodesic_lerp(mesh, polygon[2], polygon[3], d);
  else if (d != 0)
    bezier_polygon[3] = polygon[3];
  else
    bezier_polygon[3] = polygon[2];

  return bezier_polygon;
}
mesh_point de_boor(const bezier_mesh& mesh, const bezier_segment& polygon,
    const vector<float>& knot_vector, const float& t0) {
  bezier_segment old_points  = polygon;
  bezier_segment curr_points = {};
  for (auto j = 1; j <= 3; ++j) {
    for (auto i = 3; i >= j; --i) {
      auto alpha = (t0 - knot_vector[i]) /
                   (knot_vector[i - j + 4] - knot_vector[i]);
      curr_points[i] = geodesic_lerp(
          mesh, old_points[i - 1], old_points[i], alpha);
    }
    old_points = curr_points;
  }
  return curr_points.back();
}
bool spline_stop_criterion(const bezier_mesh& mesh,
    const bezier_segment& polygon, const float& treshold) {
  auto L01       = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  auto L12       = compute_geodesic_path(mesh, polygon[1], polygon[2]);
  auto L23       = compute_geodesic_path(mesh, polygon[2], polygon[3]);
  auto L03       = compute_geodesic_path(mesh, polygon[0], polygon[3]);
  auto perimeter = path_length(mesh, L01) + path_length(mesh, L12) +
                   path_length(mesh, L23);
  if (perimeter - path_length(mesh, L03) <= treshold) return true;

  return false;
}
mesh_point LR_regular(const bezier_mesh& mesh, const mesh_point& a,
    const mesh_point& b, const mesh_point& c) {
  auto Q0 = geodesic_lerp(mesh, a, b, 0.75);
  auto Q1 = geodesic_lerp(mesh, b, c, 0.25);
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}
mesh_point LR_regular(
    const bezier_mesh& mesh, const geodesic_path& l0, const geodesic_path& l1) {
  auto Q0 = eval_path_point(mesh, l0, 0.75);
  auto Q1 = eval_path_point(mesh, l1, 0.25);
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}
mesh_point LR_init(const bezier_mesh& mesh, const mesh_point& a,
    const mesh_point& b, const mesh_point& c) {
  auto Q0 = geodesic_lerp(mesh, a, b, 5 / 8.f);
  auto Q1 = geodesic_lerp(mesh, b, c, 3 / 8.f);
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}
mesh_point LR_init(
    const bezier_mesh& mesh, const geodesic_path& l0, const geodesic_path& l1) {
  auto Q0 = eval_path_point(mesh, l0, 5 / 8.f);
  auto Q1 = eval_path_point(mesh, l1, 3 / 8.f);
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}

mesh_point LR_boundary(const bezier_mesh& mesh, const mesh_point& a,
    const mesh_point& b, const mesh_point& c, const bool& left) {
  auto Q0 = mesh_point{};
  auto Q1 = mesh_point{};
  if (left) {
    Q0 = geodesic_lerp(mesh, a, b, 5 / 8.f);
    Q1 = geodesic_lerp(mesh, b, c, 0.25);

  } else {
    Q0 = geodesic_lerp(mesh, a, b, 0.75);
    Q1 = geodesic_lerp(mesh, b, c, 3 / 8.f);
  }
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}
mesh_point LR_boundary(const bezier_mesh& mesh, const geodesic_path& l0,
    const geodesic_path& l1, const bool& left) {
  auto Q0 = mesh_point{};
  auto Q1 = mesh_point{};
  if (left) {
    Q0 = eval_path_point(mesh, l0, 5 / 8.f);
    Q1 = eval_path_point(mesh, l1, 0.25);

  } else {
    Q0 = eval_path_point(mesh, l0, 0.75);
    Q1 = eval_path_point(mesh, l1, 3 / 8.f);
  }
  return geodesic_lerp(mesh, Q0, Q1, 0.5);
}
pair<bool, bezier_segment> handle_boundary(const bezier_mesh& mesh,
    const bezier_segment& control_points, const vector<int>& new_ones_entries,
    const int k) {
  bezier_segment new_ones = {};
  if (new_ones_entries[0] > 3 && new_ones_entries.back() < pow(2, k) - 1)
    return {false, {}};
  else if (new_ones_entries[0] <= 3) {
    switch (new_ones_entries[0]) {
      case 0: {
        new_ones[0] = control_points[0];
        new_ones[1] = geodesic_lerp(
            mesh, control_points[0], control_points[1], 0.5);
        new_ones[2] = geodesic_lerp(
            mesh, control_points[1], control_points[2], 0.25);
        new_ones[3] = LR_boundary(mesh, control_points[1], control_points[2],
            control_points[3], true);
      } break;
      case 1: {
        new_ones[0] = geodesic_lerp(
            mesh, control_points[0], control_points[1], 0.5);
        new_ones[1] = geodesic_lerp(
            mesh, control_points[1], control_points[2], 0.25);
        new_ones[2] = LR_boundary(mesh, control_points[1], control_points[2],
            control_points[3], true);
        new_ones[3] = geodesic_lerp(
            mesh, control_points[2], control_points[3], 0.5);

      } break;
      case 2: {
        new_ones[0] = geodesic_lerp(
            mesh, control_points[0], control_points[1], 0.25);
        new_ones[1] = LR_boundary(mesh, control_points[0], control_points[1],
            control_points[2], true);
        new_ones[2] = geodesic_lerp(
            mesh, control_points[1], control_points[2], 0.5);
        new_ones[3] = LR_regular(
            mesh, control_points[1], control_points[2], control_points[3]);

      } break;

      case 3: {
        new_ones[0] = LR_boundary(mesh, control_points[0], control_points[1],
            control_points[2], true);
        new_ones[1] = geodesic_lerp(
            mesh, control_points[1], control_points[2], 0.5);
        new_ones[2] = LR_regular(
            mesh, control_points[1], control_points[2], control_points[3]);
        new_ones[3] = geodesic_lerp(
            mesh, control_points[2], control_points[3], 0.5);

      } break;
    }
  } else {
    if (new_ones_entries.back() == pow(2, k) + 2) {
      new_ones[0] = LR_boundary(
          mesh, control_points[0], control_points[1], control_points[2], false);
      new_ones[1] = geodesic_lerp(
          mesh, control_points[1], control_points[2], 0.75);
      new_ones[2] = geodesic_lerp(
          mesh, control_points[2], control_points[3], 0.5);
      new_ones[3] = control_points[3];
    } else if (new_ones_entries.back() == pow(2, k) + 1) {
      new_ones[0] = geodesic_lerp(
          mesh, control_points[0], control_points[1], 0.5);
      new_ones[1] = LR_boundary(
          mesh, control_points[0], control_points[1], control_points[2], false);
      new_ones[2] = geodesic_lerp(
          mesh, control_points[1], control_points[2], 0.75);
      new_ones[3] = geodesic_lerp(
          mesh, control_points[2], control_points[3], 0.5);
    } else if (new_ones_entries.back() == pow(2, k)) {
      new_ones[0] = LR_regular(
          mesh, control_points[0], control_points[1], control_points[2]);
      new_ones[1] = geodesic_lerp(
          mesh, control_points[1], control_points[2], 0.5);
      new_ones[2] = LR_boundary(
          mesh, control_points[1], control_points[2], control_points[3], false);
      new_ones[3] = geodesic_lerp(
          mesh, control_points[2], control_points[3], 0.75);
    } else if (new_ones_entries.back() == pow(2, k) - 1) {
      new_ones[0] = geodesic_lerp(
          mesh, control_points[0], control_points[1], 0.5);
      new_ones[1] = LR_regular(
          mesh, control_points[0], control_points[1], control_points[2]);
      new_ones[2] = geodesic_lerp(
          mesh, control_points[1], control_points[2], 0.5);
      new_ones[3] = LR_boundary(
          mesh, control_points[1], control_points[2], control_points[3], false);
    }
  }
  return {true, new_ones};
}
bool is_bezier_straight_enough(const geodesic_path& a, const geodesic_path& b,
    const geodesic_path& c, const bezier_mesh& mesh,
    const bezier_params& params) {
  // auto a_positions = path_positions(a, triangles, positions, adjacencies),
  //      b_positions = path_positions(b, triangles, positions, adjacencies),
  //      c_positions = path_positions(c, triangles, positions, adjacencies);
  // auto angle1      = angle(
  //     a_positions.end()[-2] - b_positions[0], b_positions[0] -
  //     b_positions[1]);
  // if (length(a_positions.end()[-2] - b_positions[0]) <= 1e-6 ||
  //     length(b_positions[0] - b_positions[1]) <= 1e-6)
  //   angle1 = 0;

  // if (angle1 > pow(2, -params.precision)) return false;
  // auto angle2 = angle(
  //     b_positions.end()[-2] - c_positions[0], c_positions[0] -
  //     c_positions[1]);
  // if (length(b_positions.end()[-2] - c_positions[0]) <= 1e-6 ||
  //     length(c_positions[0] - c_positions[1]) <= 1e-6)
  //   angle2 = 0;
  // if (angle2 > pow(2, -params.precision)) return false;
  // return true;
  // {
  // On curve apex we may never reach straightess, so we check curve
  // length.
  {
    auto pos  = array<vec3f, 4>{};
    pos[0]    = eval_position(mesh, a.start);
    pos[1]    = eval_position(mesh, b.start);
    pos[2]    = eval_position(mesh, c.start);
    pos[3]    = eval_position(mesh, c.end);
    float len = 0;
    for (int i = 0; i < 3; i++) {
      len += length(pos[i] - pos[i + 1]);
    }
    if (len < params.min_curve_size) return true;
  }

  {
    auto dir0   = tangent_path_direction(mesh, a, false);  // end
    auto dir1   = tangent_path_direction(mesh, b, true);   // start
    auto angle1 = cross(dir0, dir1);
    if (fabs(angle1) > params.precision) {
      // printf("a1: %f > %f\n", angle1, params.precision);
      return false;
    }
    // if (angle1 > params.precision) return false;
  }

  {
    auto dir0   = tangent_path_direction(mesh, b, false);  // end
    auto dir1   = tangent_path_direction(mesh, c, true);   // start
    auto angle1 = angle(dir0, dir1);
    if (fabs(angle1) > params.precision) {
      // printf("a2: %f > %f\n", angle1, params.precision);
      return false;
    }
    // if (angle1 > params.precision) return false;
  }

  return true;
}

// TO DO (CLAUDIO): not optimal
std::tuple<int, vec2f, bezier_segment> find_leaf(const bezier_mesh& mesh,
    const bezier_segment& control_points, const float& t0,
    const bezier_params& params) {
  auto q = vector<mesh_point>(7);

  auto& p   = control_points;
  q[0]      = p[0];
  q[1]      = geodesic_lerp(mesh, p[0], p[1], 1 / 4.f);
  auto p0p1 = geodesic_lerp(mesh, p[0], p[1], 1 / 2.f);
  auto p1p2 = geodesic_lerp(mesh, p[1], p[2], 1 / 2.f);
  q[2]      = geodesic_lerp(mesh, p0p1, p1p2, 1 / 4.f);
  auto p2p3 = geodesic_lerp(mesh, p[2], p[3], 1 / 2.f);
  q[3]      = LR_init(mesh, p0p1, p1p2, p2p3);
  q[4]      = geodesic_lerp(mesh, p1p2, p2p3, 3 / 4.f);
  q[5]      = geodesic_lerp(mesh, p2p3, p[3], 1 / 2.f);
  q[6]      = p[3];

  bezier_segment leaf = {};
  auto           t    = zero2f;
  if (t0 <= 0.25) {
    leaf = {q[0], q[1], q[2], q[3]};
    t    = {0, 0.25};
  } else if (t0 <= 0.5) {
    leaf = {q[1], q[2], q[3], q[4]};
    t    = {0.25, 0.5};
  } else if (t0 <= 0.75) {
    leaf = {q[2], q[3], q[4], q[5]};
    t    = {0.5, 0.75};
  } else if (t0 < 1) {
    leaf = {q[3], q[4], q[5], q[6]};
    t    = {0.75, 1};
  } else
    assert(false);
  auto k = 3;
  while (true) {
    auto curr_t     = (t.x + t.y) / 2;
    int  curr_entry = (int)(pow(2, k) * curr_t);
    curr_entry += 3;
    assert(curr_entry % 2 == 0);
    if (t0 < curr_t) {
      auto curr_entries = vector<int>{
          curr_entry - 4, curr_entry - 3, curr_entry - 2, curr_entry - 1};
      auto [are_boundaries, new_ones] = handle_boundary(
          mesh, leaf, curr_entries, k);
      if (!are_boundaries) {
        new_ones[0] = geodesic_lerp(mesh, leaf[0], leaf[1], 0.5);
        new_ones[1] = LR_regular(mesh, leaf[0], leaf[1], leaf[2]);
        new_ones[2] = geodesic_lerp(mesh, leaf[1], leaf[2], 0.5);
        new_ones[3] = LR_regular(mesh, leaf[1], leaf[2], leaf[3]);
      }
      leaf = new_ones;
      t    = {t.x, curr_t};
    } else {
      auto curr_entries = vector<int>{
          curr_entry - 3, curr_entry - 2, curr_entry - 1, curr_entry};
      auto [are_boundaries, new_ones] = handle_boundary(
          mesh, leaf, curr_entries, k);
      if (!are_boundaries) {
        new_ones[0] = LR_regular(mesh, leaf[0], leaf[1], leaf[2]);
        new_ones[1] = geodesic_lerp(mesh, leaf[1], leaf[2], 0.5);
        new_ones[2] = LR_regular(mesh, leaf[1], leaf[2], leaf[3]);
        new_ones[3] = geodesic_lerp(mesh, leaf[2], leaf[3], 0.5);
      }
      leaf = new_ones;
      t    = {curr_t, t.y};
    }
    auto L0 = compute_geodesic_path(mesh, leaf[0], leaf[1]);
    auto L1 = compute_geodesic_path(mesh, leaf[1], leaf[2]);
    auto L2 = compute_geodesic_path(mesh, leaf[2], leaf[3]);
    if (is_bezier_straight_enough(L0, L1, L2, mesh, params))
      return {k, t, leaf};
    else
      ++k;
  }
}
mesh_point eval_spline_point_cheap(const bezier_mesh& mesh,
    const vector<mesh_point>& control_polygon, const vector<float>& knot_vector,
    const float& t0, const int entry) {
  assert(entry >= 3);
  auto control_points = {control_polygon[entry], control_polygon[entry - 1],
      control_polygon[entry - 2], control_polygon[entry - 3]};
  vector<mesh_point> old_points = control_points;
  vector<mesh_point> curr_points(4);
  auto               offset = entry - 3;
  for (auto j = 1; j <= 3; ++j) {
    for (auto i = entry - 3 + j; i <= entry; ++i) {
      auto alpha = (t0 - knot_vector[i]) /
                   (knot_vector[i - j + 4] - knot_vector[i]);
      curr_points[i - offset] = geodesic_lerp(
          mesh, old_points[i - 1 - offset], old_points[i - offset], alpha);
    }
    old_points = curr_points;
  }
  return curr_points.back();
}
mesh_point eval_spline_point(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    const float& t0) {
  if (t0 <= 0) return polygon[0];
  if (t0 >= 1) return polygon[3];
  auto [k, t, leaf]   = find_leaf(mesh, polygon, t0, params);
  float         step  = 1 / pow(2, k);
  vector<float> knots = {
      yocto::max(t.x - 3 * step, 0.f),
      yocto::max(t.x - 2 * step, 0.f),
      yocto::max(t.x - step, 0.f),
      t.x,
      t.y,
      yocto::min(t.y + step, 1.f),
      yocto::min(t.y + 2 * step, 1.f),
      yocto::min(t.y + 3 * step, 1.f),
  };

  return de_boor(mesh, leaf, knots, t0);
}
std::array<bezier_segment, 2> insert_point(
    const bezier_mesh& mesh, const bezier_segment& polygon, float t0) {
  auto t_start = 0.f;
  auto t_end   = 1.f;
  auto points  = polygon;
  while (true) {
    auto Q0    = geodesic_lerp(mesh, points[0], points[1], 0.5);
    auto Q1    = geodesic_lerp(mesh, points[1], points[2], 0.5);
    auto Q2    = geodesic_lerp(mesh, points[2], points[3], 0.5);
    auto R0    = geodesic_lerp(mesh, Q0, Q1, 0.5);
    auto R1    = geodesic_lerp(mesh, Q1, Q2, 0.5);
    auto S     = geodesic_lerp(mesh, R0, R1, 0.5);
    auto mid_t = (t_start + t_end) / 2.f;
    if (t0 < mid_t) {
      points[1] = Q0;
      points[2] = R0;
      points[3] = S;
      t_end     = mid_t;
    } else {
      points[0] = S;
      points[1] = R1;
      points[2] = Q2;
      t_start   = mid_t;
    }

    if (is_control_polygon_unfoldable(mesh, points)) break;
  }
  // Compute the parameter t local to the leaf control polygon.

  auto tP_local = (t0 - t_start) / (t_end - t_start);

  // Subdivide the leaf control with De Castljeau creating two new control
  // polygons. They are segment_left and segment_right.
  auto [segment_left, segment_right] = subdivide_bezier_polygon(
      mesh, points, tP_local);
  assert(check_segment(segment_left));
  assert(check_segment(segment_right));
  auto left_side = compute_geodesic_path(
      mesh, segment_left.back(), segment_left[2]);
  auto right_side = compute_geodesic_path(
      mesh, segment_right[0], segment_right[1]);
  // P is the inserted mesh point that sepraters segment_left and
  // segment_right.
  assert(segment_left[3] == segment_right[0]);
  auto P = segment_right[0];

  // left part
  {
    if (t0 == t_start) t0 = (t_end - t_start) / 100 + t_start;
    auto Pp2_len = path_length(mesh, left_side);
    auto Pp2_dir = tangent_path_direction(mesh, left_side);
    //    assert(left_leaf.start == P);
    auto delta_len = t_start * Pp2_len /
                     (t0 - t_start);  // Pp2_len * t_start / (t0 - t_start);
    auto path = straightest_path(mesh, P, Pp2_dir, delta_len + Pp2_len);
    auto Pp2  = path.end;

    auto Pp1 = geodesic_lerp(mesh, polygon[0], polygon[1], t0);

    segment_left = {polygon[0], Pp1, Pp2, P};
  }

  // right part
  {
    if (t0 == t_end) t0 = 99 * (t_end - t_start) / 100 + t_start;
    auto Pp1_len   = path_length(mesh, right_side);
    auto Pp1_dir   = tangent_path_direction(mesh, right_side);
    auto delta_len = (1 - t_end) / (t_end - t0) * Pp1_len;
    auto path      = straightest_path(mesh, P, Pp1_dir, delta_len + Pp1_len);
    auto Pp1       = path.end;
    auto Pp2       = geodesic_lerp(mesh, polygon[2], polygon[3], t0);
    segment_right  = {P, Pp1, Pp2, polygon[3]};
  }

  return {segment_left, segment_right};
}

std::array<quadratic_bezier_segment, 2> insert_point(
    const bezier_mesh& mesh, quadratic_bezier_segment& polygon, float t0) {
  quadratic_bezier_segment left = {}, right = {};
  auto                     Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t0);
  auto                     Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t0);
  auto                     Q  = geodesic_lerp(mesh, Q0, Q1, t0);
  left                        = {polygon[0], Q0, Q};
  right                       = {Q, Q1, polygon[2]};
  return {left, right};
}
std::array<bezier_segment, 2> insert_point_spline(const bezier_mesh& mesh,
    const bezier_segment& polygon, const float& t0,
    const bezier_params& params) {
  // Go down the tree and find the leaf node containing the point.
  std::array<bezier_segment, 2> result;

  if (t0 <= 0) {
    result[0] = {polygon[0], polygon[0], polygon[0], polygon[0]};
    result[1] = polygon;

  } else if (t0 >= 1) {
    result[0] = polygon;
    result[1] = {polygon[3], polygon[3], polygon[3], polygon[3]};

  } else {
    auto [depth, t, leaf]        = find_leaf(mesh, polygon, t0, params);
    float         step           = 1 / pow(2, depth);
    vector<float> knots          = {yocto::max(t.x - 3 * step, 0.f),
        yocto::max(t.x - 2 * step, 0.f), yocto::max(t.x - step, 0.f), t.x, t.y,
        yocto::min(t.y + step, 1.f), yocto::min(t.y + 2 * step, 1.f),
        yocto::min(t.y + 3 * step, 1.f)};
    auto          bezier_polygon = from_spline_to_bezier(mesh, knots, leaf);
    auto          t_rel          = (t0 - t.x) / (t.y - t.x);
    auto [left, right] = subdivide_bezier_polygon(mesh, bezier_polygon, t_rel);

    {
      auto L32       = compute_geodesic_path(mesh, left[3], left[2]);
      auto Pp2_len   = path_length(mesh, L32);
      auto Pp2_dir   = tangent_path_direction(mesh, L32);
      auto delta_len = t0 * Pp2_len / (t0 - t.x);
      auto path      = straightest_path(mesh, left[3], Pp2_dir, delta_len);
      auto Pp2       = path.end;

      auto Pp1 = geodesic_lerp(mesh, polygon[0], polygon[1], t0);

      result[0] = {polygon[0], Pp1, Pp2, left[3]};
    }

    // right part
    {
      auto L01      = compute_geodesic_path(mesh, right[0], right[1]);
      auto Pp1_len  = path_length(mesh, L01);
      auto Pp1_dir  = tangent_path_direction(mesh, L01);
      auto t0_local = (t.y - t0) / (1 - t0);
      auto path = straightest_path(mesh, right[0], Pp1_dir, Pp1_len / t0_local);
      auto Pp1  = path.end;
      auto Pp2  = geodesic_lerp(mesh, polygon[2], polygon[3], t.y);
      result[1] = {right[0], Pp1, Pp2, polygon[3]};
    }
  }
  return result;
}
std::array<bezier_segment, 2> insert_point_old(
    const bezier_mesh& mesh, bezier_tree& tree, float t0) {
  // Go down the tree and find the leaf node containing the point.
  int leaf = 0;
  while (true) {
    if (t0 < (tree.nodes[leaf].t_start + tree.nodes[leaf].t_end) / 2) {
      leaf = tree.nodes[leaf].children[0];
    } else {
      leaf = tree.nodes[leaf].children[1];
    }
    if (tree.nodes[leaf].children[0] == -1) break;
  }

  // Compute the parameter t local to the leaf control polygon.
  auto t_start  = tree.nodes[leaf].t_start;
  auto t_end    = tree.nodes[leaf].t_end;
  auto tP_local = (t0 - t_start) / (t_end - t_start);

  // Subdivide the leaf control with De Castljeau creating two new control
  // polygons. They are segment_left and segment_right.
  subdivide_bezier_node(tree, leaf, mesh, tP_local);
  auto& left_leaf     = tree.nodes[tree.nodes[leaf].children[0]];
  auto& right_leaf    = tree.nodes[tree.nodes[leaf].children[1]];
  auto  segment_left  = left_leaf.points;
  auto  segment_right = right_leaf.points;
  assert(check_segment(segment_left));
  assert(check_segment(segment_right));

  // P is the inserted mesh point that sepraters segment_left and
  // segment_right.
  assert(segment_left[3] == segment_right[0]);
  auto P = segment_right[0];

  // Process left part of the curve. Start from segment_left and go up,
  // updating segment_left when climbing up a right child.
  {
    int  node    = leaf;  // tree.nodes[leaf].children[0];
    auto Pp2_len = path_length(mesh, left_leaf.lines[2]);
    auto Pp2_dir = tangent_path_direction(mesh, left_leaf.lines[2]);

    while (true) {
      auto parent = tree.nodes[node].parent;
      if (parent == -1) break;

      if (is_left_child(tree, node)) {
        node = parent;
        continue;
      }

      auto& Nll = tree.nodes[tree.nodes[parent].children[0]];
      assert(Nll.points[3] == segment_left[0]);
      auto t = (Nll.t_end - Nll.t_start) / (t0 - Nll.t_start);

      // continue Nll.lines[0]
      assert(Nll.lines[0].start == Nll.points[0]);
      auto Pp1_dir  = tangent_path_direction(mesh, Nll.lines[0]);
      auto Pp1_len  = path_length(mesh, Nll.lines[0]) / t;
      auto Pp1_line = straightest_path(mesh, Nll.points[0], Pp1_dir, Pp1_len);
      auto Pp1      = Pp1_line.end;

      // continue left_leaf.lines[2]
      assert(left_leaf.lines[2].start == left_leaf.points[3]);
      Pp2_len /= (1 - t);
      auto Pp2_line = straightest_path(
          mesh, left_leaf.points[3], Pp2_dir, Pp2_len);
      auto Pp2 = Pp2_line.end;

      segment_left = {Nll.points[0], Pp1, Pp2, P};
      node         = parent;
    }
  }

  // Process right part of the curve. Start from segment_right and go up,
  // updating segment_right when climbing up a left child.
  {
    int  node    = leaf;  // tree.nodes[leaf].children[0];
    auto Pp1_dir = tangent_path_direction(mesh, right_leaf.lines[0]);
    auto Pp1_len = path_length(mesh, right_leaf.lines[0]);

    while (true) {
      auto parent = tree.nodes[node].parent;
      if (parent == -1) break;

      // Just go up.
      if (!is_left_child(tree, node)) {
        node = parent;
        continue;
      }

      auto& Nrr = tree.nodes[tree.nodes[parent].children[1]];
      assert(Nrr.points[0] == segment_right[3]);
      auto t = (Nrr.t_start - t0) / (Nrr.t_end - t0);

      // continue right_leaf.lines[0]
      assert(right_leaf.lines[0].start == right_leaf.points[0]);
      Pp1_len /= t;
      auto Pp1_line = straightest_path(
          mesh, right_leaf.points[0], Pp1_dir, Pp1_len);
      auto Pp1 = Pp1_line.end;

      // continue Nrr.lines[2]
      assert(Nrr.lines[2].start == Nrr.points[3]);
      auto Pp2_dir  = tangent_path_direction(mesh, Nrr.lines[2]);
      auto Pp2_len  = path_length(mesh, Nrr.lines[2]) / (1 - t);
      auto Pp2_line = straightest_path(mesh, Nrr.points[3], Pp2_dir, Pp2_len);
      auto Pp2      = Pp2_line.end;

      segment_right = {P, Pp1, Pp2, Nrr.points[3]};
      node          = parent;
    }
  }
  return {segment_left, segment_right};
}

void subdivide_bezier_tree(
    const bezier_mesh& mesh, bezier_tree& tree, bezier_params params, float t) {
  // assuming bezier_tree is empty
  auto from = 0;
  auto to   = 1;
  // auto depth = (params.subdivisions - tree.depth);
  auto depth = params.subdivisions;
  for (int i = 1; i < depth; i++) {
    for (int k = from; k < to; k++) {
      subdivide_bezier_node(tree, k, mesh, t);
    }
    auto tmp = to - from;
    from     = to;
    to += 2 * tmp;
  }
  tree.depth = params.subdivisions;
}

void subdivide_bezier_adaptive(const bezier_mesh& mesh,
    const bezier_segment& input, bezier_params params,
    vector<mesh_point>& result, int& count_path, int& count_eval,
    int depth = 0) {
  // resulting beziers: (P0, Q0, R0, S) (S, R1, Q2, P3)

  if (depth > params.max_depth) {
    printf("%s: reached max depth!\n", __FUNCTION__);
    return;
  }
  auto [P0, P1, P2, P3] = input;

  auto P0_P1 = compute_geodesic_path(mesh, P0, P1);
  auto P1_P2 = compute_geodesic_path(mesh, P1, P2);
  auto P2_P3 = compute_geodesic_path(mesh, P2, P3);
  count_path += 3;
  if (is_bezier_straight_enough(P0_P1, P1_P2, P2_P3, mesh, params)) {
    result.push_back(P0);
    result.push_back(P1);
    result.push_back(P2);
    result.push_back(P3);
    return;
  }

  auto Q0 = eval_path_midpoint(
      P0_P1, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q1 = eval_path_midpoint(
      P1_P2, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q2 = eval_path_midpoint(
      P2_P3, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q0_Q1 = compute_geodesic_path(mesh, Q0, Q1);
  auto Q1_Q2 = compute_geodesic_path(mesh, Q1, Q2);

  auto R0 = eval_path_midpoint(
      Q0_Q1, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto R1 = eval_path_midpoint(
      Q1_Q2, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto R0_R1 = compute_geodesic_path(mesh, R0, R1);

  auto S = eval_path_midpoint(
      R0_R1, mesh.triangles, mesh.positions, mesh.adjacencies);

  count_eval += 6;
  count_path += 3;
  subdivide_bezier_adaptive(
      mesh, {P0, Q0, R0, S}, params, result, count_path, count_eval, depth + 1);
  subdivide_bezier_adaptive(
      mesh, {S, R1, Q2, P3}, params, result, count_path, count_eval, depth + 1);
}
void subdivide_bezier_adaptive(const bezier_mesh& mesh,
    const bezier_segment& input, bezier_params params,
    vector<mesh_point>& result, int depth = 0) {
  // resulting beziers: (P0, Q0, R0, S) (S, R1, Q2, P3)

  // if (depth > params.max_depth) {
  //   return;
  // }
  auto [P0, P1, P2, P3] = input;

  auto P0_P1 = compute_geodesic_path(mesh, P0, P1);
  auto P1_P2 = compute_geodesic_path(mesh, P1, P2);
  auto P2_P3 = compute_geodesic_path(mesh, P2, P3);

  if (is_bezier_straight_enough(P0_P1, P1_P2, P2_P3, mesh, params)) {
    result.push_back(P0);
    result.push_back(P1);
    result.push_back(P2);
    result.push_back(P3);
    return;
  }

  auto Q0 = eval_path_midpoint(
      P0_P1, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q1 = eval_path_midpoint(
      P1_P2, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q2 = eval_path_midpoint(
      P2_P3, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto Q0_Q1 = compute_geodesic_path(mesh, Q0, Q1);
  auto Q1_Q2 = compute_geodesic_path(mesh, Q1, Q2);

  auto R0 = eval_path_midpoint(
      Q0_Q1, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto R1 = eval_path_midpoint(
      Q1_Q2, mesh.triangles, mesh.positions, mesh.adjacencies);
  auto R0_R1 = compute_geodesic_path(mesh, R0, R1);

  auto S = eval_path_midpoint(
      R0_R1, mesh.triangles, mesh.positions, mesh.adjacencies);

  subdivide_bezier_adaptive(mesh, {P0, Q0, R0, S}, params, result, depth + 1);
  subdivide_bezier_adaptive(mesh, {S, R1, Q2, P3}, params, result, depth + 1);
}
vector<mesh_point> bezier_adaptive(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params) {
  auto result     = vector<mesh_point>{};
  auto count_path = 0;
  auto count_eval = 0;
  subdivide_bezier_adaptive(mesh, control_points, params, result, 0);

  std::cout << "Number of paths";
  std::cout << count_path << std::endl;
  std::cout << "Number of eval";
  std::cout << count_eval << std::endl;
  std::cout << "Precision";
  std::cout << pow(2, -params.precision) / pif * 180 << std::endl;
  return result;
}

#include <thread>
vector<mesh_point> bezier_uniform_parallel(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params) {
  auto segments = vector<bezier_segment>{control_points};
  auto result   = vector<bezier_segment>();
  auto threads  = vector<std::thread>((int)pow(2, params.subdivisions));
  // printf("num threads: %ld\n", threads.size());

  auto f = [&](int k) {
    auto [split0, split1] = subdivide_bezier_polygon(mesh, segments[k], 0.5);
    result[k * 2]         = split0;
    result[k * 2 + 1]     = split1;
  };

  for (auto subdivision = 0; subdivision < params.subdivisions; subdivision++) {
    result.resize(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      threads[i] = std::thread(f, i);
    }
    for (auto i = 0; i < segments.size(); i++) threads[i].join();

    // for (auto i = 0; i < segments.size(); i++) {
    //   auto [split0, split1] = subdivide_bezier(mesh, segments[i],
    //   params.strip); result[i * 2]         = split0; result[i * 2 + 1] =
    //   split1;
    // }
    swap(segments, result);
  }
  // printf("num segments: %ld\n", segments.size());

  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 4};
}

vector<mesh_point> bezier_uniform(const bezier_mesh& mesh,
    const bezier_segment& control_points, const bezier_params& params) {
  time_function();
  auto segments = vector<bezier_segment>{control_points};
  auto result   = vector<bezier_segment>();
  auto count    = 0;
  for (auto subdivision = 0; subdivision < params.subdivisions; subdivision++) {
    result.clear();
    result.reserve(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      auto [split0, split1] = subdivide_bezier_polygon(mesh, segments[i], 0.5);
      count += 6;
      result.push_back(split0);
      result.push_back(split1);
    }
    swap(segments, result);
  }
  // std::cout << "Number of calls";
  // std::cout << count << std::endl;
  // std::cout << "Level";
  // std::cout << params.subdivisions << std::endl;

  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 4};
}
vector<mesh_point> bezier_uniform(const bezier_mesh& mesh,
    const quadratic_bezier_segment&                  control_points,
    const bezier_params&                             params) {
  auto segments = vector<quadratic_bezier_segment>{control_points};
  auto result   = vector<quadratic_bezier_segment>();
  auto count    = 0;
  for (auto subdivision = 0; subdivision < params.subdivisions; subdivision++) {
    result.clear();
    result.reserve(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      auto [split0, split1] = subdivide_bezier_polygon(mesh, segments[i], 0.5);
      result.push_back(split0);
      result.push_back(split1);
    }
    swap(segments, result);
  }
  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 3};
}

// vector<mesh_point> bezier(const bezier_mesh& mesh,
//     const bezier_segment& control_points, const bezier_params& params) {
//   // profile_function();
//   switch (params.algorithm) {
//     case spline_algorithm::de_casteljau_uniform: {
//       if (params.parallel)
//         return bezier_uniform_parallel(mesh, control_points, params);
//       else
//         return bezier_uniform(mesh, control_points, params);
//     }
//     case spline_algorithm::de_casteljau_adaptive: {
//       return bezier_adaptive(mesh, control_points, params);
//     }
//     case spline_algorithm::subdivision_uniform: {
//       return spline_subdivision_uniform(
//           mesh, control_points, params.subdivisions);
//     }
//     case spline_algorithm::subdivision_adaptive: {
//       return spline_subdivision_adaptive(mesh, control_points, params);
//     }
//   }
// }

// weighted averages (Note:gradients needs to be a vector such that at the
// i-th entry constains the gradient field of the squared distance field from
// the i-th control points)
vector<mesh_point> weighted_average(const bezier_mesh& mesh,
    const vector<mesh_point>& rectangle, const vector<vector<float>>& weights) {
#if HEAVY
  vector<vector<vec3f>> gradients(rectangle.size());
  vector<vector<float>> field(rectangle.size());
  vector<mesh_point>    points_on_surface(weights.size(), {-1, zero2f});
  for (auto i = 0; i < rectangle.size(); ++i) {
    auto f   = compute_distance_field(mesh.triangles, mesh.positions,
          mesh.adjacencies, mesh.v2t, mesh.solver, i, rectangle, graph);
    field[i] = f;
    std::transform(f.begin(), f.end(), f.begin(),
        [](float lambda) { return lambda * lambda; });
    gradients[i] = compute_grad(mesh.solver, mesh.triangles, mesh.positions,
        mesh.normals, mesh.Grad, f);
  }
  // set from as the control points associate to greatest weight
  vector<pair<vector<float>, int>> badones;
  for (auto i = 0; i < weights.size(); ++i) {
    auto from = optimal_seed(rectangle, weights[i]);
    auto p = gradient_descent(mesh.triangles, mesh.positions, mesh.adjacencies,
        mesh.v2t, mesh.normals, mesh.solver, mesh.angles, mesh.total_angles,
        gradients, weights[i], from);
    if (p.face != -1)
      points_on_surface[i] = p;
    else
      badones.push_back({weights[i], i});
  }
  if (points_on_surface.size() == 0) {
    std::cout << "rectangle is too big" << std::endl;
    return {};
  }
  bool stop              = badones.size() == 0;
  bool needs_subdivision = false;
  auto prev_proccessed   = badones.size();
  auto curr_processed    = 0;
  while (!stop) {
    stop = true;
    for (auto i = 0; i < badones.size(); ++i) {
      if (points_on_surface[badones[i].second].face == -1) {
        stop      = false;
        auto from = optimal_seed(points_on_surface, badones[i].first, weights);
        auto p    = gradient_descent(mesh.triangles, mesh.positions,
               mesh.adjacencies, mesh.v2t, mesh.normals, mesh.solver, mesh.angles,
               mesh.total_angles, gradients, badones[i].first, from);
        if (p.face != -1)
          points_on_surface[badones[i].second] = p;
        else
          ++curr_processed;
      }
      if (curr_processed == prev_proccessed) {
        stop              = true;
        needs_subdivision = true;

      } else {
        prev_proccessed = curr_processed;
        curr_processed  = 0;
      }
    }
  }
  // this part doesn't work yet
  if (needs_subdivision) {
    for (auto i = 0; i < badones.size(); ++i) {
      if (points_on_surface[badones[i].second].face == -1) {
        auto from = optimal_seed_with_weights(
            points_on_surface, badones[i].first, weights);
        auto curr_gradients = gradients_inside_cell(mesh.triangles,
            mesh.positions, mesh.adjacencies, mesh.v2t, mesh.solver,
            mesh.angles, mesh.normals, mesh.Grad, field, rectangle, from.first,
            from.second, badones[i].first);
        auto p              = gradient_descent(mesh.triangles, mesh.positions,
                         mesh.adjacencies, mesh.v2t, mesh.normals, mesh.solver, mesh.angles,
                         mesh.total_angles, curr_gradients, badones[i].first, from.first);
        if (p.face != -1) {
          points_on_surface[badones[i].second] = p;
          std::cout << "we subdivide" << std::endl;
        } else {
          std::cout << "rectangle is too big" << std::endl;
        }
      }
    }
  }

  return points_on_surface;
#else
  return {};
#endif
}

pair<bool, vector<mesh_point>> handle_boundary_node(const bezier_mesh& mesh,
    const spline_node& leaf, const vector<int>& new_ones_entries,
    int& count_paths, int& count_eval) {
  vector<mesh_point> new_ones(5);
  auto               k = leaf.depth + 1;
  if (new_ones_entries[0] > 3 && new_ones_entries.back() < pow(2, k) - 1)
    return {false, new_ones};
  else if (new_ones_entries[0] <= 3) {
    if (new_ones_entries[0] == 0) {
      new_ones[0] = leaf.points[0];
      new_ones[1] = eval_path_point(mesh, leaf.lines[0], 0.5);
      new_ones[2] = eval_path_point(mesh, leaf.lines[1], 0.25);
      new_ones[3] = LR_boundary(mesh, leaf.lines[1], leaf.lines[2], true);
      new_ones[4] = eval_path_point(mesh, leaf.lines[2], 0.5);
      count_paths += 1;
      count_eval += 5;
    } else if (new_ones_entries[0] == 2) {
      new_ones[0] = eval_path_point(mesh, leaf.lines[0], 0.25);
      new_ones[1] = LR_boundary(mesh, leaf.lines[0], leaf.lines[1], true);
      new_ones[2] = eval_path_point(mesh, leaf.lines[1], 0.5);
      new_ones[3] = LR_regular(mesh, leaf.lines[1], leaf.lines[2]);
      new_ones[4] = eval_path_point(mesh, leaf.lines[2], 0.5);
      count_paths += 2;
      count_eval += 7;
    } else
      assert(false);
  } else {
    if (new_ones_entries.back() == pow(2, k) + 2) {
      new_ones[0] = eval_path_point(mesh, leaf.lines[0], 0.5);
      new_ones[1] = LR_boundary(mesh, leaf.lines[0], leaf.lines[1], false);
      new_ones[2] = eval_path_point(mesh, leaf.lines[1], 0.75);
      new_ones[3] = eval_path_point(mesh, leaf.lines[2], 0.5);
      new_ones[4] = leaf.points[3];
      count_paths += 1;
      count_eval += 5;
    } else if (new_ones_entries.back() == pow(2, k)) {
      new_ones[0] = eval_path_point(mesh, leaf.lines[0], 0.5);
      new_ones[1] = LR_regular(mesh, leaf.lines[0], leaf.lines[1]);
      new_ones[2] = eval_path_point(mesh, leaf.lines[1], 0.5);
      new_ones[3] = LR_boundary(mesh, leaf.lines[1], leaf.lines[2], false);
      new_ones[4] = eval_path_point(mesh, leaf.lines[2], 0.75);
      count_paths += 2;
      count_eval += 7;
    } else
      assert(false);
  }
  return {true, new_ones};
}
pair<spline_node, spline_node> split_spline_node(const bezier_mesh& mesh,
    const spline_node& leaf, int& count_paths, int& count_eval,
    bool& max_depth_reached) {
  auto curr_t     = (leaf.t.x + leaf.t.y) / 2;
  int  curr_entry = (int)(pow(2, leaf.depth + 1) * curr_t);
  curr_entry += 3;
  if (curr_entry % 2) {
    max_depth_reached = true;
    return {leaf, {}};
  }
  auto curr_entries               = vector<int>{curr_entry - 4, curr_entry - 3,
      curr_entry - 2, curr_entry - 1, curr_entry};
  auto [are_boundaries, new_ones] = handle_boundary_node(
      mesh, leaf, curr_entries, count_paths, count_eval);
  if (!are_boundaries) {
    new_ones[0] = eval_path_point(mesh, leaf.lines[0], 0.5);
    new_ones[1] = LR_regular(mesh, leaf.lines[0], leaf.lines[1]);
    new_ones[2] = eval_path_point(mesh, leaf.lines[1], 0.5);
    new_ones[3] = LR_regular(mesh, leaf.lines[1], leaf.lines[2]);
    new_ones[4] = eval_path_point(mesh, leaf.lines[2], 0.5);
    count_eval += 7;
    count_paths += 2;
  }
  auto        L01 = compute_geodesic_path(mesh, new_ones[0], new_ones[1]);
  auto        L12 = compute_geodesic_path(mesh, new_ones[1], new_ones[2]);
  auto        L23 = compute_geodesic_path(mesh, new_ones[2], new_ones[3]);
  spline_node P0  = {{new_ones[0], new_ones[1], new_ones[2], new_ones[3]},
      {L01, L12, L23}, {leaf.t.x, curr_t}, leaf.depth + 1};
  L01             = compute_geodesic_path(mesh, new_ones[3], new_ones[4]);
  spline_node P1  = {{new_ones[1], new_ones[2], new_ones[3], new_ones[4]},
      {L12, L23, L01}, {curr_t, leaf.t.y}, leaf.depth + 1};
  count_paths += 4;
  return {P0, P1};
}

vector<mesh_point> spline_subdivision_uniform(const bezier_mesh& mesh,
    const bezier_segment& control_points, int num_subdivisions) {
  time_function();
  auto size = 7;
  struct parametric_path {
    geodesic_path path = {};
    vector<float> t    = {};
  };
  parametric_path curr_path   = {};
  parametric_path gamma01     = {};
  parametric_path gamma32     = {};
  auto            prev        = mesh_point{};
  auto            curr        = mesh_point{};
  auto            q           = vector<mesh_point>(size);
  auto            count_paths = 0;
  auto            count_eval  = 0;
  {
    auto& p      = control_points;
    gamma01.path = compute_geodesic_path(mesh, p[0], p[1]);
    gamma01.t    = path_parameters(
           gamma01.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    gamma32.path = compute_geodesic_path(mesh, p[3], p[2]);
    gamma32.t    = path_parameters(
           gamma32.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr_path.path = compute_geodesic_path(mesh, p[1], p[2]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 0.25);
    auto p0p1      = eval_geodesic_path(mesh.triangles, mesh.positions,
             mesh.adjacencies, gamma01.path, gamma01.t, 0.5);
    auto p1p2      = eval_geodesic_path(mesh.triangles, mesh.positions,
             mesh.adjacencies, curr_path.path, curr_path.t, 0.5);
    curr_path.path = compute_geodesic_path(mesh, p0p1, p1p2);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[2] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    auto p2p3 = eval_geodesic_path(mesh.triangles, mesh.positions,
        mesh.adjacencies, gamma32.path, gamma32.t, 0.5);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(mesh, p1p2, p2p3);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 3 / 8.f);
    q[3] = geodesic_lerp(mesh, prev, curr, 0.5);
    q[4] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.75);
    q[5] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma32.path, gamma32.t, 0.25);
    q[6] = p[3];
  }
  count_paths = 6;
  count_eval  = 10;
  auto p      = vector<mesh_point>{};

  for (int subdiv = 0; subdiv < num_subdivisions; subdiv++) {
    std::swap(p, q);

    auto new_size = 2 * size - 3;
    q.resize(new_size);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 1.f / pow(2, 3 + subdiv));
    curr_path.path = compute_geodesic_path(mesh, p[1], p[2]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[2] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(mesh, p[2], p[3]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    q[3] = geodesic_lerp(mesh, prev, curr, 0.5);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.5);
    count_paths += 3;
    count_eval += 6;
    for (int j = 4; j < 2 * size - 8; j += 2) {
      q[j]           = prev;
      prev           = eval_geodesic_path(mesh.triangles, mesh.positions,
                    mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(mesh, p[j / 2 + 1], p[j / 2 + 2]);
      curr_path.t    = path_parameters(
             curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
      curr     = eval_geodesic_path(mesh.triangles, mesh.positions,
              mesh.adjacencies, curr_path.path, curr_path.t, 0.25);
      q[j + 1] = geodesic_lerp(mesh, prev, curr, 1 / 2.f);
      prev     = eval_geodesic_path(mesh.triangles, mesh.positions,
              mesh.adjacencies, curr_path.path, curr_path.t, 0.5);
      count_paths += 2;
      count_eval += 4;
    }
    q[2 * size - 8] = prev;
    {
      auto qq        = &q[new_size - 4];
      auto pp        = &p[size - 4];
      prev           = eval_geodesic_path(mesh.triangles, mesh.positions,
                    mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(mesh, pp[1], pp[2]);
      curr_path.t    = path_parameters(
             curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
      curr  = eval_geodesic_path(mesh.triangles, mesh.positions,
           mesh.adjacencies, curr_path.path, curr_path.t, 3 / 8.f);
      qq[0] = geodesic_lerp(mesh, prev, curr, 0.5);
      qq[1] = eval_geodesic_path(mesh.triangles, mesh.positions,
          mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      qq[2] = eval_geodesic_path(mesh.triangles, mesh.positions,
          mesh.adjacencies, gamma32.path, gamma32.t, 1.f / pow(2, 3 + subdiv));
      qq[3] = pp[3];
      count_paths += 2;
      count_eval += 5;
    }
    size = new_size;
  }
  std::cout << "Number of paths";
  std::cout << count_paths << std::endl;
  std::cout << "Number of eval";
  std::cout << count_eval << std::endl;
  std::cout << "Level";
  std::cout << num_subdivisions + 2 << std::endl;
  return q;
}
vector<mesh_point> OLR(const bezier_mesh& mesh,
    const vector<mesh_point>& control_polygon, const int step) {
  auto  result = vector<mesh_point>{};
  auto& p      = control_polygon;
  if (step == 1) {
    result.resize(5);
    result[0] = p[0];
    for (auto i = 1; i < 4; ++i) {
      result[i] = geodesic_lerp(mesh, p[i - 1], p[i], 0.5);
    }
    result[4] = p.back();
  } else if (step == 2) {
    result.resize(7);
    result[0] = p[0];
    result[1] = geodesic_lerp(mesh, p[0], p[1], 0.5);
    result[2] = geodesic_lerp(mesh, p[1], p[2], 0.25);
    auto p0   = geodesic_lerp(mesh, p[1], p[2], 0.25);
    auto p1   = geodesic_lerp(mesh, p[2], p[3], 0.25);
    result[3] = geodesic_lerp(mesh, p0, p1, 0.75);
    result[4] = geodesic_lerp(mesh, p[2], p[3], 0.75);
    result[5] = geodesic_lerp(mesh, p[3], p[4], 0.5);
    result[6] = p.back();
  } else if (step == 3) {
    result.resize(11);
    auto size  = p.size();
    result[0]  = p[0];
    result[1]  = geodesic_lerp(mesh, p[0], p[1], 0.5);
    result[2]  = geodesic_lerp(mesh, p[1], p[2], 0.25);
    auto prev  = geodesic_lerp(mesh, p[1], p[2], 5 / 8.f);
    auto curr  = geodesic_lerp(mesh, p[2], p[3], 0.25);
    result[3]  = geodesic_lerp(mesh, prev, curr, 0.5);
    result[4]  = geodesic_lerp(mesh, p[2], p[3], 0.5);
    prev       = geodesic_lerp(mesh, p[2], p[3], 0.75);
    curr       = geodesic_lerp(mesh, p[3], p[4], 0.25);
    result[5]  = geodesic_lerp(mesh, prev, curr, 0.5);
    result[6]  = geodesic_lerp(mesh, p[3], p[4], 0.5);
    prev       = geodesic_lerp(mesh, p[3], p[4], 0.75);
    curr       = geodesic_lerp(mesh, p[4], p[5], 3 / 8.f);
    result[7]  = geodesic_lerp(mesh, prev, curr, 0.5);
    result[8]  = geodesic_lerp(mesh, p[4], p[5], 0.75);
    result[9]  = geodesic_lerp(mesh, p[5], p[6], 0.5);
    result[10] = p.back();
  }

  return result;
}
vector<mesh_point> spline_subdivision_uniform(const bezier_mesh& mesh,
    const quadratic_bezier_segment& control_points, int num_subdivisions) {
  auto size = 6;
  struct parametric_path {
    geodesic_path path = {};
    vector<float> t    = {};
  };
  parametric_path curr_path   = {};
  parametric_path gamma01     = {};
  parametric_path gamma21     = {};
  auto            q           = vector<mesh_point>(size);
  auto            count_paths = 0;
  auto            count_eval  = 0;
  {
    auto& p      = control_points;
    gamma01.path = compute_geodesic_path(mesh, p[0], p[1]);
    gamma01.t    = path_parameters(
           gamma01.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    gamma21.path = compute_geodesic_path(mesh, p[2], p[1]);
    gamma21.t    = path_parameters(
           gamma21.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 0.25);
    auto p01 = eval_geodesic_path(mesh.triangles, mesh.positions,
        mesh.adjacencies, gamma01.path, gamma01.t, 0.5);
    auto p12 = eval_geodesic_path(mesh.triangles, mesh.positions,
        mesh.adjacencies, gamma21.path, gamma21.t, 0.5);
    q[4] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma21.path, gamma21.t, 0.25);
    q[5] = p[2];

    curr_path.path = compute_geodesic_path(mesh, p01, p12);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);

    q[2] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    q[3] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.75);
  }
  count_paths    = 3;
  count_eval     = 5;
  auto p         = vector<mesh_point>{};
  curr_path.path = compute_geodesic_path(mesh, q[1], q[2]);
  curr_path.t    = path_parameters(
         curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
  for (int subdiv = 0; subdiv < num_subdivisions; subdiv++) {
    std::swap(p, q);

    auto new_size = 2 * size - 2;
    q.resize(new_size);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 1.f / pow(2, 3 + subdiv));
    curr_path.path = compute_geodesic_path(mesh, p[1], p[2]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);

    count_paths += 1;
    count_eval += 1;
    for (int j = 2; j < 2 * size - 4; j += 2) {
      q[j]           = eval_geodesic_path(mesh.triangles, mesh.positions,
                    mesh.adjacencies, curr_path.path, curr_path.t, 0.25);
      q[j + 1]       = eval_geodesic_path(mesh.triangles, mesh.positions,
                mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(mesh, p[j / 2 + 1], p[j / 2 + 2]);
      curr_path.t    = path_parameters(
             curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
      count_paths += 1;
      count_eval += 2;
    }

    q[2 * size - 4] = eval_geodesic_path(mesh.triangles, mesh.positions,
        mesh.adjacencies, gamma21.path, gamma21.t, 1.f / pow(2, 3 + subdiv));
    q[2 * size - 3] = p.back();
    count_eval += 1;
    size = new_size;
  }
  std::cout << "Number of paths";
  std::cout << count_paths << std::endl;
  std::cout << "Number of eval";
  std::cout << count_eval << std::endl;
  std::cout << "Level";
  std::cout << num_subdivisions + 2 << std::endl;
  return q;
}
vector<mesh_point> spline_subdivision_adaptive(const bezier_mesh& mesh,
    const array<mesh_point, 4>& polygon, const bezier_params& params) {
  auto q = vector<mesh_point>(7);

  auto& p         = polygon;
  q[0]            = p[0];
  q[1]            = geodesic_lerp(mesh, p[0], p[1], 1 / 4.f);
  auto p0p1       = geodesic_lerp(mesh, p[0], p[1], 1 / 2.f);
  auto p1p2       = geodesic_lerp(mesh, p[1], p[2], 1 / 2.f);
  q[2]            = geodesic_lerp(mesh, p0p1, p1p2, 1 / 4.f);
  auto p2p3       = geodesic_lerp(mesh, p[2], p[3], 1 / 2.f);
  q[3]            = LR_init(mesh, p0p1, p1p2, p2p3);
  q[4]            = geodesic_lerp(mesh, p1p2, p2p3, 3 / 4.f);
  q[5]            = geodesic_lerp(mesh, p2p3, p[3], 1 / 2.f);
  q[6]            = p[3];
  auto        L01 = compute_geodesic_path(mesh, q[0], q[1]);
  auto        L12 = compute_geodesic_path(mesh, q[1], q[2]);
  auto        L23 = compute_geodesic_path(mesh, q[2], q[3]);
  spline_node P0  = {{q[0], q[1], q[2], q[3]}, {L01, L12, L23}, {0, 0.25}, 2};
  L01             = compute_geodesic_path(mesh, q[3], q[4]);
  spline_node P1  = {{q[1], q[2], q[3], q[4]}, {L12, L23, L01}, {0.25, 0.5}, 2};
  L12             = compute_geodesic_path(mesh, q[4], q[5]);
  spline_node P2  = {{q[2], q[3], q[4], q[5]}, {L23, L01, L12}, {0.5, 0.75}, 2};
  L23             = compute_geodesic_path(mesh, q[5], q[6]);
  spline_node P3  = {{q[3], q[4], q[5], q[6]}, {L01, L12, L23}, {0.75, 1}, 2};

  P0.is_good = is_bezier_straight_enough(
      P0.lines[0], P0.lines[1], P0.lines[2], mesh, params);
  P1.is_good = is_bezier_straight_enough(
      P1.lines[0], P1.lines[1], P1.lines[2], mesh, params);
  P2.is_good = is_bezier_straight_enough(
      P2.lines[0], P2.lines[1], P2.lines[2], mesh, params);
  P3.is_good = is_bezier_straight_enough(
      P3.lines[0], P3.lines[1], P3.lines[2], mesh, params);
  auto                    count_path = 16;
  auto                    count_eval = 8;
  std::deque<spline_node> Q;
  Q.push_back(P3);
  Q.push_back(P2);
  Q.push_back(P1);
  Q.push_back(P0);
  auto P                 = vector<spline_node>{};
  bool max_depth_reached = false;
  while (!Q.empty()) {
    auto curr = Q.back();
    Q.pop_back();
    if (P.size() > 0 && P.back().depth == curr.depth && curr.is_good) {
      P.push_back(curr);

    } else {
      auto [left, right] = split_spline_node(
          mesh, curr, count_path, count_eval, max_depth_reached);
      if (max_depth_reached) {
        curr.is_good = true;
        if (P.size() > 0) {
          auto last = P.back();
          auto L = compute_geodesic_path(mesh, last.points[0], curr.points[0]);
          if (is_bezier_straight_enough(
                  L, curr.lines[0], curr.lines[1], mesh, params))
            P.push_back(curr);
          else {
            P.pop_back();
            last = {{last.points[0], curr.points[0], curr.points[1],
                        curr.points[2]},
                {L, curr.lines[0], curr.lines[1]}, last.t, last.depth, false};
            Q.push_back(curr);
            Q.push_back(last);
          }
        } else
          P.push_back(curr);

        max_depth_reached = false;
      } else {
        left.is_good = is_bezier_straight_enough(
            left.lines[0], left.lines[1], left.lines[2], mesh, params);
        right.is_good = is_bezier_straight_enough(
            right.lines[0], right.lines[1], right.lines[2], mesh, params);
        if (left.is_good && right.is_good) {
          if (P.size() == 0) {
            P.push_back(left);
            P.push_back(right);
          } else {
            auto last = P.back();
            auto L    = compute_geodesic_path(
                   mesh, last.points[0], left.points[0]);
            if (is_bezier_straight_enough(
                    L, left.lines[0], left.lines[1], mesh, params)) {
              last     = {{last.points[0], left.points[0], left.points[1],
                          left.points[2]},
                  {L, left.lines[0], left.lines[1]}, last.t, last.depth, true};
              P.back() = last;
              P.push_back(left);
              P.push_back(right);
            } else if (left.depth < last.depth) {
              left.is_good  = false;
              right.is_good = false;
              Q.push_back(right);
              Q.push_back(left);
            } else {
              P.pop_back();
              last = {{last.points[0], left.points[0], left.points[1],
                          left.points[2]},
                  {L, left.lines[0], left.lines[1]}, last.t, last.depth, false};
              Q.push_back(right);
              Q.push_back(left);
              Q.push_back(last);
            }
          }
        } else {
          Q.push_back(right);
          Q.push_back(left);
        }
      }
    }
  }
  std::cout << "Number of paths";
  std::cout << count_path << std::endl;
  std::cout << "Number of eval";
  std::cout << count_eval << std::endl;
  std::cout << "Precision";
  std::cout << pow(2, -params.precision) / pif * 180 << std::endl;
  auto polyline = vector<mesh_point>{};
  for (auto i = 0; i < P.size(); ++i) {
    if (i == 0) {
      for (auto j = 0; j < 4; ++j) {
        polyline.push_back(P[i].points[j]);
      }
    } else
      polyline.push_back(P[i].points.back());
  }

  return polyline;
}

vector<mesh_point> degree_elevation(const bezier_mesh& mesh,
    const array<mesh_point, 4>& control_points, int num_subdivisions) {
  auto p = vector<mesh_point>{control_points[0], control_points[1],
      control_points[2], control_points[3]};
  auto q = p;

  for (auto i = 0; i < num_subdivisions; ++i) {
    std::swap(p, q);
    auto new_size = p.size() + 1;
    q.resize(new_size);
    q[0] = p[0];
    for (auto j = 1; j < new_size - 1; ++j) {
      float alpha = (new_size - 1 - j) / (float)(new_size - 1);
      q[j]        = geodesic_lerp(mesh, p[j - 1], p[j], alpha);
    }
    q.back() = p.back();
  }

  return q;
}
bool samples_are_near(const bezier_mesh& mesh, const mesh_point& sample_start,
    const mesh_point& sample_to) {
  int pid_from = sample_start.face;
  int pid_to   = sample_to.face;
  int vid      = common_vertex(mesh.triangles, pid_from, pid_to);
  if (vid != -1) return true;
  if (pid_from == pid_to) return true;
  int other_prev = -1, other_next = -1;
  auto [from_is_edge, kfrom] = point_is_edge(sample_start);
  auto [to_is_edge, kto]     = point_is_edge(sample_to);
  auto [from_is_vert, hfrom] = point_is_vert(sample_start);
  auto [to_is_vert, hto]     = point_is_vert(sample_to);
  if (from_is_edge) {
    other_prev = mesh.adjacencies[pid_from][kfrom];
    if (other_prev == pid_to) return true;
    // vid = common_vertex(mesh.triangles, other_prev, pid_to);
    // if (vid != -1) return true;
  }

  if (to_is_edge) {
    other_next = mesh.adjacencies[pid_to][kto];
    if (other_next == pid_from) return true;
    // vid = common_vertex(mesh.triangles, pid_from, other_next);
    // if (vid != -1) return true;
  }

  // if (other_next != -1 && other_prev != -1) {
  //   vid = common_vertex(mesh.triangles, other_prev, other_next);
  //   if (vid != -1) return true;
  // }

  if (from_is_vert) {
    int  vid_from = mesh.triangles[pid_from][hfrom];
    auto nbr      = mesh.v2t[vid_from];
    for (int pid : nbr) {
      // vid = common_vertex(mesh.triangles, pid, pid_to);
      // if (vid != -1) return true;
      if (pid == pid_to || pid == other_next) return true;
      // if (other_next != -1) {
      //   vid = common_vertex(mesh.triangles, pid, other_next);
      //   if (vid != -1) return true;
      // }
    }
  }

  if (to_is_vert) {
    int  vid_to = mesh.triangles[pid_to][hto];
    auto nbr    = mesh.v2t[vid_to];
    for (int pid : nbr) {
      // vid = common_vertex(mesh.triangles, pid, pid_from);
      // if (vid != -1) return true;
      if (pid == pid_from || pid == other_prev) return true;
      // if (other_prev != -1) {
      //   vid = common_vertex(mesh.triangles, pid, other_prev);
      //   if (vid != -1) return true;
      // }
    }
  }

  return false;
}
mesh_point dc_point(const bezier_mesh& mesh, const geodesic_path& L0,
    const geodesic_path& L1, const geodesic_path& L2, const float& t) {
  auto p01 = eval_path_point(mesh, L0, t);
  auto p12 = eval_path_point(mesh, L1, t);
  auto p23 = eval_path_point(mesh, L2, t);

  auto L01 = compute_geodesic_path(mesh, p01, p12);
  auto L12 = compute_geodesic_path(mesh, p12, p23);

  auto a01 = eval_path_point(mesh, L01, t);
  auto a12 = eval_path_point(mesh, L12, t);

  auto L = compute_geodesic_path(mesh, a01, a12);
  return eval_path_point(mesh, L, t);
}
bool is_a_jump(const bezier_mesh& mesh, const geodesic_path& L0,
    const geodesic_path& L1, const geodesic_path& L2, const mesh_point& curr,
    const float& curr_t, const bezier_params& params) {
  auto   subdiv      = params.subdivisions;
  double step        = 1 / pow(2, 12);
  auto   prev_step   = 1 / pow(2, params.subdivisions);
  auto   prev_t      = curr_t - step;
  auto   curr_pos    = eval_position(mesh, curr);
  auto   prev_pos    = zero3f;
  auto   prev_sample = mesh_point{};
  while (curr_t - prev_t < prev_step) {
    prev_sample = dc_point(mesh, L0, L1, L2, prev_t);
    // auto path = compute_geodesic_path(mesh, curr, prev);
    // auto len  = path_length(mesh, path);
    // if (distances.size() > 0 && len < distances.back())
    //   return true;
    // else {
    prev_pos = eval_position(mesh, prev_sample);
    auto len = length(prev_pos - curr_pos);
    if (len > mesh.avg_edge_length)
      return true;
    else {
      prev_t -= step;
      curr_pos = prev_pos;
    }
  }
}
vector<mesh_point> de_casteljau_classic(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    vector<int>& badones, const bool jumps) {
  auto L0 = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  auto L1 = compute_geodesic_path(mesh, polygon[1], polygon[2]);
  auto L2 = compute_geodesic_path(mesh, polygon[2], polygon[3]);

  mesh_point prev_sample;
  badones.clear();
  double             step = 1 / pow(2, params.subdivisions);
  auto               t    = step;
  vector<mesh_point> polyline;
  polyline.push_back(polygon[0]);
  geodesic_path L01, L12, L;
  for (int i = 1; i < pow(2, params.subdivisions); ++i) {
    auto p01 = eval_path_point(mesh, L0, t);
    auto p12 = eval_path_point(mesh, L1, t);
    auto p23 = eval_path_point(mesh, L2, t);

    L01 = compute_geodesic_path(mesh, p01, p12);
    L12 = compute_geodesic_path(mesh, p12, p23);

    auto a01 = eval_path_point(mesh, L01, t);
    auto a12 = eval_path_point(mesh, L12, t);

    L                = compute_geodesic_path(mesh, a01, a12);
    auto curr_sample = eval_path_point(mesh, L, t);
    prev_sample      = polyline[i - 1];
    // if (jumps) {
    //   if (is_a_jump(mesh, L0, L1, L2, curr_sample, t, params))
    //     badones.push_back(i);
    // }

    polyline.push_back(curr_sample);
    t += step;
  }
  polyline.push_back(polygon[3]);
  return polyline;
}
std::pair<vector<geodesic_path>, vector<vector<mesh_point>>>
dc_classic_construction(
    const bezier_mesh& mesh, const bezier_segment& polygon, const float& t) {
  vector<geodesic_path>      paths(6);
  vector<vector<mesh_point>> points(4);
  points[0] = {polygon[0], polygon[1], polygon[2], polygon[3]};

  paths[0] = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  paths[1] = compute_geodesic_path(mesh, polygon[1], polygon[2]);
  paths[2] = compute_geodesic_path(mesh, polygon[2], polygon[3]);

  auto p01  = eval_path_point(mesh, paths[0], t);
  auto p12  = eval_path_point(mesh, paths[1], t);
  auto p23  = eval_path_point(mesh, paths[2], t);
  points[1] = {p01, p12, p23};

  paths[3] = compute_geodesic_path(mesh, p01, p12);
  paths[4] = compute_geodesic_path(mesh, p12, p23);

  auto p012 = eval_path_point(mesh, paths[3], t);
  auto p123 = eval_path_point(mesh, paths[4], t);
  points[2] = {p012, p123};

  paths[5] = compute_geodesic_path(mesh, p012, p123);

  points[3] = {eval_path_point(mesh, paths[5], t)};

  return {paths, points};
}
// std::pair<vector<geodesic_path>, vector<vector<mesh_point>>>
// dc_construction(
//     bezier_mesh& mesh, const bezier_tree& tree, const float& t0) {
//   auto                       leaf   = 0;
//   vector<geodesic_path>      paths  = {tree.nodes[leaf].lines[0],
//       tree.nodes[leaf].lines[1], tree.nodes[leaf].lines[1]};
//   vector<vector<mesh_point>> points = {
//       {tree.nodes[leaf].points[0], tree.nodes[leaf].points[1],
//           tree.nodes[leaf].points[2], tree.nodes[leaf].points[3]}};

//   while (true) {
//     if (t0 < (tree.nodes[leaf].t_start + tree.nodes[leaf].t_end) / 2) {
//       leaf = tree.nodes[leaf].children[0];
//     } else {
//       leaf = tree.nodes[leaf].children[1];
//     }
//     paths.push_back(tree.nodes[leaf].lines[0]);
//     paths.push_back(tree.nodes[leaf].lines[1]);
//     paths.push_back(tree.nodes[leaf].lines[2]);
//     points.push_back({tree.nodes[leaf].points[0],
//     tree.nodes[leaf].points[1],
//         tree.nodes[leaf].points[2], tree.nodes[leaf].points[3]});
//     if (tree.nodes[leaf].children[0] == -1) break;
//   }
// }
std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
dc_construction(
    bezier_mesh& mesh, const bezier_segment& polygon, const float& t0) {
  auto                          points         = polygon;
  auto                          t_start        = 0.f;
  auto                          t_end          = 1.f;
  vector<vector<geodesic_path>> paths          = {};
  vector<vector<mesh_point>>    control_points = {};
  // Recursively subdivide, shrinking the active control polygon around the
  // point of interest, until it is unfoldable.
  while (true) {
    auto L0 = compute_geodesic_path(mesh, points[0], points[1]);
    auto L1 = compute_geodesic_path(mesh, points[1], points[2]);
    auto L2 = compute_geodesic_path(mesh, points[2], points[3]);
    paths.push_back({L0, L1, L2});
    control_points.push_back({points[0], points[1], points[2], points[3]});
    auto Q0 = eval_path_midpoint(mesh, L0);
    auto Q1 = eval_path_midpoint(mesh, L1);
    auto Q2 = eval_path_midpoint(mesh, L2);

    auto R0    = geodesic_lerp(mesh, Q0, Q1, 0.5);
    auto R1    = geodesic_lerp(mesh, Q1, Q2, 0.5);
    auto S     = geodesic_lerp(mesh, R0, R1, 0.5);
    auto mid_t = (t_start + t_end) / 2;
    if (t0 < mid_t) {
      points[1] = Q0;
      points[2] = R0;
      points[3] = S;
      t_end     = mid_t;
    } else {
      points[0] = S;
      points[1] = R1;
      points[2] = Q2;
      t_start   = mid_t;
    }

    if (is_control_polygon_unfoldable(mesh, points)) break;
  }

  // Exact evaulation.
  float t = (t0 - t_start) / (t_end - t_start);
  auto  p = eval_bezier_point_cheap(mesh, points, t);
  control_points.push_back({p});
  return std::make_pair(paths, control_points);
}
std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
RDC_algorithm(bezier_mesh& mesh, const bezier_segment& control_points) {
  vector<vector<geodesic_path>> paths(3);
  vector<vector<mesh_point>>    points(3);

  auto L0   = compute_geodesic_path(mesh, control_points[0], control_points[1]);
  auto L1   = compute_geodesic_path(mesh, control_points[1], control_points[2]);
  auto L2   = compute_geodesic_path(mesh, control_points[2], control_points[3]);
  paths[0]  = {L0, L1, L2};
  points[0] = {control_points[0], control_points[1], control_points[2],
      control_points[3]};
  auto Q0   = eval_path_midpoint(mesh, L0);
  auto Q1   = eval_path_midpoint(mesh, L1);
  auto Q2   = eval_path_midpoint(mesh, L2);
  auto R0   = geodesic_lerp(mesh, Q0, Q1, 0.5);
  auto R1   = geodesic_lerp(mesh, Q1, Q2, 0.5);
  auto S    = geodesic_lerp(mesh, R0, R1, 0.5);
  paths[1].resize(3);
  paths[1][0] = compute_geodesic_path(mesh, control_points[0], Q0);
  paths[1][1] = compute_geodesic_path(mesh, Q0, R0);
  paths[1][2] = compute_geodesic_path(mesh, R0, S);
  points[1]   = {control_points[0], Q0, R0, S};

  paths[2].resize(3);
  paths[2][0] = compute_geodesic_path(mesh, S, R1);
  paths[2][1] = compute_geodesic_path(mesh, R1, Q2);
  paths[2][2] = compute_geodesic_path(mesh, Q2, control_points[3]);
  points[2]   = {S, R1, Q2, control_points[3]};

  return {paths, points};
}
std::pair<vector<vector<geodesic_path>>, vector<vector<mesh_point>>>
LR_algorithm(bezier_mesh& mesh, const bezier_segment& control_points) {
  auto size = 7;
  struct parametric_path {
    geodesic_path path = {};
    vector<float> t    = {};
  };
  parametric_path curr_path        = {};
  parametric_path gamma01          = {};
  parametric_path gamma32          = {};
  auto            prev             = mesh_point{};
  auto            curr             = mesh_point{};
  auto            q                = vector<mesh_point>(size);
  auto            num_subdivisions = 1;

  {
    auto& p      = control_points;
    gamma01.path = compute_geodesic_path(mesh, p[0], p[1]);
    gamma01.t    = path_parameters(
           gamma01.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    gamma32.path = compute_geodesic_path(mesh, p[3], p[2]);
    gamma32.t    = path_parameters(
           gamma32.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr_path.path = compute_geodesic_path(mesh, p[1], p[2]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 0.25);
    auto p0p1      = eval_geodesic_path(mesh.triangles, mesh.positions,
             mesh.adjacencies, gamma01.path, gamma01.t, 0.5);
    auto p1p2      = eval_geodesic_path(mesh.triangles, mesh.positions,
             mesh.adjacencies, curr_path.path, curr_path.t, 0.5);
    curr_path.path = compute_geodesic_path(mesh, p0p1, p1p2);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[2] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    auto p2p3 = eval_geodesic_path(mesh.triangles, mesh.positions,
        mesh.adjacencies, gamma32.path, gamma32.t, 0.5);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(mesh, p1p2, p2p3);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 3 / 8.f);
    q[3] = geodesic_lerp(mesh, prev, curr, 0.5);
    q[4] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.75);
    q[5] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma32.path, gamma32.t, 0.25);
    q[6] = p[3];
  }
  vector<mesh_point>            p;
  vector<vector<geodesic_path>> paths(2);
  vector<vector<mesh_point>>    points(2);
  points[0] = q;
  paths[0].resize(size - 1);
  for (auto i = 0; i < size - 1; ++i) {
    paths[0][i] = compute_geodesic_path(mesh, q[i], q[i + 1]);
  }
  for (int subdiv = 0; subdiv < num_subdivisions; subdiv++) {
    std::swap(p, q);

    auto new_size = 2 * size - 3;
    q.resize(new_size);
    q[0] = p[0];
    q[1] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        gamma01.path, gamma01.t, 1.f / pow(2, 3 + subdiv));
    curr_path.path = compute_geodesic_path(mesh, p[1], p[2]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    q[2] = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(mesh, p[2], p[3]);
    curr_path.t    = path_parameters(
           curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
    curr = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.25);
    q[3] = geodesic_lerp(mesh, prev, curr, 0.5);
    prev = eval_geodesic_path(mesh.triangles, mesh.positions, mesh.adjacencies,
        curr_path.path, curr_path.t, 0.5);

    for (int j = 4; j < 2 * size - 8; j += 2) {
      q[j]           = prev;
      prev           = eval_geodesic_path(mesh.triangles, mesh.positions,
                    mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(mesh, p[j / 2 + 1], p[j / 2 + 2]);
      curr_path.t    = path_parameters(
             curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
      curr     = eval_geodesic_path(mesh.triangles, mesh.positions,
              mesh.adjacencies, curr_path.path, curr_path.t, 0.25);
      q[j + 1] = geodesic_lerp(mesh, prev, curr, 1 / 2.f);
      prev     = eval_geodesic_path(mesh.triangles, mesh.positions,
              mesh.adjacencies, curr_path.path, curr_path.t, 0.5);
    }
    q[2 * size - 8] = prev;
    {
      auto qq        = &q[new_size - 4];
      auto pp        = &p[size - 4];
      prev           = eval_geodesic_path(mesh.triangles, mesh.positions,
                    mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(mesh, pp[1], pp[2]);
      curr_path.t    = path_parameters(
             curr_path.path, mesh.triangles, mesh.positions, mesh.adjacencies);
      curr  = eval_geodesic_path(mesh.triangles, mesh.positions,
           mesh.adjacencies, curr_path.path, curr_path.t, 3 / 8.f);
      qq[0] = geodesic_lerp(mesh, prev, curr, 0.5);
      qq[1] = eval_geodesic_path(mesh.triangles, mesh.positions,
          mesh.adjacencies, curr_path.path, curr_path.t, 0.75);
      qq[2] = eval_geodesic_path(mesh.triangles, mesh.positions,
          mesh.adjacencies, gamma32.path, gamma32.t, 1.f / pow(2, 3 + subdiv));
      qq[3] = pp[3];
    }
    paths[subdiv + 1].resize(new_size - 1);
    points[subdiv + 1] = q;
    for (auto i = 0; i < new_size - 1; ++i) {
      paths[subdiv + 1][i] = compute_geodesic_path(mesh, q[i], q[i + 1]);
    }
    size = new_size;
  }
  return {paths, points};
}
float u_coordinate(
    const vec2i& eid, const vec3f& pos, const vector<vec3f>& positions) {
  vec3f x           = positions[eid.x];
  float edge_length = length((positions[eid.y] - x));
  float u           = length(pos - x);

  return u / edge_length;
}

geodesic_path straightest_geodesic_biermann(const bezier_mesh& mesh,
    const mesh_point& start, const vec2f& dir, const float& path_len) {
  auto path  = geodesic_path{};
  path.start = start;
  path.strip.push_back(start.face);
  auto [is_edge, ke] = point_is_edge(start);
  auto [is_vert, kv] = point_is_vert(start);
  auto next_bary     = zero3f;
  auto curr_bary  = vec3f{1 - start.uv.x - start.uv.y, start.uv.x, start.uv.y};
  auto curr_pos   = eval_position(mesh, start);
  auto next_pos   = zero3f;
  auto curr_point = start;
  auto point_normal = zero3f;
  auto tid_normal   = zero3f;
  auto curr_tri     = start.face;
  auto next_tri     = -1;
  auto len          = 0.f;
  auto dir3d        = normalize(
             from_2d_to_3d_vector(mesh.triangles, mesh.positions, dir, start.face));

  if (is_edge) {
    auto edge = mesh.positions[mesh.triangles[curr_tri][ke]] -
                mesh.positions[mesh.triangles[curr_tri][(ke + 1) % 3]];
    if (dot(cross(edge, dir3d), tid_normal) < 0) {
      curr_point.face = mesh.adjacencies[curr_point.face][ke];
      curr_bary       = tri_bary_coords(
                mesh.positions[mesh.triangles[curr_point.face].x],
                mesh.positions[mesh.triangles[curr_point.face].y],
                mesh.positions[mesh.triangles[curr_point.face].z], curr_pos);
      tid_normal = triangle_normal(
          mesh.positions[mesh.triangles[curr_point.face].x],
          mesh.positions[mesh.triangles[curr_point.face].y],
          mesh.positions[mesh.triangles[curr_point.face].z]);
      curr_point.uv = {curr_bary.y, curr_bary.z};
      auto nextdir  = normalize(cross(point_normal, tid_normal));
      if (dot(nextdir, dir3d) < 0) nextdir *= -1;
      if (length(nextdir) > 0) dir3d = normalize(nextdir);
    }

  } else if (is_vert) {
    next_tri      = next_tid(mesh.solver, mesh.angles, mesh.positions, mesh.v2t,
             mesh.triangles, mesh.normals, mesh.triangles[curr_point.face][kv],
             dir3d);
    curr_bary     = tri_bary_coords(mesh.positions[mesh.triangles[next_tri].x],
            mesh.positions[mesh.triangles[next_tri].y],
            mesh.positions[mesh.triangles[next_tri].z], curr_pos);
    tid_normal    = triangle_normal(mesh.positions[mesh.triangles[next_tri].x],
           mesh.positions[mesh.triangles[next_tri].y],
           mesh.positions[mesh.triangles[next_tri].z]);
    curr_point.uv = {curr_bary.y, curr_bary.z};
  }
  next_tri     = curr_point.face;
  point_normal = normalize(
      curr_bary.x * mesh.normals[mesh.triangles[next_tri].x] +
      curr_bary.y * mesh.normals[mesh.triangles[next_tri].y] +
      curr_bary.z * mesh.normals[mesh.triangles[next_tri].z]);

  while (len < path_len) {
    tid_normal = triangle_normal(mesh.positions[mesh.triangles[next_tri].x],
        mesh.positions[mesh.triangles[next_tri].y],
        mesh.positions[mesh.triangles[next_tri].z]);
    curr_tri   = next_tri;

    auto nextdir = normalize(cross(point_normal, tid_normal));
    if (dot(nextdir, dir3d) < 0) nextdir *= -1;
    if (length(nextdir) > 0) dir3d = normalize(nextdir);
    trace_in_triangles(mesh.positions, mesh.triangles, dir3d, curr_bary,
        curr_tri, next_pos, next_bary);
    std::tie(is_edge, ke) = bary_is_edge(next_bary);
    std::tie(is_vert, kv) = bary_is_vert(next_bary);

    path.strip.push_back(curr_tri);
    len += length(next_pos - curr_pos);
    if (len < path_len) {
      curr_pos = next_pos;
      if (is_edge) {
        point_normal = normalize(
            next_bary[ke] * mesh.normals[mesh.triangles[curr_tri][ke]] +
            next_bary[(ke + 1) % 3] *
                mesh.normals[mesh.triangles[curr_tri][(ke + 1) % 3]]);
        next_tri = mesh.adjacencies[curr_tri][ke];
        dir3d    = project_vec(dir3d, point_normal);
        auto eid = vec2i{mesh.triangles[curr_tri][ke],
            mesh.triangles[curr_tri][(ke + 1) % 3]};
        path.lerps.push_back(u_coordinate(eid, next_pos, mesh.positions));
        next_bary = tri_bary_coords(mesh.positions[mesh.triangles[next_tri].x],
            mesh.positions[mesh.triangles[next_tri].y],
            mesh.positions[mesh.triangles[next_tri].z], next_pos);
      } else if (is_vert) {
        point_normal = mesh.normals[mesh.triangles[curr_tri][kv]];
        next_tri = next_tid(mesh.solver, mesh.angles, mesh.positions, mesh.v2t,
            mesh.triangles, mesh.normals, mesh.triangles[curr_tri][kv], dir3d);
        dir3d    = project_vec(dir3d, point_normal);
        auto vid = mesh.triangles[curr_tri][kv];
        next_pos = mesh.positions[vid];

        auto eid = common_edge(mesh.triangles, curr_tri, next_tri);
        if (eid.x != -1)
          path.lerps.push_back(u_coordinate(eid, next_pos, mesh.positions));
        else {
          path.lerps.push_back(0);
          auto poly_aux = mesh.adjacencies[curr_tri][kv];
          while (poly_aux != next_tri) {
            path.lerps.push_back(0);
            path.strip.push_back(poly_aux);
            kv = find_in_vector(mesh.triangles[poly_aux], vid);
            assert(kv != -1);
            poly_aux = mesh.adjacencies[poly_aux][kv];
          }
        }
        kv = find_in_vector(mesh.triangles[next_tri], vid);
        assert(kv != -1);
        next_bary     = zero3f;
        next_bary[kv] = 1.0;
      } else
        assert(false);
      curr_bary = next_bary;
    }
  }

  auto factor = (len - path_len);
  auto w      = normalize(curr_pos - next_pos);
  w *= factor;
  w += next_pos;
  curr_bary = tri_bary_coords(mesh.positions[mesh.triangles[curr_tri].x],
      mesh.positions[mesh.triangles[curr_tri].y],
      mesh.positions[mesh.triangles[curr_tri].z], w);
  path.lerps.pop_back();
  path.end = {curr_tri, {curr_bary.y, curr_bary.z}};
  check_point(path.end);
  return path;
}
vector<mesh_point> de_casteljau_classic_old(const bezier_mesh& mesh,
    const bezier_segment& polygon, const bezier_params& params,
    vector<int>& badones, const bool jumps) {
  auto L0 = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  auto L1 = compute_geodesic_path(mesh, polygon[1], polygon[2]);
  auto L2 = compute_geodesic_path(mesh, polygon[2], polygon[3]);

  mesh_point prev_sample;
  badones.clear();
  vector<pair<double, bool>> t = {
      std::make_pair(0, false), std::make_pair(1, true)};
  double             step  = 0;
  int                count = 0;
  vector<mesh_point> polyline;
  polyline.push_back(polygon[0]);
  polyline.push_back(polygon[3]);

  bool stop = false;
  struct added_sample {
    int        entry;
    mesh_point sample;
    double     t;
  };

  vector<added_sample> added       = {};
  int                  count_added = 0;
  int                  entry0, entry1, entry2;
  vector<mesh_point>   exit_vertices;
  geodesic_path        L01, L12, L;
  while (!stop) {
    ++count;
    step = pow(0.5, count);
    added.resize(0);
    for (int i = 0; i < t.size() - 1; ++i) {
      if (!t[i].second) {
        double curr_t = t[i].first + step;
        auto   p01    = eval_path_point(mesh, L0, curr_t);
        auto   p12    = eval_path_point(mesh, L1, curr_t);
        auto   p23    = eval_path_point(mesh, L2, curr_t);

        L01 = compute_geodesic_path(mesh, p01, p12);
        L12 = compute_geodesic_path(mesh, p12, p23);

        auto a01 = eval_path_point(mesh, L01, curr_t);
        auto a12 = eval_path_point(mesh, L12, curr_t);

        L                       = compute_geodesic_path(mesh, a01, a12);
        auto curr_sample        = eval_path_point(mesh, L, curr_t);
        prev_sample             = polyline[i];
        added_sample curr_added = {i + 1, curr_sample, curr_t};
        added.push_back(curr_added);

        if (samples_are_near(mesh, prev_sample, curr_sample))
          t[i].second = true;
        else if (step <= 1e-4) {
          t[i].second = true;
          badones.push_back(i + 1);
        }
      }
    }
    if (added.size() == 0) stop = true;
    count_added = 0;
    for (int i = 0; i < added.size(); ++i) {
      added_sample curr                 = added[i];
      int          entry_of_last_sample = curr.entry + count_added;

      polyline.insert(polyline.begin() + entry_of_last_sample, curr.sample);
      auto last_sample = polyline[entry_of_last_sample];
      auto next_sample = polyline[entry_of_last_sample + 1];

      bool next_is_connected = samples_are_near(mesh, last_sample, next_sample);
      if (!next_is_connected)
        t.insert(
            t.begin() + entry_of_last_sample, std::make_pair(curr.t, false));
      else
        t.insert(
            t.begin() + entry_of_last_sample, std::make_pair(curr.t, true));

      ++count_added;
    }
  }

  return polyline;
}
vector<bezier_segment> import_control_points(
    const bezier_mesh& mesh, const string& filename) {
  std::ifstream f;
  f.open(filename);
  if (!f.is_open()) {
    std::cout << "Anchors file not found,please check the filename"
              << std::endl;
    return {};
  }
  int                tid;
  float              uv_x, uv_y;
  vector<mesh_point> control_points;
  auto               count = 0;
  while (!f.eof()) {
    f >> tid >> uv_x >> uv_y;
    control_points.push_back(mesh_point{tid, {uv_x, uv_y}});
    ++count;
  }
  f.close();
  vector<bezier_segment> control_polygon(count / 4);
  for (auto j = 0; j < count / 4; ++j) {
    for (auto i = 0; i < 4; ++i) {
      control_polygon[j][i] = control_points[i + 4 * j];
    }
  }

  return control_polygon;
}
pair<vector<vec3f>, vector<vec3f>> import_curve(
    const bezier_mesh& mesh, const string& filename) {
  std::ifstream f;

  f.open(filename);
  if (!f.is_open()) {
    std::cout << "WA curve file not found,please check the filename"
              << std::endl;
    return {{}, {}};
  }
  int                tid;
  float              uv_x, uv_y;
  vector<mesh_point> curve;
  vector<vec3f>      polyline;
  vector<vec3f>      sampled;
  while (!f.eof()) {
    f >> tid >> uv_x >> uv_y;

    curve.push_back(mesh_point{tid, {uv_y, 1 - uv_y - uv_x}});
  }
  f.close();
  for (auto i = 0; i < curve.size() - 1; ++i) {
    if (length(eval_position(mesh, curve[i + 1]) -
               eval_position(mesh, curve[i])) > mesh.avg_edge_length * 5)
      continue;
    auto path = my_compute_geodesic_path(mesh, curve[i], curve[i + 1]);
    auto pos  = path_positions(mesh, path);
    polyline.insert(polyline.end(), pos.begin(), pos.end());
    sampled.push_back(pos[0]);
    sampled.push_back(pos.back());
  }

  return {polyline, sampled};
}
mesh_point BFS_on_tri(
    const bezier_mesh& mesh, const vec3f& pos, const int seed) {
  auto            has_been_pushed = vector<bool>(mesh.triangles.size());
  std::deque<int> Q;
  auto [is_in_tri, bary] = point_in_triangle(
      mesh.triangles, mesh.positions, seed, pos);
  if (is_in_tri) return {seed, bary};

  has_been_pushed[seed] = true;
  for (auto adj : mesh.adjacencies[seed]) Q.push_back(adj);

  while (true) {
    auto tid = Q.back();
    Q.pop_back();
    if (has_been_pushed[tid]) continue;
    std::tie(is_in_tri, bary) = point_in_triangle(
        mesh.triangles, mesh.positions, tid, pos);
    if (is_in_tri) return {tid, bary};
    has_been_pushed[tid] = true;
    for (auto adj : mesh.adjacencies[tid]) {
      if (!has_been_pushed[adj]) Q.push_front(adj);
    }
  }
}
vector<vec3f> generate_polyline_from_positions(
    const bezier_mesh& mesh, const vector<vec3f>& pos, const int seed) {
  vector<mesh_point> samples(pos.size());
  vector<vec3f>      polyline;
  auto               curr_seed = seed;
  for (auto i = 0; i < pos.size(); ++i) {
    samples[i] = BFS_on_tri(mesh, pos[i], curr_seed);
    curr_seed  = samples[i].face;
  }

  for (auto i = 0; i < pos.size() - 1; ++i) {
    auto path     = compute_geodesic_path(mesh, samples[i], samples[i + 1]);
    auto path_pos = path_positions(mesh, path);
    polyline.insert(polyline.end(), path_pos.begin(), path_pos.end());
  }

  return polyline;
}
