//
//  karcher.cpp
//  glfw
//
//  Created by Claudio Mancinelli on 21/07/2020.
//

#include "karcher.h"

#include <splinesurf/logging.h>
using namespace logging;

using namespace geometrycentral;
using namespace geometrycentral::surface;

#include <deque>
template <typename T>
inline int find(const T& vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x) return i;
  return -1;
}

// utility (geometry)
// project a vector v onto a plane having n as normal
vec3f project_vec(const vec3f& v, const vec3f& n) {
  auto proj = n * dot(v, n);

  return v - proj;
}
vec3f rot_vect(const vec3f& p, const vec3f& axis, const float angle) {
  auto M = rotation_frame(axis, angle);
  auto v = p;
  return transform_vector(M, v);
}
vec2f rot_vect(const vec2f& p, const float theta) {
  auto M = mat2f{{yocto::cos(theta), -yocto::sin(theta)},
      {yocto::sin(theta), yocto::cos(theta)}};
  auto v = M * p;
  return v;
}
inline vec3f tid_normal(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid) {
  auto p0 = positions[triangles[tid].x];
  auto p1 = positions[triangles[tid].y];
  auto p2 = positions[triangles[tid].z];

  return normalize(cross(p1 - p0, p2 - p0));
}

vec3f tid_centroid(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid) {
  vec3f p0 = positions[triangles[tid].x];
  vec3f p1 = positions[triangles[tid].y];
  vec3f p2 = positions[triangles[tid].z];

  return (p0 + p1 + p2) / 3.0;
}
int vert_from_point(const vector<vec3i>& triangles, const mesh_point& p) {
  auto bary = vector<pair<float, int>>{
      {1 - p.uv.x - p.uv.y, 0}, {p.uv.x, 1}, {p.uv.y, 2}};

  sort(bary.begin(), bary.end());
  return triangles[p.face][bary.back().second];
}
mesh_point point_from_vert(const vector<vec3i>& triangles,
    const vector<vector<int>>& v2t, const int vid, const int tid) {
  auto tr = -1;
  if (tid < 0)
    tr = v2t[vid][0];
  else {
    auto entry = find(v2t[vid], tid);
    tr         = v2t[vid][entry];
  }

  auto k    = find(triangles[tr], vid);
  auto bary = zero3f;
  bary[k]   = 1.0;
  return {tr, {bary.y, bary.z}};
}
float nbr_avg_edge_length(const geodesic_solver& G, int vid) {
  auto nbr = G.graph[vid];
  auto avg = 0.f;
  auto s   = (int)nbr.size();
  for (int i = 0; i < s; ++i) {
    avg += nbr[i].length;
  }

  return (s != 0) ? avg / s : avg;
}
vec3f polar_basis(const geodesic_solver& solver, const vector<vec3f>& positions,
    const vector<vec3f>& normals, int vid) {
  int   vid0 = solver.graph[vid][0].node;
  vec3f v    = positions[vid0] - positions[vid];
  vec3f e    = project_vec(v, normals[vid]);
  return normalize(e);
}
vec3f polar_basis(
    const vector<vec3i>& triangles, const vector<vec3f>& positions, int tid) {
  auto  c = tid_centroid(triangles, positions, tid);
  vec3f v = positions[triangles[tid].x];
  return normalize(v - c);
}
// Compute polar coordinates
float angle_in_tangent_space(const geodesic_solver& solver,
    const vector<vec3f>& positions, const vec3f& v, const int vid,
    const vec3f& n) {
  float teta;
  int   vid0 = solver.graph[vid][0].node;
  vec3f e0   = positions[vid0] - positions[vid];
  vec3f e    = normalize(project_vec(e0, n));

  teta = angle(v, e);
  if (dot(cross(e, v), n) < 0) teta *= -1;

  return teta;
}

float angle_in_tangent_space(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& v) {
  auto teta = 0.f;

  auto e = polar_basis(triangles, positions, tid);
  auto n = tid_normal(triangles, positions, tid);
  teta   = angle(v, e);
  if (dot(cross(e, v), n) < 0) teta = 2 * M_PI - teta;

  return teta;
}
inline int node_is_adjacent(const geodesic_solver& solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}
bool checking_bary(const vec3f& bary, int& entry) {
  int count      = 0;
  entry          = -1;
  float treshold = flt_max;
  for (int i = 0; i < 3; ++i) {
    if (bary[i] < 0) {
      ++count;
      if (yocto::abs(bary[i]) < treshold) {
        entry    = i;
        treshold = yocto::abs(bary[i]);
      }
    }
  }
  if (count > 1)
    return false;
  else if (yocto::abs(bary[entry] <= 1e-1))
    return true;

  return false;
}
inline bool are_barycentric_coordinates(
    const vec3f& bary, const float tol = 1e-2) {
  if (bary.x >= -tol && bary.x <= 1 + tol && bary.y >= -tol &&
      bary.y <= 1 + tol && bary.z >= -tol && bary.x <= 1 + tol)
    return true;
  return false;
}

vec3f tri_bary_coords(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& p) {
  vec3f wgts = vec3f{0.0, 0.0, 0.0};
  vec3f u = v1 - v0, v = v2 - v0, w = p - v0;
  float d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
        d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0) return zero3f;

  wgts[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(wgts[2]));
  wgts[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(wgts[1]));
  wgts[0] = 1.0 - wgts[1] - wgts[2];
  assert(!isnan(wgts[0]));

  return wgts;
}
vec3f tri_bary_coords(
    const vec2f& v0, const vec2f& v1, const vec2f& v2, const vec2f& p) {
  vec3f wgts = vec3f{0.0, 0.0, 0.0};
  vec2f u = v1 - v0, v = v2 - v0, w = p - v0;
  float d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
        d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0) return zero3f;

  wgts[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(wgts[2]));
  wgts[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(wgts[1]));
  wgts[0] = 1.0 - wgts[1] - wgts[2];
  assert(!isnan(wgts[0]));

  return wgts;
}

inline vec3f tri_bary_coords(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int pid, const vec3f& p) {
  auto px = positions[triangles[pid].x];
  auto py = positions[triangles[pid].y];
  auto pz = positions[triangles[pid].z];
  return tri_bary_coords(px, py, pz, p);
}
vec3f vector_bary_coords(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& v) {
  auto A = positions[triangles[tid].x];
  auto B = positions[triangles[tid].y];
  auto C = positions[triangles[tid].z];

  auto p    = A + v;
  auto bary = tri_bary_coords(A, B, C, p);
  return bary;
}
vec3f vector_bary_coords(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& v) {
  auto p    = p0 + v;
  auto bary = tri_bary_coords(p0, p1, p2, p);
  return bary;
}

vec3f vector_bary_coords(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec2f& v) {
  auto flat_tid = init_flat_triangle(positions, triangles[tid]);

  auto p    = flat_tid[0] + v;
  auto bary = tri_bary_coords(flat_tid[0], flat_tid[1], flat_tid[2], p);
  return bary;
}
vec3f vector_bary_coords(const unfold_triangle& flat_tid, const vec2f& v) {
  auto p    = flat_tid[0] + v;
  auto bary = tri_bary_coords(flat_tid[0], flat_tid[1], flat_tid[2], p);
  return bary;
}
pair<bool, vec2f> point_in_unfold_triangle(
    const vec2f& pos, const unfold_triangle& tr, float tol) {
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  auto bary = vec3f{0, 0, 0};
  auto v0   = tr[0];
  auto v1   = tr[1];
  auto v2   = tr[2];

  auto u = v1 - v0, v = v2 - v0, w = pos - v0;
  auto d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
       d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0) return {false, zero2f};

  bary[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(bary[2]));
  bary[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(bary[1]));
  bary[0] = 1.0 - bary[1] - bary[2];
  assert(!isnan(bary[0]));

  for (auto i = 0; i < 3; ++i) {
    if (bary[i] < -tol || bary[i] > 1.0 + tol) return {false, zero2f};
  }
  auto uv = vec2f{bary.y, bary.z};
  uv      = clamp(uv, 0.f, 1.f);
  return {true, uv};
}
mesh_point make_mesh_point(const vector<vec3i>& triangles,
    const vector<vector<int>>& v2t, const int vid, const int tid) {
  if (tid < 0) {
    auto offset  = find(triangles[v2t[vid][0]], vid);
    auto bary    = zero3f;
    bary[offset] = 1.f;
    return {v2t[vid][0], {bary.y, bary.z}};
  } else {
    auto offset  = find(triangles[tid], vid);
    auto bary    = zero3f;
    bary[offset] = 1.f;
    return {tid, {bary.y, bary.z}};
  }
}
// returns the index of the triangles in the star of from which v is pointing
// to
int next_tid(const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vector<int>>& v2t,
    const vector<vec3i>& triangles, const vector<vec3f>& normals,
    const int from, const vec3f& v) {
  auto teta = angle_in_tangent_space(solver, positions, v, from, normals[from]);
  if (teta < 0) teta += 2 * M_PI;
  auto nbr = angles[from];
  int  s   = (int)nbr.size();
  if (teta == 0) return v2t[from][0];
  for (int i = 0; i < s; ++i) {
    if (nbr[i] < teta) continue;

    if (i % 2 == 0) {
      return v2t[from][(i - 2) / 2];
    } else {
      return v2t[from][(i - 1) / 2];
    }
  }
  return v2t[from].back();
}
int next_tid_extended_graph(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<vec3f>& positions,
    const vector<vector<int>>& v2t, const vector<vec3i>& triangles,
    const vector<vec3f>& normals, const int from, const vec3f& v) {
  auto star = v2t[from];
  auto nbr  = solver.graph[from];
  auto teta = angle_in_tangent_space(solver, positions, v, from, normals[from]);
  if (teta == 0) return star[0];
  if (teta < 0) teta += 2 * M_PI;
  for (int i = 0; i < star.size(); ++i) {
    auto tid    = star[i];
    auto offset = find(triangles[tid], from);
    auto vid    = triangles[tid][(offset + 2) % 3];
    auto index  = node_is_adjacent(solver, from, vid);
    if (angles[from][index] > teta) return tid;
  }

  return star.back();
}
// get the k-ring of vid
// note: we consider the connectivity of the graph, to get the usual k- ring
// uncomment the line below
vector<int> k_ring(const geodesic_solver& solver, const int vid, const int k,
    bool mesh_connectivity = true) {
  vector<int> ring;
  vector<int> active_set = {vid};
  for (int i = 0; i < k; ++i) {
    vector<int> next_active_set;
    for (int j = 0; j < active_set.size(); ++j) {
      auto nbr = solver.graph[active_set[j]];
      for (int h = 0; h < nbr.size(); ++h) {
        if (h % 2 && mesh_connectivity) continue;
        int curr = nbr[h].node;
        if (find(ring, curr) == -1 && curr != vid) {
          next_active_set.push_back(curr);
          ring.push_back(curr);
        }
      }
    }
    active_set = next_active_set;
  }

  return ring;
}
vec3f trace_segment_vert(const vector<vec3f>& verts, const vec3f n,
    const vec3f bary, const vec3f baryV, const vec3f& dir, const int offset) {
  auto right       = verts[(offset + 1) % 3] - verts[offset];
  auto left        = verts[(offset + 2) % 3] - verts[offset];
  auto sample_bary = zero3f;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, left), n) > 0) {
      auto factor                   = bary[offset] / baryV[offset];
      sample_bary[offset]           = 0;
      sample_bary[(offset + 1) % 3] = bary[(offset + 1) % 3] -
                                      baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] = bary[(offset + 2) % 3] -
                                      baryV[(offset + 2) % 3] * factor;
    } else
      sample_bary[(offset + 2) % 3] = 1;
  } else
    sample_bary[(offset + 1) % 3] = 1;

  return sample_bary;
}

vec3f trace_segment_edge(const vector<vec3f>& verts, const vec3f n,
    const vec3f bary, const vec3f baryV, const vec3f& dir, const int offset,
    const vec3f& sample_coords) {
  auto sample_bary = zero3f;
  auto right       = verts[(offset + 1) % 3] - sample_coords;
  auto left        = verts[offset] - sample_coords;
  auto front       = verts[(offset + 2) % 3] - sample_coords;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, front), n) > 0) {
      auto factor                   = bary[offset] / baryV[offset];
      sample_bary[offset]           = 0;
      sample_bary[(offset + 1) % 3] = bary[(offset + 1) % 3] -
                                      baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] = bary[(offset + 2) % 3] -
                                      baryV[(offset + 2) % 3] * factor;

    } else {
      if (dot(cross(dir, left), n) > 0) {
        auto factor = bary[(offset + 1) % 3] / baryV[(offset + 1) % 3];
        sample_bary[(offset + 1) % 3] = 0;
        sample_bary[offset]           = bary[offset] - baryV[offset] * factor;
        sample_bary[(offset + 2) % 3] = bary[(offset + 2) % 3] -
                                        baryV[(offset + 2) % 3] * factor;

      } else {
        if (dot(left, dir) > 0) {
          sample_bary[(offset)] = 1;
        } else {
          sample_bary[(offset + 1) % 3] = 1;
        }
      }
    }
  } else {
    if (dot(right, dir) > 0) {
      sample_bary[(offset + 1) % 3] = 1;
    } else {
      sample_bary[(offset)] = 1;
    }
  }

  return sample_bary;
}
vec3f trace_segment_tri(const vector<vec3f>& verts, const vec3f n,
    const vec3f bary, const vec3f baryV, const vec3f& dir,
    const vec3f& sample_coords) {
  auto  sample_bary = zero3f;
  vec3f w0          = verts[0] - sample_coords;
  vec3f w1          = verts[1] - sample_coords;
  vec3f w2          = verts[2] - sample_coords;
  if (dot(cross(w0, dir), n) > 0 && dot(cross(dir, w1), n) > 0) {
    sample_bary = vec3f{bary[0] - bary[2] * baryV[0] / baryV[2],
        bary[1] - bary[2] * baryV[1] / baryV[2], 0};

  } else if (dot(cross(w1, dir), n) > 0 && dot(cross(dir, w2), n) > 0) {
    sample_bary = vec3f{0, bary[1] - bary[0] * baryV[1] / baryV[0],
        bary[2] - bary[0] * baryV[2] / baryV[0]};
  } else {
    sample_bary = vec3f{bary[0] - bary[1] * baryV[0] / baryV[1], 0,
        bary[2] - bary[1] * baryV[2] / baryV[1]};
  }

  return sample_bary;
}
// Identify the intersection of the polyline inside triangle pid
// tracing: https://cims.nyu.edu/gcl/papers/campen2016bms.pdf
void trace_in_triangles(const vector<vec3f>& positions,
    const vector<vec3i>& triangles, const vec3f& dir, const vec3f& bary,
    const int pid, vec3f& sample_pos, vec3f& sample_bary) {
  vec3f         baryM = zero3f, baryV = zero3f;
  vec3f         v0            = positions[triangles[pid].x];
  vec3f         v1            = positions[triangles[pid].y];
  vec3f         v2            = positions[triangles[pid].z];
  vector<vec3f> verts         = {v0, v1, v2};
  vec3f         n             = triangle_normal(v0, v1, v2);
  vec3f         sample_coords = bary.x * v0 + bary.y * v1 + bary.z * v2;
  vec3f         M             = sample_coords + dir;

  baryM = tri_bary_coords(v0, v1, v2, M);
  for (int i = 0; i < 3; ++i) {
    baryV[i] = baryM[i] - bary[i];
  }
  auto [is_vertex, k_vert]  = bary_is_vert(bary);
  auto [is_on_edge, k_edge] = bary_is_edge(bary);
  if (is_vertex) {
    sample_bary = trace_segment_vert(verts, n, bary, baryV, dir, k_vert);
    sample_pos  = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else if (is_on_edge) {
    sample_bary = trace_segment_edge(
        verts, n, bary, baryV, dir, k_edge, sample_coords);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else {
    sample_bary = trace_segment_tri(verts, n, bary, baryV, dir, sample_coords);
    sample_pos  = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  }
}
float rescaled_angle(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const int tid, const int vid, const float& lerp) {
  auto tr = triangles[tid];
  auto k  = find(tr, vid);
  assert(k != -1);
  auto v0 = positions[tr[(k + 1) % 3]], v1 = positions[tr[(k + 2) % 3]];
  auto pos          = (1 - lerp) * v0 + lerp * v1;
  auto v            = pos - positions[vid];
  auto w            = v0 - positions[vid];
  auto u            = v1 - positions[vid];
  auto theta        = angle(v, w);
  auto phi3D        = angle(w, u);
  auto s            = angles[vid].size();
  auto entry        = node_is_adjacent(solver, vid, tr[(k + 1) % 3]);
  auto theta0       = angles[vid][entry];
  auto phi2D        = (entry + 2 == s) ? 2 * pif - theta0
                                       : angles[vid][entry + 2] - theta0;
  auto scale_factor = phi3D / phi2D;

  return theta0 + theta * scale_factor;
}
bool useless_arc(const geodesic_path& arc) {
  for (auto i = 1; i < arc.lerps.size() - 1; ++i) {
    if (1 - arc.lerps[i] < 1e-2 || arc.lerps[i] < 1e-2) return true;
  }
  return false;
}

vector<pair<float, int>> connect_vert_to_neighbors(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const int k, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const int vid,
    const vector<int> k_ring, vector<pair<float, int>>& distances) {
  vector<pair<float, int>> angles_and_nodes;
  auto                     offset = find(triangles[v2t[vid][0]], vid);
  auto                     bary   = zero3f;
  bary[offset]                    = 1;
  mesh_point point                = {v2t[vid][0], vec2f{bary.y, bary.z}};
  auto       nei                  = mesh_point{};
  for (auto i = 0; i < k_ring.size(); ++i) {
    auto curr    = k_ring[i];
    offset       = find(triangles[v2t[curr][0]], curr);
    bary         = zero3f;
    bary[offset] = 1;
    nei          = {v2t[curr][0], vec2f{bary.y, bary.z}};
    vector<int> parents;
    auto strip = get_strip(solver, triangles, positions, adjacencies, v2t,
        angles, nei, point, parents);
    auto path  = shortest_path(
         triangles, positions, adjacencies, point, nei, strip);
    auto entry = node_is_adjacent(solver, vid, curr);
    if (entry >= 0) {
      angles_and_nodes.push_back({angles[vid][entry], i});
      distances[i] = {solver.graph[vid][entry].length, curr};
    } else if (!useless_arc(path)) {
      auto angle = rescaled_angle(solver, angles, triangles, positions,
          adjacencies, path.start.face, vid, path.lerps[0]);
      angles_and_nodes.push_back({angle, i});
      distances[i] = {
          path_length(path, triangles, positions, adjacencies), curr};
    }
  }
  return angles_and_nodes;
}

Eigen::SparseMatrix<double, 1> k_ring_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, geodesic_solver& solver,
    vector<vector<float>>& angles, vector<float>& total_angles) {
  auto old_solver = make_geodesic_solver(
      triangles, positions, adjacencies, v2t);

  auto old_angles = compute_angles(
      triangles, positions, adjacencies, v2t, total_angles, true);
  auto Grad = init_gradient_matrix(old_solver, old_angles, positions, normals);
  auto avg_valence = 0;
  solver.graph.resize(positions.size());
  angles.resize(positions.size());
  for (auto i = 0; i < positions.size(); ++i) {
    auto                     nbr = k_ring(old_solver, i, k, true);
    vector<pair<float, int>> dist(nbr.size());
    vector<pair<float, int>> arcs = connect_vert_to_neighbors(triangles,
        positions, adjacencies, v2t, k, old_solver, old_angles, i, nbr, dist);
    std::sort(arcs.begin(), arcs.end());
    solver.graph[i].resize(arcs.size());
    angles[i].resize(arcs.size());
    avg_valence += arcs.size();
    for (int j = 0; j < arcs.size(); ++j) {
      auto index                = arcs[j].second;
      angles[i][j]              = arcs[j].first;
      solver.graph[i][j].length = dist[index].first;
      solver.graph[i][j].node   = dist[index].second;
    }
  }
  std::cout << "avg valence is" << std::endl;
  std::cout << avg_valence / positions.size() << std::endl;
  return Grad;
}
// parallel transport

void parallel_transp(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, vec3f& v, const vector<vec3f>& normals,
    const int& from, const int& to, const int& mode) {
  switch (mode) {
    case V2V: {
      float teta = angle_in_tangent_space(
          solver, positions, v, from, normals[from]);
      auto  nbr_from = solver.graph[from];
      auto  nbr_to   = solver.graph[to];
      float phi_ij   = -1;
      float phi_ji   = -1;

      for (int i = 0; i < nbr_from.size(); ++i) {
        if (nbr_from[i].node == to) {
          phi_ij = angles[from][i];

          break;
        }
      }

      for (int j = 0; j < nbr_to.size(); ++j) {
        if (nbr_to[j].node == from) {
          phi_ji = angles[to][j];
          break;
        }
      }
      assert(phi_ij != -1);
      assert(phi_ji != -1);

      vec3f e0       = polar_basis(solver, positions, normals, to);
      float rotation = teta + phi_ji + M_PI - phi_ij;

      v = rot_vect(e0, normals[to], rotation);

    } break;

    case V2T: {
      float teta = angle_in_tangent_space(
          solver, positions, v, from, normals[from]);
      vec3i tri      = triangles[to];
      vec3f p0       = positions[tri.x];
      vec3f p1       = positions[tri.y];
      vec3f p2       = positions[tri.z];
      vec3f normal   = triangle_normal(p0, p1, p2);
      vec3f centroid = (p0 + p1 + p2) / 3.0;
      vec3f e        = normalize(p0 - centroid);

      vec3f coords = positions[from] - centroid;
      float phi_ji = angle(e, coords);
      if (dot(cross(e, coords), normal) < 0) phi_ji = 2 * M_PI - phi_ji;
      int offset = find(tri, from);
      assert(offset != -1);
      int   vid1     = tri[(offset + 1) % 3];
      int   vid2     = tri[(offset + 2) % 3];
      float factor   = 2 * M_PI / total_angles[from];
      auto  nbr_from = solver.graph[from];
      float phi_ij   = -1;
      coords *= -1;
      if (nbr_from[0].node == vid2) {
        vec3f edge       = positions[vid2] - positions[from];
        float curr_angle = angle(edge, coords);
        curr_angle *= factor;
        curr_angle = 2 * M_PI - curr_angle;
        phi_ij     = curr_angle;
      } else {
        for (int i = 0; i < nbr_from.size(); ++i) {
          if (nbr_from[i].node == vid1) {
            phi_ij = angles[from][i];
            break;
          }
        }

        vec3f edge       = positions[vid1] - positions[from];
        float curr_angle = angle(edge, coords);
        curr_angle *= factor;
        phi_ij += curr_angle;
      }

      float rot = teta + phi_ji + M_PI - phi_ij;

      e *= length(v);
      v = rot_vect(e, normal, rot);

    }

    break;

    case T2V: {
      vec3i tri      = triangles[from];
      vec3f p0       = positions[tri.x];
      vec3f p1       = positions[tri.y];
      vec3f p2       = positions[tri.z];
      vec3f n        = triangle_normal(p0, p1, p2);
      vec3f centroid = (p0 + p1 + p2) / 3.0;
      vec3f e        = normalize(p0 - centroid);
      float teta     = angle(e, v);

      if (dot(cross(e, v), n) < 0) teta = 2 * M_PI - teta;
      int offset = find(tri, to);
      assert(offset != -1);
      int vid1 = tri[(offset + 1) % 3];

      vec3f vert = positions[tri[offset]];
      vec3f v1   = positions[vid1] - vert;

      vec3f coords = vert - centroid;
      float phi_ij = angle(e, coords);
      if (dot(cross(e, coords), n) < 0) phi_ij = 2 * M_PI - phi_ij;

      coords *= -1;
      float phi_ji = angle(v1, coords);
      float factor = 2 * M_PI / total_angles[to];
      phi_ji *= factor;
      auto nbr = solver.graph[to];
      for (int i = 0; i < nbr.size(); ++i) {
        if (nbr[i].node == vid1) {
          float phi = angles[to][i];
          phi_ji += phi;
          break;
        }
      }

      float rot = teta + phi_ji + M_PI - phi_ij;
      vec3f e0  = polar_basis(solver, positions, normals, to);
      e0 *= length(v);
      v = rot_vect(e0, normals[to], rot);

    } break;

    case T2T: {
      auto flat_from = init_flat_triangle(positions, triangles[from]);
      auto k         = find(adjacencies[from], to);
      assert(k != -1);
      auto flat_to = unfold_face(
          triangles, positions, adjacencies, flat_from, from, k);
      auto bary = vec2f{0.333, 0.333};
      auto c0   = interpolate_triangle(
            flat_from[0], flat_from[1], flat_from[2], bary);
      auto c1 = interpolate_triangle(flat_to[0], flat_to[1], flat_to[2], bary);
      auto e0 = flat_from[0] - c0;
      auto e1 = flat_to[0] - c1;

      auto w      = c1 - c0;
      auto phi_ij = angle(e0, w);
      if (cross(e0, w) < 0) phi_ij = 2 * M_PI - phi_ij;
      w *= -1;
      auto phi_ji = angle(e1, w);
      if (cross(e1, w) < 0) phi_ji = 2 * M_PI - phi_ji;

      auto  n    = tid_normal(triangles, positions, from);
      auto  e    = polar_basis(triangles, positions, from);
      float teta = angle(e, v);
      if (dot(cross(e, v), n) < 0) teta = 2 * M_PI - teta;

      auto  e_to = polar_basis(triangles, positions, to);
      auto  n_to = tid_normal(triangles, positions, to);
      float rot  = teta + phi_ji + M_PI - phi_ij;
      e_to *= length(v);
      v = rot_vect(e_to, n_to, rot);

    }

    break;
  }
}
void parallel_transp(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, vec2f& v, const vector<vec3f>& normals,
    const int& from, const int& to, const int& mode) {
  // switch (mode) {
  //   case V2V: {
  //     float teta     = yocto::atan2(v.y, v.x);
  //     auto  nbr_from = solver.graph[from];
  //     auto  nbr_to   = solver.graph[to];
  //     float phi_ij   = -1;
  //     float phi_ji   = -1;

  //     for (int i = 0; i < nbr_from.size(); ++i) {
  //       if (nbr_from[i].node == to) {
  //         phi_ij = angles[from][i];

  //         break;
  //       }
  //     }

  //     for (int j = 0; j < nbr_to.size(); ++j) {
  //       if (nbr_to[j].node == from) {
  //         phi_ji = angles[to][j];
  //         break;
  //       }
  //     }
  //     assert(phi_ij != -1);
  //     assert(phi_ji != -1);

  //     auto  e0       = vec2f{1, 0};
  //     float rotation = teta + phi_ji + M_PI - phi_ij;

  //     v = rot_vect(e0, rotation);

  //   } break;

  //   case V2T: {
  //     auto teta     = yocto::atan2(v.y, v.x);
  //     auto tri      = triangles[to];
  //     auto flat_tri = init_flat_triangle(positions, tri);
  //     auto p0       = flat_tri[0];
  //     auto p1       = flat_tri[1];
  //     auto p2       = flat_tri[2];
  //     auto centroid = (p0 + p1 + p2) / 3.0;
  //     auto e        = normalize(p0 - centroid);

  //     auto  k        = find(tri, from);
  //     auto  coords   = flat_tri[k] - centroid;
  //     auto  coords3D = positions[from] - tid_centroid(positions, triangles,
  //     to); float phi_ji   = angle(e, coords); if (cross(e, coords) < 0)
  //     phi_ji = 2 * M_PI - phi_ji; int offset = find(tri, from); assert(offset
  //     != -1); int   vid1     = tri[(offset + 1) % 3]; int   vid2     =
  //     tri[(offset + 2) % 3]; float factor   = 2 * M_PI / total_angles[from];
  //     auto  nbr_from = solver.graph[from];
  //     float phi_ij   = -1;
  //     coords *= -1;
  //     if (nbr_from[0].node == vid2) {
  //       vec3f edge       = positions[vid2] - positions[from];
  //       float curr_angle = angle(edge, coords3D);
  //       curr_angle *= factor;
  //       curr_angle = 2 * M_PI - curr_angle;
  //       phi_ij     = curr_angle;
  //     } else {
  //       for (int i = 0; i < nbr_from.size(); ++i) {
  //         if (nbr_from[i].node == vid1) {
  //           phi_ij = angles[from][i];
  //           break;
  //         }
  //       }

  //       vec3f edge       = positions[vid1] - positions[from];
  //       float curr_angle = angle(edge, coords3D);
  //       curr_angle *= factor;
  //       phi_ij += curr_angle;
  //     }

  //     float rot = teta + phi_ji + M_PI - phi_ij;

  //     e *= length(v);
  //     v = rot_vect(e, rot);

  //   }

  //   break;

  //   case T2V: {
  //     auto tri      = triangles[to];
  //     auto flat_tri = init_flat_triangle(positions, tri);
  //     auto p0       = flat_tri[0];
  //     auto p1       = flat_tri[1];
  //     auto p2       = flat_tri[2];
  //     auto centroid = (p0 + p1 + p2) / 3.0;
  //     auto e        = normalize(p0 - centroid);
  //     auto teta     = angle(e, v);

  //     if (cross(e, v) < 0) teta = 2 * M_PI - teta;
  //     int k    = find(tri, to);
  //     int vid1 = tri[(k + 1) % 3];

  //     auto  coords   = flat_tri[k] - centroid;
  //     auto  coords3D = tid_centroid(positions, triangles, from) -
  //     positions[to]; vec3f v1       = positions[vid1] - vert;

  //     float phi_ij = angle(e, coords);
  //     if (cross(e, coords) < 0) phi_ij = 2 * M_PI - phi_ij;

  //     float phi_ji = angle(v1, coords3D);
  //     float factor = 2 * M_PI / total_angles[to];
  //     phi_ji *= factor;
  //     auto nbr = solver.graph[to];
  //     for (int i = 0; i < nbr.size(); ++i) {
  //       if (nbr[i].node == vid1) {
  //         float phi = angles[to][i];
  //         phi_ji += phi;
  //         break;
  //       }
  //     }

  //     float rot = teta + phi_ji + M_PI - phi_ij;
  //     vec3f e0  = polar_basis(solver, positions, normals, to);
  //     e0 *= length(v);
  //     v = rot_vect(e0, normals[to], rot);

  //   } break;

  //   case T2T: {
  //     auto flat_from = init_flat_triangle(positions, triangles[from]);
  //     auto k         = find(adjacencies[from], to);
  //     assert(k != -1);
  //     auto flat_to = unfold_face(
  //         triangles, positions, adjacencies, flat_from, from, k);
  //     auto bary = vec2f{0.333, 0.333};
  //     auto c0   = interpolate_triangle(
  //         flat_from[0], flat_from[1], flat_from[2], bary);
  //     auto c1 = interpolate_triangle(flat_to[0], flat_to[1], flat_to[2],
  //     bary); auto e0 = flat_from[0] - c0; auto e1 = flat_to[0] - c1;

  //     auto w      = c1 - c0;
  //     auto phi_ij = angle(e0, w);
  //     if (cross(e0, w) < 0) phi_ij = 2 * M_PI - phi_ij;
  //     w *= -1;
  //     auto phi_ji = angle(e1, w);
  //     if (cross(e1, w) < 0) phi_ji = 2 * M_PI - phi_ji;

  //     auto  n    = tid_normal(positions, triangles, from);
  //     auto  e    = polar_basis(triangles, positions, from);
  //     float teta = angle(e, v);
  //     if (dot(cross(e, v), n) < 0) teta = 2 * M_PI - teta;

  //     auto  e_to = polar_basis(triangles, positions, to);
  //     auto  n_to = tid_normal(positions, triangles, to);
  //     float rot  = teta + phi_ji + M_PI - phi_ij;
  //     e_to *= length(v);
  //     v = rot_vect(e_to, n_to, rot);

  //   }

  //   break;
  // }
}
vec3f transp_vec(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vec3f& v,
    const vector<vec3f>& normals, const int& from, const int& to,
    const int& mode) {
  vec3f w = v;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, w, normals, from, to, mode);
  return w;
}
vector<int> strip_to_point(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    int parent, const mesh_point& p) {
  vector<int> strip = {p.face};
  auto [is_vert, k] = bary_is_vert(get_bary(p.uv));
  if (is_vert) {
    auto vid   = triangles[p.face][k];
    auto entry = node_is_adjacent(solver, vid, parent);
    assert(entry >= 0);
    if (entry % 2) {
      auto first = (entry - 1) / 2;
      auto tid   = opposite_face(triangles, adjacencies, v2t[vid][first], vid);
      strip.push_back(tid);
    }
    return strip;

  } else {
    auto h = find(triangles[p.face], parent);

    if (h == -1) {
      for (auto i = 0; i < 3; ++i) {
        auto adj = adjacencies[p.face][i];
        h        = find(triangles[adj], parent);
        if (h != -1) {
          strip.push_back(adj);
          return strip;
        }
      }
    }
  }
  return strip;
}
vec3f tranport_vector(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vec3f& v, const vector<vec3f>& normals, const mesh_point& from,
    const mesh_point& to) {
  auto dir     = v;
  auto parents = point_to_point_geodesic_path(
      solver, triangles, positions, adjacencies, to, from);
  // handle degenerate cases
  if (parents.size() == 0) {
    auto [is_vert, k] = bary_is_vert(get_bary(from.uv));
    if (is_vert) {
      auto vid = triangles[from.face][k];
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, from.face, vid, T2V);
      auto strip = strip_to_point(
          solver, triangles, positions, adjacencies, v2t, vid, to);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, vid, strip.back(), V2T);
      if (strip.size())
        return dir;
      else {
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, dir, normals, strip.back(), strip[0], T2T);
        return dir;
      }
    }
  }
  // bring dir in the tangent space of parents[0]
  auto [from_is_vert, kf] = bary_is_vert(get_bary(from.uv));
  if (from_is_vert) {
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, from.face, triangles[from.face][kf], T2V);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, triangles[from.face][kf], parents[0], V2V);
  } else {
    auto strip = strip_to_point(
        solver, triangles, positions, adjacencies, v2t, parents[0], from);
    for (int i = 0; i < strip.size() - 1; ++i) {
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, strip[i], strip[i + 1], T2T);
    }
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, strip.back(), parents[0], T2V);
  }

  for (int i = 0; i < parents.size() - 1; ++i) {
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, parents[i], parents[i + 1], V2V);
  }
  auto [to_is_vert, kt] = bary_is_vert(get_bary(to.uv));
  if (to_is_vert) {
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, to.face, triangles[to.face][kt], T2V);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, triangles[to.face][kt], parents.back(), V2V);
    return dir;
  } else {
    auto strip = strip_to_point(
        solver, triangles, positions, adjacencies, v2t, parents.back(), to);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, parents.back(), strip.back(), V2T);
    if (strip.size())
      return dir;
    else {
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, strip.back(), strip[0], T2T);
      return dir;
    }
  }
}
std::tuple<vector<vec3f>, vector<vec3i>, unordered_map<int, int>>
add_vertices_to_mesh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<mesh_point>& points) {
  auto map      = unordered_map<int, int>{};
  auto n        = (int)points.size();
  auto Vm       = (int)positions.size();
  auto Fm       = (int)triangles.size();
  int  V        = (int)(Vm + n);
  int  F        = (int)(Fm + 2 * n);
  auto added    = 0;
  auto vertices = vector<vec3f>(V);
  auto tri      = vector<vec3i>(F);

  auto pos  = vector<vec3f>(n);
  auto tids = vector<int>(n);
  for (auto i = 0; i < n; ++i) {
    pos[i]  = eval_position(triangles, positions, points[i]);
    tids[i] = points[i].face;
  }
  for (int i = 0; i < V; ++i) {
    // if (v2t[i].size() == 0) continue;
    if (i < Vm) {
      vertices[i] = positions[i];
    } else {
      vertices[i] = pos[i - Vm];
      map[i - Vm] = i;
    }
  }
  for (int i = 0; i < Fm; ++i) {
    auto it = std::find(tids.begin(), tids.end(), i);
    if (it == tids.end()) {
      tri[i] = triangles[i];
    } else {
      tri[i].x = triangles[i].x;
      tri[i].y = triangles[i].y;
      tri[i].z = Vm + added;

      tri[Fm + 2 * added].x = triangles[i].y;
      tri[Fm + 2 * added].y = triangles[i].z;
      tri[Fm + 2 * added].z = Vm + added;

      tri[Fm + 2 * added + 1].x = triangles[i].z;
      tri[Fm + 2 * added + 1].y = triangles[i].x;
      tri[Fm + 2 * added + 1].z = Vm + added;

      ++added;
    }
  }
  return {vertices, tri, map};
}
// VTP
vector<float> exact_geodesic_distance(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vector<int>>& v2t,
    const int& source) {
  int            V = positions.size();
  int            F = triangles.size();
  vector<double> points(3 * V);
  vector<uint>   faces(3 * F);
  vector<float>  f(V);

  for (int i = 0; i < V; ++i) {
    if (v2t[i].size() == 0) continue;
    points[3 * i]     = positions[i].x;
    points[3 * i + 1] = positions[i].y;
    points[3 * i + 2] = positions[i].z;
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i]     = triangles[i].x;
    faces[3 * i + 1] = triangles[i].y;
    faces[3 * i + 2] = triangles[i].z;
  }

  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v     = verts[j];
    float                value = (float)v.geodesic_distance();
    f[j]                       = value;
  }

  return f;
}
vector<float> exact_geodesic_distance(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vector<int>>& v2t,
    const mesh_point& source) {
  auto tid               = source.face;
  auto [is_vert, offset] = point_is_vert(source);
  if (is_vert)
    return exact_geodesic_distance(
        triangles, positions, v2t, triangles[tid][offset]);
  else {
    int            V = positions.size() + 1;
    int            F = triangles.size();
    vector<double> points(3 * V);
    vector<uint>   faces(3 * (F + 2));
    vector<float>  f(V);
    auto           pos = eval_position(triangles, positions, source);

    for (int i = 0; i < V; ++i) {
      // if (v2t[i].size() == 0) continue;
      if (i != V - 1) {
        points[3 * i]     = positions[i].x;
        points[3 * i + 1] = positions[i].y;
        points[3 * i + 2] = positions[i].z;
      } else {
        points[3 * i]     = pos.x;
        points[3 * i + 1] = pos.y;
        points[3 * i + 2] = pos.z;
      }
    }
    for (int i = 0; i < F; ++i) {
      if (i != tid) {
        faces[3 * i]     = triangles[i].x;
        faces[3 * i + 1] = triangles[i].y;
        faces[3 * i + 2] = triangles[i].z;
      } else {
        faces[3 * i]     = triangles[i].x;
        faces[3 * i + 1] = triangles[i].y;
        faces[3 * i + 2] = V - 1;

        faces[3 * F]     = triangles[i].y;
        faces[3 * F + 1] = triangles[i].z;
        faces[3 * F + 2] = V - 1;

        faces[3 * (F + 1)]     = triangles[i].z;
        faces[3 * (F + 1) + 1] = triangles[i].x;
        faces[3 * (F + 1) + 2] = V - 1;
      }
    }

    geodesic_VTP::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
    algorithm.propagate(V - 1);
    vector<geodesic_VTP::Vertex> verts = mesh.vertices();
    for (int j = 0; j < V; ++j) {
      geodesic_VTP::Vertex v     = verts[j];
      float                value = (float)v.geodesic_distance();
      f[j]                       = value;
    }
    f.pop_back();
    return f;
  }
}
// utility (Compute distances)
template <typename Update, typename Stop, typename Exit>
void visit_graph(vector<float>& field, const geodesic_solver& solver,
    const vector<int>& sources, Update&& update, Stop&& stop, Exit&& exit) {
  /*
     This algortithm uses the heuristic Small Label Fisrt and Large Label Last
     https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

     Large Label Last (LLL): When extracting nodes from the queue, pick the
     front one. If it weights more than the average weight of the queue, put
     on the back and check the next node. Continue this way.
     Sometimes average_weight is less than every value due to floating point
     errors (doesn't happen with double precision).

     Small Label First (SLF): When adding a new node to queue, instead of
     always pushing it to the end of the queue, if it weights less than the
     front node of the queue, it is put on front. Otherwise the node is put at
     the end of the queue.
  */

  auto in_queue = vector<bool>(solver.graph.size(), false);
  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue = std::deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)cumulative_weight / queue.size();

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node, field[node])) continue;
    if (stop(node)) continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      // Distance of neighbor through this node
      auto new_distance = field[node] + solver.graph[node][i].length;
      auto neighbor     = solver.graph[node][i].node;

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      // Binomial Coefficient
      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      update(node, neighbor, new_distance);
    }
  }
}
static vector<pair<int, float>> check_surrounding_nodes(
    vector<pair<int, float>>& nodes) {
  sort(nodes.begin(), nodes.end());
  auto new_nodes = vector<pair<int, float>>{};
  for (auto i = 1; i < nodes.size(); ++i) {
    auto prev = nodes[i - 1];
    auto curr = nodes[i];
    if (prev.first == curr.first) {
      auto d0 = prev.second, d1 = curr.second;
      if (d0 <= d1)
        new_nodes.push_back(prev);
      else
        new_nodes.push_back(curr);
      ++i;
    } else {
      new_nodes.push_back(prev);
    }
  }
  if (nodes.back().first != new_nodes.back().first)
    new_nodes.push_back(nodes.back());

  nodes = new_nodes;
  return nodes;
}
vector<pair<int, float>> nodes_around_mesh_point(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& p) {
  auto nodes               = vector<pair<int, float>>{};
  auto [is_vertex, offset] = bary_is_vert(get_bary(p.uv));

  if (is_vertex) {
    auto vid = triangles[p.face][offset];
    nodes.push_back({vid, 0});
  } else {
    auto pid = p.face;
    auto pos = eval_position(triangles, positions, p);

    for (int i = 0; i < 3; ++i) {
      int   p0 = triangles[pid][i], p1 = triangles[pid][(i + 1) % 3];
      float d = length(positions[p0] - pos);
      nodes.push_back(std::make_pair(p0, d));

      int         CW_pid = adjacencies[pid][i];
      int         opp    = opposite_vertex(triangles, adjacencies, pid, i);
      vector<int> strip  = {CW_pid, pid};
      float       l      = length_by_flattening(
                     triangles, positions, adjacencies, p, strip);

      nodes.push_back(std::make_pair(opp, l));

      int opp_pid = opposite_face(triangles, adjacencies, CW_pid, p0);
      strip       = {opp_pid, CW_pid, pid};
      int k       = find(adjacencies[CW_pid], opp_pid);
      int q       = opposite_vertex(triangles, adjacencies, CW_pid, k);
      d = length_by_flattening(triangles, positions, adjacencies, p, strip);
      nodes.push_back(std::make_pair(q, d));

      opp_pid = opposite_face(triangles, adjacencies, CW_pid, p1);
      strip   = {opp_pid, CW_pid, pid};

      k = find(adjacencies[CW_pid], opp_pid);
      q = opposite_vertex(triangles, adjacencies, CW_pid, k);
      d = length_by_flattening(triangles, positions, adjacencies, p, strip);
      nodes.push_back(std::make_pair(q, d));
    }
  }

  return nodes;
}
unfold_triangle flat_triangle_with_origin_at_vertex(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vec3i& tr, const int vid) {
  auto k = find(tr, vid);
  assert(k != -1);
  auto tr2d         = unfold_triangle{};
  tr2d[k]           = {0, 0};
  tr2d[(k + 1) % 3] = {
      0, length(positions[tr[k]] - positions[tr[(k + 1) % 3]])};
  auto rx = length_squared(positions[tr[k]] - positions[tr[(k + 2) % 3]]);
  auto ry = length_squared(
      positions[tr[(k + 1) % 3]] - positions[tr[(k + 2) % 3]]);
  tr2d[(k + 2) % 3] = intersect_circles(tr2d[k], rx, tr2d[(k + 1) % 3], ry);
  return tr2d;
}
vector<unfold_triangle> partial_one_ring_flattening(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const int vid, const int first_tid, const int final_tid) {
  auto it    = find(v2t[vid].begin(), v2t[vid].end(), first_tid);
  auto entry = distance(v2t[vid].begin(), it);
  auto curr  = v2t[vid][entry];

  auto strip = vector<unfold_triangle>{};
  strip.push_back(flat_triangle_with_origin_at_vertex(
      triangles, positions, triangles[curr], vid));
  while (curr != final_tid) {
    auto k = find(triangles[curr], vid);
    strip.push_back(unfold_face(
        triangles, positions, adjacencies, strip.back(), curr, (k + 2) % 3));
    curr = adjacencies[curr][(k + 2) % 3];
  }

  return strip;
}

vector<pair<int, float>> nodes_around_mesh_point_extended(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const mesh_point& p) {
  auto nodes               = vector<pair<int, float>>{};
  auto [is_vertex, offset] = bary_is_vert(get_bary(p.uv));

  if (is_vertex) {
    auto vid = triangles[p.face][offset];
    nodes.push_back({vid, 0});
  } else {
    auto pid = p.face;
    auto pos = eval_position(triangles, positions, p);

    for (int i = 0; i < 3; ++i) {
      int   p0 = triangles[pid][i];
      float d  = length(positions[p0] - pos);
      nodes.push_back(std::make_pair(p0, d));
      auto k        = find(triangles[pid], p0);
      auto last_tid = adjacencies[pid][i];
      auto strip    = partial_one_ring_flattening(
             triangles, positions, adjacencies, v2t, p0, pid, last_tid);
      auto curr_tid = adjacencies[pid][(k + 2) % 3];
      auto flat_pos = interpolate_triangle(
          strip[0][0], strip[0][1], strip[0][2], p.uv);
      for (auto j = 1; j < strip.size() - 1; ++j) {
        k             = find(triangles[curr_tid], p0);
        auto curr_vid = triangles[curr_tid][(k + 2) % 3];
        d             = length(strip[j][(k + 2) % 3] - flat_pos);
        nodes.push_back(std::make_pair(curr_vid, d));
        curr_tid = adjacencies[curr_tid][(k + 2) % 3];
      }
    }
  }

  return nodes;
}

vector<float> solve_with_targets(const geodesic_solver& solver,
    const vector<pair<int, float>>&                     sources_and_dist,
    const vector<pair<int, float>>&                     targets) {
  auto update       = [](int node, int neighbor, float new_distance) {};
  auto stop         = [](int node) { return false; };
  auto max_distance = flt_min;
  auto exit_verts   = vector<int>{};
  for (auto i = 0; i < targets.size(); ++i) {
    auto it = find(exit_verts.begin(), exit_verts.end(), targets[i].first);
    if (it == exit_verts.end()) exit_verts.push_back(targets[i].first);
  }
  auto exit = [&exit_verts, &max_distance](
                  int node, const float& curr_distance) {
    auto it = find(exit_verts.begin(), exit_verts.end(), node);
    if (it != exit_verts.end()) {
      max_distance = yocto::max(curr_distance, max_distance);
      exit_verts.erase(it);
    }

    if (exit_verts.empty() && curr_distance > max_distance + max_distance / 5) {
      return true;
    }
    return false;
  };

  auto distances  = vector<float>(solver.graph.size(), flt_max);
  auto sources_id = vector<int>(sources_and_dist.size());
  for (auto i = 0; i < sources_and_dist.size(); ++i) {
    sources_id[i]                        = sources_and_dist[i].first;
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }

  visit_graph(distances, solver, sources_id, update, stop, exit);
  return distances;
}
static vector<float> solve_geodesic_distance(const geodesic_solver& solver,
    const vector<pair<int, float>>& sources_and_dist) {
  auto update = [](int node, int neighbor, float new_distance) {};
  auto stop   = [](int node) { return false; };
  auto exit   = [](int node, float& curr_distance) { return false; };

  auto distances = vector<float>{};
  distances.assign(solver.graph.size(), flt_max);
  auto sources_id = vector<int>(sources_and_dist.size());
  for (auto i = 0; i < sources_and_dist.size(); ++i) {
    sources_id[i]                        = sources_and_dist[i].first;
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }
  visit_graph(distances, solver, sources_id, update, stop, exit);

  return distances;
}

vector<float> compute_pruned_geodesic_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const mesh_point& source, const vector<mesh_point>& targets) {
  auto target_nodes = vector<pair<int, float>>{};

  auto source_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, source);

  for (auto i = 0; i < targets.size(); ++i) {
    auto curr_nodes = nodes_around_mesh_point(
        triangles, positions, adjacencies, targets[i]);
    if (curr_nodes.size() > 1) check_surrounding_nodes(curr_nodes);
    target_nodes.insert(
        target_nodes.end(), curr_nodes.begin(), curr_nodes.end());
  }
  if (source_nodes.size() > 1) check_surrounding_nodes(source_nodes);

  return solve_with_targets(solver, source_nodes, target_nodes);
}
vector<float> compute_accurate_geodesic_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<mesh_point>& sources) {
  auto source_nodes = vector<pair<int, float>>{};
  for (auto i = 0; i < sources.size(); ++i) {
    auto curr_nodes = nodes_around_mesh_point_extended(
        triangles, positions, adjacencies, v2t, sources[i]);
    source_nodes.insert(
        source_nodes.end(), curr_nodes.begin(), curr_nodes.end());
  }
  if (source_nodes.size() > 1) check_surrounding_nodes(source_nodes);
  return solve_geodesic_distance(solver, source_nodes);
}
vector<float> compute_distance_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const int entry, const vector<mesh_point>& control_points,
    const int type_of_field) {
  auto f = vector<float>{};
  if (type_of_field == exact) {
    f = exact_geodesic_distance(
        triangles, positions, v2t, control_points[entry]);
  } else {
    auto curr_targets = control_points;
    curr_targets.erase(curr_targets.begin() + entry);
    f = compute_pruned_geodesic_distances(solver, triangles, positions,
        adjacencies, v2t, control_points[entry], curr_targets);
    // f = compute_geodesic_distances(
    //     solver, triangles, positions, adjacencies, {control_points[entry]});
  }
  return f;
}
int get_parent(const vector<pair<int, float>>& nbr, const vector<float>& f) {
  auto vid    = -1;
  auto lambda = flt_max;
  for (auto i = 0; i < nbr.size(); ++i) {
    auto val = f[nbr[i].first] + nbr[i].second;
    if (val < lambda) {
      lambda = val;
      vid    = nbr[i].first;
    }
  }
  return vid;
}
int get_parent(
    const geodesic_solver& solver, const int vid, const vector<float>& f) {
  auto parent = -1;
  auto lambda = flt_max;
  auto nbr    = solver.graph[vid];

  for (auto i = 0; i < nbr.size(); ++i) {
    auto val = f[nbr[i].node] + nbr[i].length;
    if (val < lambda) {
      lambda = val;
      parent = nbr[i].node;
    }
  }

  return parent;
}
bool source_is_near(
    const vector<pair<int, float>>& source_nodes, const int vid) {
  for (auto i = 0; i < source_nodes.size(); ++i) {
    if (source_nodes[i].first == vid) return true;
  }
  return false;
}
int get_child(const geodesic_solver& solver,
    const vector<pair<int, float>>& source_nodes, const int target_parent,
    const vector<float>& scalar_field, vector<int>& path) {
  auto stop          = false;
  auto prev          = target_parent;
  auto source_parent = get_parent(solver, prev, scalar_field);
  path               = {target_parent};
  if (source_parent == -1) return target_parent;
  while (!stop) {
    prev = source_parent;
    path.push_back(prev);
    source_parent = get_parent(solver, source_parent, scalar_field);

    if (source_is_near(source_nodes, source_parent)) {
      path.push_back(source_parent);
      stop = true;
    }
  }

  return source_parent;
}
vector<int> parents_from_distance_field(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<float>& scalar_field,
    const mesh_point& source, const mesh_point& target) {
  auto source_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, source);
  auto target_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, target);
  auto parents       = vector<int>{};
  auto target_parent = -1, source_child = -1;
  if (target_nodes.size() == 1)
    target_parent = get_parent(solver, target_nodes[0].first, scalar_field);
  else
    target_parent = get_parent(target_nodes, scalar_field);
  auto path = vector<int>{};
  if (target_parent == -1) return path;
  source_child = get_child(
      solver, source_nodes, target_parent, scalar_field, path);
  if (source_nodes.size() == 1)
    path.pop_back();  // we remove "source" from the list of parents if it is
                      // a vertex
  return path;
}
vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<pair<int, float>>& sources_and_dist) {
  auto          update = [](int node, int neighbor, float new_distance) {};
  auto          stop   = [](int node) { return false; };
  auto          exit   = [](int node, float& curr_distance) { return false; };
  vector<float> distances;
  distances.assign(solver.graph.size(), flt_max);
  vector<int> sources_id(sources_and_dist.size());
  for (int i = 0; i < sources_and_dist.size(); ++i) {
    sources_id.push_back(sources_and_dist[i].first);
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }
  visit_graph(distances, solver, sources_id, update, stop, exit);

  return distances;
}
std::tuple<vector<int>, vector<float>> get_parents_and_distances(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& source, const mesh_point& target) {
  auto source_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, source);
  auto target_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, target);
  auto parents       = vector<int>{};
  auto distances     = compute_geodesic_distances(solver, source_nodes);
  auto target_parent = -1, source_child = -1;
  if (target_nodes.size() == 1)
    target_parent = parents[target_nodes[0].first];
  else
    target_parent = get_parent(target_nodes, distances);
  auto path = vector<int>{};
  if (target_parent == -1) return {path, distances};
  source_child = get_child(
      solver, source_nodes, target_parent, distances, path);
  if (source_nodes.size() == 1)
    path.pop_back();  // we remove "source" from the list of parents if it is
                      // a vertex
  return {path, distances};
}
float exact_length(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& p, const mesh_point& q, const vector<int>& strip) {
  auto coords = vector<unfold_triangle>(strip.size());
  coords[0]   = init_flat_triangle(positions, triangles[strip[0]]);
  for (auto i = 1; i < strip.size(); i++) {
    auto k = find(adjacencies[strip[i - 1]], strip[i]);
    assert(k != -1);
    auto tr = unfold_face(
        triangles, positions, adjacencies, coords[i - 1], strip[i - 1], k);
    coords[i] = tr;
  }

  auto last  = coords.back();
  auto p_pos = interpolate_triangle(
      coords[0][0], coords[0][1], coords[0][2], p.uv);
  auto q_pos = interpolate_triangle(last[0], last[1], last[2], q.uv);

  auto v = q_pos - p_pos;
  return length(v);
}
// distances is a scalar field that associates to every vertex v the geodesic
// distance from "source". This function returns the geodesic distance from
// "source" to "target" exploiting such scalar field.
float point2point_geodesic_distance(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<float>& distances, const mesh_point& source,
    const mesh_point& target) {
  auto d                   = flt_max;
  auto [is_vertex, offset] = bary_is_vert(get_bary(target.uv));

  if (is_vertex) {
    auto vid = triangles[target.face][offset];
    return distances[vid];
  } else {
    auto tid            = target.face;
    auto target_pos     = eval_position(triangles, positions, target);
    auto source_pos     = eval_position(triangles, positions, source);
    auto [inside, bary] = point_in_triangle(
        triangles, positions, source.face, target_pos);
    if (inside) return length(target_pos - source_pos);

    for (int i = 0; i < 3; ++i) {
      int p0 = triangles[tid][i];
      d      = yocto::min(d, distances[p0] + length(target_pos - p0));

      int CW_tid             = adjacencies[tid][i];
      std::tie(inside, bary) = point_in_triangle(
          triangles, positions, CW_tid, source_pos);
      if (inside)
        return (exact_length(
            triangles, positions, adjacencies, target, source, {tid, CW_tid}));

      auto k            = find(adjacencies[tid], CW_tid);
      auto q            = opposite_vertex(triangles, adjacencies, tid, k);
      k                 = find(triangles[CW_tid], q);
      auto bary3d       = zero3f;
      bary3d[k]         = 1;
      auto q_mesh_point = mesh_point{CW_tid, vec2f{bary3d.y, bary3d.z}};
      auto d_q = exact_length(triangles, positions, adjacencies, target,
          q_mesh_point, {tid, CW_tid});
      d        = yocto::min(d, distances[q] + d_q);

      int opp_tid = opposite_face(triangles, adjacencies, CW_tid, p0);
      std::tie(inside, bary) = point_in_triangle(
          triangles, positions, opp_tid, source_pos);
      if (inside)
        return (exact_length(triangles, positions, adjacencies, target, source,
            {tid, CW_tid, opp_tid}));
      k            = find(adjacencies[CW_tid], opp_tid);
      q            = opposite_vertex(triangles, adjacencies, CW_tid, k);
      k            = find(triangles[opp_tid], q);
      bary3d       = zero3f;
      bary3d[k]    = 1;
      q_mesh_point = mesh_point{opp_tid, vec2f{bary3d.y, bary3d.z}};
      d_q          = exact_length(triangles, positions, adjacencies, target,
                   q_mesh_point, {tid, CW_tid, opp_tid});
      d            = yocto::min(d, distances[q] + d_q);

      opp_tid = opposite_face(
          triangles, adjacencies, CW_tid, triangles[tid][(i + 1) % 3]);
      std::tie(inside, bary) = point_in_triangle(
          triangles, positions, opp_tid, source_pos);
      if (inside)
        return (exact_length(triangles, positions, adjacencies, target, source,
            {tid, CW_tid, opp_tid}));

      k            = find(adjacencies[CW_tid], opp_tid);
      q            = opposite_vertex(triangles, adjacencies, CW_tid, k);
      k            = find(triangles[opp_tid], q);
      bary3d       = zero3f;
      bary3d[k]    = 1;
      q_mesh_point = mesh_point{opp_tid, vec2f{bary3d.y, bary3d.z}};
      d_q          = exact_length(triangles, positions, adjacencies, target,
                   q_mesh_point, {tid, CW_tid, opp_tid});
      d            = yocto::min(d, distances[q] + d_q);
    }
  }

  return d;
}
pair<geodesic_path, vector<float>> cut_curve(const geodesic_path& path,
    const vector<float>& path_parameter_t, int entry, const mesh_point& sample,
    bool keep_the_first_portion) {
  assert(path.strip.size() == path.lerps.size() + 1);
  assert(path.strip.size() == path_parameter_t.size() - 1);
  assert(entry >= 0 && entry < path.strip.size());
  assert(sample.face == path.strip[entry]);
  // output vars
  auto cut   = geodesic_path{};
  auto cut_t = vector<float>{};

  if (keep_the_first_portion) {
    cut.start = path.start;
    cut.end   = sample;
    cut.strip.insert(
        cut.strip.end(), path.strip.begin(), path.strip.begin() + entry + 1);
    cut.lerps.insert(
        cut.lerps.end(), path.lerps.begin(), path.lerps.begin() + entry);

    cut_t.insert(cut_t.end(), path_parameter_t.begin(),
        path_parameter_t.begin() + entry + 1);
    assert(cut_t.back() <= 0.5);
    cut_t.push_back(0.5);
    for (int i = 0; i < cut_t.size(); i++) cut_t[i] *= 2;
  } else {
    cut.start = sample;
    cut.end   = path.end;
    cut.strip.insert(
        cut.strip.end(), path.strip.begin() + entry, path.strip.end());
    cut.lerps.insert(
        cut.lerps.end(), path.lerps.begin() + entry, path.lerps.end());
    cut_t.insert(
        cut_t.end(), path_parameter_t.begin() + entry, path_parameter_t.end());
    assert(cut_t[0] < 0.5);
    cut_t[0] = 0;  // we take the t corresponding to the prev sample and we
                   // put it to zero
    for (int i = 1; i < cut_t.size(); i++) cut_t[i] = (cut_t[i] - 0.5) * 2;
  }
  assert(cut.strip.size() == cut.lerps.size() + 1);
  assert(cut.strip.size() == cut_t.size() - 1);
  return std::make_pair(cut, cut_t);
}

std::tuple<geodesic_path, geodesic_path, vector<float>, vector<float>>
split_curve(geodesic_path& path, vector<float>& path_parameter_t, int entry,
    const mesh_point& sample) {
  if (path.strip.size() != path.lerps.size() + 1)
    std::cout << "HAPPENS" << std::endl;
  auto path0   = geodesic_path{};
  auto path1   = geodesic_path{};
  auto path0_t = vector<float>{};
  auto path1_t = vector<float>{};

  path0.start = path.start;
  path0.end   = sample;

  path0.strip.insert(
      path0.strip.end(), path.strip.begin(), path.strip.begin() + entry + 1);
  path0.lerps.insert(
      path0.lerps.end(), path.lerps.begin(), path.lerps.begin() + entry);

  path0_t.insert(path0_t.end(), path_parameter_t.begin(),
      path_parameter_t.begin() + entry + 1);
  assert(path0_t.back() <= 0.5);
  path0_t.push_back(0.5);
  for (int i = 0; i < path0_t.size(); i++) path0_t[i] *= 2;

  path1.start = sample;
  path1.end   = path.end;
  path1.strip.insert(
      path1.strip.end(), path.strip.begin() + entry, path.strip.end());
  path1.lerps.insert(
      path1.lerps.end(), path.lerps.begin() + entry, path.lerps.end());

  path1_t.insert(
      path1_t.end(), path_parameter_t.begin() + entry, path_parameter_t.end());
  assert(path1_t[0] < 0.5);
  path1_t[0] = 0;
  for (int i = 1; i < path1_t.size(); i++) path1_t[i] = (path1_t[i] - 0.5) * 2;

  return std::tie(path0, path1, path0_t, path1_t);
}
// utility (gradient matrix)
Eigen::VectorXd wrapper(const vector<float>& f) {
  Eigen::VectorXd F(f.size());
  for (int i = 0; i < f.size(); ++i) {
    F(i) = f[i];
  }
  return F;
}
Eigen::SparseMatrix<double, 1> PCE_grad_mat(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  Eigen::SparseMatrix<double, 1> G(triangles.size() * 3, positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T>                      entries;

  for (int i = 0; i < triangles.size(); ++i) {
    vec3f  v0   = positions[triangles[i].x];
    vec3f  v1   = positions[triangles[i].y];
    vec3f  v2   = positions[triangles[i].z];
    vec3f  n    = cross(v1 - v0, v2 - v0);
    double area = length(n);  // division by 2 is missing because
                              // there is a *2 in the computation of
                              // the gradient
    n = normalize(n);

    for (int off = 0; off < 3; ++off) {
      int   prev = triangles[i][off];
      int   curr = triangles[i][(off + 1) % 3];
      int   next = triangles[i][(off + 2) % 3];
      vec3f u    = positions[next] - positions[curr];
      vec3f v    = positions[curr] - positions[prev];
      vec3f u_90 = normalize(cross(u, n));
      vec3f v_90 = normalize(cross(v, n));

      vec3f contribute = u_90 * length(u) + v_90 * length(v);
      contribute /= area;

      int row = 3 * i;
      entries.push_back(T(row, curr, contribute.x));
      ++row;
      entries.push_back(T(row, curr, contribute.y));
      ++row;
      entries.push_back(T(row, curr, contribute.z));
    }
  }
  G.setFromTriplets(entries.begin(), entries.end());
  return G;
}
vec3f AGS_gradient(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3f>& positions, const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& gradient, const int& vid,
    const vector<vec3f>& normals) {
  vector<int> nbr = v2t[vid];

  vec3f v   = zero3f;
  float wgt = 0.0;

  for (int i = 0; i < nbr.size(); ++i) {
    int   tid      = nbr[i];
    auto  flat_tid = init_flat_triangle(positions, triangles[tid]);
    vec3f tid_grad = gradient[tid];

    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, tid_grad, normals, tid, vid, T2V);

    vec3f c0 = tid_centroid(triangles, positions, tid);

    float w0 = 1 / (1 + length(c0 - positions[vid]));

    tid_grad *= w0;

    v += tid_grad;

    wgt += w0;

    int  opp      = opposite_face(triangles, adjacencies, tid, vid);
    auto k        = find(adjacencies[tid], opp);
    auto flat_opp = unfold_face(
        triangles, positions, adjacencies, flat_tid, tid, k);
    vec3f opp_grad = gradient[opp];

    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, opp_grad, normals, opp, tid, T2T);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, opp_grad, normals, tid, vid, T2V);

    auto flat_c1 = interpolate_triangle(
        flat_opp[0], flat_opp[1], flat_opp[2], vec2f{0.33, 0.33});
    k        = find(triangles[tid], vid);
    float w1 = 1 / (1 + length(flat_c1 - flat_tid[k]));
    opp_grad *= w1;

    v += opp_grad;
    wgt += w1;
  }

  v /= wgt;
  return -v;
}
Eigen::MatrixXd rhs(int s) {
  Eigen::MatrixXd E(s, s + 1);
  Eigen::MatrixXd X     = Eigen::MatrixXd::Constant(s, 1, -1);
  E.topLeftCorner(s, 1) = X;
  Eigen::MatrixXd I(s, s);
  I.setIdentity();
  E.topRightCorner(s, s) = I;
  return E;
}
void fill_gradient_entries(vector<Eigen::Triplet<double>>& entries,
    const vector<int>& ring, const Eigen::VectorXd& a0,
    const Eigen::VectorXd& a1, const int n) {
  int vid = ring[0];
  int s   = ring.size();

  typedef Eigen::Triplet<double> T;
  for (int i = 0; i < s; ++i) {
    int entry = ring[i];
    entries.push_back(T(vid, entry, a0(i)));
    entries.push_back(T(n + vid, entry, a1(i)));
  }
}
void fill_riemannian_gradient_entries(vector<Eigen::Triplet<double>>& entries,
    const vector<int>& ring, const Eigen::VectorXd& c,
    const Eigen::VectorXd& a0, const Eigen::VectorXd& a1, const int n) {
  int             vid        = ring[0];
  int             s          = ring.size();
  double          c0_squared = pow(c[0], 2);
  double          c1_squared = pow(c[1], 2);
  Eigen::Matrix2d g_inv;
  double          det = 1 + c0_squared + c1_squared;
  g_inv << 1 + c1_squared, -c[0] * c[1], -c[0] * c[1], 1 + c0_squared;
  g_inv /= det;
  typedef Eigen::Triplet<double> T;
  for (int i = 0; i < s; ++i) {
    int entry = ring[i];
    entries.push_back(T(vid, entry, g_inv(0, 0) * a0(i) + g_inv(0, 1) * a1(i)));
    entries.push_back(
        T(n + vid, entry, g_inv(1, 0) * a0(i) + g_inv(1, 1) * a1(i)));
  }
}

// Gradient matrix
// Note: this construction allows the estimation of the gradient of a scalar
// field F through the matrix-dot-vector multiplication Grad*F, although such
// estimation it is not so accurate near sharp features.
// https://link.springer.com/content/pdf/10.1007/s40304-013-0018-2.pdf
Eigen::SparseMatrix<double, 1> init_gradient_matrix(
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vec3f>& normals) {
  typedef Eigen::Triplet<double> T;
  vector<T>                      G_entries;
  int                            V = positions.size();
  Eigen::SparseMatrix<double, 1> Grad;
  for (int i = 0; i < V; ++i) {
    auto            nbr  = solver.graph[i];
    vec3f           vert = positions[i];
    vec3f           n    = normals[i];
    int             s    = nbr.size();
    Eigen::MatrixXd Q(s, 5);
    Eigen::VectorXd h(s);
    vector<int>     ring(s + 1);
    ring[0] = i;

    for (int j = 0; j < s; ++j) {
      int curr     = nbr[j].node;
      ring[j + 1]  = curr;
      float teta   = angles[i][j];
      float d      = nbr[j].length;
      vec2f pos    = vec2f{d * std::cos(teta), d * std::sin(teta)};
      Q(j, 0)      = pos[0];
      Q(j, 1)      = pos[1];
      Q(j, 2)      = pow(pos[0], 2) / 2;
      Q(j, 3)      = pos[0] * pos[1];
      Q(j, 4)      = pow(pos[1], 2) / 2;
      vec3f coords = positions[curr] - vert;
      h(j)         = dot(coords, n);
    }
    Eigen::MatrixXd Qt = Eigen::Transpose<Eigen::MatrixXd>(Q);

    Eigen::MatrixXd                             A = Qt * Q;
    Eigen::MatrixXd                             E = rhs(s);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(A);
    Eigen::MatrixXd                             a(5, s + 1);

    if (dec.isInvertible()) {
      Eigen::MatrixXd inv = A.inverse();
      a                   = inv * Qt * E;

    } else {
      Eigen::MatrixXd Rhsa = Qt * E;
      a                    = dec.solve(Rhsa);
    }
    fill_gradient_entries(G_entries, ring, a.row(0), a.row(1), V);
  }

  Grad.resize(2 * V, V);
  Grad.setFromTriplets(G_entries.begin(), G_entries.end());

  return Grad;
}
Eigen::SparseMatrix<double, 1> init_riemannian_gradient_matrix(
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vec3f>& normals) {
  typedef Eigen::Triplet<double> T;
  vector<T>                      G_entries;
  int                            V = positions.size();
  Eigen::SparseMatrix<double, 1> Grad;
  for (int i = 0; i < V; ++i) {
    auto  nbr  = solver.graph[i];
    vec3f vert = positions[i];
    vec3f n    = normals[i];
    int   s    = nbr.size();

    Eigen::MatrixXd Q(s, 5);
    Eigen::VectorXd h(s);
    vector<int>     ring(s + 1);
    ring[0] = i;

    for (int j = 0; j < s; ++j) {
      int curr     = nbr[j].node;
      ring[j + 1]  = curr;
      float teta   = angles[i][j];
      float d      = nbr[j].length;
      vec2f pos    = vec2f{d * std::cos(teta), d * std::sin(teta)};
      Q(j, 0)      = pos[0];
      Q(j, 1)      = pos[1];
      Q(j, 2)      = pow(pos[0], 2) / 2;
      Q(j, 3)      = pos[0] * pos[1];
      Q(j, 4)      = pow(pos[1], 2) / 2;
      vec3f coords = positions[curr] - vert;
      h(j)         = dot(coords, n);
    }

    Eigen::MatrixXd Qt = Eigen::Transpose<Eigen::MatrixXd>(Q);

    Eigen::MatrixXd                             A = Qt * Q;
    Eigen::MatrixXd                             E = rhs(s);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(A);
    Eigen::VectorXd                             c(5);
    Eigen::MatrixXd                             a(5, s + 1);

    if (dec.isInvertible()) {
      Eigen::MatrixXd inv = A.inverse();
      a                   = inv * Qt * E;
      c                   = inv * Qt * h;
    } else {
      Eigen::MatrixXd Rhsa = Qt * E;
      Eigen::MatrixXd Rhsc = Qt * h;
      a                    = dec.solve(Rhsa);
      c                    = dec.solve(Rhsc);
    }

    fill_riemannian_gradient_entries(G_entries, ring, c, a.row(0), a.row(1), V);
  }

  Grad.resize(2 * V, V);
  Grad.setFromTriplets(G_entries.begin(), G_entries.end());

  return Grad;
}

vec3f polar_to_cartesian(const geodesic_solver& solver,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const double x, const double y, const int vid) {
  vec3f  g   = zero3f;
  vec2f  sol = vec2f{(float)x, (float)y};
  double phi = yocto::atan2(y, x);

  float mag = length(sol);
  vec3f e   = polar_basis(solver, positions, normals, vid);
  g         = rot_vect(e, normals[vid], phi);
  g *= mag;

  return g;
}
vec3f assembling_gradient(const geodesic_solver& solver,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const float g00_contribute, const float g01_contribute,
    const float g10_contribute, const float g11_contribute, const int vid) {
  vec3f  g   = zero3f;
  auto   x   = g00_contribute + g10_contribute;
  auto   y   = g10_contribute + g11_contribute;
  vec2f  sol = vec2f{(float)x, (float)y};
  double phi = yocto::atan2(y, x);

  float mag = length(sol);
  vec3f e   = polar_basis(solver, positions, normals, vid);
  g         = rot_vect(e, normals[vid], phi);
  g *= mag;

  return g;
}
// Compute the gradient of a scalar field f
// Note: by default we consider -grad(f);
vector<vec3f> compute_grad(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& G,
    const Eigen::VectorXd& f, bool normalized) {
  auto V = (int)positions.size();

  auto F = (int)triangles.size();

  if (G.rows() == 2 * V) {
    vector<vec3f> g(V);

    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < V; ++i) {
      g[i] = -polar_to_cartesian(
          solver, positions, normals, Grad(i), Grad(V + i), i);
      if (normalized) g[i] = normalize(g[i]);
    }
    return g;
  } else {
    vector<vec3f>   g(F);
    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < F; ++i) {
      g[i].x = Grad(3 * i);
      g[i].y = Grad(3 * i + 1);
      g[i].z = Grad(3 * i + 2);

      if (normalized) g[i] = normalize(g[i]);
    }
    return g;
  }
}
vector<vec3f> compute_grad(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& G,
    const vector<float>& f, bool normalized) {
  auto F = wrapper(f);
  return compute_grad(solver, triangles, positions, normals, G, F, normalized);
}
vector<vec3f> gradient_of_squared_distance_field(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& G,
    const vector<float>& f, bool normalized) {
  auto g = f;
  std::transform(g.begin(), g.end(), g.begin(),
      [](float lambda) { return lambda * lambda; });
  return compute_grad(solver, triangles, positions, normals, G, g, normalized);
}
vector<vector<vec3f>> gradient_of_squared_distance_field(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double, 1>& G, const vector<vector<float>>& f,
    bool normalized) {
  vector<vector<vec3f>> grad(f.size());
  for (auto i = 0; i < f.size(); ++i) {
    grad[i] = gradient_of_squared_distance_field(
        solver, triangles, positions, normals, G, f[i], normalized);
  }
  return grad;
}
bool set_order(int s, int prev_entry, int next_entry, bool nei_is_dual) {
  auto ccw_count = -1, cw_count = -1;
  if (prev_entry < next_entry) {
    ccw_count = next_entry - prev_entry;
    cw_count  = prev_entry + s - next_entry;
  } else {
    ccw_count = s - prev_entry + next_entry;
    cw_count  = prev_entry - next_entry;
  }
  if (!nei_is_dual) --ccw_count;
  if (ccw_count < cw_count)
    return true;
  else
    return false;
}
void fill_the_strip(vector<int>& strip, const vector<vector<int>>& v2t, int vid,
    int first, int last, bool nei_is_dual, bool ccw) {
  auto  start = first, end = last;
  auto& star = v2t[vid];
  auto  s    = (int)star.size();
  if (ccw && !nei_is_dual)
    end = (s - 1 + end) % s;  // I can stop one face earlier;
  if (start == end) {
    if (strip.back() != star[start]) strip.push_back(star[start]);
  } else if (ccw) {
    if (strip.back() == star[start]) start = (start + 1) % s;
    if (start > end) end += s;
    for (auto i = start; i <= end; ++i) {
      strip.push_back(star[i % s]);
    }
  } else {
    if (strip.back() == star[start % s]) start = (s - 1 + start) % s;
    if (start < end) start += s;
    for (auto i = start; i >= end; --i) {
      strip.push_back(star[i % s]);
    }
  }
}

int strip_to_point(vector<int>& strip, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, int parent, const mesh_point& p) {
  strip = {p.face};
  if (auto [is_vert, offset] = point_is_vert(p); is_vert) {
    auto vid   = triangles[p.face][offset];
    auto entry = node_is_adjacent(solver, vid, parent);
    if (entry < 0) {
      assert(vid == parent);
      return -1;
    }
    auto star        = v2t[vid];
    auto s           = (int)star.size();
    auto it          = find(star.begin(), star.end(), p.face);
    auto first       = (int)distance(star.begin(), it);
    auto last        = (entry % 2) ? (entry - 1) / 2 : (entry / 2) % s;
    auto ccw         = set_order(s, first, last, entry % 2);
    auto nei_is_dual = (bool)((entry + 2) % 2);
    fill_the_strip(strip, v2t, vid, first, last, nei_is_dual, ccw);
    if (entry % 2) {
      first    = (entry - 1) / 2;
      auto tid = opposite_face(triangles, adjacencies, v2t[vid][first], vid);
      strip.push_back(tid);
    }
    entry = node_is_adjacent(solver, parent, vid);
    return entry;
  } else {
    auto        h = find(triangles[p.face], parent);
    vector<int> adj_tri(3);
    if (h == -1) {
      for (auto i = 0; i < 3; ++i) {
        auto adj   = adjacencies[p.face][i];
        adj_tri[i] = adj;
        h          = find(triangles[adj], parent);
        if (h != -1) {
          strip.push_back(adj);
          auto entry = node_is_adjacent(
              solver, parent, triangles[adj][(h + 1) % 3]);
          assert(entry >= 0);
          return entry + 1;
        }
      }
    } else {
      auto entry = node_is_adjacent(
          solver, parent, triangles[p.face][(h + 1) % 3]);
      assert(entry >= 0);
      return entry + 1;
    }
    for (auto i = 0; i < 3; ++i) {
      auto p0 = triangles[p.face][i], p1 = triangles[p.face][(i + 1) % 3];
      auto adj = adj_tri[i];
      auto opp = opposite_face(triangles, adjacencies, adj, p0);
      h        = find(triangles[opp], parent);
      if (h != -1) {
        strip.push_back(adj);
        strip.push_back(opp);
        auto entry = node_is_adjacent(
            solver, parent, triangles[opp][(h + 1) % 3]);
        assert(entry >= 0);
        return entry + 1;
      }
      opp = opposite_face(triangles, adjacencies, adj, p1);
      h   = find(triangles[opp], parent);
      if (h != -1) {
        strip.push_back(adj);
        strip.push_back(opp);
        auto entry = node_is_adjacent(
            solver, parent, triangles[opp][(h + 1) % 3]);
        assert(entry >= 0);
        return entry + 1;
      }
    }
  }
  assert(false);
  return 0;  // TODO(fabio): cosa deve fare qui?
}

void close_the_strip(vector<int>& strip, const vector<vector<int>>& v2t,
    int vid, int prev_tri, int last_tri) {
  auto star    = v2t[vid];
  auto s       = (int)star.size();
  auto prev_it = find(star.begin(), star.end(), prev_tri);
  assert(prev_it != star.end());
  auto first   = (int)distance(star.begin(), prev_it);
  auto next_it = find(star.begin(), star.end(), last_tri);
  assert(next_it != star.end());
  auto last = (int)distance(star.begin(), next_it);
  auto ccw  = set_order(s, first, last, true);
  fill_the_strip(strip, v2t, vid, first, last, true, ccw);
}

// particular case of "get strip" when one of the two point is the parent of
// the other so the size of the strip is one or two
static vector<int> short_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  auto entry = 0;
  auto strip = vector<int>{};
  if (auto [is_vert, offset] = point_is_vert(target); is_vert) {
    auto vid = triangles[target.face][offset];
    entry    = strip_to_point(strip, solver, triangles, positions, adjacencies,
           v2t, angles, vid, source);
    if (entry < 0) return {target.face};
    if (strip.back() != target.face) {
      close_the_strip(strip, v2t, vid, strip.back(), target.face);
    }
    reverse(strip.begin(), strip.end());
  } else if (auto [is_vert, offset] = point_is_vert(source); is_vert) {
    auto vid = triangles[source.face][offset];
    entry    = strip_to_point(strip, solver, triangles, positions, adjacencies,
           v2t, angles, vid, target);
    if (strip.back() != source.face) {
      close_the_strip(strip, v2t, vid, strip.back(), source.face);
    }
  } else {
    assert(false);
  }
  return strip;
}

// returns a strip of triangles such target belongs to the first one and
// source to the last one
// TODO(fabio_): may be the names could change in order to get the call
// more consistent with the output)
vector<int> get_strip_having_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& distances,
    const mesh_point& source, const mesh_point& target, vector<int>& parents) {
  if (target.face == source.face) return {target.face};
  parents = parents_from_distance_field(
      solver, triangles, positions, adjacencies, distances, source, target);
  auto N     = (int)parents.size();
  auto first = 0, last = 0, prev_entry = 0, next_entry = 0;
  auto strip_to_mesh_point = vector<int>{}, strip = vector<int>{};
  auto ccw = false, nei_is_dual = false;
  if (N == 0) {
    return short_strip(
        solver, triangles, positions, adjacencies, v2t, angles, source, target);
  } else if (N == 1) {
    auto v      = parents[0];
    prev_entry  = strip_to_point(strip, solver, triangles, positions,
         adjacencies, v2t, angles, v, target);
    next_entry  = strip_to_point(strip_to_mesh_point, solver, triangles,
         positions, adjacencies, v2t, angles, v, source);
    first       = find(v2t[v], strip.back());
    nei_is_dual = next_entry % 2;
    last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
    ccw         = set_order((int)v2t[v].size(), first, last, nei_is_dual);
    fill_the_strip(strip, v2t, v, first, last, nei_is_dual, ccw);
    close_the_strip(strip, v2t, v, strip.back(), strip_to_mesh_point.back());
    if (strip.back() == strip_to_mesh_point.back())
      strip_to_mesh_point.pop_back();
    reverse(strip_to_mesh_point.begin(), strip_to_mesh_point.end());
    strip.insert(
        strip.end(), strip_to_mesh_point.begin(), strip_to_mesh_point.end());
    return strip;
  } else {
    prev_entry = strip_to_point(strip, solver, triangles, positions,
        adjacencies, v2t, angles, parents[0], target);
    assert(prev_entry != 0);
  }

  for (auto i = 0; i < N; ++i) {
    auto v = parents[i];
    if (i == N - 1) {
      first      = find(v2t[v], strip.back());
      next_entry = strip_to_point(strip_to_mesh_point, solver, triangles,
          positions, adjacencies, v2t, angles, v, source);
      last       = find(v2t[v], strip_to_mesh_point.back());
      ccw        = set_order((int)v2t[v].size(), first, last, next_entry % 2);
    } else {
      first = find(v2t[v], strip.back());
      assert(first != -1);
      next_entry  = node_is_adjacent(solver, v, parents[i + 1]);
      nei_is_dual = (next_entry + 2) % 2;
      last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
      ccw         = set_order((int)v2t[v].size(), first, last, nei_is_dual);
    }

    fill_the_strip(strip, v2t, v, first, last, nei_is_dual, ccw);

    if (nei_is_dual && i != N - 1) {
      auto tid = opposite_face(triangles, adjacencies, strip.back(), v);
      strip.push_back(tid);
    }
  }

  close_the_strip(
      strip, v2t, parents.back(), strip.back(), strip_to_mesh_point.back());
  if (strip.back() == strip_to_mesh_point.back())
    strip_to_mesh_point.pop_back();
  reverse(strip_to_mesh_point.begin(), strip_to_mesh_point.end());
  strip.insert(
      strip.end(), strip_to_mesh_point.begin(), strip_to_mesh_point.end());
  return strip;
}
std::tuple<vector<int>, vector<float>> get_strip_with_distances(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const mesh_point& source, const mesh_point& target) {
  if (target.face == source.face) return {{target.face}, {}};
  auto [parents, distances] = get_parents_and_distances(
      solver, triangles, positions, adjacencies, source, target);
  auto N     = (int)parents.size();
  auto first = 0, last = 0, prev_entry = 0, next_entry = 0;
  auto strip_to_mesh_point = vector<int>{}, strip = vector<int>{};
  auto ccw = false, nei_is_dual = false;
  if (N == 0) {
    return {short_strip(solver, triangles, positions, adjacencies, v2t, angles,
                source, target),
        distances};
  } else if (N == 1) {
    auto v      = parents[0];
    prev_entry  = strip_to_point(strip, solver, triangles, positions,
         adjacencies, v2t, angles, v, target);
    next_entry  = strip_to_point(strip_to_mesh_point, solver, triangles,
         positions, adjacencies, v2t, angles, v, source);
    first       = find(v2t[v], strip.back());
    nei_is_dual = next_entry % 2;
    last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
    ccw         = set_order((int)v2t[v].size(), first, last, nei_is_dual);
    fill_the_strip(strip, v2t, v, first, last, nei_is_dual, ccw);
    close_the_strip(strip, v2t, v, strip.back(), strip_to_mesh_point.back());
    if (strip.back() == strip_to_mesh_point.back())
      strip_to_mesh_point.pop_back();
    reverse(strip_to_mesh_point.begin(), strip_to_mesh_point.end());
    strip.insert(
        strip.end(), strip_to_mesh_point.begin(), strip_to_mesh_point.end());
    return {strip, distances};
  } else {
    prev_entry = strip_to_point(strip, solver, triangles, positions,
        adjacencies, v2t, angles, parents[0], target);
  }

  for (auto i = 0; i < N; ++i) {
    auto v = parents[i];
    if (i == N - 1) {
      first      = find(v2t[v], strip.back());
      next_entry = strip_to_point(strip_to_mesh_point, solver, triangles,
          positions, adjacencies, v2t, angles, v, source);
      last       = find(v2t[v], strip_to_mesh_point.back());
      ccw        = set_order((int)v2t[v].size(), first, last, next_entry % 2);
    } else {
      first = find(v2t[v], strip.back());
      assert(first != -1);
      next_entry  = node_is_adjacent(solver, v, parents[i + 1]);
      nei_is_dual = next_entry % 2;
      last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
      ccw         = set_order((int)v2t[v].size(), first, last, nei_is_dual);
    }

    fill_the_strip(strip, v2t, v, first, last, nei_is_dual, ccw);

    if (nei_is_dual && i != N - 1) {
      auto tid = opposite_face(triangles, adjacencies, strip.back(), v);
      strip.push_back(tid);
    }
  }

  close_the_strip(
      strip, v2t, parents.back(), strip.back(), strip_to_mesh_point.back());
  if (strip.back() == strip_to_mesh_point.back())
    strip_to_mesh_point.pop_back();
  reverse(strip_to_mesh_point.begin(), strip_to_mesh_point.end());
  strip.insert(
      strip.end(), strip_to_mesh_point.begin(), strip_to_mesh_point.end());
  return {strip, distances};
}
vec3f path_pos_from_entry(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, int entry) {
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
// utility (Karcher algorithm)
void correct_geodesic_distances(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const mesh_point& source,
    const mesh_point& target, vector<float>& distances) {
  auto surrounding_nodes = nodes_around_mesh_point(
      triangles, positions, adjacencies, source);
  if (surrounding_nodes.size() > 1) check_surrounding_nodes(surrounding_nodes);

  for (auto i = 0; i < surrounding_nodes.size(); ++i) {
    auto vid         = surrounding_nodes[i].first;
    auto tid         = v2t[vid][0];
    auto k           = find(triangles[tid], vid);
    auto bary        = zero3f;
    bary[k]          = 1;
    auto curr_source = mesh_point{tid, vec2f{bary.y, bary.z}};
    distances[vid]   = point2point_geodesic_distance(
          triangles, positions, adjacencies, distances, curr_source, target);
  }
}
vector<unfold_triangle> one_ring_flattening(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const int vid) {
  auto s    = v2t[vid].size();
  auto curr = v2t[vid][0];

  auto strip = vector<unfold_triangle>(s);
  strip[0]   = flat_triangle_with_origin_at_vertex(
        triangles, positions, triangles[curr], vid);
  for (auto i = 1; i < s; ++i) {
    auto k   = find(triangles[curr], vid);
    strip[i] = unfold_face(
        triangles, positions, adjacencies, strip[i - 1], curr, (k + 2) % 3);
    curr = adjacencies[curr][(k + 2) % 3];
  }
  return strip;
}

vec3f mid_point_butterfly_scheme(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<unfold_triangle>& flat_one_ring, const int vid,
    const int vid0) {
  auto control_points = vector<vec2f>(8);
  auto weights        = vector<float>(8);
  auto entry          = node_is_adjacent(solver, vid, vid0);
  entry /= 2;
  auto tid          = v2t[vid][entry];
  auto k            = find(triangles[tid], vid);
  auto s            = v2t[vid].size();
  control_points[0] = flat_one_ring[entry][k];
  weights[0] = weights[1] = 0.5;
  control_points[1]       = flat_one_ring[entry][(k + 1) % 3];
  tid                     = adjacencies[tid][(k + 2) % 3];
  k                       = find(triangles[tid], vid);
  control_points[2]       = flat_one_ring[(entry + 1) % s][(k + 1) % 3];
  weights[2]              = 1.f / 8;
  control_points[3]       = flat_one_ring[(entry + 1) % s][(k + 2) % 3];
  weights[3]              = -1.f / 16;
  tid                     = v2t[vid][(s - 2 + entry) % s];
  k                       = find(triangles[tid], vid);
  control_points[4]       = flat_one_ring[(s - 2 + entry) % s][(k + 1) % 3];
  weights[4]              = -1.f / 16;
  control_points[5]       = flat_one_ring[(s - 2 + entry) % s][(k + 2) % 3];
  weights[5]              = 1.f / 8;
  tid                     = v2t[vid][entry];
  k                       = find(triangles[tid], vid);
  auto flat0              = unfold_face(triangles, positions, adjacencies,
                   flat_one_ring[entry], tid, (k + 1) % 3);
  tid                     = adjacencies[tid][(k + 1) % 3];
  k                       = find(triangles[tid], vid0);
  control_points[6]       = flat0[(k + 1) % 3];
  weights[6]              = -1.f / 16;
  tid                     = v2t[vid][(s - 1 + entry) % s];
  k                       = find(triangles[tid], vid);
  auto flat1              = unfold_face(triangles, positions, adjacencies,
                   flat_one_ring[(s - 1 + entry) % s], tid, (k + 1) % 3);
  tid                     = adjacencies[tid][(k + 1) % 3];
  k                       = find(triangles[tid], vid0);
  control_points[7]       = flat1[(k + 2) % 3];
  weights[7]              = -1.f / 16;

  auto mid_point = zero2f;
  for (auto i = 0; i < 8; ++i) {
    mid_point += weights[i] * control_points[i];
  }

  auto alpha = length(mid_point);

  tid               = v2t[vid][entry];
  k                 = find(triangles[tid], vid0);
  auto l            = length(flat_one_ring[entry][k]);
  auto lerp         = alpha / l;
  auto bary         = zero3f;
  bary[k]           = lerp;
  bary[(k + 1) % 3] = 0;
  bary[(k + 2) % 3] = 1 - lerp;
  return bary;
}
float average_butterfly_scheme(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<float>& f, const int vid, const int vid0) {
  auto control_values = vector<float>(8);
  auto weights        = vector<float>(8);
  auto entry          = node_is_adjacent(solver, vid, vid0);
  entry /= 2;
  auto tid          = v2t[vid][entry];
  auto k            = find(triangles[tid], vid);
  auto s            = v2t[vid].size();
  control_values[0] = f[vid];
  weights[0] = weights[1] = 0.5;
  control_values[1]       = f[vid0];
  tid                     = adjacencies[tid][(k + 2) % 3];
  k                       = find(triangles[tid], vid);
  control_values[2]       = f[triangles[tid][(k + 1) % 3]];
  weights[2]              = 1.f / 8;
  control_values[3]       = f[triangles[tid][(k + 2) % 3]];
  weights[3]              = -1.f / 16;
  tid                     = v2t[vid][(s - 2 + entry) % s];
  k                       = find(triangles[tid], vid);
  control_values[4]       = f[triangles[tid][(k + 1) % 3]];
  weights[4]              = -1.f / 16;
  control_values[5]       = f[triangles[tid][(k + 2) % 3]];
  weights[5]              = 1.f / 8;
  tid                     = v2t[vid][entry];
  k                       = find(triangles[tid], vid);
  tid                     = adjacencies[tid][(k + 1) % 3];
  k                       = find(triangles[tid], vid0);
  control_values[6]       = f[triangles[tid][(k + 1) % 3]];
  weights[6]              = -1.f / 16;
  tid                     = v2t[vid][(s - 1 + entry) % s];
  k                       = find(triangles[tid], vid);
  tid                     = adjacencies[tid][(k + 1) % 3];
  k                       = find(triangles[tid], vid0);
  control_values[7]       = f[triangles[tid][(k + 2) % 3]];
  weights[7]              = -1.f / 16;

  auto avg = 0.f;
  for (auto i = 0; i < 8; ++i) {
    avg += weights[i] * control_values[i];
  }

  return avg;
}
// vid is supposed to be regular
vec3f mid_point_general_scheme(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<unfold_triangle>& flat_one_ring, const int vid, const int vid0,
    const int k_vid0) {
  auto control_points = vector<vec2f>(k_vid0);
  auto weights        = vector<float>(k_vid0);
  auto entry          = node_is_adjacent(solver, vid, vid0);
  entry /= 2;
  auto tid = v2t[vid][entry];
  auto k   = find(triangles[tid], vid);
  auto s   = v2t[vid].size();

  if (k_vid0 != 4) {
    for (auto i = 0; i < k_vid0; ++i) {
      tid               = v2t[vid][(entry + i) % s];
      k                 = find(triangles[tid], vid0);
      control_points[i] = flat_one_ring[(entry + i) % s][(k + 1) % s];
      weights[i]        = 1.f / k_vid0 *
                   (0.25 + yocto::cos(2 * pif * i / k_vid0) +
                       0.5 * yocto::cos(4 * pif * i / k_vid0));
    }
  } else {
    for (auto i = 0; i < 4; ++i) {
      tid               = v2t[vid][(entry + i) % s];
      k                 = find(triangles[tid], vid0);
      control_points[i] = flat_one_ring[(entry + i) % s][(k + 1) % s];
    }
    weights[8]  = 3.f / 8;
    weights[10] = -1.f / 8;
  }

  auto mid_point = zero2f;
  for (auto i = 0; i < weights.size(); ++i) {
    mid_point += weights[i] * control_points[i];
  }

  auto alpha = length(mid_point);

  tid               = v2t[vid][entry];
  k                 = find(triangles[tid], vid0);
  auto l            = length(flat_one_ring[entry][k]);
  auto lerp         = alpha / l;
  auto bary         = zero3f;
  bary[k]           = lerp;
  bary[(k + 1) % 3] = 0;
  bary[(k + 2) % 3] = 1 - lerp;
  return bary;
}
float average_general_scheme(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<float>& f, const int vid, const int vid0, const int k_vid0) {
  auto control_points = vector<float>(k_vid0);
  auto weights        = vector<float>(k_vid0);
  auto entry          = node_is_adjacent(solver, vid0, vid);
  entry /= 2;
  auto tid = v2t[vid0][entry];
  auto k   = find(triangles[tid], vid0);
  auto s   = v2t[vid0].size();

  if (k_vid0 != 4) {
    for (auto i = 0; i < k_vid0; ++i) {
      tid               = v2t[vid0][(entry + i) % s];
      k                 = find(triangles[tid], vid0);
      control_points[i] = f[triangles[tid][(k + 1) % 3]];
      weights[i]        = 1.f / k_vid0 *
                   (0.25 + yocto::cos(2 * pif * i / k_vid0) +
                       0.5 * yocto::cos(4 * pif * i / k_vid0));
    }
  } else {
    for (auto i = 0; i < 4; ++i) {
      tid               = v2t[vid0][(entry + i) % s];
      k                 = find(triangles[tid], vid0);
      control_points[i] = f[triangles[tid][(k + 1) % 3]];
    }
    weights[8]  = 3.f / 8;
    weights[10] = -1.f / 8;
  }

  auto avg = 0.f;
  for (auto i = 0; i < weights.size(); ++i) {
    avg += weights[i] * control_points[i];
  }

  return avg;
}
// void butterfly_control_points_and_weights(const geodesic_solver& solver,
//     const vector<vector<float>>& angles, const int vid, const int vid0,
//     vector<vec2f>& control_points, vector<float>& weights) {
//   auto nbr    = solver.graph[vid];
//   auto s      = nbr.size();
//   auto entry  = node_is_adjacent(solver, vid, vid0);
//   auto thetas = angles[vid];
//   control_points.resize(8);
//   weights.resize(8);
//   control_points[0] = vec2f{0, 0};
//   weights[0] = weights[1] = 0.5;
//   auto r1                 = solver.graph[vid][entry].length;
//   control_points[1]       = vec2f{
//       r1 * yocto::cos(thetas[entry]), r1 * yocto::sin(thetas[entry])};
//   auto r2           = solver.graph[vid][(entry + 1) % s].length;
//   control_points[2] = vec2f{r2 * yocto::cos(thetas[(entry + 1) % s]),
//       r2 * yocto::sin(thetas[(entry + 1) % s])};
//   weights[2]        = -1.f / 16;
//   auto r3           = solver.graph[vid][(entry + 2) % s].length;
//   control_points[3] = vec2f{r3 * yocto::cos(thetas[(entry + 2) % s]),
//       r3 * yocto::sin(thetas[(entry + 2) % s])};
//   weights[3]        = 1.f / 8;
//   auto r4           = solver.graph[vid][(entry + 4) % s].length;
//   control_points[4] = vec2f{r4 * yocto::cos(thetas[(entry + 4) % s]),
//       r4 * yocto::sin(thetas[(entry + 4) % s])};
//   weights[4]        = -1.f / 16;
//   auto r5           = solver.graph[vid][(s - 4 + entry) % s].length;
//   control_points[5] = vec2f{r5 * yocto::cos(thetas[(s - 4 + entry) % s]),
//       r5 * yocto::sin(thetas[(s - 4 + entry) % s])};
//   weights[5]        = -1.f / 16;
//   auto r6           = solver.graph[vid][(s - 2 + entry) % s].length;
//   control_points[6] = vec2f{r6 * yocto::cos(thetas[(s - 2 + entry) % s]),
//       r6 * yocto::sin(thetas[(s - 2 + entry) % s])};
//   weights[6]        = 1.f / 8;
//   auto r7           = solver.graph[vid][(s - 1 + entry) % s].length;
//   control_points[7] = vec2f{r7 * yocto::cos(thetas[(s - 1 + entry) % s]),
//       r7 * yocto::sin(thetas[(s - 1 + entry) % s])};
//   weights[7]        = -1.f / 16;
// }
// void mid_point_general_case(const geodesic_solver& solver,
//     const vector<vector<float>>& angles, const int vid, const int vid0,
//     vector<vec2f>& control_points, vector<float>& weights) {
//   auto nbr = solver.graph[vid];
//   auto s   = (int)nbr.size();
//   auto K   = s / 2;

//   control_points.resize(K);
//   weights.resize(K);
//   auto entry = node_is_adjacent(solver, vid, vid0);
//   if (K != 4) {
//     for (auto i = 0; i < K; ++i) {
//       auto r = nbr[(entry + 2 * i) % s].length;

//       control_points[i] = vec2f{
//           r * yocto::cos(angles[vid][(entry + 2 * i) % s]),
//           r * yocto::sin(angles[vid][(entry + 2 * i) % s])};

//       weights[i] = 1.f / K *
//                    (0.25 + yocto::cos(2 * pif * i / K) +
//                        0.5 * yocto::cos(4 * pif * i / K));
//     }
//   } else {
//     auto r            = nbr[entry].length;
//     control_points[0] = vec2f{
//         r * yocto::cos(angles[vid][entry]), r *
//         yocto::sin(angles[vid][entry])};
//     weights[0]        = 3.f / 8;
//     control_points[1] = zero2f;
//     weights[1]        = 0;
//     r                 = nbr[(entry + 2) % s].length;
//     control_points[2] = vec2f{r * yocto::cos(angles[vid][(entry + 2) % s]),
//         r * yocto::sin(angles[vid][(entry + 2) % s])};
//     weights[2]        = -1.f / 8;
//     control_points[3] = zero2f;
//     weights[3]        = 0;
//   }
// }
// vec3f mid_point_general_case(const vector<vec3i>& triangles,
//     const vector<vec3f>& positions, const vector<vector<int>>& v2t,
//     const vector<vec3f>& normals, const geodesic_solver& solver,
//     const vector<vector<float>>& angles, const int vid, const int vid0) {
//   auto nbr = solver.graph[vid];
//   auto s   = (int)nbr.size();
//   auto K   = s / 2;

//   auto control_points = vector<vec2f>(K);
//   auto weights        = vector<float>(K);
//   auto entry          = node_is_adjacent(solver, vid, vid0);
//   if (K != 4) {
//     for (auto i = 0; i < K; ++i) {
//       auto r = solver.graph[vid][(entry + 2 * i) % s].length;

//       control_points[i] = vec2f{
//           r * yocto::cos(angles[vid][(entry + 2 * i) % s]),
//           r * yocto::sin(angles[vid][(entry + 2 * i) % s])};
//       // positions[solver.graph[vid][(entry + 2 * i) % s].node];
//       weights[i] = 1.f / K *
//                    (0.25 + yocto::cos(2 * pif * i / K) +
//                        0.5 * yocto::cos(4 * pif * i / K));
//     }
//   } else {
//     // control_points[0] = positions[vid0];
//     // weights[0]        = 3.f / 8;
//     // control_points[1] = zero3f;
//     // weights[1]        = 0;
//     // control_points[2] = positions[solver.graph[vid][(entry + 2) %
//     s].node];
//     // weights[2]        = -1.f / 8;
//     // control_points[3] = zero3f;
//     // weights[3]        = 0;
//   }

//   auto mid_point = zero2f;
//   for (auto i = 0; i < K; ++i) {
//     mid_point += weights[i] * control_points[i];
//   }
//   auto d = length(mid_point);

//   // double phi = yocto::atan2(mid_point.y, mid_point.x);

//   // vec3f e       = polar_basis(solver, positions, normals, vid);
//   // auto  vec_mid = rot_vect(e, normals[vid], phi);
//   // auto  tid     = next_tid(
//   //     solver, angles, positions, v2t, triangles, normals, vid, vec_mid);
//   // auto offset = find(triangles[tid], vid);
//   // auto p0     = triangles[tid][(offset + 1) % 3];
//   // auto p1     = triangles[tid][(offset + 2) % 3];
//   // auto v      = positions[p0] - positions[vid];
//   // auto teta3D = angle(v, positions[p1] - positions[p1]);
//   // auto entry0 = node_is_adjacent(solver, vid, p0);
//   // auto teta0  = angles[vid][entry0];
//   // auto teta2D = (entry0 + 2 == s) ? 2 * pif - teta0
//   //                                 : angles[vid][entry0 + 2] - teta0;
//   // auto scale_factor = teta2D / teta3D;
//   // auto n            = tid_normal(positions, triangles, tid);
//   // auto phi3D        = (phi - teta0) * scale_factor;
//   // auto w            = normalize(rot_vect(v, n, phi3D));
//   // w *= d;
//   // w += positions[vid];

//   return (1 - d) * positions[vid] + d * positions[vid0];
// }
// note:vid0 vid1 must be in CCW order;returns the barycentric coordinates of
// midpoints ordered in CCW order
vector<vec3f> subdivide_tid(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const int vid, const int vid0, const int vid1) {
  assert(solver.graph[vid].size() % 2 == 0);
  auto k_vid  = (int)solver.graph[vid].size() / 2;
  auto k_vid0 = (int)solver.graph[vid0].size() / 2;
  auto k_vid1 = (int)solver.graph[vid1].size() / 2;

  auto m0 = zero3f, m1 = zero3f, m01 = zero3f;

  auto flat_one_ring_vid = one_ring_flattening(
      triangles, positions, adjacencies, v2t, vid);

  auto flat_one_ring_vid0 = one_ring_flattening(
      triangles, positions, adjacencies, v2t, vid0);

  auto flat_one_ring_vid1 = one_ring_flattening(
      triangles, positions, adjacencies, v2t, vid1);

  if (k_vid == 6) {
    if (k_vid0 == 6) {
      m0 = mid_point_butterfly_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid, vid, vid0);
    } else {
      m0 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid, vid, vid0, k_vid0);
    }

    if (k_vid1 == 6)
      m1 = mid_point_butterfly_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid, vid, vid1);
    else
      m1 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid, vid, vid1, k_vid1);

  } else {
    if (k_vid0 == 6)
      m0 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid0, vid0, vid, k_vid);
    else {
      auto MA0 = mid_point_general_scheme(triangles, positions, adjacencies,
          v2t, solver, flat_one_ring_vid, vid, vid0, k_vid0);
      auto MA1 = mid_point_general_scheme(triangles, positions, adjacencies,
          v2t, solver, flat_one_ring_vid0, vid0, vid, k_vid);
      m0       = (MA0 + MA1) / 2;
    }

    if (k_vid1 == 6)
      m1 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid1, vid1, vid, k_vid);
    else {
      auto MA0 = mid_point_general_scheme(triangles, positions, adjacencies,
          v2t, solver, flat_one_ring_vid, vid, vid1, k_vid1);
      auto MA1 = mid_point_general_scheme(triangles, positions, adjacencies,
          v2t, solver, flat_one_ring_vid1, vid1, vid, k_vid);
      m1       = (MA0 + MA1) / 2;
    }
  }

  if (k_vid0 == 6 && k_vid1 == 6) {
    m01 = mid_point_butterfly_scheme(triangles, positions, adjacencies, v2t,
        solver, flat_one_ring_vid0, vid0, vid1);

  } else if (k_vid0 == 6 || k_vid1 == 6) {
    if (k_vid0 == 6)
      m01 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid0, vid0, vid1, k_vid1);
    else
      m01 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
          solver, flat_one_ring_vid1, vid1, vid0, k_vid0);

  } else {
    auto MA0 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
        solver, flat_one_ring_vid0, vid0, vid1, k_vid1);
    auto MA1 = mid_point_general_scheme(triangles, positions, adjacencies, v2t,
        solver, flat_one_ring_vid1, vid1, vid0, k_vid0);
    m01      = (MA0 + MA1) / 2;
  }

  return {m0, m01, m1};
}

float average_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<float>& f, const int vid0, const int vid1) {
  auto k_vid0 = (int)solver.graph[vid0].size() / 2;
  auto k_vid1 = (int)solver.graph[vid1].size() / 2;
  auto f01    = 0.f;

  if (k_vid0 == 6 && k_vid1 == 6)
    f01 = average_butterfly_scheme(
        triangles, positions, adjacencies, v2t, solver, f, vid0, vid1);
  else if (k_vid0 == 6)
    f01 = average_general_scheme(
        triangles, positions, adjacencies, v2t, solver, f, vid0, vid1, k_vid1);
  else if (k_vid1 == 6)
    f01 = average_general_scheme(
        triangles, positions, adjacencies, v2t, solver, f, vid1, vid0, k_vid0);
  else {
    auto MA0 = average_general_scheme(
        triangles, positions, adjacencies, v2t, solver, f, vid0, vid1, k_vid1);
    auto MA1 = average_general_scheme(
        triangles, positions, adjacencies, v2t, solver, f, vid1, vid0, k_vid0);
    f01 = (MA0 + MA1) / 2;
  }

  return f01;
}
vector<bool> control_points_moved(
    const vector<mesh_point>& old_ones, const vector<mesh_point>& curr_ones) {
  auto changes = vector<bool>(4, true);
  for (auto i = 0; i < old_ones.size(); ++i) {
    if (old_ones[i].face == -1) continue;
    if (old_ones[i].face == curr_ones[i].face) {
      if (old_ones[i].uv.x == curr_ones[i].uv.x &&
          old_ones[i].uv.y == curr_ones[i].uv.y)
        changes[i] = false;
    }
  }
  return changes;
}
// Binomial Coefficient
int bin_coefficient(const int& n, const int& k) {
  int bin = 1;

  for (int i = 0; i < k; ++i) {
    bin *= (n - i);
    bin /= (i + 1);
  }

  return bin;
}
// Bernstein polynomials
void bernstein_polynomials(const int& n, const float& t, vector<float>& w) {
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    float lambda = bin_coefficient(n - 1, i) * pow(t, i) *
                   pow(1 - t, n - i - 1);
    w[i] = lambda;
    // show(w[i]);
  }
}
float delta_angles(const float& phi_i, const float& phi_j) {
  if (phi_j - phi_i < -pif) return phi_j - phi_i + 2 * pif;
  if (yocto::abs(phi_j - phi_i) < pif) return phi_j - phi_i;
  if (phi_j - phi_i > pif) return phi_j - phi_i - 2 * pif;
}
// ex, ey are the values of the vector field at the common edge while v0 and v1
// are the ones on the opposite vertices.
// note: we assume that ex,ey and v0,v1 already belong to the tangent space of
// tid0

int index_at_edge(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid0, const int tid1,
    const vec3f& w, const vec3f& ex, const vec3f& ey, const vec3f& v0,
    const vec3f& v1) {
  auto flat0    = init_flat_triangle(positions, triangles[tid0]);
  auto centroid = tid_centroid(triangles, positions, tid0);
  auto e3d      = normalize(positions[triangles[tid0].x] - centroid);
  auto c        = interpolate_triangle(
             flat0[0], flat0[1], flat0[2], vec2f{0.33, 0.33});
  auto e         = normalize(flat0[0] - c);
  auto n         = tid_normal(triangles, positions, tid0);
  auto flat_ones = vector<vec2f>(4);
  auto vectors   = vector<vec3f>{ex, ey, v0, v1};
  auto angles    = vector<float>(4);
  for (auto i = 0; i < 3; ++i) {
    auto teta = angle(e3d, vectors[i]);
    if (dot(cross(e3d, vectors[i]), n) < 0) teta = 2 * pif - teta;
    angles[i]    = teta;
    flat_ones[i] = rot_vect(e, teta);
  }
  auto index = 0;
  while (angles.size() != 0) {
    for (auto i = 0; i < angles.size(); ++i) {
      if (i == angles.size() - 1) continue;
      index += delta_angles(angles[i], angles.back());
    }
    angles.pop_back();
  }
  index /= 2 * pif;
  return index;
}
// Check if there is a singularity inside a triangle
bool bary_for_minimum(const vec3f& i, const vec3f& j, const vec3f& k,
    float& alpha, float& beta, float& gamma) {
  double d = i.x * j.y - i.x * k.y - i.y * j.x + i.y * k.x + j.x * k.y -
             j.y * k.x;
  alpha = beta = gamma = 0;
  assert(d != 0);
  alpha = (j.x * k.y - j.y * k.x) / d;
  beta  = (k.x * i.y - k.y * i.x) / d;
  gamma = (i.x * j.y - i.y * j.x) / d;
  if (alpha == 0 && beta == 0 && gamma == 0) assert(false);
  if (alpha >= 0 && beta >= 0 && gamma >= 0) return true;

  return false;
}
bool bary_for_minimum_energy(const vec3f& i, const vec3f& j, const vec3f& k,
    float& alpha, float& beta, float& gamma) {
  auto g0_norm = length_squared(i);
  auto g01     = dot(i, j);
  auto g02     = dot(i, k);
  auto g12     = dot(j, k);
  auto A       = (g0_norm + g12 - g02 - g01);
  auto det     = length_squared(i - k) * length_squared(i - j) - pow(A, 2);
  if (det == 0) assert(false);
  gamma = ((g0_norm - g02) * length_squared(i - j) + A * (g01 - g0_norm)) / det;
  beta  = (g0_norm - g01 - gamma * A) / length_squared(i - j);
  alpha = 1 - beta - gamma;
  if (alpha == 0 && beta == 0 && gamma == 0) assert(false);
  if (are_barycentric_coordinates(vec3f{alpha, beta, gamma}, 0.05)) {
    yocto::clamp(beta, 0.f, 1.f);
    yocto::clamp(gamma, 0.f, 1.f);
    alpha = 1 - beta - gamma;
    return true;
  }

  return false;
}
bool bary_for_minimum_energy(const vec2f& i, const vec2f& j, const vec2f& k,
    float& alpha, float& beta, float& gamma) {
  auto g0_norm = length_squared(i);
  auto g01     = dot(i, j);
  auto g02     = dot(i, k);
  auto g12     = dot(j, k);
  auto A       = (g0_norm + g12 - g02 - g01);
  auto det     = length_squared(i - k) * length_squared(i - j) - pow(A, 2);
  if (det == 0) return false;
  gamma = ((g0_norm - g02) * length_squared(i - j) + A * (g01 - g0_norm)) / det;
  beta  = (g0_norm - g01 - gamma * A) / length_squared(i - j);
  alpha = 1 - beta - gamma;
  if (alpha == 0 && beta == 0 && gamma == 0) return (false);
  if (alpha >= -1e-2 && beta >= -1e-2 && gamma >= -1e-2 && alpha <= 1 + 1e-2 &&
      beta <= 1 + 1e-2 && gamma <= 1 + 1e-2)
    return true;

  return false;
}
pair<bool, int> bary_for_minimum_hessian(const vec3f& i, const vec3f& j,
    const vec3f& k, float& alpha, float& beta, float& gamma) {
  auto g0_norm = length_squared(i);
  auto g1_norm = length_squared(j);
  auto g2_norm = length_squared(k);
  auto g01     = dot(i, j);
  auto g02     = dot(i, k);
  auto g12     = dot(j, k);
  auto l01     = length_squared(i - j);
  auto l02     = length_squared(i - k);

  auto            A = g0_norm + g12 - g02 - g01;
  Eigen::Matrix3d A3;
  A3 << g0_norm, g01, g02, g01, g1_norm, g12, g02, g12, g2_norm;

  auto det     = l02 * l01 - pow(A, 2);
  auto I_3     = (g01 - g0_norm) * (det - (g01 - g0_norm) * (g02 - g0_norm));
  auto det_A_3 = g0_norm * length_squared(cross(j, k)) +
                 g1_norm * length_squared(cross(i, k)) -
                 g2_norm * length_squared(cross(i, j));
  auto contribute = dot(cross(i, k), cross(j, i)) +
                    dot(cross(j, k), cross(i, j)) +
                    dot(cross(k, j), cross(i, k));
  auto det_A_2 = length_squared(cross(i, j)) + length_squared(cross(i, k)) +
                 length_squared(cross(j, k)) + 2 * contribute;
  auto simplified_det_A2 = cross(i, k) + cross(j, i) + cross(k, j);
  auto cond_A = (l01 + l02 + sqrt(pow(l01 - l02, 2) - 4 * pow(A, 2))) /
                (l01 + l02 - sqrt(pow(l01 - l02, 2) - 4 * pow(A, 2)));
  cond_A        = std::abs(cond_A);
  bool null_det = false;
  if (length_squared(cross(i, k)) <= 1e-6 ||
      length_squared(cross(i, j)) <= 1e-6 ||
      length_squared(cross(k, j)) <= 1e-6 || l01 <= 1e-4 || l02 <= 1e-4)
    null_det = true;
  if (det > 1e-10) {
    /*if (yocto::abs(j.z - k.y) >= 1e-2 || yocto::abs(k.x - i.z) >= 1e-2 ||
        yocto::abs(i.y - j.x) >= 1e-2) {*/

    gamma = ((g0_norm - g02) * l01 + A * (g01 - g0_norm)) / det;
    beta  = (g0_norm - g01 - gamma * A) / l01;
    alpha = 1 - beta - gamma;
    if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 &&
        gamma <= 1) {
      return {true, -1};

    } else
      return {false, -1};
  } else {
    /* Eigen::Matrix2f Q;
    Q << l01, A, A, l02;
    Eigen::Matrix2f I;
    I.setIdentity();

    Eigen::LLT<Eigen::Matrix2f> lltOfA(Q);
    auto                        factor = 0.00001;

    while (lltOfA.info() == Eigen::NumericalIssue) {
      Q      = Q + factor * I;
      lltOfA = Q.llt();
      factor = factor + 0.00001;
    }
    Eigen::Vector2f b;
    b << g0_norm - g01, g0_norm - g02;
    auto x = lltOfA.solve(b);
    beta   = x(0);
    gamma  = x(1);
    alpha  = 1 - beta - gamma;
    if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 &&
        gamma <= 1)
      return {true, -1};
    else*/
    return {false, -1};

    /* gamma = 0.f;
    beta  = (g0_norm - g01) / l01;
    alpha = 1 - beta;*/

  } /* else
    return {false, -1}; */
}
std::tuple<bool, vec3f, float> minimum_on_line(const vec3f& bary_in,
    const vec3f& bary_out, const vec3f& g0, const vec3f& g1, const vec3f& g2) {
  auto g0_norm = length_squared(g0);
  auto g1_norm = length_squared(g1);
  auto g2_norm = length_squared(g2);
  auto g01     = dot(g0, g1);
  auto g02     = dot(g0, g2);
  auto g12     = dot(g1, g2);
  auto l01     = length_squared(g0 - g1);
  auto l02     = length_squared(g0 - g2);

  auto A = g0_norm + g12 - g02 - g01;

  auto delta = (-bary_out.y * bary_in.x + bary_out.x * bary_in.y -
                   bary_out.y * bary_in.z + bary_out.z * bary_in.y) /
               (-bary_in.x * bary_out.z + bary_out.x * bary_in.z -
                   bary_in.y * bary_out.z + bary_out.y * bary_in.z);
  auto xi = (bary_out.y * bary_in.z - bary_out.z * bary_in.y) /
            (-bary_in.x * bary_out.z + bary_out.x * bary_in.z -
                bary_in.y * bary_out.z + bary_out.y * bary_in.z);
  if (pow(delta, 2) * l01 + l02 + 2 * delta * A >= 1e-8) {
    float gamma =
        (xi * (-delta * l01 - A) - delta * (g01 - g0_norm) - g02 + g0_norm) /
        (pow(delta, 2) * l01 + l02 + 2 * delta * A);
    float beta  = gamma * delta + xi;
    float alpha = 1 - gamma - beta;

    if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 &&
        gamma <= 1) {
      auto g = alpha * g0 + beta * g1 + gamma * g2;
      return {true, vec3f{alpha, beta, gamma}, length_squared(g)};

    } else
      return {false, zero3f, -1};
  } else
    return {false, zero3f, -1};
}
bool tri_contains_min(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const int tid, const vec3f& g0,
    const vec3f& g1, const vec3f& g2, const vector<vec3f>& normals,
    vec3f& bary) {
  int  p0 = triangles[tid].x;
  int  p1 = triangles[tid].y;
  int  p2 = triangles[tid].z;
  auto gx = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g0, normals, p0, tid, V2T);
  auto gy = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g1, normals, p1, tid, V2T);
  auto gz = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g2, normals, p2, tid, V2T);

  return bary_for_minimum_energy(gx, gy, gz, bary.x, bary.y, bary.z);
}
pair<bool, int> tri_contains_min_hessian(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const int tid, const vec3f& g0,
    const vec3f& g1, const vec3f& g2, const vector<vec3f>& normals, vec3f& bary,
    bool same_tangent_space = false) {
  //  auto G         = vector<vec3f>(3);
  //  auto grad      = vector<vec3f>(3);
  //  auto flat_grad = vector<vec2f>(3);
  //  grad           = {g0, g1, g2};
  //  auto flat_tid  = init_flat_triangle(positions, triangles[tid]);
  //  auto c         = (flat_tid[0] + flat_tid[1] + flat_tid[2]) / 3.f;
  //  for (int i = 0; i < 3; ++i) {
  //    G[i] = vector_bary_coords(triangles, positions, tid, grad[i]);
  //  }
  //  for (auto i = 0; i < 3; ++i) {
  //    flat_grad[i] = G[i].x * (flat_tid[0] - c) + G[i].y * (flat_tid[1] - c)
  //    +
  //                   G[i].z * (flat_tid[2] - c);
  //  }
  //
  //  if (yocto::abs(flat_grad[1][1] - flat_grad[0][1] - flat_grad[2][0] +
  //                 flat_grad[0][0]) < 1e-2)
  //    return {false, infty};

  if (!same_tangent_space) {
    int  p0 = triangles[tid].x;
    int  p1 = triangles[tid].y;
    int  p2 = triangles[tid].z;
    auto gx = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, g0, normals, p0, tid, V2T);
    auto gy = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, g1, normals, p1, tid, V2T);
    auto gz = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, g2, normals, p2, tid, V2T);

    return bary_for_minimum_hessian(gx, gy, gz, bary.x, bary.y, bary.z);
  } else
    return bary_for_minimum_hessian(g0, g1, g2, bary.x, bary.y, bary.z);
}
bool tri_contains_min_energy(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const int tid, const vec3f& g0,
    const vec3f& g1, const vec3f& g2, const vector<vec3f>& normals,
    vec3f& bary) {
  int  p0 = triangles[tid].x;
  int  p1 = triangles[tid].y;
  int  p2 = triangles[tid].z;
  auto gx = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g0, normals, p0, tid, V2T);
  auto gy = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g1, normals, p1, tid, V2T);
  auto gz = transp_vec(solver, angles, total_angles, triangles, positions,
      adjacencies, g2, normals, p2, tid, V2T);

  auto flat = init_flat_triangle(positions, triangles[tid]);
  auto c2D = interpolate_triangle(flat[0], flat[1], flat[2], vec2f{0.33, 0.33});
  auto c3D = interpolate_triangle(
      positions[p0], positions[p1], positions[p2], vec2f{0.33, 0.33});
  auto e3D    = positions[p0] - c3D;
  auto e2D    = flat[0] - c2D;
  auto n      = tid_normal(triangles, positions, tid);
  auto g      = vector<vec3f>{gx, gy, gz};
  auto flat_g = vector<vec2f>(3);
  for (auto i = 0; i < 3; ++i) {
    auto teta = angle(e3D, g[i]);
    if (dot(cross(e3D, g[i]), n) < 0) teta = 2 * M_PI - teta;
    auto v    = rot_vect(e2D, teta);
    flat_g[i] = v;
  }

  return bary_for_minimum_energy(
      flat_g[0], flat_g[1], flat_g[2], bary.x, bary.y, bary.z);
}
bool tri_contains_min(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const int tid, const vector<vec3f>& g0,
    const vector<vec3f>& g1, const vector<vec3f>& g2, const vector<float>& w,
    const vector<vec3f>& normals, vec3f& bary) {
  vec3f v0 = zero3f, v1 = zero3f, v2 = zero3f;

  for (int i = 0; i < g0.size(); ++i) {
    vec3f tmp0 = 2 * g0[i] * w[i];
    vec3f tmp1 = 2 * g1[i] * w[i];
    vec3f tmp2 = 2 * g2[i] * w[i];
    v0 += tmp0;
    v1 += tmp1;
    v2 += tmp2;
  }
  int p0 = triangles[tid].x;
  int p1 = triangles[tid].y;
  int p2 = triangles[tid].z;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, v0, normals, p0, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, v1, normals, p1, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, v2, normals, p2, tid, V2T);

  return bary_for_minimum(v0, v1, v2, bary.x, bary.y, bary.z);
}
vec2f from_3d_to_2d_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec3f& v, const int tid,
    const unfold_triangle& flat_tid) {
  auto baryv = vector_bary_coords(triangles, positions, tid, v);

  return baryv.y * (flat_tid[1] - flat_tid[0]) +
         baryv.z * (flat_tid[2] - flat_tid[0]);
}
vec2f from_3d_to_2d_vector(const vec3f& p0, const vec3f& p1, const vec3f& p2,
    const vec3f& v, const unfold_triangle& flat_tid) {
  auto baryv = vector_bary_coords(p0, p1, p2, v);

  return baryv.y * (flat_tid[1] - flat_tid[0]) +
         baryv.z * (flat_tid[2] - flat_tid[0]);
}
vec3f from_2d_to_3d_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec2f& v, const int tid,
    const unfold_triangle& flat_tid) {
  auto c    = (flat_tid[0] + flat_tid[1] + flat_tid[2]) / 3.f;
  auto c3d  = tid_centroid(triangles, positions, tid);
  auto e    = flat_tid[0] - c;
  auto teta = angle(e, v);
  if (cross(e, v) < 0) teta *= -1;
  auto e3d = normalize(triangles[tid].x - c3d);
  auto len = length(v);
  auto n   = tid_normal(triangles, positions, tid);
  return len * rot_vect(e3d, n, teta);
}
vec3f from_2d_to_3d_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec2f& v, const int tid) {
  auto bary = vector_bary_coords(triangles, positions, tid, v);
  return bary.y * (positions[triangles[tid].y] - positions[triangles[tid].x]) +
         bary.z * (positions[triangles[tid].z] - positions[triangles[tid].x]);
}

vector<vec2f> get_orhtonormal_basis(const unfold_triangle& tid) {
  auto e0 = normalize(tid[1] - tid[0]);
  auto e1 = rot_vect(e0, pif / 2);
  return {e0, e1};
}
std::pair<vec3f, bool> maximum_descent_direction(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& grad) {
  auto flat_tid = init_flat_triangle(positions, triangles[tid]);

  auto c = interpolate_triangle(
      flat_tid[0], flat_tid[1], flat_tid[2], vec2f{0.33, 0.33});
  auto basis = get_orhtonormal_basis(flat_tid);

  Eigen::Matrix2f A;
  A << dot(flat_tid[1] - flat_tid[0], basis[0]),
      dot(flat_tid[2] - flat_tid[0], basis[0]),
      dot(flat_tid[1] - flat_tid[0], basis[1]),
      dot(flat_tid[2] - flat_tid[0], basis[1]);
  if (A.determinant() < 1e-8) return std::make_pair(zero3f, false);
  auto A_t       = A.transpose();
  auto Q         = A_t * A;
  auto flat_grad = from_3d_to_2d_vector(
      triangles, positions, grad, tid, flat_tid);
  Eigen::Vector2f b;
  b << dot(flat_grad, basis[0]), dot(flat_grad, basis[1]);
  auto b_t = A_t * b;
  auto x   = Q.colPivHouseholderQr().solve(b_t);
  auto v   = (1 - x(0) - x(1)) * (flat_tid[0] - c) + x(0) * (flat_tid[1] - c) +
           x(1) * (flat_tid[2] - c);
  return std::make_pair(
      from_2d_to_3d_vector(triangles, positions, v, tid, flat_tid), true);
}
std::pair<vec3f, bool> maximum_descent_direction(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& grad,
    const vec3f& gx, const vec3f& gy, const vec3f& gz) {
  auto flat_tid = init_flat_triangle(positions, triangles[tid]);

  auto basis = get_orhtonormal_basis(flat_tid);
  auto baryj = vector_bary_coords(triangles, positions, tid, gy - gx);
  auto baryk = vector_bary_coords(triangles, positions, tid, gz - gx);
  auto x0    = dot(flat_tid[2] - flat_tid[0], basis[0]);
  auto x1    = dot(flat_tid[2] - flat_tid[0], basis[1]);
  Eigen::Matrix2f A;
  Eigen::Matrix2f I;
  I.setIdentity();
  auto len = length(flat_tid[1] - flat_tid[0]);

  A(0, 0) = baryj.y + baryj.z * x0 / len;
  A(0, 1) = (baryk.y * len + baryk.z * x0 - baryj.y * x0) / x1 -
            pow(x0, 2) * baryj.z / (len * x1);
  A(1, 0) = baryj.z * x1 / len;
  A(1, 1) = baryk.z - x0 * baryj.z / len;

  Eigen::LLT<Eigen::Matrix2f> lltOfA(A);

  if (lltOfA.info() == Eigen::NumericalIssue) {
    return {zero3f, false};
  }

  auto flat_grad = from_3d_to_2d_vector(
      triangles, positions, grad, tid, flat_tid);
  Eigen::Vector2f b;
  b << -flat_grad.x, -flat_grad.y;
  auto x     = lltOfA.solve(b);
  auto baryx = vector_bary_coords(flat_tid, vec2f{x(0), x(1)});

  auto xi =
      baryx.y * (positions[triangles[tid].y] - positions[triangles[tid].x]) +
      baryx.z * (positions[triangles[tid].z] - positions[triangles[tid].x]);

  return std::make_pair(xi, true);
}
unfold_triangle init_flat_triangle(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
  auto tr2d = unfold_triangle{};
  tr2d[0]   = {0, 0};
  tr2d[1]   = {0, length(p0 - p1)};
  auto rx   = length_squared(p0 - p2);
  auto ry   = length_squared(p1 - p2);
  tr2d[2]   = intersect_circles(tr2d[0], rx, tr2d[1], ry);
  return tr2d;
}
std::pair<vec3f, bool> maximum_descent_direction(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& grad, const vec3f& gx, const vec3f& gy,
    const vec3f& gz) {
  auto flat_tid = init_flat_triangle(p0, p1, p2);

  auto            basis = get_orhtonormal_basis(flat_tid);
  auto            baryj = vector_bary_coords(p0, p1, p2, gy - gx);
  auto            baryk = vector_bary_coords(p0, p1, p2, gz - gx);
  auto            x0    = dot(flat_tid[2] - flat_tid[0], basis[0]);
  auto            x1    = dot(flat_tid[2] - flat_tid[0], basis[1]);
  Eigen::Matrix2f A;
  Eigen::Matrix2f I;
  I.setIdentity();
  auto len = length(flat_tid[1] - flat_tid[0]);

  A(0, 0) = baryj.y + baryj.z * x0 / len;
  A(0, 1) = (baryk.y * len + baryk.z * x0 - baryj.y * x0) / x1 -
            pow(x0, 2) * baryj.z / (len * x1);
  A(1, 0) = baryj.z * x1 / len;
  A(1, 1) = baryk.z - x0 * baryj.z / len;

  Eigen::LLT<Eigen::Matrix2f> lltOfA(A);
  if (lltOfA.info() == Eigen::NumericalIssue) {
    return {zero3f, false};
  }

  auto            flat_grad = from_3d_to_2d_vector(p0, p1, p2, grad, flat_tid);
  Eigen::Vector2f b;
  b << -flat_grad.x, -flat_grad.y;
  auto x     = lltOfA.solve(b);
  auto baryx = vector_bary_coords(flat_tid, vec2f{x(0), x(1)});

  auto xi = baryx.y * (p1 - p0) + baryx.z * (p2 - p0);

  return std::make_pair(xi, true);
}
float weighted_average(
    const vector<vector<float>>& f, const vector<float>& w, const int vid) {
  double lambda = 0.0;
  int    n      = (int)f.size();
  for (int i = 0; i < n; ++i) {
    double value = f[i][vid] * w[i];  // pow(f[i][vid], 2) * w[i];
    lambda += value;
  }
  return (float)lambda;  // is it wrong to cast at the end?
}
inline float lagrange_polynomial_vert(const float& lambda) {
  return 2 * pow(lambda, 2) - lambda;
}
inline float lagrange_polynomial_midpoint(
    const float& lambda0, const float& lambda1) {
  return 4 * lambda0 * lambda1;
}
float field_blending(const vector<vector<vec3f>>& E,
    const vector<float>& weights, const int tid, const int eid) {
  auto lambda = 0.f;

  for (int i = 0; i < weights.size(); ++i) {
    lambda += E[i][tid][eid] * weights[i];
  }

  return lambda;
}
float field_blending(const vector<vector<float>>& f,
    const vector<float>& weights, const int vid) {
  auto lambda = 0.f;

  for (int i = 0; i < weights.size(); ++i) {
    lambda += f[i][vid] * weights[i];
  }

  return lambda;
}
std::tuple<bool, vector<pair<float, int>>, vec2f> tri_contains_min(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const geodesic_solver& solver, const vector<float>& f,
    const vector<vector<vec3f>>& E, const vector<float>& weights,
    const int tid) {
  std::tuple<bool, vector<pair<float, int>>, vec2f> out;
  auto                                              vid0 = triangles[tid].x;
  auto                                              vid1 = triangles[tid].y;
  auto                                              vid2 = triangles[tid].z;
  auto                                              f0   = f[vid0];
  auto                                              f1   = f[vid1];
  auto                                              f2   = f[vid2];

  auto f01 = field_blending(E, weights, tid, 0);
  // average_field(
  //     triangles, positions, adjacencies, v2t, solver, f, vid0, vid1);

  auto f02 = field_blending(E, weights, tid, 2);
  // average_field(
  //     triangles, positions, adjacencies, v2t, solver, f, vid0, vid2);
  // k        = find(triangles[tid], vid1);
  // auto f12 = f1 * lagrange_polynomial_vert(mid_points[1][k]) +
  //            f2 * lagrange_polynomial_vert(mid_points[2][(k + 1) % 3]);
  auto f12 = field_blending(E, weights, tid, 1);
  // average_field(
  //     triangles, positions, adjacencies, v2t, solver, f, vid1, vid2);
  Eigen::Matrix2f A;
  A << f0 + f2 - 2 * f02, f2 + f01 - f12 - f02, f2 - f12 + f01 - f02,
      f1 + f2 - 2 * f12;
  Eigen::Vector2f B;
  B << (f0 + 3 * f2) / 4 - f02, (f1 + 3 * f2) / 4 - f12;
  auto x      = A.colPivHouseholderQr().solve(B);
  auto bary   = zero2f;
  bool inside = true;
  if (x(0) == 0 && x(1) == 0) inside = false;

  for (auto i = 0; i < 2; ++i) {
    if (x(i) < -0.05 || x(i) > 1 + 0.05) inside = false;
    bary[i] = x(i);
  }
  bary        = vec2f{bary.y, 1 - bary.x - bary.y};
  auto interp = vector<pair<float, int>>{{f01, vid2}, {f02, vid1}, {f12, vid0}};
  out         = {inside, interp, bary};
  return out;
}
// return the lerp of the minum on the edge
pair<bool, float> edge_contains_min(const vector<vec3i>& triangles,
    const vector<float>& f, const vector<vector<vec3f>>& E,
    const vector<float>& weights, const int tid, const int k) {
  auto            is_sol = true;
  Eigen::Vector2f x;
  switch (k) {
    case (0): {
      auto            f01 = field_blending(E, weights, tid, 0);
      auto            f0  = f[triangles[tid].x];
      auto            f1  = f[triangles[tid].y];
      Eigen::Matrix2f A;
      A << f0, f01, f01, f1;
      Eigen::Vector2f B;
      B << f0 / 4, f1 / 4;
      x = A.colPivHouseholderQr().solve(B);

      if (x(0) == 0 && x(1) == 0) is_sol = false;
      if (x(0) > -1e-2 || x(0) < 1 + 1e-2) is_sol = false;

    } break;
    case (1): {
      auto            f12 = field_blending(E, weights, tid, 1);
      auto            f1  = f[triangles[tid].y];
      auto            f2  = f[triangles[tid].z];
      Eigen::Matrix2f A;
      A << f2, f2 - f12, f2 - f12, f2 + f1 - 2 * f12;
      Eigen::Vector2f B;
      B << 3 * f2 / 4, (f1 + 3 * f2) / 4 - f12;
      x = A.colPivHouseholderQr().solve(B);
    } break;

    case (2): {
      auto            f20 = field_blending(E, weights, tid, 1);
      auto            f0  = f[triangles[tid].y];
      auto            f2  = f[triangles[tid].z];
      Eigen::Matrix2f A;
      A << f0 - 2 * f20 + f2, f2 - f20, f2 - f20, f2;
      Eigen::Vector2f B;
      B << (f0 + 3 * f2) / 4 - f20, 3 * f2 / 4;
      x = A.colPivHouseholderQr().solve(B);
    } break;
  }
  return std::make_pair(is_sol, x(0));
}
vec3f gradient_blending(const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const int vid) {
  auto v = zero3f;

  for (int i = 0; i < weights.size(); ++i) {
    auto w = gradients[i][vid] * weights[i];
    v += w;
  }

  return v;
}
vec3f almost_gradient(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<float>& weights,
    const mesh_point& p, const vector<mesh_point>& control_points) {
  auto al_grad = zero3f;
  auto tid     = p.face;
  auto pos     = eval_position(triangles, positions, p);
  for (auto i = 0; i < control_points.size(); ++i) {
    auto strip = strip_on_dual_graph(
        dual_solver, triangles, positions, control_points[i].face, tid);
    auto path = shortest_path(
        triangles, positions, adjacencies, p, control_points[i], strip);
    auto w = normalize(
        path_pos_from_entry(triangles, positions, adjacencies, path, 0) - pos);

    auto [is_vert, k] = bary_is_vert(get_bary(p.uv));
    if (is_vert) {
      k         = find(triangles[path.strip[0]], triangles[tid][k]);
      auto bary = zero3f;
      bary[k]   = 1.f;
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, w, normals, path.strip[0], triangles[path.strip[0]][k],
          T2V);
    }

    w *= -1;          // opposite direction
    w *= weights[i];  // averaging
    al_grad += w;
  }

  return al_grad;
}
vector<vec3f> almost_gradient_blending(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<float>& weights,
    const int tid, const vector<mesh_point>& control_points) {
  auto al_grad = vector<vec3f>(3);
  for (auto i = 0; i < 3; ++i) {
    auto bary  = zero3f;
    bary[i]    = 1.f;
    al_grad[i] = almost_gradient(dual_solver, solver, triangles, positions,
        adjacencies, normals, angles, total_angles, weights,
        mesh_point{tid, {bary.y, bary.z}}, control_points);
  }
  return al_grad;
}
vec3f gradient_blending(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3f>& positions, const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& normals, const int& vid,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights) {
  auto v = zero3f;
  if (gradients[0].size() == positions.size()) {
    for (int i = 0; i < weights.size(); ++i) {
      auto w = gradients[i][vid] * weights[i];
      v += w;
    }
  } else {
    for (int i = 0; i < weights.size(); ++i) {
      auto w = AGS_gradient(solver, angles, total_angles, positions, triangles,
          adjacencies, v2t, gradients[i], vid, normals);
      v += w * weights[i];
    }
  }

  return v;
}
vec3f karcher_grad(const vector<vector<vec3f>>& gradients,
    const vector<vector<float>>& f, const vector<float>& weights,
    const int vid) {
  auto v = zero3f;
  for (auto i = 0; i < f.size(); ++i) {
    auto w = 2 * f[i][vid] * normalize(gradients[i][vid]);
    w *= weights[i];
    v += w;
  }
  return v;
}
vector<vec3f> karcher_grad(const vector<vector<vec3f>>& gradients,
    const vector<vector<float>>& f, const vector<float>& weights,
    const vector<vec3f>& normals, bool to_draw) {
  auto grad = vector<vec3f>(f[0].size());
  for (auto i = 0; i < f[0].size(); ++i) {
    grad[i] = karcher_grad(gradients, f, weights, i);
    if (to_draw) grad[i] += 0.001 * normals[i];
  }
  return grad;
}

vector<float> field_blending(const geodesic_solver& solver,
    const vector<vector<float>>& f, const vector<float>& weights, const int vid,
    const int k) {
  auto blended = vector<float>(f[0].size(), flt_max);
  auto nbr     = k_ring(solver, vid, k, true);
  for (auto j = 0; j <= nbr.size(); ++j) {
    auto lambda = 0.f;
    auto curr   = (j < nbr.size()) ? nbr[j] : vid;
    for (int i = 0; i < weights.size(); ++i) {
      lambda += f[i][curr] * weights[i];
    }
    blended[curr] = lambda;
  }

  return blended;
}

// Karcher mean in the neighborhood of a point
int minimum_on_verts(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vector<float>>& f,
    const vector<float>& w, const int source) {
  vector<pair<double, int>> means;
  vector<int>               visited = {};
  double                    lambda;
  auto                      nbr = solver.graph[source];

  int candidate;
  for (int i = 0; i <= nbr.size(); ++i) {
    int curr = (i == nbr.size()) ? source : nbr[i].node;

    lambda = weighted_average(f, w, curr);
    means.push_back(std::make_pair(lambda, i));
  }

  sort(means.begin(), means.end());

  int prev_entry = means[0].second;

  if (prev_entry == nbr.size()) {
    return source;
  } else {
    visited.push_back(source);
    bool stop = false;
    candidate = nbr[prev_entry].node;
    visited.push_back(candidate);
    nbr = solver.graph[candidate];

    while (!stop) {
      means.clear();
      for (int i = 0; i <= nbr.size(); ++i) {
        int curr = (i == nbr.size()) ? candidate : nbr[i].node;
        lambda   = weighted_average(f, w, curr);
        means.push_back(std::make_pair(lambda, i));
      }

      sort(means.begin(), means.end());

      if (means[0].second == nbr.size())
        stop = true;
      else {
        candidate = nbr[means[0].second].node;
        if (find(visited, candidate) != -1)
          return source;

        else {
          nbr = solver.graph[candidate];
          visited.push_back(candidate);
        }
      }
    }
    return candidate;
  }

  return -1;
}
// global minimum
int minimum_on_verts(const vector<vector<float>>& f, const vector<float>& w) {
  auto inf    = 0;
  auto lambda = flt_max;
  for (int i = 0; i < f[0].size(); ++i) {
    auto curr_lambda = field_blending(f, w, i);
    if (lambda > curr_lambda) {
      lambda = curr_lambda;
      inf    = i;
    }
  }
  return inf;
}
mesh_point minimum_on_verts(const vector<vec3i>& triangles,
    const vector<vector<int>>& v2t, const vector<vector<float>>& f,
    const vector<float>& w) {
  auto inf    = 0;
  auto lambda = flt_max;
  for (int i = 0; i < f[0].size(); ++i) {
    auto curr_lambda = field_blending(f, w, i);
    if (lambda > curr_lambda) {
      lambda = curr_lambda;
      inf    = i;
    }
  }
  auto tid  = v2t[inf][0];
  auto k    = find(triangles[inf], inf);
  auto bary = zero3f;
  bary[k]   = 1;
  return {tid, {bary.y, bary.z}};
}
int minimum_on_gradient_norm(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vector<float>>& f,
    const vector<vector<vec3f>>& gradients, const vector<float>& w,
    const int source) {
  vector<pair<float, int>> means;
  vector<int>              visited = {};
  auto                     nbr     = solver.graph[source];

  int candidate;
  for (int i = 0; i <= nbr.size(); ++i) {
    int  curr = (i == nbr.size()) ? source : nbr[i].node;
    auto grad = gradient_blending(gradients, w, curr);
    // auto grad = karcher_grad(gradients, f, w, curr);
    means.push_back(std::make_pair(length(grad), i));
  }

  sort(means.begin(), means.end());

  int prev_entry = means[0].second;

  if (prev_entry == nbr.size()) {
    return source;
  } else {
    visited.push_back(source);
    bool stop = false;
    candidate = nbr[prev_entry].node;
    visited.push_back(candidate);
    nbr = solver.graph[candidate];

    while (!stop) {
      means.clear();
      for (int i = 0; i <= nbr.size(); ++i) {
        int  curr = (i == nbr.size()) ? candidate : nbr[i].node;
        auto grad = gradient_blending(gradients, w, curr);
        // auto grad = karcher_grad(gradients, f, w, curr);
        means.push_back(std::make_pair(length(grad), i));
      }

      sort(means.begin(), means.end());

      if (means[0].second == nbr.size())
        stop = true;
      else {
        candidate = nbr[means[0].second].node;
        if (find(visited, candidate) != -1)
          return source;

        else {
          nbr = solver.graph[candidate];
          visited.push_back(candidate);
        }
      }
    }
    return candidate;
  }

  return -1;
}
// remember to remove tids;
int minimum_on_tri(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>& adjacencies, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vector<vec3f>>& gradients, const vector<float>& w,
    const int source, vector<int>& tids) {
  auto nbr = vector<int>{
      adjacencies[source].x, adjacencies[source].y, adjacencies[source].z};
  auto boundary = nbr;
  auto checked  = vector<bool>(triangles.size(), false);
  auto g0 = zero3f, g1 = zero3f, g2 = zero3f;
  auto arrived = false;
  for (int i = 0; i < 3; ++i) {
    int curr      = nbr[i];
    checked[curr] = true;
  }
  auto count = 0;
  while (!arrived) {
    auto aux = vector<int>{};
    for (auto i = 0; i < boundary.size(); ++i) {
      auto adj = adjacencies[boundary[i]];
      for (auto j = 0; j < 3; ++j) {
        if (!checked[adj[j]]) {
          nbr.push_back(adj[j]);
          aux.push_back(adj[j]);
        }
      }
    }
    boundary = aux;
    for (auto i = 0; i < nbr.size(); ++i) {
      if (checked[nbr[i]]) continue;
      checked[nbr[i]] = true;
      g0              = gradient_blending(gradients, w, triangles[nbr[i]].x);
      g1              = gradient_blending(gradients, w, triangles[nbr[i]].y);
      g2              = gradient_blending(gradients, w, triangles[nbr[i]].z);
      auto bary       = zero3f;
      if (tri_contains_min(solver, angles, total_angles, triangles, positions,
              adjacencies, nbr[i], g0, g1, g2, normals, bary))
        return nbr[i];
    }
    ++count;
  }
}
vec2i dir_from_bary(const vec3f& bary) {
  auto id0 = -1, id1 = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] < 0) {
      if (id0 == -1)
        id0 = i;
      else {
        id1 = (bary[(i + 1) % 3] < 0) ? (i + 2) % 3 : (i + 1) % 3;
      }
    }
  }
  return {id0, id1};  // ??
}
Eigen::VectorXf solve_quadratic_field(
    const vec3f& G0, const vec3f& G1, const vec3f& G2) {
  Eigen::MatrixXf A(6, 5);
  A << 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
      0, 0, 2, 0, 0, 1;

  auto            At = A.transpose();
  Eigen::VectorXf B(6);
  B << G0.y, G0.z, G1.y, G1.z, G2.y, G2.z;
  auto                                        rhs = At * B;
  auto                                        Q   = At * A;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(Q);
  auto                                        X = dec.solve(rhs);
  return X;
}
Eigen::VectorXf solve_quadratic_field(const vec3f& G0, const vec3f& G1,
    const vec3f& G2, const float& F1, const float& F2, const float& F3) {
  Eigen::MatrixXf A(8, 6);
  A << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0,
      0, 2, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1;

  auto            At = A.transpose();
  Eigen::VectorXf B(8);
  B << G0.y, G0.z, G1.y, G2.y, G2.z, F1, F2, F3;
  auto                                        rhs = At * B;
  auto                                        Q   = At * A;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(Q);
  auto                                        X = dec.solve(rhs);
  return X;
}

vec3f irrot_field_singularity(const Eigen::VectorXf& coeff) {
  Eigen::Matrix2f grad;
  grad << 2 * coeff(0), coeff(2), coeff(2), 2 * coeff(1);
  Eigen::Vector2f b;
  b << -coeff(3), -coeff(4);
  Eigen::ColPivHouseholderQR<Eigen::Matrix2f> dec(grad);
  auto                                        X = dec.solve(b);

  auto bary = vec2f{X(0), X(1)};
  return vec3f{1 - bary.x - bary.y, bary.x, bary.y};
}
vector<vec3f> non_irrot_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& V0,
    const vec3f& V1, const vec3f& V2, vector<vec3f>& bary) {
  vector<vec3f> G(3);
  bary.resize(3, zero3f);
  vector<vec3f> grad = {V0, V1, V2};
  auto          p0   = positions[triangles[tid].x];
  auto          p1   = positions[triangles[tid].y];
  auto          p2   = positions[triangles[tid].z];
  auto          c    = tid_centroid(triangles, positions, tid);
  for (int i = 0; i < 3; ++i) {
    G[i] = vector_bary_coords(triangles, positions, tid, grad[i]);
  }

  Eigen::VectorXf coeff = solve_quadratic_field(G[0], G[1], G[2]);

  float y = 2 * coeff(0) + coeff(3);
  float z = coeff(2) + coeff(4);
  float x = -y - z;
  bary[1] = vec3f{x, y, z};
  auto w  = x * (p0 - c) + y * (p1 - c) + z * (p2 - c);

  grad[1] = w;

  y       = coeff(2) + coeff(3);
  z       = 2 * coeff(1) + coeff(4);
  x       = -y - z;
  bary[2] = vec3f{x, y, z};
  w       = x * (p0 - c) + y * (p1 - c) + z * (p2 - c);

  grad[2] = w;

  y       = coeff(3);
  z       = coeff(4);
  x       = -y - z;
  bary[0] = vec3f{x, y, z};
  w       = x * (p0 - c) + y * (p1 - c) + z * (p2 - c);

  grad[0] = w;

  return grad;
}
vec3f handle_non_irrot_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& V0,
    const vec3f& V1, const vec3f& V2) {
  vector<vec3f> G(3);
  vector<vec3f> curr_grad = {V0, V1, V2};
  for (int i = 0; i < 3; ++i) {
    G[i] = vector_bary_coords(triangles, positions, tid, curr_grad[i]);
  }

  Eigen::VectorXf coeff = solve_quadratic_field(G[0], G[1], G[2]);
  auto            min   = irrot_field_singularity(coeff);

  return min;
}
int step_for_gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3f>& gradients,
    const vec3f& bary, const int curr_tid) {
  auto id  = dir_from_bary(bary);
  auto tid = -1;
  if (id[0] == vert) {
    auto vid = triangles[curr_tid][id[1]];
    auto v   = normalize(gradients[vid]);
    tid      = next_tid_extended_graph(
             solver, angles, positions, v2t, triangles, normals, vid, v);
  } else {
    tid = opposite_face(
        triangles, adjacencies, curr_tid, triangles[curr_tid][id[1]]);
  }

  return tid;
}
int step_for_gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const vec3f& bary, const int curr_tid) {
  auto id  = dir_from_bary(bary);
  auto tid = -1;
  if (id[0] == vert) {
    auto vid0 = triangles[curr_tid][id[1]];
    auto k    = find(triangles[curr_tid], vid0);
    auto vid1 = triangles[curr_tid][(k + 1) % 3];
    auto vid2 = triangles[curr_tid][(k + 2) % 3];
    auto g0   = gradient_blending(gradients, weights, vid0);
    auto g1   = gradient_blending(gradients, weights, vid1);
    auto g2   = gradient_blending(gradients, weights, vid2);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, g0, normals, vid0, curr_tid, V2T);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, g1, normals, vid1, curr_tid, V2T);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, g2, normals, vid2, curr_tid, V2T);
    auto t01 = length_squared(g1) / (length_squared(g0) + length_squared(g1));
    auto t02 = length_squared(g2) / (length_squared(g0) + length_squared(g2));

    auto g01 = t01 * g0 + (1 - t01) * g1;
    auto g02 = t02 * g0 + (1 - t02) * g2;
    if (length(g01) < length(g02))
      tid = adjacencies[curr_tid][(k + 2) % 3];
    else
      tid = adjacencies[curr_tid][k];
    // tid = next_tid_extended_graph(
    //     solver, angles, positions, v2t, triangles, normals, vid,
    //     normalize(v));

  } else {
    tid = opposite_face(
        triangles, adjacencies, curr_tid, triangles[curr_tid][id[1]]);
  }

  return tid;
}
pair<int, vec3f> step_for_gradient_descent_minimization(
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vec3f& bary, const vec3f& dir, const int curr_tid) {
  auto sample_bary = zero3f, sample_pos = zero3f;
  trace_in_triangles(
      positions, triangles, dir, bary, curr_tid, sample_pos, sample_bary);
  auto [is_edge, ke] = bary_is_edge(sample_bary);
  auto [is_vert, kv] = bary_is_vert(sample_bary);
  if (is_edge) {
    sample_bary = tri_bary_coords(
        triangles, positions, adjacencies[curr_tid][ke], sample_pos);
    return {adjacencies[curr_tid][ke], sample_bary};
  }

  else if (is_vert) {
    /*  auto vid = triangles[curr_tid][kv];
     auto v   = transp_vec(solver, angles, total_angles, triangles, positions,
         adjacencies, dir, normals, curr_tid, vid, T2V);
     auto next_tri = next_tid(
         solver, angles, positions, v2t, triangles, normals, vid, v);
     kv              = find(triangles[next_tri], vid);
     sample_bary     = zero3f;
     sample_bary[kv] = 1;
     return {next_tri, sample_bary}; */
    auto vid = triangles[curr_tid][kv];
    auto e0  = normalize(positions[triangles[curr_tid][(kv + 1) % 3]] - vid);
    auto e1  = normalize(positions[triangles[curr_tid][(kv + 2) % 3]] - vid);
    auto c0  = dot(e0, dir);
    auto c1  = dot(e1, dir);
    if (c0 < c1) {
      auto next_tri   = adjacencies[curr_tid][kv];
      kv              = find(triangles[next_tri], vid);
      sample_bary     = zero3f;
      sample_bary[kv] = 1;
      return {next_tri, sample_bary};

    } else {
      auto next_tri   = adjacencies[curr_tid][(kv + 2) % 3];
      kv              = find(triangles[next_tri], vid);
      sample_bary     = zero3f;
      sample_bary[kv] = 1;
      return {next_tri, sample_bary};
    }
  }
}
std::pair<mesh_point, bool> step_for_gradient_descent_minimization(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const mesh_point& from, const vec3f& gx, const vec3f& gy, const vec3f& gz,
    const vec3f& grad) {
  auto [dir, is_descent] = maximum_descent_direction(
      triangles, positions, from.face, grad, gx, gy, gz);
  if (!is_descent) return {mesh_point{}, false};
  auto step = straightest_geodesic(solver, triangles, positions, normals,
      adjacencies, v2t, angles, total_angles, from, normalize(dir),
      length(dir));
  return {step.back(), true};
}
std::pair<int, bool> step_for_gradient_descent_minimization(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& from, const vec3f& gx,
    const vec3f& gy, const vec3f& gz, const vec3f& grad, vec3f& bary) {
  auto sample_pos        = zero3f;
  auto [dir, is_descent] = maximum_descent_direction(
      triangles, positions, from.face, grad, gx, gy, gz);
  if (!is_descent) return {-1, false};
  trace_in_triangles(positions, triangles, dir,
      vec3f{1 / 3.f, 1 / 3.f, 1 / 3.f}, from.face, sample_pos, bary);
  auto [is_edge, ke] = bary_is_edge(bary);
  auto [is_vert, kv] = bary_is_vert(bary);
  if (is_edge) {
    auto next = adjacencies[from.face][ke];
    auto eid  = common_edge(triangles[from.face], triangles[next]);

    int  offset0           = find(triangles[next], eid.x);
    auto tmp               = zero3f;
    tmp[offset0]           = bary[ke];
    tmp[(offset0 + 2) % 3] = bary[(ke + 1) % 3];
    bary                   = tmp;
    return {next, true};
  } else if (is_vert) {
    auto vid = triangles[from.face][kv];
    auto e0  = normalize(positions[triangles[from.face][(kv + 1) % 3]] - vid);
    auto e1  = normalize(positions[triangles[from.face][(kv + 2) % 3]] - vid);
    auto c0  = dot(e0, dir);
    auto c1  = dot(e1, dir);

    if (c0 < c1) {
      auto next = adjacencies[from.face][kv];
      kv        = find(triangles[next], vid);
      bary      = zero3f;
      bary[kv]  = 1;
      return {next, true};
    } else {
      auto next = adjacencies[from.face][(kv + 2) % 3];
      kv        = find(triangles[next], vid);
      bary      = zero3f;
      bary[kv]  = 1;
      return {next, true};
    }
  }
}
mesh_point step_for_gradient_descent_minimization(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const mesh_point& from, vec3f& dir) {
  auto [is_vert, kv] = bary_is_vert(get_bary(from.uv));
  auto tid           = from.face;
  if (is_vert) {
    tid = next_tid(solver, angles, positions, v2t, triangles, normals,
        triangles[tid][kv], dir);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, triangles[tid][kv], tid, V2T);
  }
  auto next_bary = zero3f;
  auto next_pos  = zero3f;
  trace_in_triangles(
      positions, triangles, dir, get_bary(from.uv), tid, next_pos, next_bary);
  std::tie(is_vert, kv) = bary_is_vert(next_bary);
  auto [is_edge, ke]    = bary_is_vert(next_bary);
  if (is_edge) {
    auto next_bary = tri_bary_coords(triangles, positions, adjacencies[tid][ke],
        next_pos);  // TO DO: implement the right way to this
    return {adjacencies[tid][ke], {next_bary.y, next_bary.z}};
  } else if (is_vert) {
    auto vid = triangles[tid][kv];
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, tid, vid, T2V);
    auto next = next_tid(
        solver, angles, positions, v2t, triangles, normals, vid, dir);
    kv            = find(triangles[next], vid);
    next_bary     = zero3f;
    next_bary[kv] = 1;
    return {next, {next_bary.y, next_bary.z}};
  }

  return {-1, zero2f};
}

mesh_point gradient_descent_quadratic_interpolation(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& normals, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& E,
    const vector<float>& weights, const int tid, const int vid0) {
  auto blended                 = field_blending(solver, f, weights, vid0, 10);
  auto [arrived, interp, bary] = tri_contains_min(
      triangles, positions, adjacencies, v2t, solver, blended, E, weights, tid);
  auto prev = tid;
  auto curr = tid;
  auto next = -1;
  while (!arrived) {
    sort(interp.begin(), interp.end());
    next = opposite_face(triangles, adjacencies, curr, interp[0].second);
    if (next == prev) {
      auto vid        = interp[0].second;
      auto k          = find(triangles[curr], vid);
      auto is_on_edge = edge_contains_min(
          triangles, blended, E, weights, curr, k);
      if (is_on_edge.first) {
        auto bary         = zero2f;
        bary[k]           = is_on_edge.second;
        bary[(k + 1) % 3] = 1 - is_on_edge.second;
        return {curr, bary};
      } else {
        bary              = zero2f;
        bary[k]           = 0.5;
        bary[(k + 1) % 3] = 0.5;
        return {tid, bary};
      }

    } else {
      std::tie(arrived, interp, bary) = tri_contains_min(triangles, positions,
          adjacencies, v2t, solver, blended, E, weights, next);
      if (arrived) return {next, bary};
      prev = curr;
      curr = next;
    }
  }
  return {curr, bary};
}
// prova con AGS
mesh_point gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<float>>& f,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights,
    const mesh_point& from) {
  auto bary = zero3f;
  auto tid = from.face, prev = -1, next = -1;
  auto gx = gradient_blending(solver, angles, total_angles, positions,
      triangles, adjacencies, v2t, normals, triangles[tid].x, gradients,
      weights);
  auto gy = gradient_blending(solver, angles, total_angles, positions,
      triangles, adjacencies, v2t, normals, triangles[tid].y, gradients,
      weights);
  auto gz = gradient_blending(solver, angles, total_angles, positions,
      triangles, adjacencies, v2t, normals, triangles[tid].z, gradients,
      weights);

  /* auto gx = karcher_grad(gradients, f, weights, triangles[tid].x);
  auto gy = karcher_grad(gradients, f, weights, triangles[tid].y);
  auto gz = karcher_grad(gradients, f, weights, triangles[tid].z); */

  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary);

  while (!arrived) {
    /* next = step_for_gradient_descent(
        triangles, positions, adjacencies, bary, tid); */

    if (next == prev) {
      auto eid = common_edge(triangles[tid], triangles[prev]);
      auto ex  = gradient_blending(solver, angles, total_angles, positions,
           triangles, adjacencies, v2t, normals, eid.x, gradients, weights);
      auto ey  = gradient_blending(solver, angles, total_angles, positions,
           triangles, adjacencies, v2t, normals, eid.y, gradients, weights);
      /*  auto ex = karcher_grad(gradients, f, weights, eid.x);
       auto ey = karcher_grad(gradients, f, weights, eid.y); */

      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, ex, normals, eid.x, tid, V2T);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, ey, normals, eid.y, tid, V2T);

      auto g01   = dot(ex, ey);
      auto alpha = (length_squared(ey) - g01) / length_squared(ex - ey);

      auto k            = find(triangles[tid], eid.x);
      bary              = zero3f;
      bary[k]           = alpha;
      bary[(k + 1) % 3] = 1 - alpha;
      return {tid, vec2f{bary.y, bary.z}};
    } else {
      gx = gradient_blending(solver, angles, total_angles, positions, triangles,
          adjacencies, v2t, normals, triangles[next].x, gradients, weights);
      gy = gradient_blending(solver, angles, total_angles, positions, triangles,
          adjacencies, v2t, normals, triangles[next].y, gradients, weights);
      gz = gradient_blending(solver, angles, total_angles, positions, triangles,
          adjacencies, v2t, normals, triangles[next].z, gradients, weights);
      /* gx = karcher_grad(gradients, f, weights, triangles[next].x);
      gy = karcher_grad(gradients, f, weights, triangles[next].y);
      gz = karcher_grad(gradients, f, weights, triangles[next].z); */
      std::tie(arrived, type) = tri_contains_min_hessian(solver, angles,
          total_angles, triangles, positions, adjacencies, next, gx, gy, gz,
          normals, bary);

      if (arrived) {
        return {next, {bary.y, bary.z}};
      } else {
        prev = tid;
        tid  = next;
      }
    }
  }

  return {tid, {bary.y, bary.z}};
}

mesh_point almost_gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const dual_geodesic_solver& dual_solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<mesh_point>& control_points, const vector<float>& weights,
    const mesh_point& from, vector<int>& tids,
    vector<pair<mesh_point, float>>& badones) {
  auto tags      = vector<int>(triangles.size(), 0);
  auto bary      = zero3f;
  auto prev_bary = zero3f;

  auto tid = from.face, prev = -1, next = -1;

  auto al_grads = almost_gradient_blending(dual_solver, solver, triangles,
      positions, adjacencies, normals, angles, total_angles, weights, tid,
      control_points);
  tids.push_back(tid);

  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, al_grads[0], al_grads[1],
      al_grads[2], normals, bary);

  if (type == infty) {
    badones.push_back(
        std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
    return mesh_point{tid, vec2f{bary.y, bary.z}};
  }

  prev_bary = bary;
  tags[tid] += 1;
  auto dir = bary.x * al_grads[0] + bary.y * al_grads[1] + bary.z * al_grads[2];
  if (length(dir) < 1e-6) return {tid, {bary.y, bary.z}};
  auto curr_sample = mesh_point{tid, {bary.y, bary.z}};
  while (!arrived && tags[tid] < 2) {
    curr_sample = step_for_gradient_descent_minimization(solver, angles,
        total_angles, triangles, positions, adjacencies, normals, v2t,
        curr_sample, dir);
    next        = curr_sample.face;
    tids.push_back(next);
    tags[next] += 1;
    if (next == prev) {
      if (type == infty) {
        badones.push_back(
            std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
        return mesh_point{tid, vec2f{bary.y, bary.z}};

      } else {
        auto eid = common_edge(triangles[tid], triangles[prev]);
        auto k   = find(triangles[tid], eid.x);
        bary     = zero3f;
        bary[k]  = 1.f;
        auto ex  = almost_gradient(dual_solver, solver, triangles, positions,
             adjacencies, normals, angles, total_angles, weights,
             mesh_point{tid, {bary.y, bary.z}}, control_points);
        bary     = zero3f;
        bary[(k + 1) % 3] = 1.f;
        auto ey = almost_gradient(dual_solver, solver, triangles, positions,
            adjacencies, normals, angles, total_angles, weights,
            mesh_point{tid, {bary.y, bary.z}}, control_points);

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ex, normals, eid.x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ey, normals, eid.y, tid, V2T);

        auto g01   = dot(ex, ey);
        auto alpha = (length_squared(ey) - g01) / length_squared(ex - ey);

        bary              = zero3f;
        bary[k]           = alpha;
        bary[(k + 1) % 3] = 1 - alpha;
        return {tid, vec2f{bary.y, bary.z}};
      }

    } else {
      al_grads = almost_gradient_blending(dual_solver, solver, triangles,
          positions, adjacencies, normals, angles, total_angles, weights, next,
          control_points);

      std::tie(arrived, type) = tri_contains_min_hessian(solver, angles,
          total_angles, triangles, positions, adjacencies, next, al_grads[0],
          al_grads[1], al_grads[2], normals, bary);
      if (type == infty) {
        badones.push_back(
            std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
        return mesh_point{tid, vec2f{bary.y, bary.z}};

      } else if (arrived) {
        return {next, {bary.y, bary.z}};
      }
      dir = bary.x * al_grads[0] + bary.y * al_grads[1] + bary.z * al_grads[2];
      if (length(dir) < 1e-6) return {next, {bary.y, bary.z}};
      curr_sample = mesh_point{tid, {bary.y, bary.z}};
      prev        = tid;
      tid         = next;
      prev_bary   = bary;
    }
  }

  // if (type == saddle)
  //   badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}},
  //   -1));
  return {tid, {bary.y, bary.z}};
}
mesh_point gradient_descent_test(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const dual_geodesic_solver& dual_solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<mesh_point>&    control_points,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights,
    const mesh_point& from, vector<int>& tids,
    vector<pair<mesh_point, int>>& badones) {
  auto tags      = vector<int>(triangles.size(), 0);
  auto bary      = get_bary(from.uv);
  auto prev_bary = bary;
  auto tid = from.face, prev = -1, next = -1;
  auto curr_sample = from;
  auto is_descent  = false;
  auto gx          = gradient_blending(gradients, weights, triangles[tid].x);
  auto gy          = gradient_blending(gradients, weights, triangles[tid].y);
  auto gz          = gradient_blending(gradients, weights, triangles[tid].z);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gx, normals, triangles[tid].x, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gy, normals, triangles[tid].y, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gz, normals, triangles[tid].z, tid, V2T);
  tids.push_back(tid);

  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary, true);
  // assert(type == -1);
  if (type == infty) {
    badones.push_back(
        std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
    return mesh_point{tid, vec2f{bary.y, bary.z}};
    auto al_grad = almost_gradient(dual_solver, solver, triangles, positions,
        adjacencies, normals, angles, total_angles, weights, from,
        control_points);
    if (length(al_grad) < 1e-5) return from;

    return step_for_gradient_descent_minimization(solver, angles, total_angles,
        triangles, positions, adjacencies, normals, v2t, from, al_grad);
  }

  tags[tid] += 1;
  auto count = 0;
  while (!arrived && tags[tid] < 2) {
    //++count;
    // auto grad = (gx + gy + gz) / 3.f;
    auto grad = prev_bary.x * gx + prev_bary.y * gy + prev_bary.z * gz;
    std::pair(curr_sample, is_descent) = step_for_gradient_descent_minimization(
        solver, triangles, positions, adjacencies, v2t, normals, angles,
        total_angles, curr_sample, -gx, -gy, -gz, -grad);
    next      = curr_sample.face;
    prev_bary = get_bary(curr_sample.uv);
    if (!is_descent) {
      auto last_good_point = mesh_point{tid, vec2f{bary.y, bary.z}};
      badones.push_back(std::make_pair(last_good_point, -1));
      // if (count != 0) show(count);
      return last_good_point;
    }

    tids.push_back(next);
    tags[next] += 1;
    if (next == prev) {
      if (type == infty) {
        badones.push_back(
            std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
        return mesh_point{tid, vec2f{bary.y, bary.z}};
        auto al_grad = almost_gradient(dual_solver, solver, triangles,
            positions, adjacencies, normals, angles, total_angles, weights,
            from, control_points);
        if (length(al_grad) < 1e-5) return from;

        return step_for_gradient_descent_minimization(solver, angles,
            total_angles, triangles, positions, adjacencies, normals, v2t, from,
            al_grad);

      } else {
        auto eid = common_edge(triangles[tid], triangles[prev]);
        auto ex  = gradient_blending(gradients, weights, eid.x);
        auto ey  = gradient_blending(gradients, weights, eid.y);

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ex, normals, eid.x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ey, normals, eid.y, tid, V2T);

        auto g01   = dot(ex, ey);
        auto alpha = (length_squared(ey) - g01) / length_squared(ex - ey);

        auto k            = find(triangles[tid], eid.x);
        bary              = zero3f;
        bary[k]           = alpha;
        bary[(k + 1) % 3] = 1 - alpha;
        // if (count != 0) show(count);
        return {tid, vec2f{bary.y, bary.z}};
      }

    } else {
      gx = gradient_blending(gradients, weights, triangles[next].x);
      gy = gradient_blending(gradients, weights, triangles[next].y);
      gz = gradient_blending(gradients, weights, triangles[next].z);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gx, normals, triangles[next].x, next, V2T);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gy, normals, triangles[next].y, next, V2T);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gz, normals, triangles[next].z, next, V2T);
      std::tie(arrived, type) = tri_contains_min_hessian(solver, angles,
          total_angles, triangles, positions, adjacencies, next, gx, gy, gz,
          normals, bary, true);
      // assert(type == -1);
      if (type == infty) {
        badones.push_back(
            std::make_pair(mesh_point{tid, vec2f{bary.y, bary.z}}, -1));
        return mesh_point{tid, vec2f{bary.y, bary.z}};
        auto al_grad = almost_gradient(dual_solver, solver, triangles,
            positions, adjacencies, normals, angles, total_angles, weights,
            from, control_points);
        if (length(al_grad) < 1e-5) return from;

        return step_for_gradient_descent_minimization(solver, angles,
            total_angles, triangles, positions, adjacencies, normals, v2t, from,
            al_grad);

      } else if (arrived) {
        // if (count != 0) show(count);
        return {next, {bary.y, bary.z}};
      } else {
        prev = tid;
        tid  = next;
        bary = prev_bary;
      }
    }
  }

  // if (type == saddle)
  //   badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}},
  //   -1));
  // if (count != 0) show(count);
  return {tid, {bary.y, bary.z}};
}
std::tuple<mesh_point, bool, int, vec3f, vec3f, vec3f> find_min_along_eta(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& normals, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights,
    const vector<mesh_point>& straightest) {
  pair<int, vec3f> min_on_path = {-1, zero3f};
  auto             lambda      = flt_max;
  auto             bary        = zero3f;
  mesh_point       sample      = {};
  int              tid;
  vec3f            gx, gy, gz;
  for (auto i = 0; i < straightest.size() - 1; ++i) {
    auto bary_in  = get_bary(straightest[i].uv);
    auto bary_out = zero3f;
    if (straightest[i].face == straightest[i + 1].face) {
      bary_out = get_bary(straightest[i + 1].uv);
      tid      = straightest[i].face;
      gx       = gradient_blending(gradients, weights, triangles[tid].x);
      gy       = gradient_blending(gradients, weights, triangles[tid].y);
      gz       = gradient_blending(gradients, weights, triangles[tid].z);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gx, normals, triangles[tid].x, tid, V2T);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gy, normals, triangles[tid].y, tid, V2T);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, gz, normals, triangles[tid].z, tid, V2T);
    } else {
      bary_out = tri_bary_coords(triangles, positions, straightest[i].face,
          eval_position(triangles, positions, straightest[i + 1]));
      if (!are_barycentric_coordinates(bary_out, 1e-1)) {
        bary_out = get_bary(straightest[i + 1].uv);
        bary_in = tri_bary_coords(triangles, positions, straightest[i + 1].face,
            eval_position(triangles, positions, straightest[i]));
        assert(are_barycentric_coordinates(bary_in, 1e-1));
        tid = straightest[i + 1].face;
        gx  = gradient_blending(gradients, weights, triangles[tid].x);
        gy  = gradient_blending(gradients, weights, triangles[tid].y);
        gz  = gradient_blending(gradients, weights, triangles[tid].z);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gx, normals, triangles[tid].x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gy, normals, triangles[tid].y, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gz, normals, triangles[tid].z, tid, V2T);
      } else {
        tid = straightest[i].face;
        gx  = gradient_blending(gradients, weights, triangles[tid].x);
        gy  = gradient_blending(gradients, weights, triangles[tid].y);
        gz  = gradient_blending(gradients, weights, triangles[tid].z);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gx, normals, triangles[tid].x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gy, normals, triangles[tid].y, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, gz, normals, triangles[tid].z, tid, V2T);
      }
    }

    auto [is_in, bary_of_min, value] = minimum_on_line(
        bary_in, bary_out, gx, gy, gz);
    if (is_in && value < lambda) {
      lambda      = value;
      min_on_path = {tid, bary_of_min};
    }
    //}
  }
  if (min_on_path.first == -1)
    // return {straightest.back(), false};
    sample = straightest.back();

  else {
    // return {{min_on_path.first, {min_on_path.second.y,
    // min_on_path.second.z}},
    //     false};
    sample = {min_on_path.first, {min_on_path.second.y, min_on_path.second.z}};
  }
  tid = sample.face;
  gx  = gradient_blending(gradients, weights, triangles[tid].x);
  gy  = gradient_blending(gradients, weights, triangles[tid].y);
  gz  = gradient_blending(gradients, weights, triangles[tid].z);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gx, normals, triangles[tid].x, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gy, normals, triangles[tid].y, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gz, normals, triangles[tid].z, tid, V2T);
  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary, true);
  if (arrived) sample = {tid, {bary.y, bary.z}};
  return {sample, arrived, type, gx, gy, gz};
}
mesh_point gradient_descent_absil(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const mesh_point& from, vector<int>& tids,
    vector<vec3f>& directions, vector<mesh_point>& steps,
    vector<pair<mesh_point, int>>& badones) {
  auto       tid = from.face, prev = -1, next = -1;
  auto       tags        = vector<int>(triangles.size(), 0);
  auto       curr_sample = from;  // mesh_point{tid, {0.33, 0.33}};
  mesh_point next_sample = {-1, zero2f};
  vec3f      prev_dir    = {flt_max, 0, 0};
  auto       bary        = zero3f;
  auto       gx = gradient_blending(gradients, weights, triangles[tid].x);
  auto       gy = gradient_blending(gradients, weights, triangles[tid].y);
  auto       gz = gradient_blending(gradients, weights, triangles[tid].z);

  // tags[tid] += 1;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gx, normals, triangles[tid].x, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gy, normals, triangles[tid].y, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gz, normals, triangles[tid].z, tid, V2T);
  auto grad = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
              curr_sample.uv.x * gy + curr_sample.uv.y * gz;

  auto [dir, is_descent_direction] = maximum_descent_direction(
      triangles, positions, tid, -grad, -gx, -gy, -gz);
  if (!is_descent_direction) {
    badones.push_back(std::make_pair(from, -1));
    return from;
  }

  directions.push_back(dir);
  steps.push_back(curr_sample);
  tids.push_back(tid);
  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary, true);
  if (type == infty) {
    badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}}, -1));
    return {tid, vec2f{bary.y, bary.z}};
  }

  auto count = 0;
  while (!arrived && tags[tid] <= 2) {
    // the minus in the input is due to the fact that we compute -grad(f) by
    // default
    ++count;

    auto line = straightest_geodesic(solver, triangles, positions, normals,
        adjacencies, v2t, angles, total_angles, curr_sample, normalize(dir),
        length(dir));

    std::tie(next_sample, arrived, type, gx, gy, gz) = find_min_along_eta(
        triangles, positions, adjacencies, v2t, normals, solver, angles,
        total_angles, gradients, weights, line);
    if (arrived) return next_sample;
    if (type == infty) {
      badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}}, -1));
      return {tid, vec2f{bary.y, bary.z}};
    }
    auto step = length(eval_position(triangles, positions, next_sample) -
                       eval_position(triangles, positions, curr_sample));

    next = next_sample.face;

    /* straightest_geodesic(solver, triangles,
         positions, normals, adjacencies, v2t, angles, total_angles,
         curr_sample, normalize(dir), length(dir)) .back() .face;
         next = next_sample.face;  */

    // if (next == tid && dot(dir, prev_dir) < 0) {
    //   badones.push_back(std::make_pair(from, -1));
    //   return from;
    // }
    // tags[next] += 1;

    if (next == prev && prev != tid) {
      auto eid = common_edge(triangles[tid], triangles[prev]);
      if (eid.x != -1) {
        auto ex = gradient_blending(gradients, weights, eid.x);
        auto ey = gradient_blending(gradients, weights, eid.y);

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ex, normals, eid.x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ey, normals, eid.y, tid, V2T);

        auto g01   = dot(ex, ey);
        auto alpha = length_squared(ey - g01) / length_squared(ex - ey);

        auto k            = find(triangles[tid], eid.x);
        bary              = zero3f;
        bary[k]           = alpha;
        bary[(k + 1) % 3] = 1 - alpha;

        return {tid, vec2f{bary.y, bary.z}};
      }
      // else {
      //   badones.push_back(
      //       std::make_pair(mesh_point{next, {bary.y, bary.z}}, -1));
      //   return from;
      // }
    }  // else

    if (next == tid && step <= 1e-6) {
      badones.push_back(std::make_pair(mesh_point{next, {bary.y, bary.z}}, -1));
      return from;

    } else {
      // gx = gradient_blending(gradients, weights, triangles[next].x);
      // gy = gradient_blending(gradients, weights, triangles[next].y);
      // gz = gradient_blending(gradients, weights, triangles[next].z);

      // parallel_transp(solver, angles, total_angles, triangles, positions,
      //     adjacencies, gx, normals, triangles[next].x, next, V2T);
      // parallel_transp(solver, angles, total_angles, triangles, positions,
      //     adjacencies, gy, normals, triangles[next].y, next, V2T);
      // parallel_transp(solver, angles, total_angles, triangles, positions,
      //     adjacencies, gz, normals, triangles[next].z, next, V2T);

      // std::tie(arrived, type) = tri_contains_min_hessian(solver, angles,
      //     total_angles, triangles, positions, adjacencies, next, gx, gy, gz,
      //     normals, bary, true);
      // if (type == infty) {
      //   badones.push_back(
      //       std::make_pair(mesh_point{next, {bary.y, bary.z}}, -1));
      //   return {next, vec2f{bary.y, bary.z}};
      // }
      // if (count > 50) {
      //   badones.push_back(std::make_pair(from, -1));
      //   return from;
      // }

      prev_dir    = dir;
      prev        = tid;
      tid         = next;
      curr_sample = next_sample;
      grad        = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
             curr_sample.uv.x * gy + curr_sample.uv.y * gz;
      std::tie(dir, is_descent_direction) = maximum_descent_direction(
          triangles, positions, tid, -grad, -gx, -gy, -gz);
      if (!is_descent_direction) {
        badones.push_back(std::make_pair(curr_sample, -1));
        return from;
      }
      directions.push_back(dir);
      steps.push_back(curr_sample);
      tids.push_back(tid);

      if (length(grad) <= 1e-5) return curr_sample;
    }
  }
  /* if (tags[tid] > 2) {
    badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}},
  -1)); assert(next != -1);
  } */
  return {tid, vec2f{bary.y, bary.z}};
}
mesh_point gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const mesh_point& from) {
  auto       tid = from.face, prev = -1, next = -1;
  auto       curr_sample = from;
  mesh_point next_sample = {-1, zero2f};
  auto       bary        = zero3f;
  auto       gx = gradient_blending(gradients, weights, triangles[tid].x);
  auto       gy = gradient_blending(gradients, weights, triangles[tid].y);
  auto       gz = gradient_blending(gradients, weights, triangles[tid].z);

  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gx, normals, triangles[tid].x, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gy, normals, triangles[tid].y, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gz, normals, triangles[tid].z, tid, V2T);
  auto grad = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
              curr_sample.uv.x * gy + curr_sample.uv.y * gz;

  auto [dir, is_descent_direction] = maximum_descent_direction(
      triangles, positions, tid, -grad, -gx, -gy, -gz);
  if (!is_descent_direction) {
    return {-1, zero2f};
  }

  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary, true);
  if (type == infty) {
    return {-1, zero2f};
  }

  while (!arrived) {
    auto line = straightest_geodesic(solver, triangles, positions, normals,
        adjacencies, v2t, angles, total_angles, curr_sample, normalize(dir),
        length(dir));

    std::tie(next_sample, arrived, type, gx, gy, gz) = find_min_along_eta(
        triangles, positions, adjacencies, v2t, normals, solver, angles,
        total_angles, gradients, weights, line);
    if (arrived) return next_sample;
    if (type == infty) {
      return {-1, zero2f};
    }
    auto step = length(eval_position(triangles, positions, next_sample) -
                       eval_position(triangles, positions, curr_sample));
    next      = next_sample.face;

    if (next == prev && prev != tid) {
      auto eid = common_edge(triangles[tid], triangles[prev]);
      if (eid.x != -1) {
        auto ex = gradient_blending(gradients, weights, eid.x);
        auto ey = gradient_blending(gradients, weights, eid.y);

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ex, normals, eid.x, tid, V2T);
        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, ey, normals, eid.y, tid, V2T);

        auto g01   = dot(ex, ey);
        auto alpha = (length_squared(ey) - g01) / length_squared(ex - ey);

        auto k            = find(triangles[tid], eid.x);
        bary              = zero3f;
        bary[k]           = alpha;
        bary[(k + 1) % 3] = 1 - alpha;

        return {tid, vec2f{bary.y, bary.z}};
      }
    }

    if (next == tid && step <= 1e-6) {
      return {-1, zero2f};

    } else {
      prev        = tid;
      tid         = next;
      curr_sample = next_sample;
      grad        = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
             curr_sample.uv.x * gy + curr_sample.uv.y * gz;
      std::tie(dir, is_descent_direction) = maximum_descent_direction(
          triangles, positions, tid, -grad, -gx, -gy, -gz);
      if (!is_descent_direction) {
        return {-1, zero2f};
      }
      if (length(grad) <= 1e-5) return curr_sample;
    }
  }

  return {tid, vec2f{bary.y, bary.z}};
}
pair<int, vec3f> handle_vert(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const int vid, const int tid, const vec3f dir) {
  auto v = dir;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, v, normals, tid, vid, T2V);
  auto next = next_tid(
      solver, angles, positions, v2p_adjacencies, triangles, normals, vid, v);

  /* auto next = next_tid_extended_graph(
      solver, angles, positions, v2p_adjacencies, triangles, normals, vid,
     v);
   */

  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, v, normals, vid, next, V2T);
  return std::make_pair(next, v);
}
std::tuple<mesh_point, vector<int>> straightest_for_gradient_descent(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>&         adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const mesh_point& from, const vec3f& v, const float& l) {
  auto  vid = -1, tid = -1, next_tri = -1;
  auto  dir       = v;
  float len       = 0.0;
  auto  next_bary = zero3f, next = zero3f;
  auto  prev            = eval_position(triangles, positions, from);
  auto  bary            = get_bary(from.uv);
  auto [is_vert, kv]    = bary_is_vert(bary);
  auto [is_on_edge, ke] = bary_is_edge(bary);
  vector<int> strip;
  if (is_vert) {
    vid = triangles[from.face][kv];
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, from.face, vid, T2V);

    tid = next_tid(
        solver, angles, positions, v2p_adjacencies, triangles, normals, vid, v);

    /*  tid = next_tid_extended_graph(
         solver, angles, positions, v2p_adjacencies, triangles, normals,
       vid, v); */
    strip.push_back(tid);
    kv       = find(triangles[tid], vid);
    bary     = zero3f;
    bary[kv] = 1;
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, vid, tid, V2T);
  } else if (is_on_edge) {
    auto p0   = triangles[from.face][ke];
    auto p1   = triangles[from.face][(ke + 1) % 3];
    auto p2   = triangles[from.face][(ke + 2) % 3];
    auto n    = triangle_normal(positions[p0], positions[p1], positions[p2]);
    auto edge = normalize(positions[p1] - positions[p0]);
    if (dot(cross(edge, v), n) > 0)
      tid = from.face;
    else {
      tid  = adjacencies[from.face][ke];
      bary = tri_bary_coords(triangles, positions, tid, prev);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, from.face, tid, T2T);
      strip.push_back(tid);
    }

  } else
    tid = from.face;

  while (len < l) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    len += length(next - prev);
    if (len < l) {
      prev          = next;
      auto [V, k_v] = bary_is_vert(next_bary);
      auto [E, k_e] = bary_is_edge(next_bary);
      if (V) {
        vid      = triangles[tid][k_v];
        auto out = handle_vert(solver, triangles, positions, normals,
            adjacencies, v2p_adjacencies, angles, total_angles, vid, tid, dir);
        tid      = out.first;
        strip.push_back(tid);
        dir       = out.second;
        k_v       = find(triangles[tid], vid);
        bary      = zero3f;
        bary[k_v] = 1;
      } else if (E) {
        next_tri = adjacencies[tid][k_e];
        strip.push_back(next_tri);
        auto p0       = triangles[tid][k_e];
        auto offset0  = find(triangles[next_tri], p0);
        auto offset1  = (offset0 + 2) % 3;
        bary          = zero3f;
        bary[offset0] = next_bary[k_e];
        bary[offset1] = next_bary[(k_e + 1) % 3];

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, dir, normals, tid, next_tri, T2T);
        tid = next_tri;
      } else
        assert(false);
    }
  }

  auto factor = (len - l);
  auto w      = normalize(prev - next);
  w *= factor;
  w += next;
  bary = tri_bary_coords(triangles, positions, tid, w);

  return {{tid, vec2f{bary.y, bary.z}}, strip};
}
vector<int> strip_by_gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3f>& gradients,
    const mesh_point& from, vector<vec3f>& directions,
    vector<mesh_point>& steps, vector<pair<mesh_point, int>>& badones) {
  auto       tid = from.face, prev = -1, next = -1;
  auto       tags        = vector<int>(triangles.size(), 0);
  auto       curr_sample = from;  // mesh_point{tid, {0.33, 0.33}};
  mesh_point next_sample = {-1, zero2f};
  vec3f      prev_dir    = {flt_max, 0, 0};
  auto       bary        = zero3f;
  auto       gx          = gradients[triangles[tid].x];
  auto       gy          = gradients[triangles[tid].y];
  auto       gz          = gradients[triangles[tid].z];

  // tags[tid] += 1;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gx, normals, triangles[tid].x, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gy, normals, triangles[tid].y, tid, V2T);
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, gz, normals, triangles[tid].z, tid, V2T);
  auto grad = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
              curr_sample.uv.x * gy + curr_sample.uv.y * gz;

  auto [arrived, type] = tri_contains_min_hessian(solver, angles, total_angles,
      triangles, positions, adjacencies, tid, gx, gy, gz, normals, bary, true);
  if (type == infty) {
    badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}}, -1));
    return {tid};
  }

  auto        count   = 0;
  vector<int> strip   = {tid};
  vector<int> segment = {};
  while (!arrived) {
    // the minus in the input is due to the fact that we compute -grad(f) by
    // default
    ++count;

    auto [dir, is_descent_direction] = maximum_descent_direction(
        triangles, positions, tid, -grad, -gx, -gy, -gz);
    if (!is_descent_direction) {
      badones.push_back(std::make_pair(from, -1));
      return strip;
    }

    directions.push_back(dir);

    std::tie(next_sample, segment) = straightest_for_gradient_descent(solver,
        triangles, positions, normals, adjacencies, v2t, angles, total_angles,
        curr_sample, normalize(dir), length(dir));

    /*  find_min_along_eta(triangles, positions, adjacencies, v2t,
         normals, solver, angles, total_angles, gradients, weights,
         line);  */

    next = next_sample.face;
    strip.insert(strip.end(), segment.begin(), segment.end());
    steps.push_back(next_sample);

    gx = gradients[triangles[next].x];
    gy = gradients[triangles[next].y];
    gz = gradients[triangles[next].z];

    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, gx, normals, triangles[next].x, next, V2T);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, gy, normals, triangles[next].y, next, V2T);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, gz, normals, triangles[next].z, next, V2T);

    std::tie(arrived, type) = tri_contains_min_hessian(solver, angles,
        total_angles, triangles, positions, adjacencies, next, gx, gy, gz,
        normals, bary, true);
    if (type == infty) {
      badones.push_back(std::make_pair(mesh_point{next, {bary.y, bary.z}}, -1));
      return strip;
    }
    /*  if (count > 50) {
       badones.push_back(std::make_pair(from, -1));
       return from;
     } */
    prev_dir    = dir;
    prev        = tid;
    tid         = next;
    curr_sample = next_sample;
    grad        = (1 - curr_sample.uv.x - curr_sample.uv.y) * gx +
           curr_sample.uv.x * gy + curr_sample.uv.y * gz;
  }
  /* if (tags[tid] > 2) {
    badones.push_back(std::make_pair(mesh_point{tid, {bary.y, bary.z}},
  -1)); assert(next != -1);
  } */
  return strip;
}

mesh_point weighted_average_quadratic(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& f,
    const vector<vector<vec3f>>& gradients, const vector<vector<vec3f>>& E,
    const vector<float>& weights, int seed = 0) {
  auto local_min = minimum_on_verts(solver, triangles, f, weights, seed);
  auto v         = normalize(gradient_blending(gradients, weights, local_min));
  auto tid       = next_tid_extended_graph(
            solver, angles, positions, v2t, triangles, normals, local_min, v);

  auto min = gradient_descent_quadratic_interpolation(triangles, positions,
      adjacencies, v2t, normals, solver, angles, total_angles, f, E, weights,
      tid, local_min);
  // auto offset  = find(triangles[tid], local_min);
  // auto bary    = zero3f;
  // bary[offset] = 1;
  // mesh_point p = {tid, {bary.y, bary.z}};
  // auto min = gradient_descent(triangles, positions, adjacencies, v2t,
  // normals,
  //     solver, angles, total_angles, gradients, weights, p);

  return min;
}
mesh_point weighted_average(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& f,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights,
    int& seed) {
  auto local_min = minimum_on_verts(solver, triangles, f, weights, seed);
  seed           = local_min;
  /* local_min      = minimum_on_gradient_norm(
      solver, triangles, f, gradients, weights, local_min); */

  auto v = normalize(gradient_blending(gradients, weights, local_min));

  // auto v = normalize(karcher_grad(gradients, f, weights, local_min));

  auto tid = next_tid(
      solver, angles, positions, v2t, triangles, normals, local_min, v);

  auto offset  = find(triangles[tid], local_min);
  auto bary    = zero3f;
  bary[offset] = 1;
  mesh_point p = {tid, {bary.y, bary.z}};
  auto min = gradient_descent(triangles, positions, adjacencies, v2t, normals,
      solver, angles, total_angles, f, gradients, weights, p);

  return min;
}
mesh_point weighted_average_almost_gradients(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& f,
    const vector<float>& weights, const vector<mesh_point>& control_points,
    vector<pair<int, float>>& verts, vector<int>& tids,
    vector<pair<mesh_point, float>>& badones, int& seed) {
  auto local_min = minimum_on_verts(solver, triangles, f, weights, seed);
  /* local_min      = minimum_on_gradient_norm(
      solver, triangles, f, gradients, weights, local_min); */
  /* auto min_on_norm = minimum_on_gradient_norm(
      solver, triangles, f, gradients, weights, local_min); */
  //  if (min_on_norm != local_min) {
  //    auto entry = node_is_adjacent(solver, local_min, min_on_norm);
  //    assert(entry != -1 && entry % 2 == 0);
  //  }

  seed   = local_min;
  auto p = make_mesh_point(triangles, v2t, local_min);

  verts.push_back(std::make_pair(local_min, -1));

  auto v = normalize(almost_gradient(dual_solver, solver, triangles, positions,
      adjacencies, normals, angles, total_angles, weights, p, control_points));

  auto tid = next_tid(
      solver, angles, positions, v2t, triangles, normals, local_min, v);
  p        = make_mesh_point(triangles, v2t, local_min, tid);
  auto min = almost_gradient_descent(triangles, positions, adjacencies, v2t,
      normals, solver, dual_solver, angles, total_angles, control_points,
      weights, p, tids, badones);
  return min;
}
mesh_point weighted_average_test(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const vector<mesh_point>& control_points,
    vector<pair<int, float>>& verts, vector<int>& tids,
    vector<vec3f>& directions, vector<mesh_point>& steps,
    vector<pair<mesh_point, int>>& badones, int& seed) {
  auto local_min = minimum_on_verts(solver, triangles, f, weights, seed);
  /*auto entry     = node_is_adjacent(solver, seed, local_min);
  auto k         = 1;
  auto old_t     = t - step;
  while (entry < 0) {
     auto new_t = t + step / pow(2, k);
     bernstein_polynomials(control_points.size(), new_t, weights);
     local_min = minimum_on_verts(solver, triangles, f, weights, seed);
     entry     = node_is_adjacent(solver, seed, local_min);
     ++k;
   } */
  /*  local_min = minimum_on_gradient_norm(
       solver, triangles, f, gradients, weights, local_min); */
  /* auto min_on_norm = minimum_on_gradient_norm(
      solver, triangles, f, gradients, weights, local_min); */
  //  if (min_on_norm != local_min) {
  //    auto entry = node_is_adjacent(solver, local_min, min_on_norm);
  //    assert(entry != -1 && entry % 2 == 0);
  //  }

  seed = local_min;

  /*if (length(grad) < 1e-9) {
    auto k    = find(triangles[v2t[local_min][0]], local_min);
    auto bary = zero3f;
    bary[k]   = 1;
    return {v2t[local_min][0], vec2f{bary.y, bary.z}};
  }*/
  verts.push_back(std::make_pair(local_min, -1));
  auto v = normalize(gradient_blending(solver, angles, total_angles, positions,
      triangles, adjacencies, v2t, normals, local_min, gradients, weights));

  // auto v   = karcher_grad(gradients, f, weights, local_min);
  auto tid = next_tid(
      solver, angles, positions, v2t, triangles, normals, local_min, v);
  /* next_tid_extended_graph(
      solver, angles, positions, v2t, triangles, normals, local_min, v); */

  auto p = make_mesh_point(triangles, v2t, local_min, tid);
  auto min =
      /*auto min = gradient_descent_test(triangles, positions, adjacencies,
         v2t, normals, solver, dual_solver, angles, total_angles,
         control_points, gradients, weights, p, tids, badones); */
      gradient_descent_absil(triangles, positions, adjacencies, v2t, normals,
          solver, angles, total_angles, gradients, weights, p, tids, directions,
          steps, badones);
  return min;
}

vector<mesh_point> bezier_karcher_full_gradient(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>&             control_points) {
  //
  //  int n = control_points.size();
  //
  //  vector<mesh_point>    poly_line = {control_points[0]};
  //  vector<float>         w;
  //  vector<vector<vec3f>> gradients(n);
  //  vector<vector<float>> f(n);
  //  static auto           old_points = vector<mesh_point>(4, {-1, {-1,
  //  -1}}); static auto           cache      = vector<vector<float>>{};
  //  auto points_moved = control_points_moved(old_points, control_points);
  //  for (auto i = 0; i < n; ++i) {
  //    if (!points_moved[i]) {
  //      f[i] = cache[i];
  //    } else {
  //      f[i] = compute_distance_field(triangles, positions, adjacencies,
  //      v2t,
  //          solver, i, control_points, exact);
  //      std::transform(f[i].begin(), f[i].end(), f[i].begin(),
  //          [](float lambda) { return lambda * lambda; });
  //    }
  //  }
  //  auto E = vector<vector<vec3f>>(n, vector<vec3f>(triangles.size(),
  //  zero3f)); for (auto i = 0; i < triangles.size(); ++i) {
  //    for (auto j = 0; j < n; ++j) {
  //      auto tid0 = adjacencies[i][0];
  //      auto tid1 = adjacencies[i][1];
  //      auto tid2 = adjacencies[i][2];
  //      auto k0   = find(adjacencies[tid0], i);
  //      auto k1   = find(adjacencies[tid1], i);
  //      auto k2   = find(adjacencies[tid2], i);
  //
  //      E[j][i][0] = (E[j][tid0][k0] != 0)
  //                       ? E[j][tid0][k0]
  //                       : average_field(triangles, positions,
  //                       adjacencies, v2t,
  //                             solver, f[j], triangles[i].x,
  //                             triangles[i].y);
  //      E[j][i][1] = (E[j][tid1][k1] != 0)
  //                       ? E[j][tid1][k1]
  //                       : average_field(triangles, positions,
  //                       adjacencies, v2t,
  //                             solver, f[j], triangles[i].y,
  //                             triangles[i].z);
  //      E[j][i][2] = (E[j][tid2][k2] != 0)
  //                       ? E[j][tid2][k2]
  //                       : average_field(triangles, positions,
  //                       adjacencies, v2t,
  //                             solver, f[j], triangles[i].z,
  //                             triangles[i].x);
  //    }
  //  }
  //  cache = f;
  //
  //  for (auto i = 0; i < n; ++i) {
  //    auto F       = wrapper(f[i]);
  //    gradients[i] = compute_grad(solver, triangles, positions, normals,
  //    Grad, F);
  //  }
  //  double step = 1 / pow(2, 16);
  //  bernstein_polynomials(n, 0.0, w);
  //  // auto bary     = get_bary(control_points[1].uv);
  //  // auto max_bary = flt_min;
  //  // auto entry    = -1;
  //  // for (auto i = 0; i < 3; ++i) {
  //  //   if (bary[i] > max_bary) {
  //  //     max_bary = bary[i];
  //  //     entry    = i;
  //  //   }
  //  // }
  //  auto tid      = control_points[0].face;
  //  auto sentinel = pair<int, float>{tid, 0.0};
  //  bernstein_polynomials(n, sentinel.second, w);
  //  auto vid0 = triangles[tid].x;
  //  auto vid1 = triangles[tid].y;
  //  auto vid2 = triangles[tid].z;
  //
  //  auto gx      = gradient_blending(gradients, w, triangles[tid].x);
  //  auto gy      = gradient_blending(gradients, w, triangles[tid].y);
  //  auto gz      = gradient_blending(gradients, w, triangles[tid].z);
  //  auto blended = field_blending(solver, f, w, vid0, 5);
  //  auto [inside, interp, bary] = tri_contains_min(
  //      triangles, positions, adjacencies, v2t, solver, blended, E, w,
  //      tid);
  //  assert(inside);
  //  while (sentinel.second < 1) {
  //    // gx = gradient_blending(gradients, w, triangles[tid].x);
  //    // gy = gradient_blending(gradients, w, triangles[tid].y);
  //    // gz = gradient_blending(gradients, w, triangles[tid].z);
  //    while (inside) {
  //      std::tie(inside, interp, bary) = tri_contains_min(
  //          triangles, positions, adjacencies, v2t, solver, blended, E, w,
  //          tid);
  //      poly_line.push_back({tid, bary});
  //      sentinel.second += step;
  //      bernstein_polynomials(n, sentinel.second, w);
  //      blended = field_blending(solver, f, w, triangles[tid].x, 5);
  //      if (sentinel.second > 1) break;
  //    }
  //    sort(interp.begin(), interp.end());
  //    tid = opposite_face(triangles, adjacencies, tid, interp[0].second);
  //    std::tie(inside, interp, bary) = tri_contains_min(
  //        triangles, positions, adjacencies, v2t, solver, blended, E, w,
  //        tid);
  //    while (!inside) {
  //      sentinel.second += step;
  //      bernstein_polynomials(n, sentinel.second, w);
  //      blended = field_blending(solver, f, w, triangles[tid].x, 5);
  //      std::tie(inside, interp, bary) = tri_contains_min(
  //          triangles, positions, adjacencies, v2t, solver, blended, E, w,
  //          tid);
  //      if (sentinel.second > 1) break;
  //    }
  //    // auto arrived = tri_contains_min(solver, angles, total_angles,
  //    // triangles,
  //    //     positions, adjacencies, tid, gx, gy, gz, normals, bary);
  //
  //    //   while (tri_contains_min(solver, angles, total_angles,
  //    triangles,
  //    //   positions,
  //    //       adjacencies, sentinel.first, g0, g1, g2, w, normals, bary))
  //    {
  //    //     sentinel.second += step;
  //    //     prev_bary = bary;
  //    //     bernstein_polynomials(n, sentinel.second, w);
  //    //   }
  //    //   mesh_point p = {sentinel.first, vec2f{prev_bary.y,
  //    prev_bary.z}};
  //    //   poly_line.push_back(p);
  //    //   auto [is_vert, kv] = bary_is_vert(prev_bary);
  //    //   auto [is_edge, ke] = bary_is_edge(prev_bary);
  //    //   if (is_vert) {
  //    //     auto  vid      = triangles[p.face][kv];
  //    //     auto  nbr      = v2t[vid];
  //    //     float t        = sentinel.second;
  //    //     sentinel.first = handle_min(solver, angles, total_angles,
  //    //     triangles,
  //    //         positions, normals, adjacencies, nbr, f, gradients, t);
  //    //     if (sentinel.first != -1)
  //    //       sentinel.second = t;
  //    //     else
  //    //       return poly_line;
  //    //   } else if (is_edge) {
  //    //     sentinel.first = adjacencies[sentinel.first][ke];
  //    //     g0             = karcher_grad(f, gradients,
  //    //     triangles[sentinel.first].x); g1             =
  //    karcher_grad(f,
  //    //     gradients, triangles[sentinel.first].y); g2             =
  //    //     karcher_grad(f, gradients, triangles[sentinel.first].z);
  //    while
  //    (
  //    //         !tri_contains_min(solver, angles, total_angles,
  //    triangles,
  //    //         positions,
  //    //             adjacencies, sentinel.first, g0, g1, g2, w, normals,
  //    //             bary))
  //    //             {
  //    //       sentinel.second += step;
  //    //       bernstein_polynomials(n, sentinel.second, w);
  //    //       if (sentinel.second > 1) {
  //    //         assert(false);
  //    //         return poly_line;
  //    //       }
  //    //     }
  //    //   }
  //
  //    //   else
  //    //     assert(false);
  //    // }
  //
  //    // if (arrived)
  //
  //    // poly_line.push_back(min);
  //  }
  //
  //  return poly_line;
  return vector<mesh_point>{};
}

pair<int, float> vector_binary_search(const vector<float>& v, float t) {
  auto L = 0, R = (int)v.size(), entry = -1, m = 0;
  auto factor = 0.0;
  while (L < R) {
    m = floor((L + R) / 2);
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

pair<mesh_point, int> sample_on_curve(const vector<vec3i>& triangles,
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
    auto i = vector_binary_search(path_parameter_t, t);
    auto p = zero3f, q = zero3f;
    auto tid = -1;
    entry    = i.first;
    if (entry == 0) {
      tid = path.start.face;
      p   = eval_position(triangles, positions, path.start);
      q   = path_pos_from_entry(triangles, positions, adjacencies, path, 0);
    } else if (entry == path_parameter_t.size() - 1) {
      return std::make_pair(path.end, entry - 1);

    } else {
      p = path_pos_from_entry(
          triangles, positions, adjacencies, path, i.first - 1);
      q   = (entry == path.strip.size() - 1)
                ? eval_position(triangles, positions, path.end)
                : path_pos_from_entry(
                      triangles, positions, adjacencies, path, i.first);
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

  return std::make_pair(a, entry);
}
mesh_point eval_point_geodesic_path(const vector<vec3i>& triangles,
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
    auto i = vector_binary_search(path_parameter_t, t);
    auto p = zero3f, q = zero3f;
    auto tid = -1;
    entry    = i.first;
    if (entry == 0) {
      tid = path.start.face;
      p   = eval_position(triangles, positions, path.start);
      q   = path_pos_from_entry(triangles, positions, adjacencies, path, 0);
    } else if (entry == path_parameter_t.size() - 1) {
      return path.end;

    } else {
      p = path_pos_from_entry(
          triangles, positions, adjacencies, path, i.first - 1);
      q   = (entry == path.strip.size() - 1)
                ? eval_position(triangles, positions, path.end)
                : path_pos_from_entry(
                      triangles, positions, adjacencies, path, i.first);
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
pair<mesh_point, int> sample_on_curve(const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& start,
    const mesh_point& end, const float& t) {
  auto strip = strip_on_dual_graph(
      dual_solver, triangles, positions, end.face, start.face);
  auto path = shortest_path(
      triangles, positions, adjacencies, start, end, strip);
  auto path_t = path_parameters(path, triangles, positions, adjacencies);

  return sample_on_curve(triangles, positions, adjacencies, path, path_t, t);
}

vec3f compute_tangent_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, const int entry) {
  if (entry == 0) {
    auto p = path_pos_from_entry(triangles, positions, adjacencies, path, 0);
    return normalize(p - eval_position(triangles, positions, path.start));
  } else if (entry == path.strip.size() - 1) {
    auto p = path_pos_from_entry(
        triangles, positions, adjacencies, path, entry - 1);
    return normalize(p - eval_position(triangles, positions, path.end));
  }
  auto tid0      = path.strip[entry - 1];
  auto tid1      = path.strip[entry];
  auto tid2      = path.strip[entry + 1];
  auto flat_tid0 = init_flat_triangle(positions, triangles[tid0]);
  auto k         = find(adjacencies[tid0], tid1);
  auto pos0      = path_pos_from_entry(
           triangles, positions, adjacencies, path, entry - 1);
  auto bary0 = tri_bary_coords(triangles, positions, tid0, pos0);
  auto pos1  = path_pos_from_entry(
       triangles, positions, adjacencies, path, entry + 1);
  auto bary1     = tri_bary_coords(triangles, positions, tid2, pos1);
  auto flat_tid1 = unfold_face(
      triangles, positions, adjacencies, flat_tid0, tid0, k);
  k              = find(adjacencies[tid1], tid2);
  auto flat_tid2 = unfold_face(
      triangles, positions, adjacencies, flat_tid1, tid1, k);
  auto p_minus = interpolate_triangle(
      flat_tid0[0], flat_tid0[1], flat_tid0[2], vec2f{bary0.y, bary0.z});
  auto p_plus = interpolate_triangle(
      flat_tid2[0], flat_tid2[1], flat_tid2[2], vec2f{bary1.y, bary1.z});
  auto v    = p_plus - p_minus;
  auto bary = vector_bary_coords(flat_tid1, v);
  return bary.y *
             (positions[triangles[tid1].y] - positions[triangles[tid1].x]) +
         bary.z * (positions[triangles[tid1].z] - positions[triangles[tid1].x]);
}
// choose the optimal seed when computing the weighted averages of three
// control_points
int optimal_seed(const float& t) {
  if (t < 0.333)
    return 0;
  else if (t < 0.666)
    return 1;
  else
    return 2;
}
// general version
float length_of_vector(const vector<float>& v) {
  auto len = 0.f;
  for (auto i = 0; i < v.size(); ++i) {
    len += v[i] * v[i];
  }
  return yocto::sqrt(len);
}
vector<float> subtract(const vector<float>& v, const vector<float>& w) {
  assert(v.size() == w.size());
  vector<float> sub(v.size());
  for (auto i = 0; i < v.size(); ++i) {
    sub[i] = v[i] - w[i];
  }
  return sub;
}
mesh_point optimal_seed(
    const vector<mesh_point>& control_points, const vector<float>& weights) {
  vector<pair<float, int>> ordering(weights.size());
  for (auto i = 0; i < weights.size(); ++i) {
    ordering[i] = {weights[i], i};
  }
  sort(ordering.begin(), ordering.end());

  return control_points[ordering.back().second];
}
mesh_point optimal_seed(const vector<mesh_point>& points,
    const vector<float>& point_weights, const vector<vector<float>>& weights) {
  auto seed_entry = 0;
  auto lambda     = flt_max;
  for (auto i = 0; i < points.size(); ++i) {
    auto magnitude = length_of_vector(subtract(weights[i], point_weights));
    if (magnitude < lambda) {
      lambda     = magnitude;
      seed_entry = i;
    }
  }
  return points[seed_entry];
}
std::pair<mesh_point, vector<float>> optimal_seed_with_weights(
    const vector<mesh_point>& points, const vector<float>& point_weights,
    const vector<vector<float>>& weights) {
  auto seed_entry = 0;
  auto lambda     = flt_max;
  for (auto i = 0; i < points.size(); ++i) {
    auto magnitude = length_of_vector(subtract(weights[i], point_weights));
    if (magnitude < lambda) {
      lambda     = magnitude;
      seed_entry = i;
    }
  }
  return {points[seed_entry], weights[seed_entry]};
}
vec2f get_cartesian_coordinates(const vector<float>& weights) {
  vec2f coords = zero2f;
  if (weights[0] + weights[2] != 0) {
    coords.x = weights[2] / (weights[0] + weights[2]);
  } else {
    coords.x = weights[3] / (weights[1] + weights[3]);
  }

  if (weights[1] + weights[0] != 0) {
    coords.y = weights[1] / (weights[0] + weights[1]);
  } else {
    coords.y = weights[3] / (weights[2] + weights[3]);
  }
  return coords;
}
vector<vector<vec3f>> gradients_inside_cell(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<vector<float>>&          fields,
    const vector<mesh_point>& control_points, const mesh_point& new_source,
    const vector<float>& source_weights, vector<float>& target_weights) {
  mesh_point            p0 = {}, p1 = {}, p2 = {};
  float                 w01 = 0.f, w02 = 0.f;
  vector<float>         f;
  vector<vector<vec3f>> new_grads(control_points.size());
  auto                  coords = get_cartesian_coordinates(source_weights);
  auto target_coords           = get_cartesian_coordinates(target_weights);
  auto uv                      = zero2f;
  if (target_weights[0] > source_weights[0]) {
    p0  = control_points[0];
    p1  = control_points[1];
    p2  = control_points[2];
    w01 = 1 - coords.y;
    w02 = 1 - coords.x;
    f   = fields[0];
    uv  = {target_coords.x / coords.x, target_coords.y / coords.y};
    target_weights[2] = (1 - uv.x) * (1 - uv.y);
    target_weights[1] = uv.y * (1 - uv.x);
    target_weights[3] = uv.x * (1 - uv.y);
    target_weights[0] = uv.x * uv.y;

  } else if (target_weights[1] > source_weights[1]) {
    p0  = control_points[1];
    p1  = control_points[0];
    p2  = control_points[3];
    w01 = coords.y;
    w02 = coords.x;
    f   = fields[1];
    uv  = {target_coords.x / coords.x, (target_coords.y - coords.y) / coords.y};
    target_weights[1] = (1 - uv.x) * (1 - uv.y);
    target_weights[2] = uv.y * (1 - uv.x);
    target_weights[0] = uv.x * (1 - uv.y);
    target_weights[3] = uv.x * uv.y;
  } else if (target_weights[2] > source_weights[2]) {
    p0  = control_points[2];
    p1  = control_points[0];
    p2  = control_points[3];
    f   = fields[2];
    w01 = coords.x;
    w02 = 1 - coords.y;
    uv  = {(target_coords.x - coords.x) / coords.x, target_coords.y / coords.y};
    target_weights[1] = (1 - uv.x) * (1 - uv.y);
    target_weights[0] = uv.y * (1 - uv.x);
    target_weights[2] = uv.x * (1 - uv.y);
    target_weights[3] = uv.x * uv.y;
  } else {
    p0                = control_points[3];
    p1                = control_points[1];
    p2                = control_points[2];
    f                 = fields[3];
    w01               = coords.x;
    w02               = coords.y;
    uv                = {(target_coords.x - coords.x) / coords.x,
        (target_coords.y - coords.y) / coords.y};
    target_weights[0] = (1 - uv.x) * (1 - uv.y);
    target_weights[1] = uv.y * (1 - uv.x);
    target_weights[3] = uv.x * (1 - uv.y);
    target_weights[2] = uv.x * uv.y;
  }

  vector<int> parents;
  auto        strip = get_strip_having_distances(solver, triangles, positions,
             adjacencies, v2t, angles, f, p0, p1, parents);
  auto gamma  = shortest_path(triangles, positions, adjacencies, p1, p0, strip);
  auto gammat = path_parameters(gamma, triangles, positions, adjacencies);
  auto q01    = sample_on_curve(
         triangles, positions, adjacencies, gamma, gammat, w01)
                 .first;

  strip  = get_strip_having_distances(solver, triangles, positions, adjacencies,
       v2t, angles, f, p0, p2, parents);
  gamma  = shortest_path(triangles, positions, adjacencies, p2, p0, strip);
  gammat = path_parameters(gamma, triangles, positions, adjacencies);
  auto q02 = sample_on_curve(
      triangles, positions, adjacencies, gamma, gammat, w02)
                 .first;

  vector<mesh_point> new_vertices{new_source, q01, p0, q02};
  auto               curr_field = compute_distance_field(
                    triangles, positions, adjacencies, v2t, solver, 0, new_vertices, graph);
  std::transform(curr_field.begin(), curr_field.end(), curr_field.begin(),
      [](float lambda) { return lambda * lambda; });
  new_grads[0] = compute_grad(
      solver, triangles, positions, normals, Grad, curr_field);
  for (auto i = 1; i < new_vertices.size(); ++i) {
    if (i % 2) {
      curr_field = compute_distance_field(triangles, positions, adjacencies,
          v2t, solver, 0, new_vertices, graph);
      std::transform(curr_field.begin(), curr_field.end(), curr_field.begin(),
          [](float lambda) { return lambda * lambda; });
      new_grads[i] = compute_grad(
          solver, triangles, positions, normals, Grad, curr_field);
    } else {
      std::transform(f.begin(), f.end(), f.begin(),
          [](float lambda) { return lambda * lambda; });
      new_grads[i] = compute_grad(
          solver, triangles, positions, normals, Grad, f);
    }
  }

  return new_grads;
}

vector<mesh_point> split_control_poligon(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& f,
    const vector<vector<vec3f>>& gradients, const float& t,
    const mesh_point& knot, const vector<mesh_point>& control_points,
    bool& decast) {
  vector<pair<int, float>>      verts{};
  vector<int>                   tids{};
  vector<mesh_point>            steps{};
  vector<vec3f>                 directions = {};
  vector<pair<mesh_point, int>> badones{};
  decast                 = false;
  auto               p0  = control_points[0];
  auto               p1  = control_points[1];
  auto               p2  = control_points[2];
  auto               p3  = control_points[3];
  mesh_point         p12 = {-1, zero2f};
  vector<mesh_point> new_control_points(7);
  vector<float>      weights(4);
  new_control_points[0] = p0;
  new_control_points[6] = p3;
  // compute b1 as weighted average of p0,p1,p2
  auto entry_of_seed = optimal_seed(t);
  auto seed          = triangles[control_points[entry_of_seed].face].x;
  // we should take the vertex with greatest barycentric coordinate
  weights[0]            = pow(1 - t, 2);
  weights[1]            = 2 * t * (1 - t);
  weights[2]            = pow(t, 2);
  weights[3]            = 0.f;
  auto start            = make_mesh_point(triangles, v2t, seed);
  new_control_points[2] = gradient_descent_absil(triangles, positions,
      adjacencies, v2t, normals, solver, angles, total_angles, gradients,
      weights, start, tids, directions, steps, badones);
  /*   weighted_average_test(dual_solver, solver, angles, total_angles,
        triangles, positions, normals, adjacencies, v2t, f, gradients,
        weights, control_points, verts, tids, directions, badones, seed); */
  if (badones.size() != 0) {
    /* strip = strip_on_dual_graph(
        dual_solver, triangles, positions, p2.face, p1.face);
    auto l12 = shortest_path(triangles, positions, adjacencies, p1, p2,
    strip); auto l12_t = path_parameters(l12, triangles, positions,
    adjacencies); p12 = sample_on_curve(triangles, positions, adjacencies,
    l12, l12_t, t).first; new_control_points[2] =
    sample_on_curve(dual_solver, triangles, positions, adjacencies,
    new_control_points[1], p12, t) .first; */
    return {};
  }
  badones.clear();
  // compute b2 as weighted average of p1,p2,p3
  seed                  = triangles[control_points[entry_of_seed + 1].face].x;
  start                 = make_mesh_point(triangles, v2t, seed);
  weights[0]            = 0.f;
  weights[1]            = pow(1 - t, 2);
  weights[2]            = 2 * t * (1 - t);
  weights[3]            = pow(t, 2);
  new_control_points[4] = gradient_descent_absil(triangles, positions,
      adjacencies, v2t, normals, solver, angles, total_angles, gradients,
      weights, start, tids, directions, steps, badones);
  /* weighted_average_test(dual_solver, solver, angles, total_angles,
      triangles, positions, normals, adjacencies, v2t, f, gradients,
      weights, control_points, verts, tids, directions, badones, seed); */
  if (badones.size() != 0) {
    /* if (p12.face == -1) {
      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, p2.face, p1.face);
      auto l12 = shortest_path(
          triangles, positions, adjacencies, p1, p2, strip);
      auto l12_t = path_parameters(l12, triangles, positions, adjacencies);
      p12 = sample_on_curve(triangles, positions, adjacencies, l12, l12_t,
    t) .first;
    }
    new_control_points[4] = sample_on_curve(dual_solver, triangles,
    positions, adjacencies, p12, new_control_points[5], t) .first; */
    return {};
  }
  // compute b0
  auto strip = strip_on_dual_graph(
      dual_solver, triangles, positions, p1.face, p0.face);
  auto l01   = shortest_path(triangles, positions, adjacencies, p0, p1, strip);
  auto l01_t = path_parameters(l01, triangles, positions, adjacencies);
  new_control_points[1] =
      sample_on_curve(triangles, positions, adjacencies, l01, l01_t, t).first;
  // compute b3
  strip = strip_on_dual_graph(
      dual_solver, triangles, positions, p3.face, p2.face);
  auto l23   = shortest_path(triangles, positions, adjacencies, p2, p3, strip);
  auto l23_t = path_parameters(l23, triangles, positions, adjacencies);
  new_control_points[5] =
      sample_on_curve(triangles, positions, adjacencies, l23, l23_t, t).first;

  // check C^1 continuty at knot
  strip               = strip_on_dual_graph(dual_solver, triangles, positions,
                    new_control_points[4].face, new_control_points[2].face);
  auto l              = shortest_path(triangles, positions, adjacencies,
                   new_control_points[2], new_control_points[4], strip);
  auto l_t            = path_parameters(l, triangles, positions, adjacencies);
  auto candidate_knot = sample_on_curve(
      triangles, positions, adjacencies, l, l_t, t);
  if (candidate_knot.first.face == knot.face) {
    auto tangent_vector = normalize(compute_tangent_vector(
        triangles, positions, adjacencies, l, candidate_knot.second));
    auto v              = normalize(
                     eval_position(triangles, positions, knot) -
                     eval_position(triangles, positions, candidate_knot.first));
    auto teta = yocto::acos(dot(tangent_vector, v));
    if (teta <= 1e-2) {
      new_control_points[3] = knot;
      decast                = false;
    } else
      new_control_points[3] = candidate_knot.first;
  } else
    new_control_points[3] = candidate_knot.first;

  // new_control_points[3] = knot;

  return new_control_points;
}

std::tuple<vector<mesh_point>, vector<vector<float>>, vector<vector<vec3f>>>
split_control_poligon(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const Eigen::SparseMatrix<double, 1>& Grad,
    const bezier_splits& curve, const std::unordered_map<int, int>& subdivision,
    const int entry) {
  vector<pair<int, float>>      verts{};
  vector<int>                   tids{};
  vector<mesh_point>            steps{};
  vector<vec3f>                 directions = {};
  vector<pair<mesh_point, int>> badones    = {{{-1, zero2f}, -1}};
  vector<int>                   strip;
  auto                          p0     = curve.control_polygons[entry][0];
  auto                          p1     = curve.control_polygons[entry][1];
  auto                          p2     = curve.control_polygons[entry][2];
  auto                          p3     = curve.control_polygons[entry][3];
  auto                          offset = subdivision.at(entry);
  mesh_point                    p12    = {-1, zero2f};
  vector<mesh_point>            new_control_points(5);
  vector<float>                 weights(4);
  vector<float>                 f;
  vector<vector<float>>         scalar_field(5);
  vector<vector<vec3f>>         gradient_field(5);
  vector<int>                   parents;
  // compute b0
  strip = get_strip_having_distances(solver, triangles, positions, adjacencies,
      v2t, angles, curve.distance_fields[offset][1], p1, p0, parents);

  auto l01   = shortest_path(triangles, positions, adjacencies, p0, p1, strip);
  auto l01_t = path_parameters(l01, triangles, positions, adjacencies);
  new_control_points[0] =
      sample_on_curve(triangles, positions, adjacencies, l01, l01_t, 0.5).first;
  f = compute_geodesic_distances(
      solver, triangles, positions, adjacencies, {new_control_points[0]});
  scalar_field[0] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  auto F            = wrapper(f);
  gradient_field[0] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // compute b3
  strip = get_strip_having_distances(solver, triangles, positions, adjacencies,
      v2t, angles, curve.distance_fields[offset][3], p3, p2, parents);
  auto l23   = shortest_path(triangles, positions, adjacencies, p2, p3, strip);
  auto l23_t = path_parameters(l23, triangles, positions, adjacencies);
  new_control_points[4] =
      sample_on_curve(triangles, positions, adjacencies, l23, l23_t, 0.5).first;
  f = compute_geodesic_distances(
      solver, triangles, positions, adjacencies, {new_control_points[4]});
  scalar_field[4] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[4] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // compute b1 as weighted average of p0,p1,p2

  auto seed = triangles[p1.face].x;
  // we should take the vertex with greatest barycentric coordinate
  weights[0] = 1 / 4.f;
  weights[1] = 0.5;
  weights[2] = 1 / 4.f;
  weights[3] = 0.f;
  auto start = make_mesh_point(triangles, v2t, seed);
  // new_control_points[1] = gradient_descent_absil(triangles, positions,
  //     adjacencies, v2t, normals, solver, angles, total_angles,
  //     curve.gradient_fields[offset], weights, start, tids, directions,
  //     steps, badones);
  if (badones.size() != 0) {
    // info("happens");
    auto Q = sample_on_curve(
        triangles, positions, adjacencies, l01, l01_t, 0.66)
                 .first;
    strip      = get_strip_having_distances(solver, triangles, positions,
             adjacencies, v2t, angles, curve.distance_fields[offset][2], p2, Q,
             parents);
    auto l12   = shortest_path(triangles, positions, adjacencies, Q, p2, strip);
    auto l12_t = path_parameters(l12, triangles, positions, adjacencies);
    // p12 = sample_on_curve(triangles, positions, adjacencies, l12, l12_t,
    // 0.5)
    //           .first;
    // strip = get_strip_having_distances(solver, triangles, positions,
    //     adjacencies, v2t, angles, scalar_field[0], new_control_points[0],
    //     p12, parents);
    // l12   = shortest_path(
    //     triangles, positions, adjacencies, p12, new_control_points[0],
    //     strip);
    // l12_t = path_parameters(l12, triangles, positions, adjacencies);
    new_control_points[1] = sample_on_curve(
        triangles, positions, adjacencies, l12, l12_t, 0.25)
                                .first;
  }
  f = compute_geodesic_distances(
      solver, triangles, positions, adjacencies, {new_control_points[1]});
  scalar_field[1] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[1] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // badones.clear();
  // compute b2 as weighted average of p1,p2,p3
  seed       = triangles[p2.face].x;
  start      = make_mesh_point(triangles, v2t, seed);
  weights[0] = 0.f;
  weights[1] = 1 / 4.f;
  weights[2] = 0.5;
  weights[3] = 1 / 4.f;
  // new_control_points[3] = gradient_descent_absil(triangles, positions,
  //     adjacencies, v2t, normals, solver, angles, total_angles,
  //     curve.gradient_fields[offset], weights, start, tids, directions,
  //     steps, badones);

  if (badones.size() != 0) {
    // if (p12.face == -1) {
    //   strip    = get_strip_having_distances(solver, triangles, positions,
    //       adjacencies, v2t, angles, curve.distance_fields[offset][2], p2,
    //       p1, parents);
    //   auto l12 = shortest_path(
    //       triangles, positions, adjacencies, p1, p2, strip);
    //   auto l12_t = path_parameters(l12, triangles, positions,
    //   adjacencies);

    //   p12 = sample_on_curve(triangles, positions, adjacencies, l12,
    //   l12_t, 0.5)
    //             .first;
    // }
    // info("happens");
    strip    = get_strip_having_distances(solver, triangles, positions,
           adjacencies, v2t, angles, curve.distance_fields[offset][2], p2, p1,
           parents);
    auto l12 = shortest_path(triangles, positions, adjacencies, p1, p2, strip);
    auto l12_t = path_parameters(l12, triangles, positions, adjacencies);
    auto Q     = sample_on_curve(
            triangles, positions, adjacencies, l12, l12_t, 0.33)
                 .first;
    // strip    = get_strip_having_distances(solver, triangles, positions,
    //     adjacencies, v2t, angles, scalar_field[4], new_control_points[4],
    //     p12, parents);
    // auto l23 = shortest_path(
    //     triangles, positions, adjacencies, p12, new_control_points[4],
    //     strip);
    // auto l23_t = path_parameters(l23, triangles, positions, adjacencies);
    strip      = get_strip_having_distances(solver, triangles, positions,
             adjacencies, v2t, angles, curve.distance_fields[offset][3], p3, Q,
             parents);
    auto l23   = shortest_path(triangles, positions, adjacencies, Q, p3, strip);
    auto l23_t = path_parameters(l23, triangles, positions, adjacencies);
    new_control_points[3] = sample_on_curve(
        triangles, positions, adjacencies, l23, l23_t, 0.75)
                                .first;
  }
  f = compute_geodesic_distances(
      solver, triangles, positions, adjacencies, {new_control_points[3]});
  scalar_field[3] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[3] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // check C^1 continuty at knot
  strip = get_strip_having_distances(solver, triangles, positions, adjacencies,
      v2t, angles, scalar_field[3], new_control_points[3],
      new_control_points[1], parents);

  auto l   = shortest_path(triangles, positions, adjacencies,
        new_control_points[1], new_control_points[3], strip);
  auto l_t = path_parameters(l, triangles, positions, adjacencies);
  new_control_points[2] =
      sample_on_curve(triangles, positions, adjacencies, l, l_t, 0.5).first;
  f = compute_geodesic_distances(
      solver, triangles, positions, adjacencies, {new_control_points[2]});
  scalar_field[2] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[2] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  return {new_control_points, scalar_field, gradient_field};
}
vector<mesh_point> split_control_poligon(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<vec3f>>& gradients,
    const vector<mesh_point>& control_points) {
  vector<pair<int, float>>      verts{};
  vector<int>                   tids{};
  vector<mesh_point>            steps{};
  vector<vec3f>                 directions = {};
  vector<pair<mesh_point, int>> badones{};
  vector<int>                   strip;
  auto                          p0  = control_points[0];
  auto                          p1  = control_points[1];
  auto                          p2  = control_points[2];
  auto                          p3  = control_points[3];
  mesh_point                    p12 = {-1, zero2f};
  vector<mesh_point>            new_control_points(7);
  vector<float>                 weights(4);
  new_control_points[0] = p0;
  new_control_points[6] = p3;

  // compute b0
  strip = strip_on_dual_graph(
      dual_solver, triangles, positions, p1.face, p0.face);
  auto l01   = shortest_path(triangles, positions, adjacencies, p0, p1, strip);
  auto l01_t = path_parameters(l01, triangles, positions, adjacencies);
  new_control_points[1] =
      sample_on_curve(triangles, positions, adjacencies, l01, l01_t, 0.5).first;
  // compute b3
  strip = strip_on_dual_graph(
      dual_solver, triangles, positions, p3.face, p2.face);
  auto l23   = shortest_path(triangles, positions, adjacencies, p2, p3, strip);
  auto l23_t = path_parameters(l23, triangles, positions, adjacencies);
  new_control_points[5] =
      sample_on_curve(triangles, positions, adjacencies, l23, l23_t, 0.5).first;
  // compute b1 as weighted average of p0,p1,p2

  auto seed = triangles[control_points[1].face].x;
  // we should take the vertex with greatest barycentric coordinate
  weights[0]            = 1 / 4.f;
  weights[1]            = 0.5;
  weights[2]            = 1 / 4.f;
  weights[3]            = 0.f;
  auto start            = make_mesh_point(triangles, v2t, seed);
  new_control_points[2] = gradient_descent_absil(triangles, positions,
      adjacencies, v2t, normals, solver, angles, total_angles, gradients,
      weights, start, tids, directions, steps, badones);
  /*   weighted_average_test(dual_solver, solver, angles, total_angles,
        triangles, positions, normals, adjacencies, v2t, f, gradients,
        weights, control_points, verts, tids, directions, badones, seed); */
  if (badones.size() != 0) {
    strip = strip_on_dual_graph(
        dual_solver, triangles, positions, p2.face, p1.face);
    auto l12 = shortest_path(triangles, positions, adjacencies, p1, p2, strip);
    auto l12_t = path_parameters(l12, triangles, positions, adjacencies);
    p12 = sample_on_curve(triangles, positions, adjacencies, l12, l12_t, 0.5)
              .first;
    new_control_points[2] = sample_on_curve(dual_solver, triangles, positions,
        adjacencies, new_control_points[1], p12, 0.5)
                                .first;
  }
  badones.clear();
  // compute b2 as weighted average of p1,p2,p3
  seed                  = triangles[control_points[2].face].x;
  start                 = make_mesh_point(triangles, v2t, seed);
  weights[0]            = 0.f;
  weights[1]            = 1 / 4.f;
  weights[2]            = 0.5;
  weights[3]            = 1 / 4.f;
  new_control_points[4] = gradient_descent_absil(triangles, positions,
      adjacencies, v2t, normals, solver, angles, total_angles, gradients,
      weights, start, tids, directions, steps, badones);
  /* weighted_average_test(dual_solver, solver, angles, total_angles,
      triangles, positions, normals, adjacencies, v2t, f, gradients,
      weights, control_points, verts, tids, directions, badones, seed); */
  if (badones.size() != 0) {
    if (p12.face == -1) {
      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, p2.face, p1.face);
      auto l12 = shortest_path(
          triangles, positions, adjacencies, p1, p2, strip);
      auto l12_t = path_parameters(l12, triangles, positions, adjacencies);
      p12 = sample_on_curve(triangles, positions, adjacencies, l12, l12_t, 0.5)
                .first;
    }
    new_control_points[4] = sample_on_curve(dual_solver, triangles, positions,
        adjacencies, p12, new_control_points[5], 0.5)
                                .first;
  }

  // check C^1 continuty at knot
  strip               = strip_on_dual_graph(dual_solver, triangles, positions,
                    new_control_points[4].face, new_control_points[2].face);
  auto l              = shortest_path(triangles, positions, adjacencies,
                   new_control_points[2], new_control_points[4], strip);
  auto l_t            = path_parameters(l, triangles, positions, adjacencies);
  auto candidate_knot = sample_on_curve(
      triangles, positions, adjacencies, l, l_t, 0.5);
  new_control_points[3] = candidate_knot.first;

  return new_control_points;
}

vector<mesh_point> bezier_karcher_almost_gradients(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const int number_of_subdivision, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<pair<mesh_point, float>>& badones) {
  int n = control_points.size();

  vector<mesh_point> poly_line /*=   {control_points[0]} */;
  verts.clear();
  vector<float>         w;
  vector<vector<float>> f(n);

  for (auto i = 0; i < n; ++i) {
    f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
        solver, i, control_points, exact);

    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }

  double step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);
  auto bary     = get_bary(control_points[0].uv);
  auto max_bary = flt_min;
  auto entry    = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] > max_bary) {
      max_bary = bary[i];
      entry    = i;
    }
  }
  int seed = triangles[control_points[0].face][entry];
  tids.resize(pow(2, number_of_subdivision) + 1);
  auto t           = 0.f;
  auto prev_sample = control_points[0];
  for (int i = 0; i <= pow(2, number_of_subdivision); ++i) {
    auto curr_tids = vector<int>{};
    auto size      = badones.size();
    auto curr =
        /* gradient_descent_absil(triangles, positions, adjacencies, v2t,
            normals, solver, angles, total_angles, gradients, w,
           prev_sample, badones); */
        weighted_average_almost_gradients(dual_solver, solver, angles,
            total_angles, triangles, positions, normals, adjacencies, v2t, f, w,
            control_points, verts, curr_tids, badones, seed);
    tids[i]             = curr_tids;
    verts.back().second = t;
    t += step;
    bernstein_polynomials(n, t, w);
    if (badones.size() != size) {
      // badones.back().first  = poly_line.back();
      badones.back().second = t - step;
      /* size                  = badones.size();
      curr = gradient_descent_test(triangles, positions, adjacencies, v2t,
          normals, solver, dual_solver, angles, total_angles,
      control_points, gradients, w, curr, curr_tids, badones); */
    } else {
      poly_line.push_back(curr);
      prev_sample = curr;
    }
  }

  return poly_line;
}
vector<mesh_point> bezier_karcher_test(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<vector<vec3f>>& directions,
    vector<vector<mesh_point>>& steps, vector<pair<mesh_point, int>>& badones) {
  time_function();
  int n = control_points.size();

  vector<mesh_point> poly_line = {control_points[0]};
  verts.clear();
  vector<float>         w;
  vector<vector<vec3f>> gradients(n);
  vector<vector<float>> f(n);

  for (auto i = 0; i < n; ++i) {
    f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
        solver, i, control_points, type_of_field);

    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }

  for (auto i = 0; i < n; ++i) {
    auto F       = wrapper(f[i]);
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }
  double step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);
  auto bary     = get_bary(control_points[0].uv);
  auto max_bary = flt_min;
  auto entry    = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] > max_bary) {
      max_bary = bary[i];
      entry    = i;
    }
  }
  int seed = triangles[control_points[0].face][entry];

  tids.resize(pow(2, number_of_subdivision) + 1);
  directions.resize(pow(2, number_of_subdivision) + 1);
  steps.resize(pow(2, number_of_subdivision) + 1);
  auto t           = 0.f;
  auto prev_sample = control_points[0];
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // show(i);
    auto curr_tids  = vector<int>{};
    auto curr_dir   = vector<vec3f>{};
    auto curr_steps = vector<mesh_point>{};
    auto size       = badones.size();
    /*   auto curr = weighted_average_test(dual_solver, solver, angles,
       total_angles, triangles, positions, normals, adjacencies, v2t, f,
       gradients, w,
          control_points, verts, curr_tids, curr_dir, badones, seed); */
    auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
        normals, solver, angles, total_angles, gradients, w, prev_sample,
        curr_tids, curr_dir, curr_steps, badones);
    tids[i]   = curr_tids;
    directions[i] = curr_dir;
    steps[i]      = curr_steps;
    verts.push_back({-1, t});
    t += step;
    bernstein_polynomials(n, t, w);
    if (badones.size() != size) {
      badones.back().first = poly_line.back();

      badones.back().second = i - 1;

    } else {
      poly_line.push_back(curr);
      prev_sample = curr;
    }
  }
  poly_line.push_back(control_points.back());
  return poly_line;
}
vector<mesh_point> karcher_test(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    vector<int>& jumps) {
  jumps.clear();
  int                           n = control_points.size();
  vector<pair<int, float>>      verts;
  vector<vector<int>>           tids;
  vector<vector<vec3f>>         directions;
  vector<vector<mesh_point>>    steps;
  vector<pair<mesh_point, int>> badones;
  vector<mesh_point>            poly_line = {control_points[0]};
  verts.clear();
  vector<float>         w;
  vector<vector<vec3f>> gradients(n);
  vector<vector<float>> f(n);

  for (auto i = 0; i < n; ++i) {
    f[i] = exact_geodesic_distance(
        triangles, positions, v2t, control_points[i]);

    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }

  for (auto i = 0; i < n; ++i) {
    auto F       = wrapper(f[i]);
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }
  double step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);

  tids.resize(pow(2, number_of_subdivision) + 1);
  directions.resize(pow(2, number_of_subdivision) + 1);
  steps.resize(pow(2, number_of_subdivision) + 1);
  auto t           = 0.f;
  auto prev_sample = control_points[0];

  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // show(i);
    auto curr_tids  = vector<int>{};
    auto curr_dir   = vector<vec3f>{};
    auto curr_steps = vector<mesh_point>{};
    auto size       = badones.size();
    auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
        normals, solver, angles, total_angles, gradients, w, prev_sample,
        curr_tids, curr_dir, curr_steps, badones);

    t += step;
    bernstein_polynomials(n, t, w);
    if (badones.size() == size) {
      poly_line.push_back(curr);
      prev_sample = curr;
    } else {
      jumps.push_back(i);
      poly_line.push_back(badones.back().first);
    }
  }

  return poly_line;
}
vector<mesh_point> karcher_test(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<int>& control_points, const int number_of_subdivision) {
  int                           n = control_points.size();
  vector<pair<int, float>>      verts;
  vector<vector<int>>           tids;
  vector<vector<vec3f>>         directions;
  vector<vector<mesh_point>>    steps;
  vector<pair<mesh_point, int>> badones;
  vector<mesh_point>            poly_line = {
      point_from_vert(triangles, v2t, control_points[0])};
  verts.clear();
  vector<float> w;
  double        step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);

  tids.resize(pow(2, number_of_subdivision) + 1);
  directions.resize(pow(2, number_of_subdivision) + 1);
  steps.resize(pow(2, number_of_subdivision) + 1);
  auto t           = 0.f;
  auto prev_sample = point_from_vert(triangles, v2t, control_points[0]);

  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // show(i);
    auto curr_tids  = vector<int>{};
    auto curr_dir   = vector<vec3f>{};
    auto curr_steps = vector<mesh_point>{};
    auto size       = badones.size();
    auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
        normals, solver, angles, total_angles, gradients, w, prev_sample,
        curr_tids, curr_dir, curr_steps, badones);

    if (badones.size() == size) {
      poly_line.push_back(curr);
      prev_sample = curr;
    } else {
      poly_line.push_back(badones.back().first);
    }
    t += step;
    bernstein_polynomials(n, t, w);
  }

  return poly_line;
}
SurfacePoint exp_map(const std::unique_ptr<SurfaceMesh>& mesh,
    const std::unique_ptr<VertexPositionGeometry>&       geometry,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>&         adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const SurfacePoint& start, const Vector2& dir, const float& len) {
  auto    basisx    = geometry->faceTangentBasis[start.face];
  Vector3 direction = basisx[0] * dir.x + basisx[1] * dir.y;
  auto    yocto_dir = vec3f{
      (float)direction.x, (float)direction.y, (float)direction.z};
  auto tid      = (int)start.face.getIndex();
  auto bary     = start.faceCoords;
  auto p        = mesh_point{(int)tid, vec2f{(float)bary.y, (float)bary.z}};
  auto path     = straightest_geodesic(solver, triangles, positions, normals,
          adjacencies, v2p_adjacencies, angles, total_angles, p, yocto_dir, len);
  auto end      = path.back();
  auto bary_end = get_bary(end.uv);
  return SurfacePoint(
      mesh->face(end.face), Vector3{bary_end.x, bary_end.y, bary_end.z});
}

mesh_point weighted_average(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const float& t) {
  int                n         = control_points.size();
  vector<mesh_point> poly_line = {control_points[0]};
  vector<float>      w;
  bernstein_polynomials(n, t, w);
  vector<vector<vec3f>>         gradients(n);
  vector<vector<float>>         f(n);
  vector<pair<mesh_point, int>> badones;

  for (auto i = 0; i < n; ++i) {
    f[i] = exact_geodesic_distance(
        triangles, positions, v2t, control_points[i]);

    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }

  for (auto i = 0; i < n; ++i) {
    auto F       = wrapper(f[i]);
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }

  auto prev_sample = minimum_on_verts(triangles, v2t, f, w);

  auto curr_tids  = vector<int>{};
  auto curr_dir   = vector<vec3f>{};
  auto curr_steps = vector<mesh_point>{};

  auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
      normals, solver, angles, total_angles, gradients, w, prev_sample,
      curr_tids, curr_dir, curr_steps, badones);

  if (badones.size() != 0)
    info(
        "Control points are too far from each other, please change the configuration and retry");
  return curr;
}
// good one
vector<mesh_point> bezier_karcher_test(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, pair<mesh_point, int>& badone) {
  int n = control_points.size();

  vector<mesh_point>    poly_line = {control_points[0]};
  vector<float>         w;
  vector<vector<vec3f>> gradients(n);
  vector<vector<float>> f(n);

  for (auto i = 0; i < n; ++i) {
    f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
        solver, i, control_points, type_of_field);

    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }

  for (auto i = 0; i < n; ++i) {
    auto F       = wrapper(f[i]);
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }
  double step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);
  auto bary     = get_bary(control_points[0].uv);
  auto max_bary = flt_min;
  auto entry    = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] > max_bary) {
      max_bary = bary[i];
      entry    = i;
    }
  }
  int                           seed = triangles[control_points[0].face][entry];
  auto                          t    = 0.f;
  auto                          prev_sample = control_points[0];
  vector<pair<mesh_point, int>> badones;
  vector<pair<int, float>>      verts;
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // show(i);
    auto curr_tids  = vector<int>{};
    auto curr_dir   = vector<vec3f>{};
    auto curr_steps = vector<mesh_point>{};
    // auto curr = weighted_average_test(dual_solver, solver, angles,
    // total_angles,
    //     triangles, positions, normals, adjacencies, v2t, f, gradients, w,
    //     control_points, verts, curr_tids, curr_dir, badones, seed);
    auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
        normals, solver, angles, total_angles, gradients, w, prev_sample,
        curr_tids, curr_dir, curr_steps, badones);
    t += step;
    bernstein_polynomials(n, t, w);
    if (badones.size() != 0) {
      badone.first = poly_line.back();

      badone.second = i - 1;
      return poly_line;

    } else {
      poly_line.push_back(curr);
      prev_sample = curr;
    }
  }
  poly_line.push_back(control_points.back());
  return poly_line;
}
vector<mesh_point> bezier_karcher_test(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    pair<mesh_point, int>& badone) {
  int n = control_points.size();

  vector<mesh_point> poly_line = {control_points[0]};
  vector<float>      w;

  double step = 1 / pow(2, number_of_subdivision);
  bernstein_polynomials(n, 0.0, w);
  auto bary     = get_bary(control_points[0].uv);
  auto max_bary = flt_min;
  auto entry    = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] > max_bary) {
      max_bary = bary[i];
      entry    = i;
    }
  }
  int                           seed = triangles[control_points[0].face][entry];
  auto                          t    = 0.f;
  auto                          prev_sample = control_points[0];
  vector<pair<mesh_point, int>> badones;
  vector<pair<int, float>>      verts;
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // show(i);
    auto curr_tids  = vector<int>{};
    auto curr_dir   = vector<vec3f>{};
    auto curr_steps = vector<mesh_point>{};
    auto curr = gradient_descent_absil(triangles, positions, adjacencies, v2t,
        normals, solver, angles, total_angles, gradients, w, prev_sample,
        curr_tids, curr_dir, curr_steps, badones);
    t += step;
    bernstein_polynomials(n, t, w);
    if (badones.size() != 0) {
      badone.first = poly_line.back();

      badone.second = i - 1;
      return poly_line;

    } else {
      poly_line.push_back(curr);
      prev_sample = curr;
    }
  }
  poly_line.push_back(control_points.back());
  return poly_line;
}

// knot.x= entry along the subdivion knot.y=entry in the subdivided polyline
vector<mesh_point> bezier_karcher(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, splits& composite_bezier, const bool all_in_once,
    const vec2i& knot) {
  time_function();
  auto                          n = (int)control_points.size();
  vector<mesh_point>            polyline;
  vector<vector<vec3f>>         gradients(n);
  vector<vector<float>>         f(n);
  vector<pair<int, float>>      verts;
  vector<vector<int>>           tids;
  vector<vector<mesh_point>>    steps;
  vector<vector<vec3f>>         directions;
  vector<pair<mesh_point, int>> first_badones;
  vector<mesh_point>            curr_control_points;
  int                           entry    = 0;
  pair<mesh_point, int>         curr_bad = {{}, -1};

  if (composite_bezier.control_polygons.size() == 0) {
    // initialization
    if (all_in_once) {
      vector<vector<float>> f(n);

      for (auto i = 0; i < n; ++i) {
        f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
            solver, i, control_points, type_of_field);

        std::transform(f[i].begin(), f[i].end(), f[i].begin(),
            [](float lambda) { return lambda * lambda; });
        auto F       = wrapper(f[i]);
        gradients[i] = compute_grad(
            solver, triangles, positions, normals, Grad, F);
      }
      auto samples = bezier_karcher_test(dual_solver, solver, triangles,
          positions, adjacencies, normals, v2t, angles, total_angles, Grad,
          control_points, number_of_subdivision, type_of_field, curr_bad);
      composite_bezier.polyline.push_back(samples);
      composite_bezier.control_polygons.push_back(control_points);
      composite_bezier.badones.push_back({curr_bad});
      if (curr_bad.second == -1)
        composite_bezier.goodones.push_back(true);
      else
        composite_bezier.goodones.push_back(false);
    } else {
      auto samples = bezier_karcher_test(dual_solver, solver, triangles,
          positions, adjacencies, normals, v2t, angles, total_angles, Grad,
          control_points, number_of_subdivision, type_of_field, verts, tids,
          directions, steps, first_badones);
      composite_bezier.polyline.push_back(samples);
      composite_bezier.control_polygons.push_back(control_points);
      composite_bezier.badones.push_back(first_badones);
      if (first_badones.size() == 0)
        composite_bezier.goodones.push_back(true);
      else
        composite_bezier.goodones.push_back(false);
    }
  }
  if (all_in_once) {
    //   auto fine_curve = false;
    //   while (!fine_curve) {
    //     entry = -1;
    //     for (auto i = 0; i < composite_bezier.goodones.size(); ++i) {
    //       if (composite_bezier.goodones[i]) continue;

    //       curr_bad            = composite_bezier.badones[i][0];
    //       curr_control_points = composite_bezier.control_polygons[i];
    //       entry               = i;
    //       for (auto j = 0; j < 4; ++j) {
    //         f[j] = compute_distance_field(triangles, positions,
    //         adjacencies, v2t,
    //             solver, j, curr_control_points, type_of_field);

    //         std::transform(f[j].begin(), f[j].end(), f[j].begin(),
    //             [](float lambda) { return lambda * lambda; });
    //         auto F       = wrapper(f[j]);
    //         gradients[j] = compute_grad(
    //             solver, triangles, positions, normals, Grad, F);
    //       }

    //       break;
    //     }
    //     if (entry >= 0) {

    //       tids.clear();
    //       directions.clear();
    //       // verts.clear();
    //       auto [poly0, poly1, new_control_points, badones] =
    //           smoothing_bezier_curve(dual_solver, solver, triangles,
    //           positions,
    //               adjacencies, normals, v2t, angles, total_angles, Grad,
    //               gradients, curr_control_points, number_of_subdivision,
    //               type_of_field);
    //       if (poly0.size() == 0) {
    //         show(curr_bad.second);
    //         auto bad_entry = (curr_bad.second <= 1) ? 1 : curr_bad.second
    //         / 2; auto new_bad   = std::make_pair(
    //             composite_bezier.polyline[entry][bad_entry], bad_entry);
    //         composite_bezier.badones[entry].insert(
    //             composite_bezier.badones[entry].begin(), new_bad);
    //       } else {
    //         composite_bezier.control_polygons[entry] =
    //         {new_control_points[0],
    //             new_control_points[1], new_control_points[2],
    //             new_control_points[3]};
    //         composite_bezier.control_polygons.insert(
    //             composite_bezier.control_polygons.begin() + entry + 1,
    //             {new_control_points[3], new_control_points[4],
    //                 new_control_points[5], new_control_points[6]});

    //         composite_bezier.badones[entry] = {badones[0]};
    //         composite_bezier.badones.insert(
    //             composite_bezier.badones.begin() + entry + 1,
    //             {badones[1]});
    //         if (badones[0].second == -1)
    //           composite_bezier.goodones[entry] = true;
    //         else
    //           composite_bezier.goodones[entry] = false;

    //         if (badones[1].second == -1)
    //           composite_bezier.goodones.insert(
    //               composite_bezier.goodones.begin() + entry + 1, true);
    //         else
    //           composite_bezier.goodones.insert(
    //               composite_bezier.goodones.begin() + entry + 1, false);

    //         composite_bezier.polyline[entry] = poly0;
    //         composite_bezier.polyline.insert(
    //             composite_bezier.polyline.begin() + entry + 1, poly1);
    //       }
    //     } else
    //       fine_curve = true;
    //}
  } else {
    bool de_cast;
    auto curr_bad       = composite_bezier.badones[knot.x][knot.y];
    curr_control_points = composite_bezier.control_polygons[knot.x];
    for (auto j = 0; j < 4; ++j) {
      f[j] = compute_distance_field(triangles, positions, adjacencies, v2t,
          solver, j, curr_control_points, type_of_field);

      std::transform(f[j].begin(), f[j].end(), f[j].begin(),
          [](float lambda) { return lambda * lambda; });
      auto F       = wrapper(f[j]);
      gradients[j] = compute_grad(
          solver, triangles, positions, normals, Grad, F);
    }

    auto [poly0, poly1, new_control_points, badones] = smoothing_bezier_curve(
        dual_solver, solver, triangles, positions, adjacencies, normals, v2t,
        angles, total_angles, Grad, gradients, curr_control_points,
        number_of_subdivision, type_of_field);

    if (poly0.size() == 0) {
      info("you have to split earlier");
      return {};

    } else {
      composite_bezier.control_polygons[knot.x] = {new_control_points[0],
          new_control_points[1], new_control_points[2], new_control_points[3]};
      composite_bezier.control_polygons.insert(
          composite_bezier.control_polygons.begin() + knot.x + 1,
          {new_control_points[3], new_control_points[4], new_control_points[5],
              new_control_points[6]});

      composite_bezier.badones[knot.x] = badones[0];
      composite_bezier.badones.insert(
          composite_bezier.badones.begin() + knot.x + 1, badones[1]);
      if (badones[0].size() == 0)
        composite_bezier.goodones[knot.x] = true;
      else
        composite_bezier.goodones[knot.x] = false;

      if (badones[1].size() == 0)
        composite_bezier.goodones.insert(
            composite_bezier.goodones.begin() + knot.x + 1, true);
      else
        composite_bezier.goodones.insert(
            composite_bezier.goodones.begin() + knot.x + 1, false);

      composite_bezier.polyline[knot.x] = poly0;
      composite_bezier.polyline.insert(
          composite_bezier.polyline.begin() + knot.x + 1, poly1);
    }
  }

  for (auto i = 0; i < composite_bezier.polyline.size(); ++i) {
    polyline.insert(polyline.end(), composite_bezier.polyline[i].begin(),
        composite_bezier.polyline[i].end());
  }
  return polyline;
}
std::tuple<vector<mesh_point>, vector<float>, vector<int>>
sort_points_and_weigts(const vector<mesh_point>& control_points,
    const vector<float>& weights, const float& t) {
  if (t <= 0.25)
    return {control_points, weights, {0, 1, 2, 3}};
  else if (t <= (yocto::sqrt(3) - 1) / 2)
    return {{control_points[1], control_points[0], control_points[2],
                control_points[3]},
        {weights[1], weights[0], weights[2], weights[3]}, {1, 0, 2, 3}};
  else if (t <= 0.5)
    return {{control_points[1], control_points[2], control_points[0],
                control_points[3]},
        {weights[1], weights[2], weights[0], weights[3]}, {1, 2, 0, 3}};
  else if (t <= (3 - yocto::sqrt(3) / 2))
    return {{control_points[2], control_points[1], control_points[3],
                control_points[0]},
        {weights[2], weights[1], weights[3], weights[0]}, {2, 1, 3, 0}};
  else if (t <= 0.75)
    return {{control_points[2], control_points[3], control_points[1],
                control_points[0]},
        {weights[2], weights[3], weights[1], weights[0]}, {2, 3, 1, 0}};
  else
    return {{control_points[3], control_points[2], control_points[1],
                control_points[0]},
        {weights[3], weights[2], weights[1], weights[0]}, {3, 2, 1, 0}};
}
std::tuple<vector<mesh_point>, vector<float>, vector<int>>
sort_points_and_weigts_quadric(const vector<mesh_point>& control_points,
    const vector<float>& weights, const float& t) {
  if (t <= 1 / 3.f)
    return {control_points, weights, {0, 1, 2}};
  else if (t <= 0.5)
    return {{control_points[1], control_points[0], control_points[2]},
        {weights[1], weights[0], weights[2]}, {1, 0, 2}};
  else if (t <= 2 / 3.f)
    return {{control_points[1], control_points[2], control_points[0]},
        {weights[1], weights[2], weights[0]}, {1, 2, 0}};
  else
    return {{control_points[2], control_points[1], control_points[0]},
        {weights[2], weights[1], weights[0]}, {2, 1, 0}};
}
std::tuple<geodesic_path, vector<float>> set_gamma01(
    const vector<geodesic_path>& paths, const vector<vector<float>>& paths_t,
    const vector<int>& order) {
  switch (order[0]) {
    case 0: {
      return {paths[0], paths_t[0]};

    } break;
    case 1: {
      if (order[1] == 0)
        return {paths[0], paths_t[0]};
      else
        return {paths[1], paths_t[1]};

    } break;

    case 2: {
      if (order[1] == 1)
        return {paths[1], paths_t[1]};
      else
        return {paths[2], paths_t[2]};
    } break;

    case 3: {
      return {paths[2], paths_t[2]};
    } break;
  }
}
std::tuple<geodesic_path, vector<float>> set_gamma01_quadric(
    const vector<geodesic_path>& paths, const vector<vector<float>>& paths_t,
    const vector<int>& order) {
  switch (order[0]) {
    case 0: {
      return {paths[0], paths_t[0]};

    } break;
    case 1: {
      if (order[1] == 0)
        return {paths[0], paths_t[0]};
      else
        return {paths[1], paths_t[1]};

    } break;

    case 2: {
      return {paths[1], paths_t[1]};
    } break;
  }
}
mesh_point weighted_average_dyn_quadric(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<vector<float>>& field,
    const vector<geodesic_path>& control_polygon,
    const vector<vector<float>>& control_polygon_t,
    const vector<mesh_point>& control_points, const vector<float>& weights,
    const vector<int>& order) {
  auto        w1 = weights[1] / (1 - weights[2]);
  vector<int> parents;

  auto [gamma01, gamma01t] = set_gamma01_quadric(
      control_polygon, control_polygon_t, order);

  auto P01 = (order[0] < order[1]) ? sample_on_curve(triangles, positions,
                                         adjacencies, gamma01, gamma01t, w1)
                                         .first
                                   : sample_on_curve(triangles, positions,
                                         adjacencies, gamma01, gamma01t, 1 - w1)
                                         .first;

  auto strip = get_strip_having_distances(solver, triangles, positions,
      adjacencies, v2t, angles, field[order[2]], control_points[2], P01,
      parents);
  // get_strip(solver, triangles, positions, adjacencies, v2t, angles,
  //     control_points[2], P01, parents);
  auto gamma02 = shortest_path(
      triangles, positions, adjacencies, P01, control_points[2], strip);
  auto gamma02_t = path_parameters(gamma02, triangles, positions, adjacencies);

  return sample_on_curve(
      triangles, positions, adjacencies, gamma02, gamma02_t, weights[2])
      .first;
}
mesh_point weighted_average_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<vector<float>>& field,
    const vector<geodesic_path>& control_polygon,
    const vector<vector<float>>& control_polygon_t,
    const vector<mesh_point>& control_points, const vector<float>& weights,
    const vector<int>& order) {
  auto        w1 = weights[1] / (1 - weights[2] - weights[3]);
  auto        w2 = weights[2] / (1 - weights[3]);
  vector<int> parents;

  auto [gamma01, gamma01t] = set_gamma01(
      control_polygon, control_polygon_t, order);

  auto P01 = (order[0] < order[1]) ? sample_on_curve(triangles, positions,
                                         adjacencies, gamma01, gamma01t, w1)
                                         .first
                                   : sample_on_curve(triangles, positions,
                                         adjacencies, gamma01, gamma01t, 1 - w1)
                                         .first;

  auto strip = get_strip_having_distances(solver, triangles, positions,
      adjacencies, v2t, angles, field[order[2]], control_points[2], P01,
      parents);
  // strip = get_strip(solver, triangles, positions, adjacencies, v2t,
  // angles,
  //     control_points[2], P01, parents);
  auto gamma02 = shortest_path(
      triangles, positions, adjacencies, P01, control_points[2], strip);
  auto gamma02_t = path_parameters(gamma02, triangles, positions, adjacencies);
  auto P012      = sample_on_curve(
           triangles, positions, adjacencies, gamma02, gamma02_t, w2)
                  .first;
  strip = strip_on_dual_graph(
      dual_solver, triangles, positions, control_points[3].face, P012.face);
  // strip = get_strip_having_distances(solver, triangles, positions,
  // adjacencies,
  //     v2t, angles, field[order[3]], control_points[3], P012, parents);
  // strip   = get_strip(solver, triangles, positions, adjacencies, v2t,
  // angles,
  //     control_points[3], P012, parents);
  gamma02 = shortest_path(
      triangles, positions, adjacencies, P012, control_points[3], strip);
  gamma02_t = path_parameters(gamma02, triangles, positions, adjacencies);

  return sample_on_curve(
      triangles, positions, adjacencies, gamma02, gamma02_t, weights[3])
      .first;
}
mesh_point weighted_average_dyn(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles,
    const vector<mesh_point>& control_points, const vector<float>& weights) {
  vector<pair<float, int>> order(4);
  for (auto i = 0; i < 4; ++i) {
    order[i] = {weights[i], i};
  }

  sort(order.begin(), order.end(), std::greater<pair<float, int>>());
  auto w1 = order[1].first / (1 - order[2].first - order[3].first);
  auto w2 = order[2].first / (1 - order[3].first);

  auto        p0 = control_points[order[0].second];
  auto        p1 = control_points[order[1].second];
  auto        p2 = control_points[order[2].second];
  auto        p3 = control_points[order[3].second];
  vector<int> parents;
  auto        strip = get_strip(
             solver, triangles, positions, adjacencies, v2t, angles, p1, p0, parents);
  auto gamma  = shortest_path(triangles, positions, adjacencies, p0, p1, strip);
  auto gammat = path_parameters(gamma, triangles, positions, adjacencies);
  auto P01    = sample_on_curve(
         triangles, positions, adjacencies, gamma, gammat, w1)
                 .first;
  strip = get_strip(
      solver, triangles, positions, adjacencies, v2t, angles, p2, P01, parents);
  gamma     = shortest_path(triangles, positions, adjacencies, P01, p2, strip);
  gammat    = path_parameters(gamma, triangles, positions, adjacencies);
  auto P012 = sample_on_curve(
      triangles, positions, adjacencies, gamma, gammat, w2)
                  .first;
  strip  = get_strip(solver, triangles, positions, adjacencies, v2t, angles, p3,
       P012, parents);
  gamma  = shortest_path(triangles, positions, adjacencies, P012, p3, strip);
  gammat = path_parameters(gamma, triangles, positions, adjacencies);
  return sample_on_curve(
      triangles, positions, adjacencies, gamma, gammat, order[3].first)
      .first;
}
std::tuple<vector<pair<mesh_point, int>>, vector<vector<float>>>
b_spline_subdivision(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles,
    const vector<geodesic_path>& control_polygon,
    const vector<vector<float>>& control_polygon_t,
    const vector<mesh_point>&    control_points) {
  vector<pair<mesh_point, int>> new_points(5);
  vector<vector<float>>         fields(5);
  vector<int>                   parents;
  vector<float>                 f;
  // first edge
  new_points[0] = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[0], control_polygon_t[0], 0.5);
  // first vertex
  auto q01 = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[0], control_polygon_t[0], 0.75)
                 .first;
  auto q12 = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[1], control_polygon_t[1], 0.25)
                 .first;
  auto strip = get_strip(solver, triangles, positions, adjacencies, v2t, angles,
      q12, q01, parents);
  auto gamma = shortest_path(
      triangles, positions, adjacencies, q01, q12, strip);
  auto gammat   = path_parameters(gamma, triangles, positions, adjacencies);
  new_points[1] = sample_on_curve(
      triangles, positions, adjacencies, gamma, gammat, 0.5);
  // second edge
  new_points[2] = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[1], control_polygon_t[1], 0.5);
  // second vertex
  q01 = sample_on_curve(triangles, positions, adjacencies, control_polygon[1],
      control_polygon_t[1], 0.75)
            .first;
  q12 = sample_on_curve(triangles, positions, adjacencies, control_polygon[2],
      control_polygon_t[2], 0.25)
            .first;
  strip = get_strip(solver, triangles, positions, adjacencies, v2t, angles, q12,
      q01, parents);
  gamma = shortest_path(triangles, positions, adjacencies, q01, q12, strip);
  gammat        = path_parameters(gamma, triangles, positions, adjacencies);
  new_points[3] = sample_on_curve(
      triangles, positions, adjacencies, gamma, gammat, 0.5);
  // third edge
  new_points[4] = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[2], control_polygon_t[2], 0.5);

  vector<mesh_point> left  = {control_points[0], new_points[0].first,
      new_points[1].first, new_points[2].first};
  vector<mesh_point> right = {new_points[2].first, new_points[3].first,
      new_points[4].first, control_points[3]};

  for (auto i = 1; i < 3; ++i) {
    f = compute_distance_field(
        triangles, positions, adjacencies, v2t, solver, i, left, graph);
    fields[i - 1] = f;
  }
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 3,
      {control_points[0], new_points[0].first, new_points[1].first,
          new_points[2].first, new_points[3].first, new_points[4].first,
          control_points[3]},
      graph);
  fields[2] = f;
  for (auto i = 1; i < 3; ++i) {
    f = compute_distance_field(
        triangles, positions, adjacencies, v2t, solver, i, right, graph);
    fields[i + 2] = f;
  }

  return {new_points, fields};
}
std::tuple<std::pair<mesh_point, int>, vector<float>> b_spline_subdivision_even(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const geodesic_path& gamma, const vector<float>& gamma_t,
    const vector<mesh_point>& control_points) {
  auto p = sample_on_curve(
      triangles, positions, adjacencies, gamma, gamma_t, 0.5);
  auto f = compute_distance_field(triangles, positions, adjacencies, v2t,
      solver, 0, {p.first, control_points[0], control_points[1]}, graph);
  return {p, f};
}
std::tuple<std::pair<mesh_point, int>, vector<float>> b_spline_subdivision_odd(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<geodesic_path>& control_polygon,
    const vector<vector<float>>& control_polygon_t,
    const vector<mesh_point>&    control_points) {
  vector<int> parents;
  auto        q01 = sample_on_curve(triangles, positions, adjacencies,
             control_polygon[0], control_polygon_t[0], 0.75)
                 .first;
  auto q12 = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[1], control_polygon_t[1], 0.25)
                 .first;
  auto strip = get_strip(solver, triangles, positions, adjacencies, v2t, angles,
      q12, q01, parents);
  auto gamma = shortest_path(
      triangles, positions, adjacencies, q01, q12, strip);
  auto gammat = path_parameters(gamma, triangles, positions, adjacencies);
  auto p      = sample_on_curve(
           triangles, positions, adjacencies, gamma, gammat, 0.5);
  auto f = compute_distance_field(triangles, positions, adjacencies, v2t,
      solver, 0, {p.first, control_points[0], control_points[1]}, graph);
  return {p, f};
}
bool stopping_criterion(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& p0, const mesh_point& p1,
    const mesh_point& p2, const mesh_point& p3, double treshold) {
  auto A = eval_position(triangles, positions, p0);
  auto B = eval_position(triangles, positions, p1);
  auto C = eval_position(triangles, positions, p2);
  auto D = eval_position(triangles, positions, p3);

  auto ab = A - B;
  auto ac = C - B;
  auto cd = D - C;

  auto teta01 = angle(ab, ac);
  ac *= -1;
  auto teta12 = angle(ac, cd);

  auto Sur01 = yocto::sin(teta01) * length(ab) * length(ac) / 2;
  auto Sur12 = yocto::sin(teta12) * length(ac) * length(cd) / 2;

  if (Sur01 <= treshold && Sur12 <= treshold) return true;

  return false;
}
// tenere conto dell'entrata
vector<mesh_point> bezier_through_spline_subdivisions(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<mesh_point>& control_points) {
  // auto poly_line = vector<mesh_point>{};
  // struct parametric_path {
  //   geodesic_path path;
  //   vector<float> t;
  // };
  // struct control_polygon {
  //   vector<parametric_path> geodesics      = {};
  //   vector<mesh_point>      control_points = {};
  // };

  // auto        treshold = pow(10.0, (int)-3);
  // vector<int> parents;
  // auto        a = control_points[0];
  // auto        b = control_points[1];
  // auto        c = control_points[2];
  // auto        d = control_points[3];
  // auto L0 = parametric_path{}, L1 = parametric_path{}, L2 =
  // parametric_path{},
  //      L01 = parametric_path{}, L12 = parametric_path{}, L =
  //      parametric_path{};
  // auto strip = get_strip(
  //     solver, triangles, positions, adjacencies, v2t, angles, b, a, parents);
  // L0.path = shortest_path(triangles, positions, adjacencies, a, b, strip);
  // L0.t    = path_parameters(L0.path, triangles, positions, adjacencies);
  // strip   = get_strip(
  //     solver, triangles, positions, adjacencies, v2t, angles, c, b, parents);
  // L1.path = shortest_path(triangles, positions, adjacencies, b, c, strip);
  // L1.t    = path_parameters(L0.path, triangles, positions, adjacencies);
  // strip   = get_strip(
  //     solver, triangles, positions, adjacencies, v2t, angles, d, c, parents);
  // L2.path = shortest_path(triangles, positions, adjacencies, c, d, strip);
  // L2.t    = path_parameters(L2.path, triangles, positions, adjacencies);

  // // struct bezier_split {
  // //   vector<parametric_path> control_polygon;
  // // };

  // // auto s            = std::deque<bezier_split>();
  // // auto first_bezier = vector<parametric_path>{L0, L1, L2};

  // // auto first_one = bezier_split{first_bezier};
  // // s.push_back(first_one);
  // control_polygon curr_polygon = {{L0, L1, L2}, {a, b, c, d}};
  // poly_line.push_back(a);
  // auto count   = 0;
  // bool arrived = false;
  // while (count < 50) {
  //   ++count;
  //   // auto curr_bezier = s.back();
  //   // s.pop_back();
  //   // L0 = curr_bezier.control_polygon[0];

  //   // L1 = curr_bezier.control_polygon[1];

  //   // L2 = curr_bezier.control_polygon[2];

  //   // if (stopping_criterion(triangles, positions, L0.path.start,
  //   L0.path.end,
  //   //         L2.path.start, L2.path.end, treshold)) {
  //   //   poly_line.push_back(L2.path.end);
  //   // } else {
  //   vector<std::tuple<vector<pair<mesh_point, int>>, vector<vector<float>>>>
  //                   new_points;
  //   control_polygon four_polygon = {
  //       {curr_polygon.geodesics[0], curr_polygon.geodesics[1],
  //           curr_polygon.geodesics[2]},
  //       {curr_polygon.control_points[0], curr_polygon.control_points[1],
  //           curr_polygon.control_points[2], curr_polygon.control_points[3]}};
  //   for (auto i = 0; i < curr_polygon.geodesics.size(); ++i) {
  //     //   if (i % 2) {
  //     //     auto [p, f] = b_spline_subdivision_odd(solver, triangles,
  //     //     positions,
  //     //         adjacencies, v2t, angles,
  //     //         {curr_polygon.geodesics[i - 1].path,
  //     //         curr_polygon.geodesics[i].path,
  //     //             curr_polygon.geodesics[i + 1].path},
  //     //         {curr_polygon.geodesics[i - 1].t,
  //     curr_polygon.geodesics[i].t,
  //     //             curr_polygon.geodesics[i + 1].t},
  //     //         {a, d});

  //     //     new_points.push_back({p, f});
  //     //   } else
  //     //     new_points.push_back(
  //     //         b_spline_subdivision_even(solver, triangles, positions,
  //     //         adjacencies,
  //     //             v2t, angles, {curr_polygon.geodesics[i].path},
  //     //             {curr_polygon.geodesics[i].t}, {a, d}));
  //     // }
  //     new_points.push_back(b_spline_subdivision(solver, triangles, positions,
  //         adjacencies, v2t, angles,
  //         {four_polygon.geodesics[0].path, four_polygon.geodesics[1].path,
  //             four_polygon.geodesics[2].path},
  //         {four_polygon.geodesics[0].t, four_polygon.geodesics[1].t,
  //             four_polygon.geodesics[2].t},
  //         {four_polygon.geodesics[0].path.start,
  //             four_polygon.geodesics[1].path.start,
  //             four_polygon.geodesics[2].path.start,
  //             four_polygon.geodesics[2].path.end}));

  //     // strip   = get_strip_having_distances(solver, triangles, positions,
  //     //     adjacencies, v2t, angles, scalar_fields[4],
  //     //     new_control_points[4].first, new_control_points[3].first,
  //     parents);
  //     // l2.path = shortest_path(triangles, positions, adjacencies,
  //     //     new_control_points[3].first, new_control_points[4].first,
  //     strip);
  //     // l2.t    = path_parameters(l2.path, triangles, positions,
  //     adjacencies);
  //     // std::tie(L2.path, L2.t) = cut_curve(L2.path, L2.t,
  //     //     new_control_points[4].second, new_control_points[4].first,
  //     false);
  //     // s.push_back({{l1, l2, L2}});

  //     // std::tie(L0.path, L0.t) = cut_curve(L0.path, L0.t,
  //     //     new_control_points[0].second, new_control_points[0].first,
  //     true);
  //     // strip   = get_strip_having_distances(solver, triangles, positions,
  //     //     adjacencies, v2t, angles, scalar_fields[1],
  //     //     new_control_points[1].first, new_control_points[0].first,
  //     parents);
  //     // l1.path = shortest_path(triangles, positions, adjacencies,
  //     //     new_control_points[0].first, new_control_points[1].first,
  //     strip);
  //     // l1.t    = path_parameters(l1.path, triangles, positions,
  //     adjacencies);
  //     // strip   = get_strip_having_distances(solver, triangles, positions,
  //     //     adjacencies, v2t, angles, scalar_fields[2],
  //     //     new_control_points[2].first, new_control_points[1].first,
  //     parents);
  //     // l2.path = shortest_path(triangles, positions, adjacencies,
  //     //     new_control_points[1].first, new_control_points[2].first,
  //     strip);
  //     // l2.t    = path_parameters(l2.path, triangles, positions,
  //     adjacencies);

  //     // s.push_back({{L0, l1, l2}});
  //     //}
  //   }
  //   //  for (auto i = 0; i < curr_polygon.control_points.size(); ++i) {
  //   //    poly_line.push_back(curr_polygon.control_points[i]);
  //   //  }

  //   return curr_polygon.control_points;
}
std::tuple<vector<mesh_point>, vector<vector<float>>, vector<vector<vec3f>>>
split_control_poligon_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const Eigen::SparseMatrix<double, 1>& Grad, const bezier_splits& curve,
    const std::unordered_map<int, int>& subdivision, const int entry,
    const float& t) {
  vector<int>           strip;
  auto                  p0     = curve.control_polygons[entry][0];
  auto                  p1     = curve.control_polygons[entry][1];
  auto                  p2     = curve.control_polygons[entry][2];
  auto                  p3     = curve.control_polygons[entry][3];
  auto                  offset = subdivision.at(entry);
  mesh_point            p12    = {-1, zero2f};
  vector<mesh_point>    new_control_points(5);
  vector<float>         weights(3);
  vector<float>         f;
  vector<vector<float>> scalar_field(5);
  vector<vector<vec3f>> gradient_field(5);
  vector<int>           parents;
  vector<geodesic_path> control_polygon(3);
  vector<vector<float>> control_polygon_t(3);

  for (auto i = 0; i < 3; ++i) {
    strip = get_strip_having_distances(solver, triangles, positions,
        adjacencies, v2t, angles, curve.distance_fields[offset][i + 1],
        curve.control_polygons[entry][i + 1], curve.control_polygons[entry][i],
        parents);

    control_polygon[i]   = shortest_path(triangles, positions, adjacencies,
          curve.control_polygons[entry][i], curve.control_polygons[entry][i + 1],
          strip);
    control_polygon_t[i] = path_parameters(
        control_polygon[i], triangles, positions, adjacencies);
  }
  // compute p01
  new_control_points[0] = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[0], control_polygon_t[0], t)
                              .first;

  // compute p23
  new_control_points[4] = sample_on_curve(triangles, positions, adjacencies,
      control_polygon[2], control_polygon_t[2], t)
                              .first;

  weights = {(float)pow(1 - t, 2), 2 * t * (1 - t), (float)pow(t, 2)};
  vector<int>        order              = {};
  vector<mesh_point> curr_points        = {p0, p1, p2};
  std::tie(curr_points, weights, order) = sort_points_and_weigts_quadric(
      curr_points, weights, t);
  new_control_points[1] = weighted_average_dyn_quadric(dual_solver, solver,
      triangles, positions, adjacencies, v2t, angles,
      {curve.distance_fields[offset][0], curve.distance_fields[offset][1],
          curve.distance_fields[offset][2]},
      {control_polygon[0], control_polygon[1]},
      {control_polygon_t[0], control_polygon_t[1]}, curr_points, weights,
      order);

  weights     = {(float)pow(1 - t, 2), 2 * t * (1 - t), (float)pow(t, 2)};
  curr_points = {p1, p2, p3};
  std::tie(curr_points, weights, order) = sort_points_and_weigts_quadric(
      curr_points, weights, t);
  new_control_points[3] = weighted_average_dyn_quadric(dual_solver, solver,
      triangles, positions, adjacencies, v2t, angles,
      {curve.distance_fields[offset][1], curve.distance_fields[offset][2],
          curve.distance_fields[offset][3]},
      {control_polygon[1], control_polygon[2]},
      {control_polygon_t[1], control_polygon_t[2]}, curr_points, weights,
      order);
  // compute p (we are computing two times a distance field)
  strip = get_strip(solver, triangles, positions, adjacencies, v2t, angles,
      new_control_points[3], new_control_points[1], parents);

  auto l   = shortest_path(triangles, positions, adjacencies,
        new_control_points[1], new_control_points[3], strip);
  auto l_t = path_parameters(l, triangles, positions, adjacencies);
  new_control_points[2] =
      sample_on_curve(triangles, positions, adjacencies, l, l_t, 0.5).first;

  // udpate distance and gradient fields
  // p01
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 1,
      {p0, new_control_points[0], new_control_points[1], new_control_points[2]},
      graph);
  scalar_field[0] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  auto F            = wrapper(f);
  gradient_field[0] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // p012
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 2,
      {p0, new_control_points[0], new_control_points[1], new_control_points[2]},
      graph);
  scalar_field[1] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[1] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // P (pruned so that all the control points are reached)
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 3,
      {p0, new_control_points[0], new_control_points[1], new_control_points[2],
          new_control_points[3], new_control_points[4], p3},
      graph);
  scalar_field[2] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[2] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // p123
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 1,
      {new_control_points[2], new_control_points[3], new_control_points[4], p3},
      graph);
  scalar_field[3] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[3] = compute_grad(
      solver, triangles, positions, normals, Grad, F);
  // p23
  f = compute_distance_field(triangles, positions, adjacencies, v2t, solver, 2,
      {new_control_points[2], new_control_points[3], new_control_points[4], p3},
      graph);
  scalar_field[4] = f;
  std::transform(f.begin(), f.end(), f.begin(),
      [](float lambda) { return lambda * lambda; });
  F                 = wrapper(f);
  gradient_field[4] = compute_grad(
      solver, triangles, positions, normals, Grad, F);

  return {new_control_points, scalar_field, gradient_field};
}
std::tuple<vector<mesh_point>, vector<mesh_point>, vector<mesh_point>,
    vector<vector<pair<mesh_point, int>>>>
smoothing_bezier_curve(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<vector<vec3f>>&          gradients,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field) {
  vector<vector<vec3f>>                 dir   = {};
  vector<vector<int>>                   tri   = {};
  vector<pair<int, float>>              verts = {};
  vector<vector<mesh_point>>            steps = {};
  vector<vector<pair<mesh_point, int>>> badones(2);
  vector<mesh_point>                    new_control_points = {};

  /* info("first number of subdivision $",
  first_number_of_subdivision); info("first number of subdivision $",
  second_number_of_subdivision);*/
  new_control_points = split_control_poligon(dual_solver, solver, angles,
      total_angles, triangles, positions, normals, adjacencies, v2t, gradients,
      control_points);

  auto first_split  = bezier_karcher_test(dual_solver, solver, triangles,
       positions, adjacencies, normals, v2t, angles, total_angles, Grad,
       {new_control_points[0], new_control_points[1], new_control_points[2],
          new_control_points[3]},
       number_of_subdivision / 2, type_of_field, verts, tri, dir, steps,
       badones[0]);
  auto second_split = bezier_karcher_test(dual_solver, solver, triangles,
      positions, adjacencies, normals, v2t, angles, total_angles, Grad,
      {new_control_points[3], new_control_points[4], new_control_points[5],
          new_control_points[6]},
      number_of_subdivision / 2, type_of_field, verts, tri, dir, steps,
      badones[1]);

  return {first_split, second_split, new_control_points, badones};
}
vector<mesh_point> bezier_karcher(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field) {
  time_function();
  bezier_splits bezier_curve = {};
  bezier_curve.distance_fields.resize(
      20, vector<vector<float>>(4, vector<float>(positions.size())));
  bezier_curve.gradient_fields.resize(
      20, vector<vector<vec3f>>(4, vector<vec3f>(positions.size())));

  vector<int> slots(20, 0);
  bezier_curve.control_polygons.push_back(control_points);
  vector<float>                f;
  std::unordered_map<int, int> subdivision;
  subdivision[0] = 0;
  for (auto j = 0; j < 4; ++j) {
    f = compute_distance_field(triangles, positions, adjacencies, v2t, solver,
        j, control_points, type_of_field);
    bezier_curve.distance_fields[0][j] = f;
    std::transform(f.begin(), f.end(), f.begin(),
        [](float lambda) { return lambda * lambda; });
    auto F                             = wrapper(f);
    bezier_curve.gradient_fields[0][j] = compute_grad(
        solver, triangles, positions, normals, Grad, F);
  }
  bezier_curve.subidivisions            = {number_of_subdivision};
  slots[0]                              = 1;
  auto                       fine_curve = false;
  auto                       count      = 0;
  auto                       entry = -1, next_entry = -1;
  std::pair<mesh_point, int> bad_sample;
  vector<mesh_point>         polyline = {};
  vector<vector<vec3f>>      gradients(4);
  while (/*!fine_curve &&*/ count < 20) {
    polyline.clear();
    fine_curve = true;

    for (auto i = 0; i < bezier_curve.control_polygons.size(); ++i) {
      bad_sample = {{-1, zero2f}, -1};

      auto poly = bezier_karcher_test(solver, triangles, positions, adjacencies,
          normals, v2t, angles, total_angles,
          bezier_curve.gradient_fields[subdivision.at(i)],
          bezier_curve.control_polygons[i], bezier_curve.subidivisions[i],
          bad_sample);

      if (bad_sample.first.face >= 0) {
        info("happens");
        ++count;

        entry            = subdivision.at(i);
        fine_curve       = false;
        auto curr_subdiv = bezier_curve.subidivisions[i];
        auto t           = bad_sample.second / pow(2, curr_subdiv);
        auto p0          = bezier_curve.control_polygons[i][0];
        auto p3          = bezier_curve.control_polygons[i][3];

        auto f3 = bezier_curve.distance_fields[entry][3];
        auto G3 = bezier_curve.gradient_fields[entry][3];
        auto [new_control_points, new_scalar_fields, new_gradient_fields] =
            split_control_poligon_dyn(dual_solver, solver, angles, total_angles,
                triangles, positions, normals, adjacencies, v2t, Grad,
                bezier_curve, subdivision, i, t);
        // split_control_poligon(solver, angles, total_angles, triangles,
        //     positions, normals, adjacencies, v2t, Grad, bezier_curve,
        //     subdivision, i);
        bezier_curve.subidivisions[i] = (int)std::log2f(bad_sample.second);
        bezier_curve.subidivisions.insert(
            bezier_curve.subidivisions.begin() + i + 1,
            (int)std::log2f(pow(2, curr_subdiv) - bad_sample.second));
        bezier_curve.control_polygons[i] = {p0, new_control_points[0],
            new_control_points[1], new_control_points[2]};
        bezier_curve.control_polygons.insert(
            bezier_curve.control_polygons.begin() + i + 1,
            {new_control_points[2], new_control_points[3],
                new_control_points[4], p3});

        bezier_curve.distance_fields[entry] = {
            bezier_curve.distance_fields[entry][0], new_scalar_fields[0],
            new_scalar_fields[1], new_scalar_fields[2]};
        for (auto j = 0; j < 20; ++j) {
          if (!slots[j]) {
            next_entry = j;
            slots[j]   = 1;
            break;
          }
        }
        if (subdivision.find(i + 1) == subdivision.end())
          subdivision[i + 1] = next_entry;
        else {
          for (auto j = i + 1; j < bezier_curve.control_polygons.size() - 1;
               ++j) {
            auto old_entry = subdivision.at(j);
            subdivision[j] = next_entry;
            if (subdivision.find(j + 1) == subdivision.end())
              subdivision[j + 1] = old_entry;
            else
              next_entry = subdivision.at(j + 1);
          }
        }

        bezier_curve.distance_fields[next_entry] = {new_scalar_fields[2],
            new_scalar_fields[3], new_scalar_fields[4], f3};
        bezier_curve.gradient_fields[entry]      = {
            bezier_curve.gradient_fields[entry][0], new_gradient_fields[0],
            new_gradient_fields[1], new_gradient_fields[2]};
        bezier_curve.gradient_fields[next_entry] = {new_gradient_fields[2],
            new_gradient_fields[3], new_gradient_fields[4], G3};

        break;
      }

      polyline.insert(polyline.end(), poly.begin(), poly.end());
      if (i == bezier_curve.control_polygons.size() - 1) return polyline;
      if (count >= 20) return polyline;
    }
  }
}

std::tuple<geodesic_path, geodesic_path, geodesic_path, vector<float>,
    vector<float>, vector<float>>
compute_tet(const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>&      adjacencies,
    const vector<mesh_point>& control_points) {
  auto strip = strip_on_dual_graph(dual_solver, triangles, positions,
      control_points[1].face, control_points[0].face);
  auto L01 = shortest_path(triangles, positions, adjacencies, control_points[0],
      control_points[1], strip);
  auto L01t = path_parameters(L01, triangles, positions, adjacencies);
  strip     = strip_on_dual_graph(dual_solver, triangles, positions,
          control_points[2].face, control_points[0].face);
  auto L02 = shortest_path(triangles, positions, adjacencies, control_points[0],
      control_points[2], strip);
  auto L02t = path_parameters(L02, triangles, positions, adjacencies);
  strip     = strip_on_dual_graph(dual_solver, triangles, positions,
          control_points[3].face, control_points[0].face);
  auto L03 = shortest_path(triangles, positions, adjacencies, control_points[0],
      control_points[3], strip);
  auto L03t = path_parameters(L03, triangles, positions, adjacencies);

  return {L01, L02, L03, L01t, L02t, L03t};
}
mesh_point weighted_average_tet(const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& L01,
    const geodesic_path& L02, const geodesic_path& L03,
    const vector<float>& L01t, const vector<float>& L02t,
    const vector<float>& L03t, const float& t) {
  vector<float> w;
  bernstein_polynomials(4, t, w);
  auto d               = w[1] + w[2] + w[3];
  auto [P01_d, entry0] = sample_on_curve(
      triangles, positions, adjacencies, L01, L01t, d);
  auto [P02_d, entry1] = sample_on_curve(
      triangles, positions, adjacencies, L02, L02t, d);
  auto [P03_d, entry2] = sample_on_curve(
      triangles, positions, adjacencies, L03, L03t, d);
  auto strip = strip_on_dual_graph(
      dual_solver, triangles, positions, P02_d.face, P01_d.face);
  auto path = shortest_path(
      triangles, positions, adjacencies, P01_d, P02_d, strip);
  auto path_t      = path_parameters(path, triangles, positions, adjacencies);
  auto [Q, entry3] = sample_on_curve(
      triangles, positions, adjacencies, path, path_t, w[2] / (w[2] + w[1]));
  strip = strip_on_dual_graph(
      dual_solver, triangles, positions, P03_d.face, Q.face);
  path   = shortest_path(triangles, positions, adjacencies, Q, P03_d, strip);
  path_t = path_parameters(path, triangles, positions, adjacencies);
  auto [P, entry4] = sample_on_curve(
      triangles, positions, adjacencies, path, path_t, w[3] / d);

  return P;
}
mesh_point weighted_average_tet(const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& L01,
    const geodesic_path& L02, const geodesic_path& L03,
    const vector<float>& L01t, const vector<float>& L02t,
    const vector<float>& L03t, const vector<float>& weights) {
  auto d               = weights[1] + weights[2] + weights[3];
  auto [P01_d, entry0] = sample_on_curve(
      triangles, positions, adjacencies, L01, L01t, d);
  auto [P02_d, entry1] = sample_on_curve(
      triangles, positions, adjacencies, L02, L02t, d);
  auto [P03_d, entry2] = sample_on_curve(
      triangles, positions, adjacencies, L03, L03t, d);

  auto strip = strip_on_dual_graph(
      dual_solver, triangles, positions, P02_d.face, P01_d.face);
  auto path = shortest_path(
      triangles, positions, adjacencies, P01_d, P02_d, strip);
  auto path_t      = path_parameters(path, triangles, positions, adjacencies);
  auto [Q, entry3] = sample_on_curve(triangles, positions, adjacencies, path,
      path_t, weights[2] / (weights[2] + weights[1]));
  strip            = strip_on_dual_graph(
                 dual_solver, triangles, positions, P03_d.face, Q.face);
  path   = shortest_path(triangles, positions, adjacencies, Q, P03_d, strip);
  path_t = path_parameters(path, triangles, positions, adjacencies);
  auto [P, entry4] = sample_on_curve(
      triangles, positions, adjacencies, path, path_t, weights[3] / d);

  return P;
}
vector<mesh_point> bezier_karcher_tet(const geodesic_solver& solver,
    const dual_geodesic_solver& dual_solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision) {
  int n = control_points.size();

  vector<mesh_point> poly_line;  // = {control_points[0]};
  vector<float>      w;
  double             step = 1 / pow(2, number_of_subdivision);

  auto [L01, L02, L03, L01t, L02t, L03t] = compute_tet(
      dual_solver, triangles, positions, adjacencies, control_points);

  auto t = step;
  bernstein_polynomials(n, t, w);
  auto Q = mesh_point{};
  auto P = mesh_point{};

  auto P01_d = mesh_point{};
  auto P02_d = mesh_point{};
  auto P03_d = mesh_point{};
  int  entry;
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    // auto d = w[0] + w[1] + w[3];

    /* if (t <= 0.5) {
      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P21_d.face, P20_d.face);
      auto path = shortest_path(
          triangles, positions, adjacencies, P20_d, P21_d, strip);
      auto path_t = path_parameters(path, triangles, positions,
    adjacencies);

      std::tie(Q, entry) = sample_on_curve(triangles, positions,
    adjacencies, path, path_t, 1 - w[0] / (w[0] + w[1]));

      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P23_d.face, Q.face);
      path = shortest_path(triangles, positions, adjacencies, Q, P23_d,
    strip); path_t = path_parameters(path, triangles, positions,
    adjacencies); std::tie(P, entry) = sample_on_curve(triangles, positions,
    adjacencies, path, path_t, 1 - (w[0] + w[1]) / d);
      poly_line.push_back(P);
    } else {
      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P23_d.face, P21_d.face);
      auto path = shortest_path(
          triangles, positions, adjacencies, P21_d, P23_d, strip);
      auto path_t = path_parameters(path, triangles, positions,
    adjacencies);

      std::tie(Q, entry) = sample_on_curve(triangles, positions,
    adjacencies, path, path_t, 1 - w[1] / (w[3] + w[1]));

      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P20_d.face, Q.face);
      path = shortest_path(triangles, positions, adjacencies, Q, P20_d,
    strip); path_t = path_parameters(path, triangles, positions,
    adjacencies); std::tie(P, entry) = sample_on_curve(triangles, positions,
    adjacencies, path, path_t, 1 - (w[3] + w[1]) / d);
      poly_line.push_back(P);
    } */
    // if (t <= 0.25) {
    auto d                 = w[1] + w[2] + w[3];
    std::tie(P01_d, entry) = sample_on_curve(
        triangles, positions, adjacencies, L01, L01t, d);
    std::tie(P02_d, entry) = sample_on_curve(
        triangles, positions, adjacencies, L02, L02t, d);
    std::tie(P03_d, entry) = sample_on_curve(
        triangles, positions, adjacencies, L03, L03t, d);
    auto strip = strip_on_dual_graph(
        dual_solver, triangles, positions, P02_d.face, P01_d.face);
    auto path = shortest_path(
        triangles, positions, adjacencies, P01_d, P02_d, strip);
    auto path_t = path_parameters(path, triangles, positions, adjacencies);

    std::tie(Q, entry) = sample_on_curve(
        triangles, positions, adjacencies, path, path_t, w[2] / (w[2] + w[1]));

    strip = strip_on_dual_graph(
        dual_solver, triangles, positions, P03_d.face, Q.face);
    path   = shortest_path(triangles, positions, adjacencies, Q, P03_d, strip);
    path_t = path_parameters(path, triangles, positions, adjacencies);
    std::tie(P, entry) = sample_on_curve(
        triangles, positions, adjacencies, path, path_t, w[3] / d);
    poly_line.push_back(P);
    /* } else {
      auto d                 = w[0] + w[1] + w[3];
      std::tie(P20_d, entry) = sample_on_curve(
          triangles, positions, adjacencies, L20, L20t, t);
      std::tie(P21_d, entry) = sample_on_curve(
          triangles, positions, adjacencies, L21, L21t, t);
      std::tie(P23_d, entry) = sample_on_curve(
          triangles, positions, adjacencies, L23, L23t, t);
      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P21_d.face, P20_d.face);
      auto path = shortest_path(
          triangles, positions, adjacencies, P20_d, P21_d, strip);
      auto path_t = path_parameters(path, triangles, positions,
    adjacencies);

      std::tie(Q, entry) = sample_on_curve(triangles, positions,
    adjacencies, path, path_t, w[1] / (w[0] + w[1]));

      strip = strip_on_dual_graph(
          dual_solver, triangles, positions, P23_d.face, Q.face);
      path = shortest_path(triangles, positions, adjacencies, Q, P23_d,
    strip); path_t = path_parameters(path, triangles, positions,
    adjacencies); std::tie(P, entry) = sample_on_curve( triangles,
    positions, adjacencies, path, path_t, w[3] / d); poly_line.push_back(P);
    }
 */
    t += step;
    bernstein_polynomials(n, t, w);
  }
  // poly_line.push_back(control_points.back());
  return poly_line;
}

mesh_point de_casteljau_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& L0,
    const vector<float>& L0t, const geodesic_path& L1, const vector<float>& L1t,
    const geodesic_path& L2, const vector<float>& L2t, const float& t) {
  auto p01 = sample_on_curve(triangles, positions, adjacencies, L0, L0t, t);
  auto p12 = sample_on_curve(triangles, positions, adjacencies, L1, L1t, t);
  auto p23 = sample_on_curve(triangles, positions, adjacencies, L2, L2t, t);

  auto strip = strip_on_dual_graph(
      solver, triangles, positions, p12.first.face, p01.first.face);
  auto L01 = shortest_path(
      triangles, positions, adjacencies, p01.first, p12.first, strip);
  auto L01t = path_parameters(L01, triangles, positions, adjacencies);

  strip = strip_on_dual_graph(
      solver, triangles, positions, p23.first.face, p12.first.face);
  auto L12 = shortest_path(
      triangles, positions, adjacencies, p12.first, p23.first, strip);
  auto L12t = path_parameters(L12, triangles, positions, adjacencies);

  auto p012 = sample_on_curve(triangles, positions, adjacencies, L01, L01t, t);
  auto p123 = sample_on_curve(triangles, positions, adjacencies, L12, L12t, t);

  strip = strip_on_dual_graph(
      solver, triangles, positions, p123.first.face, p012.first.face);
  auto L = shortest_path(
      triangles, positions, adjacencies, p012.first, p123.first, strip);
  auto Lt = path_parameters(L, triangles, positions, adjacencies);
  auto p  = sample_on_curve(triangles, positions, adjacencies, L, Lt, t);

  return p.first;
}
std::tuple<mesh_point, geodesic_path, vector<float>, geodesic_path,
    vector<float>, geodesic_path, vector<float>, geodesic_path, vector<float>>
de_casteljau_point_with_subdivision(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, geodesic_path& L0, vector<float>& L0t,
    geodesic_path& L1, vector<float>& L1t, geodesic_path& L2,
    vector<float>& L2t, const float& t) {
  auto p01 = sample_on_curve(triangles, positions, adjacencies, L0, L0t, t);
  auto p12 = sample_on_curve(triangles, positions, adjacencies, L1, L1t, t);
  auto p23 = sample_on_curve(triangles, positions, adjacencies, L2, L2t, t);

  auto strip = strip_on_dual_graph(
      solver, triangles, positions, p12.first.face, p01.first.face);
  auto L01 = shortest_path(
      triangles, positions, adjacencies, p01.first, p12.first, strip);
  auto L01t = path_parameters(L01, triangles, positions, adjacencies);

  strip = strip_on_dual_graph(
      solver, triangles, positions, p23.first.face, p12.first.face);
  auto L12 = shortest_path(
      triangles, positions, adjacencies, p12.first, p23.first, strip);
  auto L12t = path_parameters(L12, triangles, positions, adjacencies);

  auto p012 = sample_on_curve(triangles, positions, adjacencies, L01, L01t, t);
  auto p123 = sample_on_curve(triangles, positions, adjacencies, L12, L12t, t);

  strip = strip_on_dual_graph(
      solver, triangles, positions, p123.first.face, p012.first.face);
  auto L = shortest_path(
      triangles, positions, adjacencies, p012.first, p123.first, strip);
  auto Lt = path_parameters(L, triangles, positions, adjacencies);
  auto p  = sample_on_curve(triangles, positions, adjacencies, L, Lt, t);

  std::tie(L0, L0t) = cut_curve(L0, L0t, p01.second, p01.first, true);
  auto [gamma0, gamma1, gamma0t, gamma1t] = split_curve(
      L, Lt, p.second, p.first);
  std::tie(L01, L01t) = cut_curve(L01, L01t, p012.second, p012.first, true);
  std::tie(L12, L12t) = cut_curve(L12, L12t, p123.second, p123.first, false);
  std::tie(L2, L2t)   = cut_curve(L2, L2t, p23.second, p23.first, false);
  return std::tie(
      p.first, L01, L01t, gamma0, gamma0t, gamma1, gamma1t, L12, L12t);
}
vector<mesh_point> bezier_karcher_tet_subdivision(const geodesic_solver& solver,
    const dual_geodesic_solver& dual_solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision) {
  vector<mesh_point> poly_line;  // = {control_points[0]};
  vector<float>      w;
  struct parametric_path {
    geodesic_path path;
    vector<float> t;
  };
  struct bezier_split {
    vector<parametric_path> control_polygon;
    vector<parametric_path> previous_polyon;
  };
  auto Q                                 = std::deque<bezier_split>();
  auto [L01, L02, L03, L01t, L02t, L03t] = compute_tet(
      dual_solver, triangles, positions, adjacencies, control_points);
  auto L1 = parametric_path{L01, L01t};
  auto L2 = parametric_path{L02, L02t};
  auto L3 = parametric_path{L03, L03t};

  auto p0 = L1.path.start;
  auto p1 = L1.path.end;
  auto p2 = L2.path.end;
  auto p3 = L3.path.end;
  auto q0 = weighted_average_tet(dual_solver, triangles, positions, adjacencies,
      L1.path, L2.path, L3.path, L1.t, L2.t, L3.t, 1 / 3.f);
  auto q1 = weighted_average_tet(dual_solver, triangles, positions, adjacencies,
      L1.path, L2.path, L3.path, L1.t, L2.t, L3.t, 2 / 3.f);
  /* auto a  = weighted_average_tet(dual_solver, triangles, positions,
      adjacencies, L1.path, L2.path,
      L3.path, L1.t, L2.t,
      L3.t, 1 / 2.f); */
  auto b0 = weighted_average_tet(dual_solver, triangles, positions, adjacencies,
      L1.path, L2.path, L3.path, L1.t, L2.t, L3.t, 1 / 6.f);
  auto b1 = weighted_average_tet(dual_solver, triangles, positions, adjacencies,
      L1.path, L2.path, L3.path, L1.t, L2.t, L3.t, 5 / 6.f);
  auto strip = strip_on_dual_graph(
      dual_solver, triangles, positions, p2.face, p1.face);
  auto L12  = shortest_path(triangles, positions, adjacencies, p1, p2, strip);
  auto L12t = path_parameters(L12, triangles, positions, adjacencies);
  strip     = strip_on_dual_graph(
          dual_solver, triangles, positions, p3.face, p2.face);
  auto L23    = shortest_path(triangles, positions, adjacencies, p2, p3, strip);
  auto L23t   = path_parameters(L23, triangles, positions, adjacencies);
  auto a      = mesh_point{};
  auto gamma0 = geodesic_path{};
  auto gamma1 = geodesic_path{};
  auto gamma2 = geodesic_path{};
  auto gamma3 = geodesic_path{};

  auto gamma0t = vector<float>{};
  auto gamma1t = vector<float>{};
  auto gamma2t = vector<float>{};
  auto gamma3t = vector<float>{};

  std::tie(a, gamma0, gamma0t, gamma1, gamma1t, gamma2, gamma2t, gamma3,
      gamma3t) = de_casteljau_point_with_subdivision(dual_solver, triangles,
      positions, adjacencies, L1.path, L1.t, L12, L12t, L23, L23t, 0.5);

  std::tie(L01, L02, L03, L01t, L02t, L03t) = compute_tet(
      dual_solver, triangles, positions, adjacencies, {a, b1, q1, p3});
  auto split = bezier_split{
      vector<parametric_path>{parametric_path{L01, L01t},
          parametric_path{L02, L02t}, parametric_path{L03, L03t}},
      vector<parametric_path>{parametric_path{gamma2, gamma2t},
          parametric_path{gamma3, gamma3t}, parametric_path{L23, L23t}}};
  Q.push_back(split);
  std::tie(L01, L02, L03, L01t, L02t, L03t) = compute_tet(
      dual_solver, triangles, positions, adjacencies, {p0, b0, q0, a});
  split = bezier_split{
      vector<parametric_path>{parametric_path{L01, L01t},
          parametric_path{L02, L02t}, parametric_path{L03, L03t}},
      vector<parametric_path>{L1, parametric_path{gamma0, gamma0t},
          parametric_path{gamma1, gamma1t}}};

  Q.push_back(split);

  while (!Q.empty()) {
    auto curr_polygon = Q.back().control_polygon;
    auto prev_polygon = Q.back().previous_polyon;
    Q.pop_back();
    p0 = curr_polygon[0].path.start;
    p1 = curr_polygon[0].path.end;
    p2 = curr_polygon[1].path.end;
    p3 = curr_polygon[2].path.end;

    if (stopping_criterion(triangles, positions, p0, p1, p2, p3, 1e-6))
      poly_line.push_back(p3);
    else {
      auto q0 = weighted_average_tet(dual_solver, triangles, positions,
          adjacencies, curr_polygon[0].path, curr_polygon[1].path,
          curr_polygon[2].path, curr_polygon[0].t, curr_polygon[1].t,
          curr_polygon[2].t, 1 / 3.f);
      auto q1 = weighted_average_tet(dual_solver, triangles, positions,
          adjacencies, curr_polygon[0].path, curr_polygon[1].path,
          curr_polygon[2].path, curr_polygon[0].t, curr_polygon[1].t,
          curr_polygon[2].t, 2 / 3.f);
      /* auto a  = weighted_average_tet(dual_solver, triangles, positions,
          adjacencies, curr_polygon[0].path, curr_polygon[1].path,
          curr_polygon[2].path, curr_polygon[0].t, curr_polygon[1].t,
          curr_polygon[2].t, 1 / 2.f); */
      auto b0    = weighted_average_tet(dual_solver, triangles, positions,
             adjacencies, curr_polygon[0].path, curr_polygon[1].path,
             curr_polygon[2].path, curr_polygon[0].t, curr_polygon[1].t,
             curr_polygon[2].t, 1 / 6.f);
      auto b1    = weighted_average_tet(dual_solver, triangles, positions,
             adjacencies, curr_polygon[0].path, curr_polygon[1].path,
             curr_polygon[2].path, curr_polygon[0].t, curr_polygon[1].t,
             curr_polygon[2].t, 5 / 6.f);
      auto strip = strip_on_dual_graph(
          dual_solver, triangles, positions, p2.face, p1.face);
      auto L12 = shortest_path(
          triangles, positions, adjacencies, p1, p2, strip);
      auto L12t = path_parameters(L12, triangles, positions, adjacencies);
      strip     = strip_on_dual_graph(
              dual_solver, triangles, positions, p3.face, p2.face);
      auto L23 = shortest_path(
          triangles, positions, adjacencies, p2, p3, strip);
      auto L23t = path_parameters(L23, triangles, positions, adjacencies);

      std::tie(a, gamma0, gamma0t, gamma1, gamma1t, gamma2, gamma2t, gamma3,
          gamma3t) = de_casteljau_point_with_subdivision(dual_solver, triangles,
          positions, adjacencies, prev_polygon[0].path, prev_polygon[0].t,
          prev_polygon[1].path, prev_polygon[1].t, prev_polygon[2].path,
          prev_polygon[2].t, 0.5);

      std::tie(L01, L02, L03, L01t, L02t, L03t) = compute_tet(
          dual_solver, triangles, positions, adjacencies, {a, b1, q1, p3});
      split = bezier_split{
          vector<parametric_path>{parametric_path{L01, L01t},
              parametric_path{L02, L02t}, parametric_path{L03, L03t}},
          vector<parametric_path>{parametric_path{gamma2, gamma2t},
              parametric_path{gamma3, gamma3t}, prev_polygon[2]}};
      Q.push_back(split);
      std::tie(L01, L02, L03, L01t, L02t, L03t) = compute_tet(
          dual_solver, triangles, positions, adjacencies, {p0, b0, q0, a});
      split = bezier_split{
          vector<parametric_path>{parametric_path{L01, L01t},
              parametric_path{L02, L02t}, parametric_path{L03, L03t}},
          vector<parametric_path>{prev_polygon[0],
              parametric_path{gamma0, gamma0t},
              parametric_path{gamma1, gamma1t}}};

      Q.push_back(split);
    }
  }
  return poly_line;
}

vector<vector<float>> control_points_distance_field(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<mesh_point>& control_points,
    const int type_of_field) {
  vector<vector<float>> f(control_points.size());

  for (auto i = 0; i < control_points.size(); ++i) {
    f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
        solver, i, control_points, type_of_field);
    std::transform(f[i].begin(), f[i].end(), f[i].begin(),
        [](float lambda) { return lambda * lambda; });
  }
  return f;
}

vector<vector<vec3f>> control_points_gradients(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<vector<float>>& f) {
  vector<vector<vec3f>> grad(f.size());

  for (auto i = 0; i < f.size(); ++i) {
    auto F  = wrapper(f[i]);
    grad[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }

  return grad;
}
int get_seed(const vector<vec3i>& triangles, const mesh_point& control_point) {
  auto bary     = get_bary(control_point.uv);
  auto max_bary = flt_min;
  auto entry    = -1;
  for (auto i = 0; i < 3; ++i) {
    if (bary[i] > max_bary) {
      max_bary = bary[i];
      entry    = i;
    }
  }
  return triangles[control_point.face][entry];
}
vector<mesh_point> solve_for_control_points(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& gradients,
    const vector<float>& t, const float& step,
    vector<pair<mesh_point, int>>& badones) {
  auto n       = control_points.size();
  auto seed    = -1;
  auto w       = vector<float>{};
  auto samples = vector<mesh_point>(t.size());
  for (auto i = 0; i < t.size(); ++i) {
    if (t[i] == 1 / 6.f)
      seed = get_seed(triangles, control_points[0]);
    else if (t[i] <= 0.5)
      seed = get_seed(triangles, control_points[1]);
    else if (t[i] <= 2 / 3.f)
      seed = get_seed(triangles, control_points[2]);
    else
      seed = get_seed(triangles, control_points[3]);

    bernstein_polynomials(n, t[i], w);
    auto tids  = vector<int>{};
    auto dir   = vector<vec3f>{};
    auto steps = vector<mesh_point>{};
    auto verts = vector<pair<int, float>>{};
    auto size  = badones.size();
    samples[i] = weighted_average_test(dual_solver, solver, angles,
        total_angles, triangles, positions, normals, adjacencies, v2t, f,
        gradients, w, control_points, verts, tids, dir, steps, badones, seed);
    if (size != badones.size()) badones.back().second = t[i];
  }

  return samples;
}
vector<mesh_point> solve_for_intermediate_points(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& gradients,
    const vector<float>& t, const float& step,
    vector<pair<mesh_point, int>>& badones) {
  auto seed    = -1;
  auto w       = vector<float>{};
  auto samples = vector<mesh_point>(3);
  for (auto i = 0; i < t.size(); ++i) {
    if (t[i] <= 0.5)
      seed = get_seed(triangles, control_points[0]);
    else
      seed = get_seed(triangles, control_points[1]);

    switch (i) {
      case 0: {
        w = {135 / 432.f, 405 / 432.f, -135 / 432.f, 27 / 432.f};
      } break;
      case 1: {
        w = {-1 / 16.f, 27 / 48.f, 27 / 48.f, -1 / 16.f};
      } break;
      case 2: {
        w = {27 / 432.f, -135 / 432.f, 405 / 432.f, 135 / 432.f};
      } break;
    }
    auto size  = badones.size();
    auto tids  = vector<int>{};
    auto dir   = vector<vec3f>{};
    auto steps = vector<mesh_point>{};
    auto verts = vector<pair<int, float>>{};
    samples[i] = weighted_average_test(dual_solver, solver, angles,
        total_angles, triangles, positions, normals, adjacencies, v2t, f,
        gradients, w, control_points, verts, tids, dir, steps, badones, seed);
    if (size != badones.size()) badones.back().second = t[i];
  }

  return samples;
}

// we assume that v is defined in the tangent space of from.face
vector<mesh_point> straightest_geodesic(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const mesh_point& from, const vec3f& v, const float& l) {
  auto  vid = -1, tid = -1, next_tri = -1;
  auto  dir       = v;
  float len       = 0.0;
  auto  next_bary = zero3f, next = zero3f;
  auto  prev            = eval_position(triangles, positions, from);
  auto  samples         = vector<mesh_point>{from};
  auto  bary            = get_bary(from.uv);
  auto [is_vert, kv]    = bary_is_vert(bary);
  auto [is_on_edge, ke] = bary_is_edge(bary);
  if (is_vert) {
    vid = triangles[from.face][kv];
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, from.face, vid, T2V);

    tid = next_tid(
        solver, angles, positions, v2p_adjacencies, triangles, normals, vid, v);

    /* tid = next_tid_extended_graph(
        solver, angles, positions, v2p_adjacencies, triangles, normals, vid,
       v);
     */
    kv       = find(triangles[tid], vid);
    bary     = zero3f;
    bary[kv] = 1;
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, dir, normals, vid, tid, V2T);
  } else if (is_on_edge) {
    auto p0   = triangles[from.face][ke];
    auto p1   = triangles[from.face][(ke + 1) % 3];
    auto p2   = triangles[from.face][(ke + 2) % 3];
    auto n    = triangle_normal(positions[p0], positions[p1], positions[p2]);
    auto edge = normalize(positions[p1] - positions[p0]);
    if (dot(cross(edge, v), n) > 0)
      tid = from.face;
    else {
      tid  = adjacencies[from.face][ke];
      bary = tri_bary_coords(triangles, positions, tid, prev);
      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, from.face, tid, T2T);
    }

  } else
    tid = from.face;

  while (len < l) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    samples.push_back({tid, vec2f{next_bary.y, next_bary.z}});
    len += length(next - prev);
    if (len < l) {
      prev          = next;
      auto [V, k_v] = bary_is_vert(next_bary);
      auto [E, k_e] = bary_is_edge(next_bary);
      if (V) {
        vid       = triangles[tid][k_v];
        auto out  = handle_vert(solver, triangles, positions, normals,
             adjacencies, v2p_adjacencies, angles, total_angles, vid, tid, dir);
        tid       = out.first;
        dir       = out.second;
        k_v       = find(triangles[tid], vid);
        bary      = zero3f;
        bary[k_v] = 1;
      } else if (E) {
        next_tri      = adjacencies[tid][k_e];
        auto p0       = triangles[tid][k_e];
        auto offset0  = find(triangles[next_tri], p0);
        auto offset1  = (offset0 + 2) % 3;
        bary          = zero3f;
        bary[offset0] = next_bary[k_e];
        bary[offset1] = next_bary[(k_e + 1) % 3];

        parallel_transp(solver, angles, total_angles, triangles, positions,
            adjacencies, dir, normals, tid, next_tri, T2T);
        tid = next_tri;
      } else
        assert(false);
    }
  }

  auto factor = (len - l);
  auto w      = normalize(prev - next);
  w *= factor;
  w += next;
  bary = tri_bary_coords(triangles, positions, tid, w);
  samples.pop_back();
  samples.push_back({tid, vec2f{bary.y, bary.z}});

  return samples;
}

pair<int, float> straightest_to_vert(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const int& from, const int k,
    const vec3f& v) {
  auto tid = -1, next_tri = -1, crossed = 0;
  auto dir       = v;
  auto len       = 0.f;
  auto next_bary = zero3f, bary = zero3f, next = zero3f, prev = positions[from];
  tid = next_tid(solver, angles, positions, v2t, triangles, normals, from, v);
  auto offset  = find(triangles[tid], from);
  bary[offset] = 1;
  parallel_transp(solver, angles, total_angles, triangles, positions,
      adjacencies, dir, normals, from, tid, V2T);

  while (crossed < k) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    len += length(next - prev);

    prev          = next;
    auto [V, k_v] = bary_is_vert(next_bary, 1e-1);
    auto [E, k_e] = bary_is_edge(next_bary);
    if (V) {
      auto vid = triangles[tid][k_v];
      len -= length(next - prev);
      len += length(positions[vid] - prev);
      return std::make_pair(vid, len);
    } else if (E) {
      ++crossed;
      next_tri      = adjacencies[tid][k_e];
      auto p0       = triangles[tid][k_e];
      auto offset0  = find(triangles[next_tri], p0);
      auto offset1  = (offset0 + 2) % 3;
      bary          = zero3f;
      bary[offset0] = next_bary[k_e];
      bary[offset1] = next_bary[(k_e + 1) % 3];

      parallel_transp(solver, angles, total_angles, triangles, positions,
          adjacencies, dir, normals, tid, next_tri, T2T);
      tid = next_tri;
    } else
      assert(false);
  }

  return std::make_pair(-1, len);
}
mesh_point transport_vector(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const geodesic_path& path,
    const mesh_point& point) {
  auto from = path.start;

  auto pos   = zero3f;
  auto start = eval_position(triangles, positions, path.start);
  if (path.strip.size() <= 1)
    pos = eval_position(triangles, positions, path.end) - start;
  else {
    auto eid = get_edge(
        triangles, positions, adjacencies, path.strip[0], path.strip[1]);
    assert(eid.x != -1);
    auto lerp = path.lerps[0];
    pos       = (1 - lerp) * positions[eid.x] + lerp * positions[eid.y];
  }
  auto displacement = normalize(pos - start);
  auto v = tranport_vector(solver, angles, total_angles, triangles, positions,
      adjacencies, v2t, displacement, normals, from, point);
  auto [is_vert, kv] = bary_is_vert(get_bary(point.uv));
  auto tid           = -1;
  if (is_vert) {
    auto vid = triangles[point.face][kv];
    tid = next_tid(solver, angles, positions, v2t, triangles, normals, vid, v);
    parallel_transp(solver, angles, total_angles, triangles, positions,
        adjacencies, v, normals, vid, tid, V2T);
  } else
    tid = point.face;

  auto px = triangles[point.face].x, py = triangles[point.face].y,
       pz        = triangles[point.face].z;
  auto bary      = get_bary(point.uv);
  auto next_bary = zero3f;
  trace_in_triangles(positions, triangles, v, bary, point.face, pos, next_bary);
  auto point_pos = eval_position(triangles, positions, point);
  auto l         = length(pos - point_pos);
  auto path_len  = path_length(path, triangles, positions, adjacencies);
  if (l > path_len) {
    auto w = normalize(pos - point_pos);
    w *= path_len;
    w += point_pos;
    bary = tri_bary_coords(positions[px], positions[py], positions[pz], w);
    return {point.face, vec2f{bary.x, bary.y}};
  }

  auto straightest = straightest_geodesic(solver, triangles, positions, normals,
      adjacencies, v2t, angles, total_angles, point, v, path_len);
  return straightest.back();
}
inline vector<float> subdivide_angles(const int number_of_subdivision) {
  auto tetas = vector<float>(number_of_subdivision);
  auto step  = 2 * M_PI / number_of_subdivision;
  auto phi   = 0.f;
  for (auto i = 0; i < number_of_subdivision; ++i) {
    tetas[i] = phi;
    phi += step;
  }
  return tetas;
}
vector<pair<vec2f, int>> connect_vert_to_neighbors(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const int vid, const vector<float> tetas) {
  vector<pair<vec2f, int>> angles_and_nodes;
  auto                     e = polar_basis(solver, positions, normals, vid);
  auto                     n = normals[vid];
  for (auto i = 0; i < tetas.size(); ++i) {
    auto v         = rot_vect(e, n, tetas[i]);
    auto neighbors = straightest_to_vert(solver, triangles, positions, normals,
        adjacencies, v2t, angles, total_angles, vid, k + 1, v);
    if (neighbors.first != -1) {
      auto last = (angles_and_nodes.size() == 0) ? std::make_pair(vec2f{}, -1)
                                                 : angles_and_nodes.back();
      if (neighbors.first != last.second)
        angles_and_nodes.push_back(
            std::make_pair(vec2f{tetas[i], neighbors.second}, neighbors.first));
      else if (neighbors.second < last.first.y) {
        angles_and_nodes.pop_back();
        angles_and_nodes.push_back(
            std::make_pair(vec2f{tetas[i], neighbors.second}, neighbors.first));
      }
    }
  }
  return angles_and_nodes;
}
Eigen::SparseMatrix<double, 1> spanning_direction_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, geodesic_solver& solver,
    vector<vector<float>>& angles, vector<float>& total_angles) {
  auto old_solver = make_geodesic_solver(
      triangles, positions, adjacencies, v2t);

  auto old_angles = compute_angles(
      triangles, positions, adjacencies, v2t, total_angles, true);

  auto avg_valence = 0;
  solver.graph.resize(positions.size());
  angles.resize(positions.size());
  auto tetas = subdivide_angles(1000);

  for (auto i = 0; i < positions.size(); ++i) {
    auto arcs = connect_vert_to_neighbors(triangles, positions, adjacencies,
        normals, v2t, k, old_solver, old_angles, total_angles, i, tetas);
    solver.graph[i].resize(arcs.size());
    angles[i].resize(arcs.size());
    avg_valence += arcs.size();
    for (int j = 0; j < arcs.size(); ++j) {
      angles[i][j]              = arcs[j].first.x;
      solver.graph[i][j].length = arcs[j].first.y;
      solver.graph[i][j].node   = arcs[j].second;
    }
  }
  auto Grad = init_gradient_matrix(solver, angles, positions, normals);
  std::cout << "avg valence is" << std::endl;
  std::cout << avg_valence / positions.size() << std::endl;
  return Grad;
}
geodesic_solver spanning_direction_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, vector<vector<float>>& angles,
    vector<float>& total_angles) {
  auto old_solver = make_geodesic_solver(
      triangles, positions, adjacencies, v2t);

  auto old_angles = compute_angles(
      triangles, positions, adjacencies, v2t, total_angles, true);
  auto avg_valence = 0;
  auto solver      = geodesic_solver{};
  solver.graph.resize(positions.size());
  angles.resize(positions.size());
  auto tetas = subdivide_angles(1000);

  for (auto i = 0; i < positions.size(); ++i) {
    auto arcs = connect_vert_to_neighbors(triangles, positions, adjacencies,
        normals, v2t, k, old_solver, old_angles, total_angles, i, tetas);
    solver.graph[i].resize(arcs.size());
    angles[i].resize(arcs.size());
    avg_valence += arcs.size();
    for (int j = 0; j < arcs.size(); ++j) {
      angles[i][j]              = arcs[j].first.x;
      solver.graph[i][j].length = arcs[j].first.y;
      solver.graph[i][j].node   = arcs[j].second;
    }
  }
  std::cout << "avg valence is" << std::endl;
  std::cout << avg_valence / positions.size() << std::endl;
  return solver;
}
// Inverse Problem
inline vec3f interp_grad(
    const vec3f& gx, const vec3f& gy, const vec3f& gz, const vec3f& bary) {
  return gx * bary.x + gy * bary.y + gz + bary.z;
}
inline float interp_field(
    const float& fx, const float& fy, const float& fz, const vec3f& bary) {
  return fx * bary.x + fy * bary.y + fz + bary.z;
}
vector<float> coefficients_from_control_points(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vector<vec3f>>& gradients, const vector<vector<float>>& f,
    const vector<mesh_point>& control_points, const mesh_point& point) {
  auto n    = control_points.size();
  auto g    = vector<vec3f>(n);
  auto bary = get_bary(point.uv);
  auto p0 = triangles[point.face].x, p1 = triangles[point.face].y,
       p2     = triangles[point.face].z;
  auto normal = triangle_normal(positions[p0], positions[p1], positions[p2]);
  for (auto i = 0; i < n; ++i) {
    auto v0 = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, gradients[i][p0], normals, p0, point.face, V2T);
    auto v1 = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, gradients[i][p1], normals, p1, point.face, V2T);
    auto v2 = transp_vec(solver, angles, total_angles, triangles, positions,
        adjacencies, gradients[i][p2], normals, p2, point.face, V2T);
    g[i]    = normalize(interp_grad(v0, v1, v2, bary));

    g[i] *= 2 * interp_field(f[i][p0], f[i][p1], f[i][p2], bary);
  }

  auto weights = vector<float>(n);
  auto sum     = 0.f;
  // https://pdfs.semanticscholar.org/b59e/7cfd344c6d4562ec58761bdf51236410867a.pdf
  for (auto i = 0; i < n; ++i) {
    if (length(g[i]) <= 1e-5) {
      weights.resize(n, 0.f);
      weights[i] = 1;
      return weights;
    }
    auto prev_teta = angle(g[(n + i - 1) % n], g[i]);
    auto curr_teta = angle(g[i], g[(i + 1) % n]);
    if (dot(cross(g[(n + i - 1) % n], g[i]), normal) < 0) prev_teta *= -1;
    if (dot(cross(g[i], g[(i + 1) % n]), normal) < 0) curr_teta *= -1;

    weights[i] = (yocto::tan(prev_teta / 2) + yocto::tan(curr_teta / 2)) /
                 length(g[i]);
    sum += weights[i];
  }
  for (auto& w : weights) {
    w /= sum;
  }
  return weights;
}

vector<mesh_point> bezier_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field) {
  time_function();
  auto                  n        = (int)control_points.size();
  vector<mesh_point>    polyline = {control_points[0]};
  vector<vector<float>> f(n);
  vector<geodesic_path> control_polygon(n - 1);
  vector<vector<float>> control_polygon_t(n - 1);

  for (auto i = 0; i < n; ++i) {
    f[i] = compute_distance_field(triangles, positions, adjacencies, v2t,
        solver, i, control_points, type_of_field);
  }
  vector<int> parents;
  for (auto i = 0; i < n - 1; ++i) {
    auto strip = strip_on_dual_graph(dual_solver, triangles, positions,
        control_points[i + 1].face, control_points[i].face);
    // get_strip_having_distances(solver, triangles, positions,
    //     adjacencies, v2t, angles, f[i + 1], control_points[i + 1],
    //     control_points[i], parents);
    // auto strip = get_strip(solver, triangles, positions, adjacencies,
    // v2t,
    //     angles, control_points[i + 1], control_points[i], parents);
    control_polygon[i]   = shortest_path(triangles, positions, adjacencies,
          control_points[i], control_points[i + 1], strip);
    control_polygon_t[i] = path_parameters(
        control_polygon[i], triangles, positions, adjacencies);
  }

  auto          step = 1 / pow(2, number_of_subdivision);
  auto          t    = step;
  vector<float> w;
  bernstein_polynomials(n, t, w);
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    auto [curr_points, weights, order] = sort_points_and_weigts(
        control_points, w, t);
    auto p = weighted_average_dyn(dual_solver, solver, triangles, positions,
        adjacencies, v2t, angles, f, control_polygon, control_polygon_t,
        curr_points, weights, order);
    polyline.push_back(p);
    t += step;
    bernstein_polynomials(n, t, w);
  }

  polyline.push_back(control_points[3]);
  return polyline;
}

vector<mesh_point> bezier_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<vector<float>>& distances,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field) {
  time_function();
  auto                  n        = (int)control_points.size();
  vector<mesh_point>    polyline = {control_points[0]};
  vector<vector<float>> f(n);
  vector<geodesic_path> control_polygon(n - 1);
  vector<vector<float>> control_polygon_t(n - 1);
  vector<int>           parents;
  for (auto i = 0; i < n - 1; ++i) {
    auto strip = get_strip_having_distances(solver, triangles, positions,
        adjacencies, v2t, angles, f[i + 1], control_points[i + 1],
        control_points[i], parents);
    // auto strip = get_strip(solver, triangles, positions, adjacencies,
    // v2t,
    //     angles, control_points[i + 1], control_points[i], parents);
    control_polygon[i]   = shortest_path(triangles, positions, adjacencies,
          control_points[i], control_points[i + 1], strip);
    control_polygon_t[i] = path_parameters(
        control_polygon[i], triangles, positions, adjacencies);
  }

  auto          step = 1 / pow(2, number_of_subdivision);
  auto          t    = step;
  vector<float> w;
  bernstein_polynomials(n, t, w);
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    auto [curr_points, weights, order] = sort_points_and_weigts(
        control_points, w, t);
    auto p = weighted_average_dyn(dual_solver, solver, triangles, positions,
        adjacencies, v2t, angles, f, control_polygon, control_polygon_t,
        curr_points, weights, order);
    polyline.push_back(p);
    t += step;
    bernstein_polynomials(n, t, w);
  }

  polyline.push_back(control_points[3]);
  return polyline;
}
vector<vector<float>> knots_insertion(const vector<float>& t,
    const vector<float>& u, const int k, const int m, const int n) {
  auto h = 3;

  vector<vector<float>> a(4, vector<float>(4));
  for (auto j = 0; j < 4; ++j) {
    for (auto i = 0; i < 4; ++i) {
      a[j][i] = (t[j + i] - t[j] != 0) ? (0.5 - t[j]) / (t[j + i] - t[j]) : 0;
    }
  }

  // vector<float> d(4);
  // for (auto r = 0; r < 4; ++r) {
  //   for (auto j = r; j < 4; ++j) {

  //   }
  // }

  // for (auto j = 0; j < m + 1; ++j) {
  //   auto sum = 0;
  //   for (auto i = 0; i < n; ++i) {
  //     if (t[i + k - 1] - t[i] == 0 && t[i + k] - t[i + 1] == 0)
  //       a_2[i][j] = 0.f;
  //     else if (t[i + k - 1] - t[i] == 0) {
  //       a_2[i][j] = ((t[i + k] - u[j + k - 1]) * a_1[i + 1][j]) /
  //                   (t[i + k] - t[i + 1]);
  //     } else if (t[i + k] - t[i + 1] == 0) {
  //       a_2[i][j] = ((u[j + k - 1] - t[i]) * a_1[i][j]) / (t[i + k - 1] -
  //       t[i]);
  //     } else {
  //       a_2[i][j] =
  //           ((u[j + k - 1] - t[i]) * a_1[i][j]) / (t[i + k - 1] - t[i]) +
  //           ((t[i + k] - u[j + k - 1]) * a_1[i + 1][j]) / (t[i + k] - t[i +
  //           1]);
  //     }

  //     sum += a_2[i][j];
  //     if (sum == 1) break;
  //   }
  // }
  // for (auto j = 0; j < m + 1; ++j) {
  //   auto sum = 0;
  //   for (auto i = 0; i < n; ++i) {
  //     if (t[i + k - 1] - t[i] == 0 && t[i + k] - t[i + 1] == 0)
  //       a_3[i][j] = 0.f;
  //     else if (t[i + k - 1] - t[i] == 0) {
  //       a_3[i][j] = ((t[i + k] - u[j + k - 1]) * a_2[i + 1][j]) /
  //                   (t[i + k] - t[i + 1]);
  //     } else if (t[i + k] - t[i + 1] == 0) {
  //       a_3[i][j] = ((u[j + k - 1] - t[i]) * a_2[i][j]) / (t[i + k - 1] -
  //       t[i]);
  //     } else {
  //       a_3[i][j] =
  //           ((u[j + k - 1] - t[i]) * a_2[i][j]) / (t[i + k - 1] - t[i]) +
  //           ((t[i + k] - u[j + k - 1]) * a_2[i + 1][j]) / (t[i + k] - t[i +
  //           1]);
  //     }

  //     sum += a_3[i][j];
  //     if (sum == 1) break;
  //   }
  // }
  // for (auto j = 0; j < m + 1; ++j) {
  //   auto sum = 0;
  //   for (auto i = 0; i < 4; ++i) {
  //     if (t[i + k - 1] - t[i] == 0 && t[i + k] - t[i + 1] == 0)
  //       a_4[i][j] = 0.f;
  //     else if (t[i + k - 1] - t[i] == 0) {
  //       a_4[i][j] = ((t[i + k] - u[j + k - 1]) * a_3[i + 1][j]) /
  //                   (t[i + k] - t[i + 1]);
  //     } else if (t[i + k] - t[i + 1] == 0) {
  //       a_4[i][j] = ((u[j + k - 1] - t[i]) * a_3[i][j]) / (t[i + k - 1] -
  //       t[i]);
  //     } else {
  //       a_4[i][j] =
  //           ((u[j + k - 1] - t[i]) * a_3[i][j]) / (t[i + k - 1] - t[i]) +
  //           ((t[i + k] - u[j + k - 1]) * a_3[i + 1][j]) / (t[i + k] - t[i +
  //           1]);
  //     }

  //     sum += a_4[i][j];
  //     if (sum == 1) break;
  //   }
  // }
  // return a_3;
}
std::tuple<vector<int>, mesh_point, mesh_point> handle_short_strips(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<int>& strip, const mesh_point& start, const mesh_point& end) {
  if (strip.size() == 1) {
    return {strip, start, end};
  } else if (strip.size() == 2) {
    auto [inside, b2f] = point_in_triangle(triangles, positions, start.face,
        eval_position(triangles, positions, end));
    if (inside) {
      auto new_end = mesh_point{start.face, b2f};
      return {{start.face}, start, new_end};
    }
    std::tie(inside, b2f) = point_in_triangle(triangles, positions, end.face,
        eval_position(triangles, positions, start));

    if (inside) {
      auto new_start = mesh_point{end.face, b2f};
      return {{end.face}, new_start, end};
    }

    return {strip, start, end};
  }
  return {{-1}, {}, {}};
}
std::tuple<vector<int>, mesh_point, mesh_point> cleaned_strip(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<int>& strip,
    const mesh_point& start, const mesh_point& end) {
  vector<int> cleaned = strip;

  auto start_entry = 0, end_entry = (int)strip.size() - 1;
  auto b3f           = zero3f;
  auto new_start     = start;
  auto new_end       = end;
  auto [is_vert, kv] = point_is_vert(end);
  auto [is_edge, ke] = point_is_edge(end);
  if (strip.size() <= 2)
    return handle_short_strips(triangles, positions, strip, start, end);
  // Erasing from the bottom
  if (is_vert) {
    auto vid      = triangles[end.face][kv];
    auto curr_tid = strip[end_entry - 1];
    kv            = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.pop_back();
      --end_entry;
      if (end_entry == 1) break;
      // see comment below
      auto curr_tid = strip[end_entry - 1];
      kv            = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned.back()], vid);
    assert(kv != -1);
    b3f[kv] = 1;
    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}};  // updating end
  } else if (is_edge) {
    if (end.face != strip.back()) {
      assert(adjacencies[end.face][ke] == strip.back());

      if (end.face == strip[end_entry - 1]) cleaned.pop_back();
    } else if (adjacencies[end.face][ke] == strip[end_entry - 1]) {
      cleaned.pop_back();
    }
    auto p             = eval_position(triangles, positions, end);
    auto [inside, b2f] = point_in_triangle(
        triangles, positions, cleaned.back(), p);
    assert(inside);
    new_end = mesh_point{cleaned.back(), b2f};  // updating end
  }
  std::tie(is_vert, kv) = point_is_vert(start);
  std::tie(is_edge, ke) = point_is_vert(start);

  if (is_vert) {
    auto vid      = triangles[start.face][kv];
    auto curr_tid = strip[start_entry + 1];
    kv            = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.erase(cleaned.begin());
      ++start_entry;
      if (start_entry == end_entry - 1) break;
      // size of the strip must be at least two(see
      // handle_degenerate_case_for_tracing or get_strip)
      auto curr_tid = strip[start_entry + 1];
      kv            = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned[0]], vid);
    assert(kv != -1);
    b3f       = zero3f;
    b3f[kv]   = 1;
    new_start = mesh_point{cleaned[0], vec2f{b3f.y, b3f.z}};  // udpdating start

  } else if (is_edge) {
    if (start.face != strip[0]) {
      assert(adjacencies[start.face][ke] == strip[0]);
      if (start.face == strip[1]) cleaned.erase(cleaned.begin());
    } else if (adjacencies[start.face][ke] == strip[1]) {
      cleaned.erase(cleaned.begin());
    }
    auto p             = eval_position(triangles, positions, start);
    auto [inside, b2f] = point_in_triangle(triangles, positions, cleaned[0], p);
    assert(inside);
    new_start = {cleaned[0], b2f};  // updating start
  }
  return {cleaned, new_start, new_end};
}
