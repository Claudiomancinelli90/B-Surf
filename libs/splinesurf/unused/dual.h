
#pragma once

#include <deque>
#include <iostream>
using namespace std;

#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

#include "algorithms.h"

// Local data structures

// Dual solver -- using graph consisting of centroids of tris and dual edges
struct dual_geodesic_solver {
  vector<vec3f> centroids = {};  // for each tri, coords of centroid (nodes)
  vector<vec3i> arcs =
      {};  // for each tri, its adjacent tris (clone of mesh.adjacency)
  vector<vec3f> distances = {};  // for each tri, distances of centroid from
                                 // neighbors' centroids (arc length)
};

// Initialize solver by computing centroids and lengths of dual edges
void init_dual_solver(const Mesh& mesh, dual_geodesic_solver& solver) {
  int        ntri = mesh.triangles.size();
  vec3f      v0, v1, v2;
  Triangle2D t, tn;
  solver.arcs = mesh.adjacency;
  solver.centroids.resize(ntri);
  solver.distances.resize(ntri);
  for (auto i = 0; i < ntri; i++) {
    v0                 = mesh.positions[mesh.triangles[i][0]];
    v1                 = mesh.positions[mesh.triangles[i][1]];
    v2                 = mesh.positions[mesh.triangles[i][2]];
    solver.centroid[i] = (v0 + v1 + v2) / 3.0;
    t                  = init_flat_triangle(mesh.positions, mesh.triangles[i]);
    for (auto j = 0; i < 3; i++) {
      if (solver.arcs[i][j] == -1) continue;  // skip false arcs at boundary
      tn = unfold_face(mesh.triangles, mesh.positions, mesh.adjacency, t,
          mesh.triangles[i], j);
      solver.distances[i][j] = length(
          (t.x + t.y + t.z - tn.x - tn.y - tn.z) / 3.0);
    }
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Basic solver (private)
// sources are centroids, field reports distances from sources per triangle (per
// centroid) Needs an initialized field (I/O parameter)

template <typename Update, typename Stop, typename Exit>
void visit_dual_graph(vector<float>& field, const dual_geodesic_solver& solver,
    const vector<int>& sources, Update&& update, Stop&& stop, Exit&& exit) {
  auto in_queue = vector<bool>(solver.centroids.size(), false);

  // Setup queue
  // Cumulative weights of elements in queue used to keep track of the average
  // weight of the queue.
  double cumulative_weight = 0.0;
  auto   queue             = deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  // Main loop
  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)cumulative_weight / queue.size();

    // Large Label Last heuristics
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

    // Check early exit and stop propagation conditions.
    if (exit(node)) break;
    if (stop(node)) continue;

    // Process neighbors of current node
    for (auto i = 0; i < 3; i++) {
      if (solver.arcs[node][i] == -1) continue;  // skip false arcs at boundary
      // Distance of neighbor through this node
      auto new_distance = field[node] + solver.distances[node][i];
      auto neighbor     = solver.arcs[node][i];

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label First
        // heuristics
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
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Public solvers

// From a set of sources (centroids) to all centroids -- optional max_distance
// to stop propagation return distance field

void compute_dual_distances(const dual_geodesic_solver& solver,
    const vector<int>& sources, vector<float>& distances,
    float max_distance = flt_max) {
  distances = vector<float>(solver.centroids.size(), flt_max);
  for (auto source : sources) distances[source] = 0.0f;
  auto update = [](int node, int neighbor, float new_distance) {};
  auto stop   = [&](int node) { return distances[node] > max_distance; };
  auto exit   = [](int node) { return false; };
  visit_dual_graph(distances, solver, sources, update, stop, exit);
}

// From one sources to one target (both centroids) --- with early exit
// return discrete path from target to source (including extremes)
vector<int> point_to_point_geodesic_path(
    const dual_geodesic_solver& solver, int source, int target) {
  vector<float> distances.assign(solver.centroids.size(), flt_max);
  vector<int>   sources = {source};
  vector<int>   parents.assign(solver.centroids.size(), -1);
  auto update = [&parents](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
  };
  auto stop = [](int node) { return false; };
  auto exit = [&](int node) { return node == target };
  visit_dual_graph(distances, solver, sources, update, stop, exit);
  vector<int> discrete_path;
  int         p = parents[target];
  discrete_path.push_back(target);
  while (p != source) discrete_path.push_back(p);
  discrete_path.push_back(target);
  return discrete_path;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
