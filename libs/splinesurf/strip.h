#include <deque>

#include "spline.h"

template <typename T>
inline int find_in_vector(const T& vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x) return i;
  return -1;
}

struct strip_arena {
  vector<float> field    = {};
  vector<int>   parents  = {};
  vector<bool>  in_queue = {};
};

void init_arena(strip_arena& arena, size_t size) {
  arena.parents  = vector<int>(size, -1);
  arena.field    = vector<float>(size, flt_max);
  arena.in_queue = vector<bool>(size, false);
}

template <typename Update, typename Stop, typename Exit>
void search_strip(vector<float>& field, vector<bool>& in_queue,
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int start, int end, Update&& update,
    Stop&& stop, Exit&& exit) {
  auto destination_pos = eval_position(
      triangles, positions, {end, {1.0f / 3, 1.0f / 3}});

  auto estimate_dist = [&](int face) {
    auto p = eval_position(triangles, positions, {face, {1.0f / 3, 1.0f / 3}});
    return length(p - destination_pos);
  };
  field[start] = estimate_dist(start);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue      = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

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
    if (exit(node)) break;
    if (stop(node)) continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1) continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

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
      if (update(node, neighbor, new_distance)) return;
    }
  }
}

vector<int> compute_strip(
    const bezier_mesh& mesh, strip_arena& arena, int start, int end) {
  if (start == end) return {start};

  // initialize once for all and sparsely cleanup at the end of every solve
  // auto parents  = vector<int>(mesh.dual_solver.graph.size(), -1);
  // auto field    = vector<float>(mesh.dual_solver.graph.size(), flt_max);
  // auto in_queue = vector<bool>(mesh.dual_solver.graph.size(), false);
  auto visited = vector<int>{start};

  auto sources = vector<int>{start};
  auto update  = [&arena, &visited, end](
                    int node, int neighbor, float new_distance) {
    arena.parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(arena.field, arena.in_queue, mesh.dual_solver, mesh.triangles,
      mesh.positions, start, end, update, stop, exit);
  // visit_dual_graph(field, mesh.dual_solver, mesh.triangles, mesh.positions,
  //     in_queue, start, end, update, visited);

  // extract_strip
  auto strip = vector<int>{};
  auto node  = end;
  assert(arena.parents[end] != -1);
  strip.reserve((int)sqrt(arena.parents.size()));
  while (node != -1) {
    assert(find_in_vector(strip, node) != 1);
    strip.push_back(node);
    node = arena.parents[node];
  }

  // cleanup buffers
  for (auto& v : visited) {
    arena.parents[v]  = -1;
    arena.field[v]    = flt_max;
    arena.in_queue[v] = false;
  }
  // assert(check_strip(mesh.adjacencies, strip));
  return strip;
}

vector<int> compute_strip_tlv(const bezier_mesh& mesh, int start, int end) {
  if (start == end) return {start};

  thread_local static auto parents  = vector<int>{};
  thread_local static auto field    = vector<float>{};
  thread_local static auto in_queue = vector<bool>{};

  if (parents.size() != mesh.dual_solver.graph.size()) {
    parents.assign(mesh.dual_solver.graph.size(), -1);
    field.assign(mesh.dual_solver.graph.size(), flt_max);
    in_queue.assign(mesh.dual_solver.graph.size(), false);
  }

  // initialize once for all and sparsely cleanup at the end of every solve
  auto visited = vector<int>{start};
  auto sources = vector<int>{start};
  auto update  = [&visited, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(field, in_queue, mesh.dual_solver, mesh.triangles,
      mesh.positions, start, end, update, stop, exit);
  // visit_dual_graph(field, mesh.dual_solver, mesh.triangles, mesh.positions,
  //     in_queue, start, end, update, visited);

  // extract_strip
  auto strip = vector<int>{};
  auto node  = end;
  strip.reserve((int)sqrt(parents.size()));
  while (node != -1) {
    assert(find_in_vector(strip, node) != 1);
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  for (auto& v : visited) {
    parents[v]  = -1;
    field[v]    = flt_max;
    in_queue[v] = false;
  }
  // assert(check_strip(mesh.adjacencies, strip));
  return strip;
}
