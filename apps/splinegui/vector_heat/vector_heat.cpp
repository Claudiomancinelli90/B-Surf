#include "vector_heat.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

SurfacePoint exp_map(const flipout::flipout_mesh& mesh,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>&         adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const SurfacePoint& start, const Vector2& dir, const float& len) {
  auto    basisx    = mesh.geometry->faceTangentBasis[start.face];
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
  return SurfacePoint(mesh.topology->face(end.face),
      Vector3{bary_end.x, bary_end.y, bary_end.z});
}
Vector2 karcher_logmap(const VertexData<Vector2>& logmap,
    const vector<Vertex>& points, const vector<float>& weights) {
  Vector2 average = Vector2::zero();
  for (auto i = 0; i < points.size(); ++i) {
    average = average + logmap[points[i]] * weights[i];
  }
  return average;
}
std::vector<mesh_point> weighted_average(const flipout::flipout_mesh& mesh,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>&         adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const std::vector<int>& control_points, const int number_of_subdivision) {
  auto                   n = (int)control_points.size();
  VectorHeatMethodSolver heat_solver(*mesh.geometry);
  auto                   points = std::vector<Vertex>(n);
  vector<mesh_point>     result = {
      point_from_vert(triangles, v2p_adjacencies, control_points[0])};
  mesh.geometry->requireFaceTangentBasis();

  for (auto i = 0; i < n; ++i) {
    points[i] = mesh.topology->vertex(control_points[i]);
  }

  vector<float> w;
  double        step             = 1 / pow(2, number_of_subdivision);
  auto          prev_sample      = SurfacePoint(points[0]).inSomeFace();
  auto          t                = step;
  auto          prev_prev_sample = prev_sample;
  bernstein_polynomials(n, step, w);
  for (int i = 1; i < pow(2, number_of_subdivision); ++i) {
    auto grad  = Vector2{1, 1};
    auto count = 0;
    while (grad.norm() > 1e-5 && count < 100) {
      ++count;
      auto logmap = heat_solver.computeLogMap(prev_sample);
      grad        = karcher_logmap(logmap, points, w);
      prev_sample = exp_map(mesh, solver, triangles, positions, normals,
          adjacencies, v2p_adjacencies, angles, total_angles, prev_sample, grad,
          grad.norm());
    }
    if (count != 100) {
      auto coords = prev_sample.faceCoords;
      result.push_back(mesh_point{(int)prev_sample.face.getIndex(),
          vec2f{(float)coords.y, (float)coords.z}});
      prev_prev_sample = prev_sample;
    } else {
      auto damping = 0.5;
      while (damping > 1e-6 && grad.norm() > 1e-5) {
        auto logmap = heat_solver.computeLogMap(prev_sample);
        grad        = karcher_logmap(logmap, points, w);
        prev_sample = exp_map(mesh, solver, triangles, positions, normals,
            adjacencies, v2p_adjacencies, angles, total_angles, prev_sample,
            grad, damping * grad.norm());

        damping /= 2;
      }
      auto coords = prev_sample.faceCoords;
      result.push_back(mesh_point{(int)prev_sample.face.getIndex(),
          vec2f{(float)coords.y, (float)coords.z}});
      prev_sample = prev_prev_sample;
    }
    t += step;
    bernstein_polynomials(n, t, w);
  }

  return result;
}
Eigen::VectorXd geometry_central_wrapper(
    const VertexData<double>& f, const bool squared = false) {
  Eigen::VectorXd field(f.size());
  for (auto i = 0; i < f.size(); ++i) {
    field(i) = (squared) ? f[i] * f[i] : f[i];
  }
  return field;
}
vector<float> geometry_central_wrapper_yocto(
    const VertexData<double>& f, const bool squared = false) {
  vector<float> field(f.size());
  for (auto i = 0; i < f.size(); ++i) {
    field[i] = (squared) ? f[i] * f[i] : f[i];
  }
  return field;
}
std::tuple<vector<vector<vec3f>>, vector<vector<float>>> compute_karcher_grad(
    const flipout::flipout_mesh& mesh, const vector<int>& control_points,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double, 1>& Grad) {
  HeatMethodDistanceSolver heat_solver(*mesh.geometry);

  vector<vector<vec3f>> gradients(control_points.size());
  vector<vector<float>> fields(control_points.size());
  for (auto i = 0; i < control_points.size(); ++i) {
    auto dist = heat_solver.computeDistance(
        mesh.topology->vertex(control_points[i]));
    fields[i]    = geometry_central_wrapper_yocto(dist);
    auto F       = geometry_central_wrapper(dist, true);
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, F);
  }

  return {gradients, fields};
}
vector<mesh_point> bezier_karcher_final(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const flipout::flipout_mesh& mesh, bool use_vector_heat) {
  vector<int> vertices(control_points.size());
  for (int i = 0; i < 4; i++) {
    vertices[i] = vert_from_point(triangles, control_points[i]);
  }
  if (use_vector_heat)
    return weighted_average(mesh, solver, triangles, positions, normals,
        adjacencies, v2t, angles, total_angles, vertices,
        number_of_subdivision);
  auto [gradients, fields] = compute_karcher_grad(
      mesh, vertices, solver, triangles, positions, normals, Grad);
  return karcher_test(solver, triangles, positions, adjacencies, normals, v2t,
      angles, total_angles, gradients, vertices, number_of_subdivision);
}