#pragma once
#include <VTP/geodesic_algorithm_exact.h>
#include <VTP/geodesic_mesh.h>
#include <geometry-central/include/geometrycentral/surface/meshio.h>
#include <geometry-central/include/geometrycentral/surface/vector_heat_method.h>
#include <stdio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
using namespace yocto;

enum transport_mode { V2V, V2T, T2V, T2T };
enum simplex { vert, edge };
enum solutions { non_irrot, infty };
enum field { exact, graph };
enum badones { control, intermediate };

struct splits {
  vector<vector<mesh_point>>            control_polygons = {};
  vector<vector<pair<mesh_point, int>>> badones          = {};
  vector<bool>                          goodones         = {};
  vector<vector<mesh_point>>            polyline         = {};
};

struct bezier_splits {
  vector<vector<mesh_point>>    control_polygons = {};
  vector<vector<vector<float>>> distance_fields  = {};
  vector<vector<vector<vec3f>>> gradient_fields  = {};
  vector<int>                   subidivisions    = {};
};

// function used in core.pp
vec3f tid_centroid(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int pid);

vec3f tri_bary_coords(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& p);

vec3f tri_bary_coords(
    const vec2f& v0, const vec2f& v1, const vec2f& v2, const vec2f& p);

vec3f vector_bary_coords(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& v);

vec3f vector_bary_coords(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& v);

vec3f vector_bary_coords(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec2f& v);

vec3f vector_bary_coords(const unfold_triangle& flat_tid, const vec2f& v);

vec3f project_vec(const vec3f& v, const vec3f& n);

int vert_from_point(const vector<vec3i>& triangles, const mesh_point& p);

mesh_point        point_from_vert(const vector<vec3i>& triangles,
           const vector<vector<int>>& v2t, const int vid, const int tid = -1);
pair<bool, vec2f> point_in_unfold_triangle(
    const vec2f& pos, const unfold_triangle& tr, float tol = 1e-2);

vec3f from_2d_to_3d_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec2f& v, const int tid);

vec2f rot_vect(const vec2f& p, const float theta);

mesh_point make_mesh_point(const vector<vec3i>& triangles,
    const vector<vector<int>>& v2t, const int vid, const int tid = -1);

std::tuple<vector<vec3f>, vector<vec3i>, unordered_map<int, int>>
add_vertices_to_mesh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<mesh_point>& points);

Eigen::SparseMatrix<double, 1> init_gradient_matrix(
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vec3f>& normals);

vec3f tranport_vector(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vec3f& v, const vector<vec3f>& normals, const mesh_point& from,
    const mesh_point& to);

std::tuple<vector<int>, mesh_point, mesh_point> cleaned_strip(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<int>& strip,
    const mesh_point& start, const mesh_point& end);

void trace_in_triangles(const vector<vec3f>& positions,
    const vector<vec3i>& triangles, const vec3f& dir, const vec3f& bary,
    const int pid, vec3f& sample_pos, vec3f& sample_bary);

int next_tid(const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vector<int>>& v2t,
    const vector<vec3i>& triangles, const vector<vec3f>& normals,
    const int from, const vec3f& v);

vector<mesh_point> straightest_geodesic(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec3i>& adjacencies,
    const vector<vector<int>>&   v2p_adjacencies,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const mesh_point& from, const vec3f& v, const float& l);

void bernstein_polynomials(const int& n, const float& t, vector<float>& w);
vector<vec3f>      karcher_grad(const vector<vector<vec3f>>& gradients,
         const vector<vector<float>>& f, const vector<float>& weights,
         const vector<vec3f>& normals, bool to_draw = false);
vector<mesh_point> bezier_karcher(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, splits& composite_bezier,
    const bool all_in_once = true, const vec2i& knot = {-1, -1});

vector<mesh_point> bezier_karcher(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision);

vector<mesh_point> karcher_test(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    vector<int>& jumps);

mesh_point weighted_average(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const float& t);

vector<mesh_point> bezier_karcher_test(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<vector<vec3f>>& directions,
    vector<vector<mesh_point>>& steps, vector<pair<mesh_point, int>>& badones);

vector<mesh_point> bezier_karcher_test(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    pair<mesh_point, int>& badone);

vector<mesh_point> bezier_karcher_tet(const geodesic_solver& solver,
    const dual_geodesic_solver& dual_solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision);

vector<mesh_point> karcher_test(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<int>& control_points, const int number_of_subdivision);

vector<mesh_point> bezier_karcher_tet_subdivision(const geodesic_solver& solver,
    const dual_geodesic_solver& dual_solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision);

vector<mesh_point> bezier_karcher_almost_gradients(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const int number_of_subdivision, vector<pair<int, float>>& verts,
    vector<vector<int>>& tids, vector<pair<mesh_point, float>>& badones);

std::tuple<vector<mesh_point>, vector<vector<float>>, vector<vector<vec3f>>>
split_control_poligon_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const Eigen::SparseMatrix<double, 1>& Grad, const bezier_splits& curve,
    const std::unordered_map<int, int>& subdivision, const int entry,
    const float& t = 0.5);

mesh_point weighted_average_tet(const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& L01,
    const geodesic_path& L02, const geodesic_path& L03,
    const vector<float>& L01t, const vector<float>& L02t,
    const vector<float>& L03t, const vector<float>& weights);

vector<vector<float>> control_points_distance_field(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<mesh_point>& control_points,
    const int type_of_field);

vector<vector<vec3f>> control_points_gradients(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<vector<float>>& f);

vector<mesh_point> solve_for_control_points(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& gradients,
    const vector<float>& t, const float& step,
    vector<pair<mesh_point, int>>& badones);

vector<mesh_point> solve_for_intermediate_points(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<mesh_point>& control_points,
    const vector<vector<float>>& f, const vector<vector<vec3f>>& gradients,
    const vector<float>& t, const float& step,
    vector<pair<mesh_point, int>>& badones);

vector<mesh_point> bezier_karcher_bisection(
    const dual_geodesic_solver& dual_solver, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<vector<mesh_point>>& final_control_polygons,
    vector<vector<mesh_point>>& bad_control_polygons,
    vector<int>&                type_of_badones);

vector<mesh_point> bezier_karcher_bisection_hybrid(
    const geodesic_solver& solver, const dual_geodesic_solver& dual_solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field, vector<vector<mesh_point>>& final_control_polygons,
    vector<vector<mesh_point>>& bad_control_polygons,
    vector<int>&                type_of_badones);

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
    const int type_of_field);

vector<mesh_point> bezier_karcher_full_gradient(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>&             control_points);

pair<mesh_point, int> sample_on_curve(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, const vector<float>& path_parameter_t,
    const float& t);

mesh_point eval_geodesic_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, const vector<float>& path_parameter_t,
    const float& t);

Eigen::SparseMatrix<double, 1> k_ring_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, geodesic_solver& solver,
    vector<vector<float>>& angles, vector<float>& total_angles);

Eigen::SparseMatrix<double, 1> spanning_direction_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, geodesic_solver& solver,
    vector<vector<float>>& angles, vector<float>& total_angles);
geodesic_solver spanning_direction_geodesic_solver(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int k, vector<vector<float>>& angles,
    vector<float>& total_angles);
int next_tid_for_tracing(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const int tid, const vec3f& last_bary,
    vec3f& v, vec3f& updated_bary);

mesh_point transport_vector(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const geodesic_path& path,
    const mesh_point& point);

vector<float> coefficients_from_control_points(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vector<vec3f>>& gradients, const vector<vector<float>>& f,
    const vector<mesh_point>& control_points, const mesh_point& point);

vector<float> exact_geodesic_distance(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vector<int>>& v2t,
    const mesh_point& source);

vector<float> compute_distance_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const int entry, const vector<mesh_point>& control_points,
    const int type_of_field);

vec3f gradient_blending(const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const int vid);

vec3f gradient_blending(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3f>& positions, const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& normals, const int& vid,
    const vector<vector<vec3f>>& gradients, const vector<float>& weights);

float field_blending(const vector<vector<float>>& f,
    const vector<float>& weights, const int vid);

vector<vec3f> subdivide_tid(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const int vid, const int vid0, const int vid1);

Eigen::SparseMatrix<double, 1> init_riemannian_gradient_matrix(
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<vec3f>& positions, const vector<vec3f>& normals);
Eigen::SparseMatrix<double, 1> PCE_grad_mat(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
vector<vec3f> compute_riemannian_gradient(const geodesic_solver& solver,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double>& G, const Eigen::VectorXd& f,
    bool normalized = false);

vec3f AGS_gradient(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3f>& positions, const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vec3f>& gradient, const int& vid,
    const vector<vec3f>& normals);

vector<vec3f> compute_grad(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& G,
    const Eigen::VectorXd& f, bool normalized = false);

vector<vec3f> compute_grad(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const Eigen::SparseMatrix<double, 1>& G,
    const vector<float>& f, bool normalized = false);
bool bary_for_minimum_energy(const vec3f& i, const vec3f& j, const vec3f& k,
    float& alpha, float& beta, float& gamma);

bool bary_for_minimum_energy(const vec2f& i, const vec2f& j, const vec2f& k,
    float& alpha, float& beta, float& gamma);
pair<bool, int> bary_for_minimum_hessian(const vec3f& i, const vec3f& j,
    const vec3f& k, float& alpha, float& beta, float& gamma);
void            bumped_sphere(const float& treshold, vector<vec3f>& positions);

Eigen::VectorXf solve_quadratic_field(
    const vec3f& G0, const vec3f& G1, const vec3f& G2);

Eigen::VectorXf solve_quadratic_field(const vec3f& G0, const vec3f& G1,
    const vec3f& G2, const float& F1, const float& F2, const float& F3);
vec3f           irrot_field_singularity(const Eigen::VectorXf& coeff);
vector<vec3f>   non_irrot_field(const vector<vec3i>& triangles,
      const vector<vec3f>& positions, const int tid, const vec3f& V0,
      const vec3f& V1, const vec3f& V2, vector<vec3f>& bary);

std::pair<vec3f, bool> maximum_descent_direction(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const int tid, const vec3f& grad,
    const vec3f& gx, const vec3f& gy, const vec3f& gz);

std::pair<vec3f, bool> maximum_descent_direction(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& grad, const vec3f& gx, const vec3f& gy,
    const vec3f& gz);

vec3f almost_gradient(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vec3f>& normals, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<float>& weights,
    const mesh_point& p, const vector<mesh_point>& control_points);

mesh_point step_for_gradient_descent_minimization(const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<float>& total_angles,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const mesh_point& from, vec3f& dir);

vec3f compute_tangent_vector(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path, const int entry);

vector<int> strip_by_gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vec3f>& gradients,
    const mesh_point& from, vector<vec3f>& directions,
    vector<mesh_point>& steps, vector<pair<mesh_point, int>>& badones);

vector<int> get_strip_having_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const vector<float>& distances,
    const mesh_point& source, const mesh_point& target, vector<int>& parents);
std::tuple<vector<int>, vector<float>> get_strip_with_distances(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const mesh_point& source, const mesh_point& target);
std::tuple<vector<mesh_point>, vector<float>, vector<int>>
sort_points_and_weigts(const vector<mesh_point>& control_points,
    const vector<float>& weights, const float& t);

std::tuple<geodesic_path, vector<float>> set_gamma01(
    const vector<geodesic_path>& paths, const vector<vector<float>>& paths_t,
    const vector<int>& order);

vector<mesh_point> bezier_dyn(const dual_geodesic_solver& dual_solver,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const int type_of_field);

vector<float> compute_pruned_geodesic_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const mesh_point& source, const vector<mesh_point>& targets);

vector<mesh_point> spline_control_points(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles,
    const vector<mesh_point>&    control_points);

vector<vector<float>> knots_insertion(const vector<float>& t,
    const vector<float>& u, const int k, const int m, const int n);

mesh_point optimal_seed(
    const vector<mesh_point>& control_points, const vector<float>& weights);

mesh_point optimal_seed(const vector<mesh_point>& points,
    const vector<float>& point_weights, const vector<vector<float>>& weights);

mesh_point gradient_descent(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vec3f>& normals,
    const geodesic_solver& solver, const vector<vector<float>>& angles,
    const vector<float>& total_angles, const vector<vector<vec3f>>& gradients,
    const vector<float>& weights, const mesh_point& from);

Eigen::VectorXd wrapper(const vector<float>& f);

vector<mesh_point> bezier_through_spline_subdivisions(
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<mesh_point>& control_points);

vector<vector<vec3f>> gradients_inside_cell(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, const geodesic_solver& solver,
    const vector<vector<float>>& angles, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<vector<float>>&          fields,
    const vector<mesh_point>& control_points, const mesh_point& new_source,
    const vector<float>& source_weights, vector<float>& target_weights);

std::pair<mesh_point, vector<float>> optimal_seed_with_weights(
    const vector<mesh_point>& points, const vector<float>& point_weights,
    const vector<vector<float>>& weights);

inline vec3f get_bary(const vec2f& uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}