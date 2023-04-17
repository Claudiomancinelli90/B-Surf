#pragma once
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <geometrycentral/surface/vector_heat_method.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>
#include <splinesurf/karcher.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_mesh.h>

#include "../flipout/flipout.h"
using namespace yocto;
using namespace geometrycentral;
using namespace geometrycentral::surface;

vector<mesh_point> bezier_karcher_final(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vec3f>& normals,
    const vector<vector<int>>& v2t, const vector<vector<float>>& angles,
    const vector<float>&                  total_angles,
    const Eigen::SparseMatrix<double, 1>& Grad,
    const vector<mesh_point>& control_points, const int number_of_subdivision,
    const flipout::flipout_mesh& mesh, bool use_vector_heat);

std::tuple<vector<vector<vec3f>>, vector<vector<float>>> compute_karcher_grad(
    const flipout::flipout_mesh& mesh, const vector<int>& control_points,
    const geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const Eigen::SparseMatrix<double, 1>& Grad);