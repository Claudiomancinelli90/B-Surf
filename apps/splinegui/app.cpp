#include "app.h"

#include <splinesurf/logging.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_modelio.h>
#include <yocto/yocto_shape.h>

void update_control_points_shape(Gui_Spline& spline, const bezier_mesh& mesh) {
  auto positions = vector<vec3f>(spline.control_points.size());
  // We are going to draw only anchor points, not tangents (hence i += 3).
  for (int i = 0; i < positions.size(); i += 3) {
    assert(check_point(spline.control_points[i]));
    positions[i] = eval_position(mesh, spline.control_points[i]);
  }
  set_points_shape(spline.anchor_points_shape, positions);
}
void update_glvector_field(App& app, const vector<vec3f>& vector_field,
    const float& scale, const string& name) {
  gpu::delete_shape(app.gpu_shapes[name]);
  auto& triangles = app.mesh.triangles;
  auto& positions = app.mesh.positions;
  if (vector_field.size() == triangles.size()) {
    app.gpu_shapes[name] = gpu::make_vector_field_shape(
        vector_field, triangles, positions, scale);
  } else {
    app.gpu_shapes[name] = gpu::make_vector_field_shape(
        vector_field, positions, scale);
  }
}
void de_casteljau(Gui_Curve& curve, const array<mesh_point, 4>& control_points,
    const bezier_mesh& mesh, const bezier_params params) {
  auto& tree = curve.tree;
  tree.nodes.clear();
  tree.nodes.push_back({});
  tree.nodes[0].points = control_points;

  // set lines
  tree.nodes[0].lines[0] = compute_geodesic_path(
      mesh, control_points[0], control_points[1]);
  tree.nodes[0].lines[1] = compute_geodesic_path(
      mesh, control_points[1], control_points[2]);
  // This is reversed!!!
  tree.nodes[0].lines[2] = compute_geodesic_path(
      mesh, control_points[3], control_points[2]);

  subdivide_bezier_tree(mesh, tree, params);
  int  depth         = params.subdivisions;
  auto bezier_points = bezier_tree_points(mesh, tree, depth);
  curve.positions    = make_polyline_positions(mesh, bezier_points);

  auto& t0     = curve.tangents[0];
  t0.path      = tree.nodes[0].lines[0];  // avoidable copy
  t0.positions = path_positions(mesh, t0.path);
  t0.length    = path_length(t0.positions);

  auto& t1     = curve.tangents[1];
  t1.path      = tree.nodes[0].lines[2];  // avoidable copy
  t1.positions = path_positions(mesh, t1.path);
  t1.length    = path_length(t1.positions);
}

void update_curve_shape(Gui_Spline& spline, const int curve_id,
    const bezier_mesh& mesh, const bezier_params& params) {
  if (curve_id < 0 || curve_id >= spline.curves.size()) return;
  auto& curve   = spline.curves[curve_id];
  auto  polygon = get_control_polygon(spline, curve_id);

  auto& t0     = curve.tangents[0];
  t0.path      = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  t0.positions = path_positions(mesh, t0.path);
  t0.length    = path_length(t0.positions);

  auto& t1     = curve.tangents[1];
  t1.path      = compute_geodesic_path(mesh, polygon[3], polygon[2]);
  t1.positions = path_positions(mesh, t1.path);
  t1.length    = path_length(t1.positions);
  auto points  = vector<mesh_point>{};

  if (params.algorithm == spline_algorithm::de_casteljau_uniform) {
    points = bezier_uniform(mesh, polygon, params);
    printf("num points: %d\n", points.size());
  } else if (params.algorithm == spline_algorithm::de_casteljau_adaptive) {
    points = bezier_adaptive(mesh, polygon, params);
    printf("num points: %d\n", points.size());
  } else if (params.algorithm == spline_algorithm::subdivision_uniform) {
    points = spline_subdivision_uniform(mesh, polygon, params.subdivisions);
    printf("num points: %d\n", points.size());
  } else if (params.algorithm == spline_algorithm::subdivision_adaptive) {
    points = spline_subdivision_adaptive(mesh, polygon, params);
    printf("num points: %d\n", points.size());
  }

  curve.positions = make_polyline_positions(mesh, points);

  if (curve.positions.size() == 0)
    curve.positions = {eval_position(mesh, points[0])};
}
void update_curve_shape(Gui_Spline& spline, int curve_id,
    const bezier_mesh& mesh, const bezier_params& params, const bool quadric,
    const bool use_vector_heat) {
  if (curve_id < 0 || curve_id >= spline.curves.size()) return;
  auto& curve           = spline.curves[curve_id];
  auto  polygon         = get_control_polygon(spline, curve_id);
  auto  quadric_polygon = quadratic_bezier_segment{};
  if (quadric) quadric_polygon = {polygon[0], polygon[1], polygon[3]};

  auto& t0     = curve.tangents[0];
  t0.path      = compute_geodesic_path(mesh, polygon[0], polygon[1]);
  t0.positions = path_positions(mesh, t0.path);
  t0.length    = path_length(t0.positions);

  auto& t1     = curve.tangents[1];
  t1.path      = compute_geodesic_path(mesh, polygon[3], polygon[2]);
  t1.positions = path_positions(mesh, t1.path);
  t1.length    = path_length(t1.positions);

  if (params.algorithm == spline_algorithm::flipout) {
    auto vertices = vector<int>(4);
    for (int i = 0; i < 4; i++) {
      vertices[i] = vert_from_point(
          mesh.triangles, polygon[i]);  // mesh.triangles[polygon[i].face].x;
    }
    try {
      auto bezier = flipout::compute_bezier_curve(
          *mesh.flipout, vertices, params.subdivisions);
      curve.positions = flipout::path_positions(bezier.get());
    } catch (std::exception& e) {
      printf("[flipout exception]: %s\n", e.what());
    }
    return;
  }

  // bezier_segment polygon;
  // polygon[0] = {1404, {0.094313, 0.87842}};
  // polygon[1] = {809, {0.079695, 0.13082}};
  // polygon[2] = {1370, {0.457860, 0.18341}};
  // polygon[3] = {1493, {0.002444, 0.58035}};
  // ellipsoid broken curve
  // polygon[0] = {16169, {0.189351, 0.558249}};
  // polygon[1] = {52623, {0.645337, 0.112507}};
  // polygon[2] = {47488, {0.490417, 0.230988}};
  // polygon[3] = {11104, {0.074401, 0.270324}};
  // bugged subdivision adaptive, bunny_70k
  // polygon[0] = {29322, {0.517030, 0.058723}};
  // polygon[1] = {9519, {0.001291, 0.779000}};
  // polygon[2] = {51630, {0.252770, 0.27111}};
  // polygon[3] = {23768, {0.60169, 0.15421}};
  // fig4 pole
  // polygon[0] = {22504, {0.056851, 0.083361}};
  // polygon[1] = {20227, {0.220981, 0.383081}};
  // polygon[2] = {17510, {0.706649, 0.23103}};
  // polygon[3] = {19778, {0.082177, 0.446810}};
  // bad karcher bunny 70k
  // polygon[0] = {20567, {0.271632, 0.18985}};
  // polygon[1] = {36458, {0.412200, 0.535899}};
  // polygon[2] = {18009, {0.627066, 0.145965}};
  // polygon[3] = {19284, {0.320517, 0.506506}};
  // if (params.algorithm <= spline_algorithm::de_casteljau_adaptive) {
  //   de_casteljau(curve, polygon, mesh, params);
  // }
  curve.samples.clear();
  if (params.algorithm == spline_algorithm::de_casteljau_uniform) {
    curve.samples = (quadric) ? bezier_uniform(mesh, quadric_polygon, params)
                              : bezier_uniform(mesh, polygon, params);
    printf("num curve.samples: %d\n", curve.samples.size());
  } else if (params.algorithm == spline_algorithm::de_casteljau_adaptive) {
    curve.samples = bezier_adaptive(mesh, polygon, params);
    printf("num curve.samples: %d\n", curve.samples.size());
  } else if (params.algorithm == spline_algorithm::subdivision_uniform) {
    curve.samples = (quadric) ? spline_subdivision_uniform(
                                    mesh, quadric_polygon, params.subdivisions)
                              : spline_subdivision_uniform(
                                    mesh, polygon, params.subdivisions);
    printf("num curve.samples: %d\n", curve.samples.size());
  } else if (params.algorithm == spline_algorithm::subdivision_adaptive) {
    curve.samples = spline_subdivision_adaptive(mesh, polygon, params);
    printf("num curve.samples: %d\n", curve.samples.size());
  } else if (params.algorithm == spline_algorithm::de_casteljau_classic) {
    vector<int> badones;
    curve.samples = de_casteljau_classic(mesh, polygon, params, badones);

    // printf("num curve.samples: %d\n", curve.samples.size());
    curve.positions.clear();
    // printf("num badones: %d\n", badones.size());
    bool prev_is_good = true;
    curve.positions.clear();
    curve.sampling_rate.resize(curve.samples.size());
    curve.offsets.resize(curve.samples.size());
    for (auto i = 1; i < curve.samples.size(); ++i) {
      auto pos = make_polyline_positions(
          mesh, {curve.samples[i - 1], curve.samples[i]});
      curve.sampling_rate[i - 1] = path_length(pos);
      curve.offsets[i - 1].x     = (int)curve.positions.size();
      curve.positions.insert(curve.positions.end(), pos.begin(), pos.end());
      curve.offsets[i - 1].y = (int)curve.positions.size();
      //   auto curr_is_good =
      //       (std::find(badones.begin(), badones.end(), i) == badones.end());
      //   if (curr_is_good && prev_is_good) {
      //     auto pos = make_polyline_positions(
      //         mesh, {curve.samples[i - 1], curve.samples[i]});
      //     curve.positions.insert(curve.positions.end(), pos.begin(),
      //     pos.end()); prev_is_good = true;
      //   } else if (!curr_is_good) {
      //     prev_is_good = false;
      //     // auto n       = triangle_normal(
      //     //     mesh.positions[mesh.triangles[curve.samples[i - 1].face].x],
      //     //     mesh.positions[mesh.triangles[curve.samples[i - 1].face].y],
      //     //     mesh.positions[mesh.triangles[curve.samples[i -
      //       1].face].z]);
      //       // auto curr_pos = eval_position(mesh, curve.samples[i - 1]) -
      //       0.01 *n;
      //       // curve.positions.push_back(curr_pos);
      //       // n =
      //       triangle_normal(mesh.positions[mesh.triangles[curve.samples[i].face].x],
      //       //     mesh.positions[mesh.triangles[curve.samples[i].face].y],
      //       //     mesh.positions[mesh.triangles[curve.samples[i].face].z]);
      //       auto curr_pos = eval_position(mesh, curve.samples[i]);
      //       // curve.positions.push_back(zero3f);
      //       curve.positions.push_back(curr_pos);

      //   } else {
      //     prev_is_good = true;
      //     // auto n       = triangle_normal(
      //     //     mesh.positions[mesh.triangles[curve.samples[i].face].x],
      //     //     mesh.positions[mesh.triangles[curve.samples[i].face].y],
      //     //     mesh.positions[mesh.triangles[curve.samples[i].face].z]);
      //     auto curr_pos = eval_position(mesh, curve.samples[i]);
      //     curve.positions.push_back(curr_pos);
      //   }
    }

  } else if (params.algorithm == spline_algorithm::karcher) {
    // curve.samples = bezier_karcher_test(mesh, polygon, params.subdivisions,
    // badones); curve.positions.clear(); auto vertices = vector<int>(4); for
    // (int i = 0; i < 4; i++) {
    //   vertices[i] = vert_from_point(
    //       mesh.triangles, polygon[i]);  // mesh.triangles[polygon[i].face].x;
    // }
    curve.samples = karcher_average(mesh,
        {polygon[0], polygon[1], polygon[2], polygon[3]}, params.subdivisions,
        use_vector_heat);
    curve.positions.clear();
    curve.sampling_rate.resize(curve.samples.size());
    curve.offsets.resize(curve.samples.size());
    for (auto i = 1; i < curve.samples.size(); ++i) {
      auto pos = make_polyline_positions(
          mesh, {curve.samples[i - 1], curve.samples[i]});
      curve.sampling_rate[i - 1] = path_length(pos);
      curve.offsets[i - 1].x     = (int)curve.positions.size();
      curve.positions.insert(curve.positions.end(), pos.begin(), pos.end());
      curve.offsets[i - 1].y = (int)curve.positions.size();
    }

    auto pos = make_polyline_positions(
        mesh, {curve.samples.back(), polygon[3]});
    curve.sampling_rate.back() = path_length(pos);
    curve.offsets.back().x     = (int)curve.positions.size();
    curve.positions.insert(curve.positions.end(), pos.begin(), pos.end());
    curve.offsets.back().y = (int)curve.positions.size();

    printf("num curve.samples: %d\n", curve.samples.size());
  }
  if (params.algorithm != spline_algorithm::karcher &&
      params.algorithm != spline_algorithm::de_casteljau_classic &&
      params.algorithm != spline_algorithm::flipout)
    curve.positions = make_polyline_positions(mesh, curve.samples);

  if (curve.positions.size() == 0)
    curve.positions = {eval_position(mesh, curve.samples[0])};
}

void update_all_splines(App& app) {
  for (auto& spline : app.splines()) {
    spline.curves.resize(spline.control_points.size() / 3);
    spline.curves_to_update = make_set(spline.curves.size());
  }
}

#include <yocto/yocto_sampling.h>
vector<vec4f> bake_ambient_occlusion(App& app, int num_samples) {
  auto result = vector<vec4f>(app.mesh.positions.size(), zero4f);
  auto rng    = rng_state{};

  auto f = [&](int i) {
    auto frame = basis_fromz(app.mesh.normals[i]);
    for (int sample = 0; sample < num_samples; sample++) {
      auto dir  = sample_hemisphere_cos(rand2f(rng));
      dir       = transform_direction(frame, dir);
      auto ray  = ray3f{app.mesh.positions[i], dir};
      auto isec = intersect_triangles_bvh(
          app.bvh, app.mesh.triangles, app.mesh.positions, ray);
      if (!isec.hit && dir.y > 0) {
        result[i] += vec4f{1, 1, 1, 1};
      } else {
        result[i] += vec4f{0, 0, 0, 1};
      }
    }
  };
  parallel_for((int)app.mesh.positions.size(), f);

  for (auto& r : result) r /= r.w;
  return result;
}

// TODO(giacomo): rename after removing realtime lib
// TODO(giacomo): consider adding these functions to yocto_gui
shade_camera _make_lookat_camera(
    const vec3f& from, const vec3f& to, const vec3f& up = {0, 1, 0}) {
  auto camera  = shade_camera{};
  camera.frame = lookat_frame(from, to, {0, 1, 0});
  camera.focus = length(from - to);
  return camera;
}

shade_camera _make_framing_camera(const vector<vec3f>& positions) {
  auto direction = vec3f{0, 1, 2};
  auto box       = bbox3f{};
  for (auto& p : positions) {
    expand(box, p);
  }
  auto box_center = center(box);
  auto box_size   = max(size(box));
  return _make_lookat_camera(direction * box_size + box_center, box_center);
}

void init_camera(App& app, const vec3f& from, const vec3f& to) {
  *app.camera      = _make_framing_camera(app.mesh.positions);
  app.camera_focus = app.camera->focus;
}

void set_points_shape(ogl_shape* shape, const vector<vec3f>& positions) {
  auto sphere = make_uvsphere({16, 16}, 1, {1, 1});
  if (positions.empty()) return;

  set_vertex_buffer(shape, sphere.positions, 0);
  set_index_buffer(shape, quads_to_triangles(sphere.quads));

  set_vertex_buffer(shape, positions, 1);
  set_instance_buffer(shape, 1, true);
}
void set_points_shape(ogl_shape* shape, const bezier_mesh& mesh,
    const vector<mesh_point>& points) {
  auto pos = vector<vec3f>(points.size());
  for (int i = 0; i < points.size(); i++) {
    pos[i] = eval_position(mesh, points[i]);
  }
  set_points_shape(shape, pos);
}

void set_polyline_shape(ogl_shape* shape, const vector<vec3f>& positions) {
  if (positions.empty()) return;
  auto cylinder = make_uvcylinder({16, 1, 1}, {1, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_vertex_buffer(shape, cylinder.positions, 0);
  set_vertex_buffer(shape, cylinder.normals, 1);
  set_index_buffer(shape, quads_to_triangles(cylinder.quads));

  auto froms  = vector<vec3f>();
  auto tos    = vector<vec3f>();
  auto colors = vector<vec3f>();
  froms.reserve(positions.size() - 1);
  tos.reserve(positions.size() - 1);
  // colors.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size() - 1; i++) {
    if (positions[i] == positions[i + 1]) continue;
    froms.push_back(positions[i]);
    tos.push_back(positions[i + 1]);
    // colors.push_back({sinf(i) * 0.5f + 0.5f, cosf(i) * 0.5f + 0.5f, 0});
  }

  // TODO(giacomo): solve rendering bug with degenerate polyline
  if (froms.empty()) {
    shape->num_instances = 0;
    set_vertex_buffer(shape, {}, 0);
    set_vertex_buffer(shape, {}, 1);
    set_vertex_buffer(shape, {}, 2);
    set_vertex_buffer(shape, {}, 3);
    set_index_buffer(shape, vector<vec3i>{});
  } else {
    set_vertex_buffer(shape, froms, 2);
    set_instance_buffer(shape, 2, true);
    set_vertex_buffer(shape, tos, 3);
    set_instance_buffer(shape, 3, true);
  }
}

void set_mesh_shape(ogl_shape* shape, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals) {
  set_vertex_buffer(shape, positions, 0);
  set_vertex_buffer(shape, normals, 1);
  set_index_buffer(shape, triangles);
}

void _set_polyline_shape(ogl_shape* shape, const vector<vec3f>& positions,
    const vector<vec3f>& normals) {
  set_vertex_buffer(shape, positions, 0);
  if (normals.size()) {
    set_vertex_buffer(shape, normals, 1);
  }
  shape->elements = ogl_element_type::lines;
}

vector<vec3f> make_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = zero3f;
  for (auto& t : triangles) {
    auto normal = cross(
        positions[t.y] - positions[t.x], positions[t.z] - positions[t.x]);
    normals[t.x] += normal;
    normals[t.y] += normal;
    normals[t.z] += normal;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

void init_bvh(App& app) {
  app.bvh = make_triangles_bvh(app.mesh.triangles, app.mesh.positions, {});
}

ray3f camera_ray(const App& app, vec2f mouse) {
  auto camera_ray = [](const frame3f& frame, float lens, const vec2f& film,
                        const vec2f& image_uv) {
    auto e = zero3f;
    auto q = vec3f{
        film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), lens};
    auto q1  = -q;
    auto d   = normalize(q1 - e);
    auto ray = ray3f{transform_point(frame, e), transform_direction(frame, d)};
    return ray;
  };

  mouse += 1;
  mouse /= 2;
  mouse.y      = 1 - mouse.y;
  auto& camera = *app.camera;
  return camera_ray(camera.frame, camera.lens,
      {camera.film, camera.film / camera.aspect}, mouse);
}

vec2f screenspace_from_worldspace(App& app, const vec3f& position) {
  auto [x, y, z] = position;
  auto uv4f      = app.matrices.projection_view * vec4f{x, y, z, 1};
  return vec2f{uv4f.x / uv4f.w, uv4f.y / uv4f.w};
};

mesh_point intersect_mesh(const App& app, vec2f mouse) {
  auto ray  = camera_ray(app, mouse);
  auto isec = intersect_triangles_bvh(
      app.bvh, app.mesh.triangles, app.mesh.positions, ray);

  if (isec.hit) {
    return mesh_point{isec.element, isec.uv};
  } else {
    return {-1, {0, 0}};
  }
}

bool load_program(ogl_program* program, const string& vertex_filename,
    const string& fragment_filename) {
  auto error           = ""s;
  auto vertex_source   = ""s;
  auto fragment_source = ""s;

  if (!load_text(vertex_filename, vertex_source, error)) {
    printf("error loading vertex shader (%s): \n%s\n", vertex_filename.c_str(),
        error.c_str());
    return false;
  }
  if (!load_text(fragment_filename, fragment_source, error)) {
    printf("error loading fragment shader (%s): \n%s\n",
        fragment_filename.c_str(), error.c_str());
    return false;
  }

  auto error_buf = ""s;
  if (!set_program(program, vertex_source, fragment_source, error, error_buf)) {
    printf("\nerror: %s\n", error.c_str());
    printf("    %s\n", error_buf.c_str());
    return false;
  }
  return true;
}

void init_gpu(App& app, bool envlight) {
  string error;
  init_ogl(error);

  app.scene = new shade_scene{};

  app.camera = add_camera(app.scene);
  // app.scene->cameras.push_back(&app.camera);
  app.spline_material = add_material(
      app.scene, {0, 0, 0}, {1, 0, 0}, 1, 0, 0.4);

  init_camera(app);

  // Init opengl mesh
  auto& mesh = app.mesh;
  // auto  colors     = bake_ambient_occlusion(app, 1024);
  app.mesh_shape = add_shape(app.scene, {}, {}, mesh.triangles, {},
      mesh.positions, mesh.normals, {}, {});
  // app.mesh_material = add_material(
  //     app.scene, {0, 0, 0}, {.113, .309, .164}, 0.04, 0, 0.5);
  // app.mesh_material = add_material(
  //     app.scene, {0, 0, 0}, {0.183, 0.5, 0.265}, 0.04, 0, 0.4);
  app.mesh_material = add_material(
      app.scene, {0, 0, 0}, {0.9, 0.9, 0.9}, 0.04, 0, 0.4);
  add_instance(app.scene, identity3x4f, app.mesh_shape, app.mesh_material);

  app.shade_params.hide_environment = true;
  app.shade_params.exposure         = -0.5;
  app.shade_params.background       = {1, 1, 1, 1};
  app.shade_params.lighting         = envlight ? shade_lighting_type::envlight
                                               : shade_lighting_type::eyelight;

  init_scene(app.scene, true);

  // setup IBL
  if (envlight) {
    auto img = image<vec4f>{};
    load_image("data/env.hdr", img, error);
    auto texture = new shade_texture{};
    set_texture(texture, img, true, true, true);
    auto environment = add_environment(app.scene);
    set_emission(environment, {1, 1, 1}, texture);
    init_environments(app.scene);
  }

  // Init shaders.
  auto base = string(SHADERS_PATH);  // defined in parent CMakeLists.txt

  load_program(&app.shaders["mesh"], base + "mesh.vert", base + "mesh.frag");
  load_program(
      &app.shaders["lines"], base + "points.vert", base + "points.frag");

  if (envlight) {
    load_program(&app.shaders["points"], base + "sphere.vert",
        base + "sphere-envlight.frag");
    load_program(&app.shaders["polyline"], base + "polyline.vert",
        base + "polyline-envlight.frag");
  } else {
    load_program(&app.shaders["points"], base + "sphere.vert",
        base + "sphere-eyelight.frag");
    load_program(&app.shaders["polyline"], base + "polyline.vert",
        base + "polyline-eyelight.frag");
  }
  load_program(&app.shaders["flat"], base + "flat.vert", base + "flat.frag");

  app.gpu_shaders["points"] = gpu::make_shader_from_file(
      base + "points.vert", base + "points.frag");
  // init edge shape
  auto surface_offset = 0.001f;
  auto positions      = mesh.positions;
  for (int i = 0; i < mesh.positions.size(); i++) {
    positions[i] += mesh.normals[i] * surface_offset;
  }
  auto edges = vector<vec2i>();
  edges.reserve(mesh.positions.size() * 3);
  for (auto& t : mesh.triangles) {
    if (t.x < t.y) edges.push_back({t.x, t.y});
    if (t.y < t.z) edges.push_back({t.y, t.z});
    if (t.z < t.x) edges.push_back({t.z, t.x});
  }
  set_vertex_buffer(&app.edges_shape, positions, 0);
  set_index_buffer(&app.edges_shape, edges);
}

void delete_app(App& app) {
  for (auto& [name, shape] : app.shapes) {
    clear_shape(&shape);
  }
  for (auto& [name, shader] : app.shaders) {
    clear_program(&shader);
  }
}

#include <yocto/yocto_sceneio.h>

#include "scene_exporter.h"
#ifdef _WIN32
#undef near
#undef far
#endif

void export_yocto_scene(const App& app, const string& name) {
  auto error  = ""s;
  auto output = name;

  // make a directory if needed
  auto dirname       = name;
  auto shapes_path   = path_join(dirname, "shapes");
  auto textures_path = path_join(dirname, "textures");
  if (!make_directory(dirname, error)) print_fatal(error);
  if (!make_directory(shapes_path, error)) print_fatal(error);
  if (!make_directory(textures_path, error)) print_fatal(error);

  auto scene  = new sceneio_scene{};
  scene->name = "name";
  auto camera = add_camera(scene);

  camera->frame    = app.scene->cameras[0]->frame;
  camera->lens     = app.scene->cameras[0]->lens;
  camera->aspect   = app.scene->cameras[0]->aspect;
  camera->film     = app.scene->cameras[0]->film;
  camera->aperture = app.scene->cameras[0]->aperture;
  camera->focus    = app.scene->cameras[0]->focus;

  // add mesh
  {
    auto shape       = add_shape(scene);
    shape->name      = "mesh";
    shape->triangles = app.mesh.triangles;
    shape->positions = app.mesh.positions;
    shape->normals   = app.mesh.normals;

    auto material       = add_material(scene);
    material->color     = {0.5, 0.5, 0.9};
    material->metallic  = 0.04;
    material->roughness = 0.4;

    auto instance      = add_instance(scene);
    instance->shape    = shape;
    instance->material = material;

    auto success = save_shape(path_join(shapes_path, "mesh.obj"), {}, {},
        app.mesh.triangles, {}, {}, {}, {}, app.mesh.positions,
        app.mesh.normals, {}, {}, {}, error);
    if (!success) {
      printf("error saving mesh: %s \n", error.c_str());
    }
  }

  {  // add splines
    auto shape  = add_shape(scene);
    shape->name = "splines";
    for (auto& spline : app.splines()) {
      for (int i = 0; i < spline.curves.size(); i++) {
        auto& curve = spline.curves[i];
        auto  lines = vector<vec2i>(curve.positions.size() - 1);
        for (int k = 0; k < lines.size(); k++) {
          lines[k] = {k, k + 1};
          lines[k] += shape->positions.size();
          append(shape->lines, lines);
        }
        append(shape->positions, curve.positions);
      }
    }
    shape->radius = vector<float>(shape->positions.size(), 0.0015);

    auto material      = add_material(scene);
    material->color    = {0.8, 0.1, 0.1};
    auto instance      = add_instance(scene);
    instance->shape    = shape;
    instance->material = material;
    auto success       = save_shape(path_join(shapes_path, "splines.obj"), {},
              shape->lines, {}, {}, {}, {}, {}, shape->positions, {}, {}, {}, {},
              error);
    if (!success) {
      printf("error saving splines: %s \n", error.c_str());
    }
  }

  {  // add paths
    for (int i = 0; i < app.added_paths.size(); i++) {
      auto& path       = app.added_paths[i]->path;
      auto& color      = app.added_paths[i]->color;
      auto& radius     = app.added_paths[i]->radius;
      auto  shape      = add_shape(scene);
      shape->positions = path_positions(app.mesh, path);
      shape->lines.resize(shape->positions.size() - 1);
      for (int i = 0; i < shape->lines.size(); i++) {
        shape->lines[i] = {i, i + 1};
      }
      shape->radius = vector<float>(shape->positions.size(), radius);

      auto material      = add_material(scene);
      material->color    = color;
      auto instance      = add_instance(scene);
      instance->shape    = shape;
      instance->material = material;
      auto success       = save_shape(
                path_join(shapes_path, std::to_string(i) + "path.obj"), {},
                shape->lines, {}, {}, {}, {}, {}, shape->positions, {}, {}, {}, {},
                error);
      if (!success) {
        printf("error saving path: %s \n", error.c_str());
      }
    }
  }

  {  // add points
    for (int i = 0; i < app.added_points.size(); i++) {
      auto& points = app.added_points[i]->points;
      auto& color  = app.added_points[i]->color;
      auto& radius = app.added_points[i]->radius;
      auto  sphere = make_sphere(32, radius);

      auto material   = add_material(scene);
      material->color = color;

      auto shape       = add_shape(scene);
      shape->positions = sphere.positions;
      shape->triangles = quads_to_triangles(sphere.quads);
      shape->normals   = sphere.normals;

      for (int i = 0; i < points.size(); i++) {
        auto instance      = add_instance(scene);
        instance->shape    = shape;
        instance->material = material;
        instance->frame.o  = eval_position(app.mesh, points[i]);
      }

      auto success = save_shape(
          path_join(shapes_path, std::to_string(i) + "point.obj"),
          shape->points, {}, {}, {}, {}, {}, {}, shape->positions, {}, {}, {},
          {}, error);
      if (!success) {
        printf("error saving point: %s \n", error.c_str());
      }
    }
  }

  auto environment_tex = add_texture(scene);
  load_image("data/env.png", environment_tex->hdr, error);
  environment_tex->name     = "env";
  auto environment          = add_environment(scene);
  environment->emission     = {0, 0, 0};
  environment->emission_tex = environment_tex;

  add_instance(scene, "arealight1",
      lookat_frame({-1, 1, 1}, {0, 0.1, 0}, {0, 1, 0}, true),
      add_shape(scene, "arealight1", make_rect({1, 1}, {0.2, 0.2})),
      add_emission_material(scene, "arealight1", {10, 10, 10}, nullptr));
  add_instance(scene, "arealight2",
      lookat_frame({1, 1, 0.5}, {0, 0.1, 0}, {0, 1, 0}, true),
      add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
      add_emission_material(scene, "arealight2", {10, 10, 10}, nullptr));
  add_instance(scene, "arealight3",
      lookat_frame({0, 1, -1}, {0, 0.1, 0}, {0, 1, 0}, true),
      add_shape(scene, "arealight3", make_rect({1, 1}, {0.2, 0.2})),
      add_emission_material(scene, "arealight3", {10, 10, 10}, nullptr));

  auto floor_material   = add_material(scene, "floor_material");
  floor_material->color = {1, 1, 1};
  auto floor_frame      = identity3x4f;
  std::swap(floor_frame.z, floor_frame.y);
  floor_frame.o.y = -0.5;

  add_instance(scene, "floor", floor_frame,
      add_shape(scene, "floor", make_rect({1, 1}, {1000, 1000})),
      floor_material);

  save_scene(path_join(output, "scene.json"), scene, error);
}
vector<vector<vec3f>> split_broken_curve(
    const vector<vec3f>& positions, const float& threshold) {
  auto curves     = vector<vector<vec3f>>{};
  auto curr_curve = vector<vec3f>{};
  for (auto i = 0; i < positions.size(); ++i) {
    auto from = positions[i];
    auto to   = positions[i + 1];
    if (from == to) continue;
    if (length(to - from) >= threshold) {
      curves.push_back(curr_curve);
      curr_curve.clear();
      continue;
    }
    curr_curve.push_back(from);
  }
  return curves;
}
void update_path_shape(shade_shape* shape, const bezier_mesh& mesh,
    const vector<vec3f>& positions, const float& radius,
    const float& threshold) {
  auto frames = vector<frame3f>();
  frames.reserve(positions.size() - 1);

  for (int i = 0; i < positions.size() - 1; i++) {
    auto from = positions[i];
    auto to   = positions[i + 1];
    if (from == to || length(to - from) >= threshold) continue;

    auto frame = frame_fromz(from, normalize(to - from));
    frame.z *= length(to - from);
    frames.push_back(frame);
  }

  auto cylinder = make_uvcylinder({8, 1, 1}, {radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}

void update_path_shape(Gui_Curve& curve, const bezier_mesh& mesh,
    const float& radius, const float& treshold, const bool cut_after_first) {
  auto frames = vector<frame3f>();
  frames.reserve(curve.positions.size() - 1);
  auto len = path_length(curve.positions);
  for (auto i = 0; i < curve.sampling_rate.size(); ++i) {
    if (curve.sampling_rate[i] > treshold * len && !cut_after_first) continue;
    if (curve.sampling_rate[i] > treshold * len && cut_after_first) break;
    auto j0 = curve.offsets[i].x;
    auto j1 = curve.offsets[i].y;
    for (auto j = j0; j < j1 - 1; ++j) {
      auto from = curve.positions[j];
      auto to   = curve.positions[j + 1];

      auto frame = frame_fromz(from, normalize(to - from));
      frame.z *= length(to - from);
      frames.push_back(frame);
    }
  }

  auto cylinder = make_uvcylinder({8, 1, 1}, {radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }

  auto& shape = curve.shape;

  set_quads(shape, cylinder.quads);
  set_positions(shape, cylinder.positions);
  set_normals(shape, cylinder.normals);
  set_texcoords(shape, cylinder.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}
void update_path_shape(shade_shape* shape, const bezier_mesh& mesh,
    const geodesic_path& path, const float& radius, const float& threshold) {
  auto positions = path_positions(
      path, mesh.triangles, mesh.positions, mesh.adjacencies);

  update_path_shape(shape, mesh, positions, radius);
}

Added_Path* add_path_shape(
    App& app, const geodesic_path& path, float radius, const vec3f& color) {
  auto added_path    = app.added_paths.emplace_back(new Added_Path{});
  added_path->path   = path;
  added_path->color  = color;
  added_path->radius = radius;

  auto shape     = add_shape(app.scene);
  auto material  = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  auto positions = path_positions(
      path, app.mesh.triangles, app.mesh.positions, app.mesh.adjacencies);
  added_path->positions = positions;
  update_path_shape(shape, app.mesh, positions, radius);
  added_path->instance = add_instance(
      app.scene, identity3x4f, shape, material, false);

  return added_path;
}
// change this function in order to create more than one path accordingly to the
// jumps
Added_Path* add_path_shape(App& app, const vector<vec3f>& positions,
    float radius, const vec3f& color, const bool jumps) {
  auto added_path       = app.added_paths.emplace_back(new Added_Path{});
  added_path->color     = color;
  added_path->radius    = radius;
  added_path->positions = positions;
  auto shape            = add_shape(app.scene);
  auto material         = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  if (!jumps)
    update_path_shape(shape, app.mesh, positions, radius);
  else
    update_path_shape(
        shape, app.mesh, positions, radius, app.mesh.avg_edge_length * 10);

  added_path->instance = add_instance(
      app.scene, identity3x4f, shape, material, false);

  return added_path;
}
Added_Points* add_points_shape(App& app, const vector<mesh_point>& points,
    float radius, const vec3f& color) {
  auto added_points    = app.added_points.emplace_back(new Added_Points{});
  added_points->points = points;
  added_points->color  = color;
  added_points->radius = radius;

  auto shape    = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_points_shape(shape, app.mesh, points, radius);
  added_points->instance = add_instance(
      app.scene, identity3x4f, shape, material, false);
  return added_points;
}
Added_Points* add_points_shape(
    App& app, const vector<vec3f>& points, float radius, const vec3f& color) {
  auto added_points = app.added_points.emplace_back(new Added_Points{});
  // added_points->points = points; TODO(giacomo): fix this
  added_points->color  = color;
  added_points->radius = radius;

  auto shape    = add_shape(app.scene);
  auto material = add_material(app.scene, {0, 0, 0}, color, 1, 0, 0.4);
  update_points_shape(shape, points, radius);
  added_points->instance = add_instance(
      app.scene, identity3x4f, shape, material, false);
  return added_points;
}

void update_points_shape(
    shade_shape* shape, const vector<vec3f>& positions, float radius) {
  auto frames = vector<frame3f>(positions.size(), identity3x4f);
  frames.reserve(positions.size() - 1);
  for (int i = 0; i < positions.size(); i++) {
    frames[i].o = positions[i];
  }

  auto sphere = make_sphere(32, radius);

  set_quads(shape, sphere.quads);
  set_positions(shape, sphere.positions);
  set_normals(shape, sphere.normals);
  set_texcoords(shape, sphere.texcoords);
  set_colors(shape, {});
  set_instances(shape, frames);
}

void update_points_shape(shade_shape* shape, const bezier_mesh& mesh,
    const vector<mesh_point>& points, float radius) {
  auto positions = vector<vec3f>(points.size());
  for (int i = 0; i < positions.size(); i++) {
    positions[i] = eval_position(mesh, points[i]);
  }
  update_points_shape(shape, positions, radius);
}

using json = nlohmann::json;
inline void to_json(json& js, const Gui_Spline& spline) {
  js["control_points"] = spline.control_points;
  js["is_smooth"]      = spline.is_smooth;
}
inline void from_json(const json& j, Gui_Spline& spline) {
  j.at("control_points").get_to(spline.control_points);
  j.at("is_smooth").get_to(spline.is_smooth);
}
inline void to_json(json& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}
inline void from_json(const json& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}
inline void to_json(json& js, const shade_camera& camera) {
  //  js["frame"]    = camera.frame;
  js["lens"]     = camera.lens;
  js["aspect"]   = camera.aspect;
  js["film"]     = camera.film;
  js["near"]     = camera.near;
  js["far"]      = camera.far;
  js["aperture"] = camera.aperture;
  js["focus"]    = camera.focus;
}

inline void from_json(const json& js, shade_camera& camera) {
  //  js.at("frame").get_to(camera.frame);
  js.at("lens").get_to(camera.lens);
  js.at("aspect").get_to(camera.aspect);
  js.at("film").get_to(camera.film);
  js.at("near").get_to(camera.near);
  js.at("far").get_to(camera.far);
  js.at("aperture").get_to(camera.aperture);
  js.at("focus").get_to(camera.focus);
}

void save_editing(const App& app, const string& filename) {
  auto js             = json{};
  js["splines"]       = app.splines();
  shade_camera camera = *app.camera;
  js["camera"]        = json{{"lens", camera.lens}, {"aspect", camera.aspect},
      {"film", camera.film}, {"near", camera.near}, {"far", camera.far},
      {"aperture", camera.aperture}, {"focus", camera.focus}};
  nlohmann::to_json(
      js["camera"]["frame"], (const array<float, 12>&)camera.frame);
  string error;
  save_text(filename, js.dump(2), error);
}

bool load_editing(App& app, const string& filename) {
  auto error = ""s;
  auto js    = json{};

  // error helpers
  auto parse_error = [filename, &error]() -> json {
    error = filename + ": parse error in json";
    printf("error loading json %s\n:  %s\n", filename.c_str(), error.c_str());
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    js = json::parse(text);
  } catch (std::exception& e) {
    return parse_error();
  }

  try {
    app.commit_state();
    app.splines()          = js["splines"].get<vector<Gui_Spline>>();
    app.camera->lens       = js["camera"]["lens"];
    app.camera->aspect     = js["camera"]["aspect"];
    app.camera->film       = js["camera"]["film"];
    app.camera->near       = js["camera"]["near"];
    app.camera->far        = js["camera"]["far"];
    app.camera->aperture   = js["camera"]["aperture"];
    app.camera->focus      = js["camera"]["focus"];
    array<float, 12> frame = js["camera"]["frame"];
    app.camera->frame      = *(frame3f*)frame.data();

    //    nlohmann::from_json(
    //        js["camera"]["frame"], (const array<float,
    //        12>&)app.camera->frame);
    app.input() = {};
    update_all_splines(app);
  } catch (std::exception& e) {
    error = "error parsing json: ("s + e.what() + ")"s;
    return false;
  }

  return true;
}
pair<vector<vec3f>, vector<vec3f>> points_pos_and_colors(
    const bezier_mesh& mesh, const vector<Added_Points*>& points) {
  auto pos    = vector<vec3f>{};
  auto colors = vector<vec3f>{};
  for (auto point : points) {
    auto curr_pos   = vector<vec3f>(point->points.size());
    auto curr_color = vector<vec3f>(
        point->points.size(), point->instance->material->color);
    for (auto i = 0; i < point->points.size(); ++i) {
      curr_pos[i] = eval_position(
          mesh.triangles, mesh.positions, point->points[i]);
    }
    pos.insert(pos.end(), curr_pos.begin(), curr_pos.end());
    colors.insert(colors.end(), curr_color.begin(), curr_color.end());
  }
  return {pos, colors};
}
void export_scene(App& app, const string& filename, const bool export_edges,
    const vector<vec3f>& spline_colors) {
  auto final_pos     = app.mesh.positions;
  auto final_normals = app.mesh.normals;
  auto final_tri     = app.mesh.triangles;
  auto final_colors = vector<vec3f>(final_pos.size(), app.mesh_material->color);
  auto size         = (int)app.mesh.positions.size();
  auto point_pos    = vector<vec3f>{};
  if (export_edges) {
    for (auto tr : app.mesh.triangles) {
      for (int k = 0; k < 3; ++k) {
        auto a = tr[k];
        auto b = tr[(k + 1) % 3];
        if (a > b) continue;
        auto quads   = vector<vec4i>{};
        auto pos     = vector<vec3f>{};
        auto normals = vector<vec3f>{};
        auto tex     = vector<vec2f>{};
        polyline_to_cylinders(quads, pos, normals, tex,
            {app.mesh.positions[a], app.mesh.positions[b]}, 8, 0.0004);
        auto tri = quads_to_triangles(quads);
        merge_triangles(final_tri, tri, size);
        final_pos.insert(final_pos.end(), pos.begin(), pos.end());
        final_normals.insert(
            final_normals.end(), normals.begin(), normals.end());
        auto curr_colors = vector<vec3f>(pos.size(), {0, 0, 0});
        final_colors.insert(
            final_colors.end(), curr_colors.begin(), curr_colors.end());
        size += (int)pos.size();
      }
    }
  }
  auto spline_c = spline_colors;
  if (spline_c.size() != app.splines().size()) {
    spline_c = vector<vec3f>(app.splines().size(), vec3f{1, 0, 0});
    // std::cout << "wrong number of splines" << std::endl;
  }

  for (auto i = 0; i < app.splines().size(); ++i) {
    auto& spline = app.splines()[i];
    for (auto curve : spline.curves) {
      auto quads   = vector<vec4i>{};
      auto pos     = vector<vec3f>{};
      auto normals = vector<vec3f>{};
      auto tex     = vector<vec2f>{};
      if (app.enable_jumps) {
        auto len = path_length(curve.positions);
        for (auto j = 0; j < curve.sampling_rate.size(); ++j) {
          if (curve.sampling_rate[j] > app.scale_factor * len) continue;
          auto j0       = curve.offsets[j].x;
          auto j1       = curve.offsets[j].y;
          auto curr_pos = vector<vec3f>{};
          for (auto h = j0; h < j1; ++h) {
            curr_pos.push_back(curve.positions[h]);
          }
          if (curr_pos.size() <= 1) continue;
          polyline_to_cylinders(
              quads, pos, normals, tex, curr_pos, 4, 0.003f * app.curve_size);

          auto tri = quads_to_triangles(quads);
          merge_triangles(final_tri, tri, size);
          final_pos.insert(final_pos.end(), pos.begin(), pos.end());
          final_normals.insert(
              final_normals.end(), normals.begin(), normals.end());
          auto curr_colors = vector<vec3f>(pos.size(), spline_c[i]);
          final_colors.insert(
              final_colors.end(), curr_colors.begin(), curr_colors.end());
          size += (int)pos.size();
        }
      } else {
        polyline_to_cylinders(quads, pos, normals, tex, curve.positions, 12,
            0.003f * app.curve_size);

        auto tri = quads_to_triangles(quads);
        merge_triangles(final_tri, tri, size);
        final_pos.insert(final_pos.end(), pos.begin(), pos.end());
        final_normals.insert(
            final_normals.end(), normals.begin(), normals.end());
        auto curr_colors = vector<vec3f>(pos.size(), spline_c[i]);
        final_colors.insert(
            final_colors.end(), curr_colors.begin(), curr_colors.end());
        size += (int)pos.size();
      }
    }
    for (auto j = 0; j < spline.control_points.size() - 1; ++j) {
      // auto path = compute_geodesic_path(
      //     app.mesh, spline.control_points[j], spline.control_points[j +
      //     1]);
      // auto quads   = vector<vec4i>{};
      // auto pos     = vector<vec3f>{};
      // auto normals = vector<vec3f>{};
      // auto tex     = vector<vec2f>{};

      // polyline_to_cylinders(
      //     quads, pos, normals, tex, path_positions(app.mesh, path), 8,
      //     0.0020);

      // auto tri = quads_to_triangles(quads);
      // merge_triangles(final_tri, tri, size);
      // final_pos.insert(final_pos.end(), pos.begin(), pos.end());
      // final_normals.insert(final_normals.end(), normals.begin(),
      // normals.end()); auto curr_colors = vector<vec3f>(pos.size(), {0, 0,
      // 1}); final_colors.insert(
      //     final_colors.end(), curr_colors.begin(), curr_colors.end());
      // size += (int)pos.size();
      if (i > 0 && j > 0) continue;
      point_pos.push_back(eval_position(app.mesh, spline.control_points[j]));
    }
    point_pos.push_back(eval_position(app.mesh, spline.control_points.back()));
  }

  for (auto& path : app.added_paths) {
    auto quads   = vector<vec4i>{};
    auto pos     = vector<vec3f>{};
    auto normals = vector<vec3f>{};
    auto tex     = vector<vec2f>{};

    polyline_to_cylinders(
        quads, pos, normals, tex, path->positions, 12, path->radius);
    auto tri = quads_to_triangles(quads);
    merge_triangles(final_tri, tri, size);
    final_pos.insert(final_pos.end(), pos.begin(), pos.end());
    final_normals.insert(final_normals.end(), normals.begin(), normals.end());
    auto curr_colors = vector<vec3f>(
        pos.size(), path->instance->material->color);
    final_colors.insert(
        final_colors.end(), curr_colors.begin(), curr_colors.end());
    size += (int)pos.size();
  }
  // for (auto points : app.added_points) {
  //   auto quads   = vector<vec4i>{};
  //   auto pos     = vector<vec3f>{};
  //   auto normals = vector<vec3f>{};
  //   auto tex     = vector<vec2f>{};
  //   point_pos.resize(points->points.size());
  //   for (auto j = 0; j < points->points.size(); ++j) {
  //     point_pos[j] = eval_position(app.mesh, points->points[j]);
  //   }

  //   points_to_spheres(quads, pos, normals, tex, point_pos, 4,
  //   points->radius); for (auto& n : normals) n *= -1; auto tri =
  //   quads_to_triangles(quads); merge_triangles(final_tri, tri, size);
  //   final_pos.insert(final_pos.end(), pos.begin(), pos.end());
  //   final_normals.insert(final_normals.end(), normals.begin(),
  //   normals.end()); auto curr_colors = vector<vec3f>(
  //       pos.size(), points->instance->material->color);
  //   final_colors.insert(
  //       final_colors.end(), curr_colors.begin(), curr_colors.end());
  //   size += (int)pos.size();
  // }

  auto quads   = vector<vec4i>{};
  auto pos     = vector<vec3f>{};
  auto normals = vector<vec3f>{};
  auto tex     = vector<vec2f>{};

  points_to_spheres(quads, pos, normals, tex, point_pos, 8, 0.0075);
  auto point_color = vector<vec3f>(pos.size(), vec3f{0, 0, 1});
  for (auto& n : normals) n *= -1;
  auto tri = quads_to_triangles(quads);
  merge_triangles(final_tri, tri, size);
  final_pos.insert(final_pos.end(), pos.begin(), pos.end());
  final_normals.insert(final_normals.end(), normals.begin(), normals.end());
  final_colors.insert(
      final_colors.end(), point_color.begin(), point_color.end());

  // auto colors8 = vector<vec3b>(final_colors.size());
  // for (auto idx = 0; idx < final_colors.size(); idx++)
  //   colors8[idx] = float_to_byte(final_colors[idx]);

  string err = "";
  // save_mesh(filename, final_tri, final_pos, final_normals, tex, colors8,
  // err,
  //           true);
  save_mesh(filename, final_tri, final_pos, final_normals, tex, final_colors,
      err, true);
}
