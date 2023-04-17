#pragma once

using namespace yocto;

inline bool is_anchor_control_point(int i) {
  if (i == -1) return false;
  return i % 3 == 0;
}
inline bool is_handle_control_point(int i) {
  if (i == -1) return false;
  return !is_anchor_control_point(i);
}
inline void set_selected_spline(App& app, int selected) {
  app.input().selected_spline = selected;
}
inline int add_control_point(
    App& app, Gui_Spline& spline, const mesh_point& point) {
  auto id = (int)spline.control_points.size();
  spline.control_points.push_back(point);
  // Same type as last point.
  // auto is_smooth = app.newspline.is_smooth.back();
  spline.is_smooth.push_back(true);

  return id;
}

inline void update_local_coordinates(
    Gui_Spline& spline, const mesh_point& center, const bezier_mesh& mesh);

inline void update_local_coordinates(
    Gui_Spline& spline, int selected, const bezier_mesh& mesh) {
  update_local_coordinates(spline, spline.control_points[selected], mesh);
}

inline void set_selected_control_point(
    App& app, int selected, const mesh_point& pos, const vec2f& mouse) {
  auto& input = app.input();

  if (selected == -1) {
    input.selected_control_point = -1;
    input.active_control_point   = -1;
    return;
  }

  update_local_coordinates(
      app.spline(), app.spline().control_points[selected], app.mesh);

  auto center = screenspace_from_worldspace(app, eval_position(app.mesh, pos));

  input.active_control_point_offset = mouse - center;
  input.selected_control_point      = selected;
  input.active_control_point        = selected;
}

inline vec2i intersect_control_points(App& app, const vec2f& mouse) {
  // Find index of clicked control point.
  float min_dist        = flt_max;
  int   selected_spline = -1;
  int   selected_point  = -1;
  float threshold       = 0.1;
  for (int spline_id = 0; spline_id < app.splines().size(); spline_id++) {
    auto& spline = app.splines()[spline_id];
    for (int i = 0; i < spline.control_points.size(); ++i) {
      // Skip handle points of non-selected anchors.
      if (is_handle_control_point(i)) {
        if (app.input().selected_spline != spline_id) {
          continue;
        }
        if (yocto::abs(i - app.input().selected_control_point) > 1) {
          continue;
        }
      }
      auto& point  = spline.control_points[i];
      auto  pos    = eval_position(app.mesh, point);
      auto  pos_ss = screenspace_from_worldspace(app, pos);
      float dist   = length(pos_ss - mouse);
      if (dist < threshold && dist < min_dist) {
        selected_point  = i;
        selected_spline = spline_id;
        min_dist        = dist;
      }
    }
  }
  return {selected_spline, selected_point};
}

bool move_handle_point(
    App& app, Gui_Spline& spline, const mesh_point& point, int selected) {
  // Update curve with edited tangent.
  auto& curve = spline.curves[selected / 3];
  spline.curves_to_update.insert(selected / 3);

  // Adjacent curve. Its the the previous or the next.

  int curve_id;
  int offset;
  int t;       // index of edited handle
  int anchor;  // index of edited handle
  if (selected % 3 == 1) {
    curve_id = (selected / 3) - 1;
    offset   = -2;
    t        = 0;
    anchor   = selected - 1;
  } else {
    curve_id = (selected / 3) + 1;
    offset   = 2;
    t        = 1;
    anchor   = selected + 1;
  }
  if (!spline.is_smooth[anchor]) return true;

  if (curve_id >= 0 && curve_id < spline.curves.size()) {
    auto& other_curve = spline.curves[curve_id];

    auto& tangent       = curve.tangents[t];
    auto& other_tangent = other_curve.tangents[1 - t];
    auto  len           = other_tangent.length;
    // path_length(app.mesh, other_tangent.path);
    other_tangent.path = continue_path(app.mesh, tangent.path, -len);
    spline.control_points[selected + offset] = other_tangent.path.end;
    spline.curves_to_update.insert(curve_id);
  }

  return true;
}

inline mat2f parallel_transport_rotation(
    const bezier_mesh& mesh, const mesh_point& start, const mesh_point& end) {
  auto path = compute_geodesic_path(mesh, start, end);
  return parallel_transport_rotation(
      mesh.triangles, mesh.positions, mesh.adjacencies, path);
}

bool move_anchor_point(
    App& app, Gui_Spline& spline, const mesh_point& point, int selected) {
  auto start = spline.control_points[selected];
  if (start == point) {
    return false;
  }

  auto has_point_double_tangent = true;

  auto transport_rotation = parallel_transport_rotation(app.mesh, start, point);
  auto direction          = zero2f;

  for (int i = 0; i < 2; i++) {
    auto offset = 2 * i - 1;
    auto curve  = (selected / 3) + i - 1;
    if (selected + offset < 0 ||
        selected + offset >= spline.control_points.size()) {
      has_point_double_tangent = false;
      continue;
    }
    if (curve < 0 || curve >= spline.curves.size()) {
      has_point_double_tangent = false;
      continue;
    }

    auto& handle  = spline.control_points[selected + offset];
    auto& tangent = spline.curves[curve].tangents[1 - i].path;
    auto  len     = path_length(app.mesh, tangent);
    if (direction == zero2f) {
      direction = tangent_path_direction(app.mesh, tangent);
      direction = transport_rotation * direction;
    } else {
      if (spline.is_smooth[selected]) {
        direction = -direction;
      } else {
        // TODO(giacomo): non-smooth tangent are not handeld properly. Hence,
        // the angle between non-smooth tangents could slightly change.
        direction = tangent_path_direction(app.mesh, tangent);
        direction = transport_rotation * direction;
      }
    }

    auto new_tangent =
        // straightest_geodesic_biermann(
        //     app.mesh, point, direction, len);
        straightest_path(app.mesh.triangles, app.mesh.positions,
            app.mesh.adjacencies, point, direction, len);
    // handle = new_tangent.end;
    spline.curves_to_update.insert(curve);
  }

  spline.control_points[selected] = point;
  return true;
}

void move_active_control_point(
    App& app, Gui_Spline& spline, const vec2f& mouse) {
  // Move selected point.
  auto selected = app.input().active_control_point;
  auto point    = intersect_mesh(
         app, mouse - app.input().active_control_point_offset);

  if (spline.control_points[selected] == point) return;

  if (point.face == -1) return;

  bool success = false;
  if (is_handle_control_point(selected)) {
    spline.control_points[selected] = point;
    success = move_handle_point(app, spline, point, selected);
  } else {
    success = move_anchor_point(app, spline, point, selected);
  }
}

void add_new_curve_to_spline(
    App& app, Gui_Spline& spline, const mesh_point& point, const vec2f& mouse) {
  // Add new curve to spline
  auto curve_id = (int)spline.curves.size();
  spline.curves.push_back({});
  if (spline.curves.size() > 1) {
    if (spline.is_smooth.back()) {
      auto& prev_tangent = spline.curves[spline.curves.size() - 2].tangents[1];
      auto  len          = prev_tangent.length;
      auto  tangent      = continue_path(app.mesh, prev_tangent.path, -len);
      add_control_point(app, spline, tangent.end);
    } else {
      add_control_point(app, spline, spline.control_points.back());
    }
  }

  auto selected = (int)add_control_point(app, spline, point);
  add_control_point(app, spline, point);
  set_selected_control_point(app, selected, point, mouse);
  spline.curves_to_update.insert(curve_id);
}

inline void update_local_coordinates(
    Gui_Spline& spline, const mesh_point& center, const bezier_mesh& mesh) {
  if (spline.local_coords.size() == spline.control_points.size()) return;
  spline.local_coords.resize(spline.control_points.size());
  spline.center = center;
  for (int i = 0; i < spline.control_points.size(); i++) {
    if (spline.control_points[i] == center) {
      spline.local_coords[i] = {0, 0};
      continue;
    }
    auto path = compute_geodesic_path(mesh, center, spline.control_points[i]);
    auto direction         = tangent_path_direction(mesh, path);
    auto length            = path_length(mesh, path);
    spline.local_coords[i] = direction * length;
  }
}

inline void regenerate_from_local_coordinate(
    Gui_Spline& spline, const bezier_mesh& mesh, const mesh_point& center) {
  auto f = [&](int i) {
    if (spline.local_coords[i] == vec2f{0, 0}) {
      spline.control_points[i] = center;
      return;
    }
    auto coord = spline.frame * spline.local_coords[i];
    auto path  = straightest_path(mesh, center, coord);
    // straightest_geodesic_biermann(
    //     mesh, center, coord);
    check_point(path.end);
    spline.control_points[i] = path.end;
  };
  parallel_for((int)spline.control_points.size(), f);
  // serial_for((int)spline.control_points.size(), f);
}

inline bool rotate_spline(
    const bezier_mesh& mesh, Gui_Spline& spline, float angle) {
  if (angle == 0) return false;

  auto v        = vec2f{yocto::cos(angle), yocto::sin(angle)};
  auto rotation = mat2f{{v.x, -v.y}, {v.y, v.x}};
  spline.frame  = rotation * spline.frame;
  regenerate_from_local_coordinate(spline, mesh, spline.center);
  return true;
}

inline bool translate_spline(
    const bezier_mesh& mesh, Gui_Spline& spline, const mesh_point& point) {
  if (point == spline.center) return false;

  auto rotation = parallel_transport_rotation(mesh, spline.center, point);
  spline.frame  = rotation * spline.frame;
  spline.center = point;
  regenerate_from_local_coordinate(spline, mesh, spline.center);
  return true;
}

inline bool scale_spline(
    const bezier_mesh& mesh, Gui_Spline& spline, float scaling) {
  if (scaling == 1) return false;
  spline.frame = spline.frame * scaling;
  regenerate_from_local_coordinate(spline, mesh, spline.center);
  return true;
}
