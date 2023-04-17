#include <stdio.h>
#include <yocto/yocto_common.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>

#include <thread>
#include <vector>
using namespace std;

#include <splinesurf/spline.h>
#include <splinesurf/splineio.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_window.h>

#include "app.h"
using namespace flipout;
using namespace yocto;

//
#include "editing.h"
#include "playback.h"

void set_common_uniforms(const App& app, const ogl_program* program) {
  auto& view       = app.matrices.view;
  auto& projection = app.matrices.projection;
  set_uniform(program, "frame", identity4x4f);
  set_uniform(program, "view", view);
  set_uniform(program, "projection", projection);
  set_uniform(program, "eye", app.camera->frame.o);
  set_uniform(program, "envlight", (int)app.envlight);
  set_uniform(program, "gamma", app.shade_params.gamma);
  set_uniform(program, "exposure", app.shade_params.exposure);
  // set_uniform(program, "size", app.line_size);
  if (app.scene->environments.size()) {
    auto& env = app.scene->environments.front();
    if (env->envlight_diffuse)
      set_uniform(program, "envlight_irradiance", env->envlight_diffuse, 6);
    if (env->envlight_specular)
      set_uniform(program, "envlight_reflection", env->envlight_specular, 7);
    if (env->envlight_brdflut)
      set_uniform(program, "envlight_brdflut", env->envlight_brdflut, 8);
  }
}

void draw_scene(const App& app, const vec4i& viewport) {
  clear_ogl_framebuffer(vec4f{0, 0, 0, 1});

  // Draw mesh and environment.
  draw_scene(app.scene, app.camera, viewport, app.shade_params);

  if (app.show_points) {
    auto program = &app.shaders.at("points");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "size", 3.0f * 0.0015f * app.line_size);

    // TODO(giacomo): delete this
    set_uniform(program, "color", vec3f{0, 1, 0});

    draw_shape(app.input().temp_points);
    // Draw all anchor points of all splines.
    set_uniform(program, "color", vec3f{1, 0, 0});
    if (app.input().selected_spline != -1) {
      draw_shape(app.spline().anchor_points_shape);
    }

    auto draw_mesh_point = [&](const mesh_point& point, const vec3f& color) {
      if (point.face == -1) return;
      auto        p     = eval_position(app.mesh, point);
      static auto shape = ogl_shape{};
      set_points_shape(&shape, {p});
      set_uniform(program, "color", color);
      draw_shape(&shape);
    };

    // Draw new point.
    draw_mesh_point(app.input().new_point, {1, 0.5, 0});
    draw_mesh_point(app.input().eval_point, {0, 1, 1});
    draw_mesh_point(app.input().eval_point_cheap, {1, 1, 0});

    for (int i = 0; i < app.eval_points.size(); i++) {
      draw_mesh_point(app.eval_points[i], {1, 1, 1});
    }

    // Draw fixed handle of first curve of spline.
    if (app.spline().curves.empty() &&
        app.spline().control_points.size() == 2) {
      draw_mesh_point(app.spline().control_points[1], {1, 0, 0});
    }
  }
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(
                     app.camera->film / (camera_aspect * 2 * app.camera->lens)))
          : (2 * yocto::atan(app.camera->film / (2 * app.camera->lens)));
  auto view       = frame_to_mat(inverse(app.camera->frame));
  auto projection = perspective_mat(
      camera_yfov, camera_aspect, app.shade_params.near, app.shade_params.far);
  if (app.temp_levels > 0) draw_shape(app.temp_points[app.temp_levels]);
  if (app.gpu_shapes.find("vector_field") != app.gpu_shapes.end())
    gpu::draw_shape(app.gpu_shapes.at("vector_field"),
        app.gpu_shaders.at("points"), gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));
  // Draw Spline
  {
    auto program = &app.shaders.at("polyline");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "size", 0.0015f * app.line_size);

    // for (auto& spline : app.splines()) {
    //   for (int i = 0; i < spline.curves.size(); i++) {
    //     auto& curve = spline.curves[i];
    //     set_uniform(program, "color", spline.color);
    //     draw_shape(curve.shape->shape);
    //   }
    // }

    // TODO: move somewhere
    auto handle_points = vector<vec3f>();
    if (app.show_points && app.input().selected_control_point != -1) {
      auto selected = app.input().selected_control_point;
      if (is_handle_control_point(selected)) {
        if (selected % 3 == 1) selected -= 1;
        if (selected % 3 == 2) selected += 1;
      }
      set_uniform(program, "color", vec3f{0, 0, 1});
      auto prev_id = (selected / 3) - 1;
      if (prev_id >= 0) {
        auto& prev_curve = app.spline().curves[prev_id];
        draw_shape(prev_curve.tangents[1].shape);
        handle_points.push_back(
            eval_position(app.mesh, app.spline().control_points[selected - 1]));
      }
      auto next_id = (selected / 3);
      if (next_id < app.spline().curves.size()) {
        auto& next_curve = app.spline().curves[next_id];
        draw_shape(next_curve.tangents[0].shape);
        handle_points.push_back(
            eval_position(app.mesh, app.spline().control_points[selected + 1]));
      }

      auto program = &app.shaders.at("points");
      bind_program(program);
      set_common_uniforms(app, program);
      set_uniform(program, "size", 3 * 0.0015f * app.line_size);
      if (app.envlight) {
        auto& env = app.scene->environments.front();
        set_uniform(program, "envlight_irradiance", env->envlight_diffuse, 6);
        set_uniform(program, "envlight_reflection", env->envlight_specular, 7);
        set_uniform(program, "envlight_brdflut", env->envlight_brdflut, 8);
      }
      static auto points_shape = ogl_shape{};
      set_points_shape(&points_shape, handle_points);
      draw_shape(&points_shape);
    }
  }

  if (app.show_edges) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{0, 0, 0});
    draw_shape(&app.edges_shape);
  }

  // for (auto& command : app.render_commands) {
  //   command();
  // }
  // app.render_commands.clear();
}

void clear_input(App& app) {
  //  app.spline().control_points.clear();
  //  clear_shape(&app.shapes["control_points"]);
  //  clear_shape(&app.shapes["tangent0"]);
  //  clear_shape(&app.shapes["tangent1"]);
  //  clear_shape(&app.shapes["extension0"]);
  //  clear_shape(&app.shapes["extension1"]);
  // clear_spline(app.spline());
}

inline void sleep(int ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

void process_gui_input(App& app, gui_window* win) {
  // if (win->gui_buttons["Clear"]) {
  //   clear_input(app);
  // }
  // if (win->gui_buttons["Spline [enter]"]) {
  //   if (app.spline().control_points.size() < 4) {
  //     //      if (!load_bezier_params(app.testname,
  //     app.spline().control_points,
  //     //              app.bezier_params, app.error)) {
  //     //        print_fatal(app.error);
  //     //      }
  //     // TODO(giacomo): restore loading tests
  //   }
  // }
}

inline bool is_pressing(gui_button button) {
  return button.state == gui_button::state::pressing;
}
inline bool is_releasing(gui_button button) {
  return button.state == gui_button::state::releasing;
}
inline bool is_down(gui_button button) {
  return button.state == gui_button::state::down ||
         button.state == gui_button::state::pressing;
}

inline bool is_pressing(const gui_input& input, gui_key key) {
  return is_pressing(input.key_buttons[(int)key]);
}

void process_key_input(App& app, const gui_input& input) {
  for (int i = 0; i < input.key_buttons.size(); i++) {
    auto key = gui_key(i);
    if (!is_pressing(input, key)) continue;

    printf("%c pressed!\n", (char)key);

    // if (key == gui_key::enter) compute_bezier(app, app.spline.back());
    if (key == gui_key::enter) {
      app.commit_state();
      app.input().selected_spline = (int)app.splines().size();
      app.state.splines.push_back({});
      app.input().selected_control_point = -1;
      app.input().active_control_point   = -1;
    }
    if (key == gui_key::escape) {
      app.show_points = !app.show_points;
    }

    // Affine transforms
    if (key == gui_key('T')) {
      app.input().rotating    = false;
      app.input().scaling     = false;
      app.input().translating = true;
      update_local_coordinates(
          app.spline(), app.input().selected_control_point, app.mesh);
    }
    if (key == gui_key('R')) {
      app.input().rotating    = true;
      app.input().scaling     = false;
      app.input().translating = false;
      update_local_coordinates(
          app.spline(), app.input().selected_control_point, app.mesh);
    }
    if (key == gui_key('S')) {
      app.input().rotating    = false;
      app.input().scaling     = true;
      app.input().translating = false;
      update_local_coordinates(
          app.spline(), app.input().selected_control_point, app.mesh);
    }

    if (key == gui_key('C')) {
      app.commit_state();
      app.spline().control_points.back() = app.spline().control_points[0];
      app.spline().curves_to_update.insert(app.spline().curves.size() - 1);
    }

    if (key == gui_key('D')) {
      app.commit_state();
      app.input().selected_spline = (int)app.splines().size();
      app.state.splines.push_back(app.splines().back());
    }
    if (key == gui_key(' ')) {
      auto id = get_selected_anchor_point(app);
      if (id != -1) {
        app.spline().is_smooth[id].flip();
      }
    }
    if (key == gui_key('Z')) {
      if (app.editing_history.size()) {
        clear_all_splines(app);
        app.undo();
        update_all_splines(app);
      }
    }
    if (key == gui_key('Y')) {
      if (app.editing_history.size()) {
        clear_all_splines(app);
        app.redo();
        update_all_splines(app);
      }
    }

    if (key == gui_key('P')) {
      app.commit_state();
      auto&                    spline = app.spline();
      auto&                    curve  = spline.curves[0];
      quadratic_bezier_segment left_q = {}, right_q = {};
      bezier_segment           left = {}, right = {};
      if (app.quadric_curve) {
        auto                     polygon      = get_control_polygon(spline, 0);
        quadratic_bezier_segment quad_polygon = {
            polygon[0], polygon[1], polygon[3]};

        std::tie(left_q, right_q) = insert_point(
            app.mesh, quad_polygon, app.input().parameter_t);

        spline        = {};
        spline.curves = {{}, {}};
        for (int i = 0; i < 4; i++) {
          if (i < 2)
            add_control_point(app, spline, left_q[i]);
          else if (i == 2)
            add_control_point(app, spline, polygon.back());
          else
            add_control_point(app, spline, left_q.back());
        }
        for (int i = 1; i < 4; i++) {
          if (i < 2)
            add_control_point(app, spline, right_q[i]);
          else if (i == 2)
            add_control_point(app, spline, polygon.back());
          else
            add_control_point(app, spline, right_q.back());
        }
        update_control_points_shape(spline, app.mesh);
        spline.curves_to_update = {0, 1};

      } else {
        auto curve_id = (app.input().selected_control_point / 3);
        curve_id      = clamp(curve_id, 0, (int)app.spline().curves.size());
        auto polygon  = get_control_polygon(app.spline(), curve_id);
        if (app.bezier_params.algorithm ==
                spline_algorithm::de_casteljau_uniform ||
            app.bezier_params.algorithm ==
                spline_algorithm::de_casteljau_adaptive)
          std::tie(left, right) = insert_point(app.mesh, polygon, 0.5);
        else
          std::tie(left, right) = insert_point_spline(app.mesh,
              {app.spline().control_points[0], app.spline().control_points[1],
                  app.spline().control_points[2],
                  app.spline().control_points[3]},
              app.input().parameter_t, app.bezier_params);
        for (auto& curve : spline.curves) {
          clear_shape(curve.shape);
          clear_shape(curve.tangents[0].shape);
          clear_shape(curve.tangents[1].shape);
        }

        spline        = {};
        spline.curves = {{}, {}};

        for (int i = 0; i < 4; i++) {
          add_control_point(app, spline, left[i]);
        }
        for (int i = 1; i < 4; i++) {
          add_control_point(app, spline, right[i]);
        }
        update_control_points_shape(spline, app.mesh);
        spline.curves_to_update = {0, 1};
      }

      // auto& new_spline = app.state.splines.emplace_back();
      // for (int i = 0; i < 3; i++) {
      //   new_spline.control_points.push_back(left[i]);
      // }
      // for (int i = 0; i < 4; i++) {
      //   new_spline.control_points.push_back(right[i]);
      // }
      // new_spline.curves.resize(2);
      // new_spline.curves_to_update.insert(0);
      // new_spline.curves_to_update.insert(1);

      // auto temp_result =
    }

    if (key == gui_key('C')) clear_input(app);
    if (key == gui_key('I')) {  // Record input
      if (!app.recording_input) {
        app.recording_input = !app.recording_input;
      } else {
        app.recording_input = !app.recording_input;
        save_input_record(app.input_record, "apps/splinegui/test.recording");
      }
    }
    if (key == gui_key('R') && input.modifier_shift) {
      for (auto& [name, shader] : app.shaders) {
        auto base = string(SHADERS_PATH);  // defined in parent CMakeLists.txt
        printf("compiling %s shader\n", name.c_str());
        auto vert = base + name + ".vert";
        auto frag = base + name + ".frag";
        load_program(&shader, vert, frag);
      }
    }
  }
}

bool process_camera_move(App& app, const gui_input& input) {
  // Rotate or pan camera by dragging left mouse click.
  auto  rotating = input.modifier_shift;
  auto  panning  = input.modifier_alt;
  auto& camera   = *app.camera;

  auto update_camera_frame = [&](frame3f& frame, float& focus, bool rotating,
                                 bool panning, bool zooming) {
    auto last_pos    = input.mouse_last;
    auto mouse_pos   = input.mouse_pos;
    auto mouse_left  = is_down(input.mouse_left);
    auto mouse_right = is_down(input.mouse_right);
    // handle mouse and keyboard for navigation
    if (mouse_left) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left) rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right) dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_left) pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x    = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    }
  };

  if (is_down(input.mouse_left) && (rotating || panning)) {
    update_camera_frame(
        camera.frame, app.camera_focus, rotating, panning, false);
    return true;
  }

  // Zoom-in/out by scrolling;
  float zoom = input.scroll.y * 0.1;
  if (zoom != 0) {
    update_turntable(camera.frame, app.camera_focus, zero2f, zoom, zero2f);
    return true;
  }

  return false;
}

bool process_user_input(App& app, const gui_input& input) {
  //  static bool yyy = false;

  if (process_camera_move(app, input)) {
    update_camera_info(app, input);
    return false;
  }

  auto mouse = input.mouse_pos;
  auto size  = vec2f{(float)input.window_size.x, (float)input.window_size.y};
  mouse      = vec2f{2 * (mouse.x / size.x) - 1, 1 - 2 * (mouse.y / size.y)};

  // Process left click:
  // - Selecting existing control point?
  // - Adding control point?
  //    - First curve? Or just a new curve?
  //
  auto editing = app.input().rotating || app.input().scaling ||
                 app.input().translating;

  if (is_pressing(input.mouse_left)) {
    auto [spline_id, point_id] = intersect_control_points(app, mouse);
    // app.input().drag_start_uv  = mouse;
    // set_selected_spline(app, spline_id);
    if (point_id != -1) {
      set_selected_spline(app, spline_id);
      auto p = app.spline().control_points[point_id];
      set_selected_control_point(app, point_id, p, mouse);
      return true;
    } else {
      set_selected_control_point(app, point_id, {}, mouse);
    }
  }

  if (is_pressing(input.mouse_left) && !editing) {
    // Here if pressing, but not clicked on an existing control point.
    auto point = intersect_mesh(app, mouse);
    if (!app.weighted_average) {
      if (point.face != -1) {
        app.input().new_point = point;
        app.xxx_point         = point;

        // If first two points of first curve, do nothing.
        // When releasing, this point will be added to the control points of the
        // spline.
        if (app.spline().curves.empty() &&
            app.spline().control_points.size() <= 1) {
          return true;
        }

        // Create new curve and enter in editing mode;
        add_new_curve_to_spline(app, app.spline(), point, mouse);
        return true;
      }
    } else {
      app.control_points.push_back(point);
      if (app.added_points.size() != 0) {
        for (auto& points : app.added_points) {
          clear_shape(points->instance->shape);
        }
        app.added_points.clear();
      }
      auto new_control_points = add_points_shape(
          app, app.control_points, 0.003, {0, 0, 1});
      return true;
    }
  }
  auto drag = input.mouse_pos - input.mouse_last;
  if (is_down(input.mouse_left) && length(drag) > 0.05) {
    // Rotate
    if (app.input().rotating) {
      for (auto& spline : app.splines()) {
        auto angle = drag.x * 0.1;
        rotate_spline(app.mesh, spline, angle);
        auto params = app.bezier_params;
        params.subdivisions += 1;
        params.precision += 1;
      }
      update_all_splines(app);
      return true;
    }

    // Scale
    if (app.input().scaling) {
      for (auto& spline : app.splines()) {
        auto scaling = clamp(1 + drag.y * 0.1f, 0.9f, 1.1f);
        scale_spline(app.mesh, spline, scaling);
        auto params = app.bezier_params;
        params.subdivisions += 1;
        params.precision += 1;
      }
      update_all_splines(app);
      return true;
    }

    auto point = intersect_mesh(app, mouse);

    // Translate or move control point
    if (point.face != -1) {
      if (app.input().translating) {
        for (auto& spline : app.splines()) {
          translate_spline(app.mesh, spline, point);

          auto params = app.bezier_params;
          params.subdivisions += 1;
          params.precision += 1;
        }
        update_all_splines(app);
        return true;
      }

      if (app.input().active_control_point == -1) {
        app.input().new_point = point;
      } else {
        move_active_control_point(app, app.spline(), mouse);
      }
      // return?
    }
  }

  if (is_releasing(input.mouse_left)) {
    app.commit_state();  // TODO(giacomo): Bug?
    if (app.input().new_point.face != -1) {
      if (app.spline().curves.empty()) {
        add_control_point(app, app.spline(), app.input().new_point);
      }
    }

    if (app.input().active_control_point != -1) {
      auto selected = app.input().active_control_point;
      auto params   = app.bezier_params;
      params.subdivisions += 1;
      params.precision += 1;
      auto& spline = app.spline();
      if (is_handle_control_point(selected)) {
        spline.curves_to_update.insert(selected / 3);
        if (app.input().enforce_tangents) {
          auto id = selected / 3;
          if (selected % 3 == 1) id -= 1;
          if (selected % 3 == 2) id += 1;
          if (id >= 0 && id < spline.curves.size()) {
            spline.curves_to_update.insert(id);
          }
        }
      }
      if (is_anchor_control_point(selected)) {
        auto prev = (selected / 3) - 1;
        if (prev >= 0 && prev < spline.curves.size()) {
          spline.curves_to_update.insert(prev);
        }
        auto next = (selected / 3);
        if (next >= 0 && next < spline.curves.size()) {
          spline.curves_to_update.insert(next);
        }
      }
    }
    auto& input                = app.input();
    input.active_control_point = -1;
    input.new_point            = {};
    app.context                = editing_context::is_doing_nothing;
  }

  return false;
}

void update_app(App& app, const gui_input& input) {
  // process_gui_input(app, input); TODO(giacomo)

  if (is_active(app.widget)) return;

  app.window_size = input.window_size;
  // TODO(giacomo): make this optional, curve evaluation viz
  /*if (app.input().selected_spline != -1 && app.state.splines.size()) {
    double i = 0;
    auto   t = modf(input.time_now * 0.1, &i);
    if (app.spline().curves.size()) {
      auto segment = bezier_segment{
          app.spline().control_points[0],
          app.spline().control_points[1],
          app.spline().control_points[2],
          app.spline().control_points[3],
      };
      app.input().eval_point = eval_bezier_point(app.mesh, segment, t);
    }
  }*/

  process_key_input(app, input);

  process_user_input(app, input);

  auto tasks = vector<vec2i>{};

  for (int i = 0; i < app.splines().size(); i++) {
    if (app.splines()[i].curves_to_update.size()) {
      update_control_points_shape(app.splines()[i], app.mesh);
    }
    for (auto& k : app.splines()[i].curves_to_update) {
      tasks.push_back({i, k});
    }
  }
  auto f = [&](int i) {
    auto& spline   = app.splines()[tasks[i].x];
    auto  curve_id = tasks[i].y;
    // update_curve_shape(spline, curve_id, app.mesh, app.bezier_params);
    update_curve_shape(spline, curve_id, app.mesh, app.bezier_params,
        app.quadric_curve, app.use_vector_heat);
  };

  if (tasks.size()) {
#if PROFILE
    auto timer = print_timed("compute_bezier");
#endif
    // serial_for(tasks.size(), f);
    parallel_for(tasks.size(), f);
  }

  if (tasks.size()) {
#if PROFILE
    auto timer = print_timed("update_gpu");
#endif
    for (auto& spline : app.splines()) {
      for (auto& k : spline.curves_to_update) {
        auto& curve = spline.curves[k];
        // set_polyline_shape(curve.shape->shape, curve.positions);
        if (!curve.shape) {
          curve.shape = add_shape(app.scene);
          add_instance(
              app.scene, identity3x4f, curve.shape, app.spline_material, false);
        }
        if (app.enable_jumps)
          update_path_shape(
              curve, app.mesh, 0.003f * app.curve_size, app.scale_factor);
        else
          update_path_shape(
              curve.shape, app.mesh, curve.positions, 0.003f * app.curve_size);

        set_polyline_shape(
            curve.tangents[0].shape, curve.tangents[0].positions);
        set_polyline_shape(
            curve.tangents[1].shape, curve.tangents[1].positions);
      }

      spline.curves_to_update.clear();
    }
  }
}

void init_from_svg(App& app, const Svg& svg) {
  auto timer  = print_timed("load_svg");
  auto center = app.xxx_point;

  auto p0    = eval_position(app.mesh, {center.face, {0, 0}});
  auto p1    = eval_position(app.mesh, {center.face, {0, 1}});
  auto v     = normalize(p1 - p0);
  auto frame = basis_fromz(eval_normal(app.mesh, {center.face, {0, 0}}));
  auto rot   = vec2f{dot(v, frame.x), dot(v, frame.y)};
  //
  app.commit_state();
  app.splines() = {};

  for (auto& shape : svg) {
    for (auto& path : shape.paths) {
      auto& spline  = app.splines().emplace_back();
      spline.center = center;
      spline.frame  = mat2f{rot, vec2f{-rot.y, rot.x}};
      spline.color  = shape.color;
      for (auto& segment : path) {
        spline.curves.push_back({});
        for (int i = 0; i < 3; i++) {
          vec2f uv = clamp(segment[i], 0.0f, 1.0f);
          uv -= vec2f{0.5, 0.5};
          uv *= app.svg_size;
          // auto line = straightest_path(app.mesh, center, uv);
          add_control_point(app, spline, {});
          spline.local_coords.push_back(uv);
        }
      }
      auto& segment = path.back();
      vec2f uv      = clamp(segment[3], 0.0f, 1.0f);
      uv -= vec2f{0.5, 0.5};
      uv *= app.svg_size;
      // auto line = straightest_path(app.mesh, center, uv);
      add_control_point(app, spline, {});
      spline.local_coords.push_back(uv);
    }
  }
  update_all_splines(app);

  {
    auto timer = print_timed("compute_svg");
    for (auto& spline : app.splines()) {
      regenerate_from_local_coordinate(spline, app.mesh, spline.center);
    }
  }
}

void init_from_svg_weights(App& app, const Svg& svg) {
  auto start = app.xxx_point;

  app.commit_state();
  app.splines() = {};
  vector<vector<float>> weights;

  for (auto& shape : svg) {
    for (auto& path : shape.paths) {
      auto& spline = app.splines().emplace_back();
      for (auto& segment : path) {
        spline.curves.push_back({});
        for (int i = 0; i < 3; i++) {
          // vec2f uv = clamp(segment[i], 0.0f, 1.0f);
          vec2f uv = segment[i];
          if (uv.x < 0) {
            printf("skip!\n");
            continue;
          }
          if (uv.y < 0) {
            printf("skip!\n");
            continue;
          }
          if (uv.x > 1) {
            printf("skip!\n");
            continue;
          }
          if (uv.y > 1) {
            printf("skip!\n");
            continue;
          }
          auto w0 = (1 - uv.x) * (1 - uv.y);
          auto w1 = uv.y * (1 - uv.x);
          auto w2 = uv.x * (1 - uv.y);
          auto w3 = uv.x * uv.y;
          weights.push_back({w0, w1, w2, w3});
        }
      }
      auto& segment = path.back();
      vec2f uv      = clamp(segment[3], 0.0f, 1.0f);
      auto  w0      = (1 - uv.x) * (1 - uv.y);
      auto  w1      = uv.y * (1 - uv.x);
      auto  w2      = uv.x * (1 - uv.y);
      auto  w3      = uv.x * uv.y;
      weights.push_back({w0, w1, w2, w3});
    }
  }

  auto rectangle = vector<mesh_point>{};
  for (vec2f dir : vector<vec2f>{{1, 1}, {-1, 1}, {1, -1}, {-1, -1}}) {
    auto line = straightest_path(app.mesh, app.xxx_point, dir * 0.23f);
    rectangle.push_back(line.end);
  }
  // set_points_shape(app.input().temp_points, app.mesh, rectangle);
  auto result = weighted_average(app.mesh, rectangle, weights);
  set_points_shape(app.input().temp_points, app.mesh, result);
}

void draw(const gui_input& input, void* data) {
  auto& app = *(App*)data;
  // TODO(giacom): Make the cpu rest when window's app is not focused.
  //    if (started && !win->input.is_window_focused) {
  //      sleep(30);
  //      return;
  //    }
  app.started = true;

  update_camera_info(app, input);

  if (app.recording_input) {
    app.input_record.push_back(input);
  }

  // Do everything
  auto& t = app.playback_tick;
  if (app.playback && t < app.input_record.size()) {
    update_app(app, app.input_record[t]);
    t += 1;
  } else {
    update_app(app, input);
  }

  // Draw mesh, lines and stuff.
  draw_scene(app, input.framebuffer_viewport);

  if (is_pressing(input, gui_key('M'))) {
    auto img = capture_screenshot();
    for (auto& c : img) {
      c.x = yocto::pow(c.x, 2.2);
      c.y = yocto::pow(c.y, 2.2);
      c.z = yocto::pow(c.z, 2.2);
    }

    auto error   = string{};
    auto success = save_image(
        path_join("figures/", path_basename(app.testname)) + ".png", img,
        error);
    if (!success) {
      printf("[%s]: %s\n", __FUNCTION__, error.c_str());
    }
  }

  // Draw gui.
  // gui_begin(win, "Splinegui");
  auto widget = app.widget;
  begin_widget(widget, "splinegui");

  draw_textinput(widget, "mesh name", app.filename);
  draw_textinput(widget, "test name", app.testname);
  draw_textinput(widget, "scene_name", app.scene_name);
  draw_textinput(widget, "Import Curve", app.comparison);
  if (draw_button(widget, "Save test")) {
    save_editing(app, app.testname);
    // if (!save_bezier_params(app.testname, app.spline().control_points,
    //         app.bezier_params, app.error))
    //   print_fatal(app.error);
  }
  continue_line(widget);
  if (draw_button(widget, "Load test")) {
    load_editing(app, app.testname);
    //      if (!load_bezier_params(app.testname,
    //      app.spline().control_points,
    //              app.bezier_params, app.error))
    //        print_fatal(app.error);
    // TODO(giacomo): restore loading tests
  }
  if (draw_checkbox(widget, "Show Gradient", app.show_gradient)) {
    app.vector_field.resize(app.mesh.positions.size(), zero3f);
    if (app.show_gradient) {
      auto t       = app.input().parameter_t;
      auto polygon = get_control_polygon(app.spline(), 0);
      auto w       = vector<float>(polygon.size());
      bernstein_polynomials((int)polygon.size(), t, w);
      vector<int> vertices(polygon.size());
      for (int i = 0; i < 4; i++) {
        vertices[i] = vert_from_point(app.mesh.triangles, polygon[i]);
      }
      if (app.karcher_fields.size() == 0)
        std::tie(app.karcher_grads, app.karcher_fields) = compute_karcher_grad(
            *app.mesh.flipout, vertices, app.mesh.solver, app.mesh.triangles,
            app.mesh.positions, app.mesh.normals, app.mesh.Grad);

      app.vector_field = karcher_grad(
          app.karcher_grads, app.karcher_fields, w, app.mesh.normals, true);
      update_glvector_field(app, app.vector_field,
          app.vector_size_factor * app.vector_size, "vector_field");
    } else {
      gpu::delete_shape(app.gpu_shapes.at("vector_field"));
      app.vector_field.clear();
    }
  }
  draw_separator(widget);
  draw_textinput(widget, "Scene Name", app.exported_scene_name);
  if (draw_button(widget, "Export Scene")) {
    auto spline_colors = vector<vec3f>{};
    export_scene(app, app.exported_scene_name, app.export_edges, spline_colors);
  }
  continue_line(widget);
  draw_checkbox(widget, "Export Edges", app.export_edges);
  if (draw_slider(widget, "Vector Size", app.vector_size_factor, 0, 10)) {
    if (app.vector_field.size() > 0)
      update_glvector_field(app, app.vector_field,
          app.vector_size_factor * app.vector_size, "vector_field");
  }
  if (draw_button(widget, "Import Anchors")) {
    auto control_polygon = import_control_points(
        app.mesh, "WA_anchors/anchors_" + app.comparison + ".txt");
    if (control_polygon.size() == 0) return;
    if (app.show_olr) {
      app.input().selected_spline = (int)app.splines().size();
      app.state.splines.push_back({});
    }

    vector<vector<int>>        vertices(control_polygon.size(), vector<int>(4));
    vector<vector<vec3f>>      flipout_pos(control_polygon.size());
    vector<vector<mesh_point>> points(
        control_polygon.size(), vector<mesh_point>(4));

    // WA
    if (app.show_wa_curve) {
      std::tie(app.WA_curve, app.WA_samples) = import_curve(
          app.mesh, "WA_curves/result_" + app.comparison + ".txt");
      if (app.WA_curve.size() == 0) return;
      auto curves = split_broken_curve(
          app.WA_curve, app.mesh.avg_edge_length * 10);
      for (auto& curve : curves) {
        add_path_shape(
            app, curve, 0.003f * app.curve_size, vec3f{0, 0.470, 0}, true);
      }

      // if (app.show_WA_samples)
      //   app.WA_samples_shape = add_points_shape(
      //       app, app.WA_samples, 0.003, vec3f{0, 0, 0});
    }
    // OLR
    if (app.show_olr) {
      auto& spline = app.splines().back();
      spline.curves.resize(control_polygon.size());
      app.show_points = false;
      for (auto j = 0; j < control_polygon.size(); ++j) {
        for (int i = 0; i < 4; i++) {
          vertices[j][i] = vert_from_point(
              app.mesh.triangles, control_polygon[j][i]);
          points[j][i] = control_polygon[j][i];
          // make_point_from_vert(app.mesh, vertices[j][i]);
          add_control_point(app, spline, points[j][i]);
        }

        add_points_shape(app,
            {points[j][0], points[j][1], points[j][2], points[j][3]},
            0.003f * app.line_size, {0, 0, 1});
        auto L0 = my_compute_geodesic_path(
            app.mesh, points[j][0], points[j][1]);
        auto L1 = my_compute_geodesic_path(
            app.mesh, points[j][1], points[j][2]);
        auto L2 = my_compute_geodesic_path(
            app.mesh, points[j][2], points[j][3]);
        if (!point_are_close(L0.start, L0.end))
          add_path_shape(app, L0, 0.0010f * app.line_size, {0, 0, 1});
        if (!point_are_close(L1.start, L1.end))
          add_path_shape(app, L1, 0.0010f * app.line_size, {0, 0, 1});
        if (!point_are_close(L2.start, L2.end))
          add_path_shape(app, L2, 0.0010f * app.line_size, {0, 0, 1});

        spline.curves_to_update.insert(j);
        update_curve_shape(spline, j, app.mesh, app.bezier_params, false);
      }
    }

    // flipout
    if (app.show_flipout_curve) {
      for (auto j = 0; j < control_polygon.size(); ++j)
        for (int i = 0; i < 4; i++) {
          vertices[j][i] = vert_from_point(
              app.mesh.triangles, control_polygon[j][i]);
          // points[j][i] = make_point_from_vert(app.mesh, vertices[j][i]);
        }
      for (auto j = 0; j < vertices.size(); ++j) {
        try {
          auto bezier = flipout::compute_bezier_curve(
              *app.mesh.flipout, vertices[j], app.bezier_params.subdivisions);
          auto curr_pos = generate_polyline_from_positions(app.mesh,
              flipout::path_positions(bezier.get()), points[j][0].face);

          add_path_shape(
              app, curr_pos, 0.005f * app.curve_size, vec3f{1, 0, 0}, true);
          // if (app.show_WA_samples)
          //   app.WA_samples_shape = add_points_shape(
          //       app, app.WA_samples, 0.002, vec3f{0, 0, 0});
        } catch (std::exception& e) {
          printf("[flipout exception]: %s\n", e.what());
        }
      }
    }

    // for (auto path : app.added_paths) {
    //   auto pos = path_positions(app.mesh, path->path);
    //   update_path_shape(path->instance->shape, app.mesh, pos,
    //       0.005f * app.curve_size,
    //       app.avg_edge_length_factor * app.mesh.avg_edge_length);
    // }
  }
  continue_line(widget);
  draw_checkbox(widget, "OLR", app.show_olr);
  continue_line(widget);
  draw_checkbox(widget, "WA", app.show_wa_curve);
  // continue_line(widget);
  if (draw_checkbox(widget, "samples", app.show_WA_samples)) {
    //   auto it = std::find(app.added_points.begin(), app.added_points.end(),
    //       app.good_samples_shape);
    //   if (app.show_WA_samples && app.good_samples.size() > 0) {
    //     if (it == app.added_points.end())
    //       app.good_samples_shape = add_points_shape(
    //           app, app.good_samples, 0.003 * app.line_size, vec3f{0, 0, 0});
    //     else
    //       update_points_shape(app.good_samples_shape->instance->shape,
    //           app.good_samples, 0.003 * app.line_size);
    //   } else if (app.good_samples.size() > 0 && it != app.added_points.end())
    //     clear_shape(app.good_samples_shape->instance->shape);
    //   it = std::find(app.added_points.begin(), app.added_points.end(),
    //       app.bad_samples_shape);
    //   if (app.show_WA_samples && app.bad_samples.size() > 0) {
    //     if (it == app.added_points.end())
    //       app.bad_samples_shape = add_points_shape(
    //           app, app.bad_samples, 0.003 * app.line_size, vec3f{0, 0, 0});
    //     else
    //       update_points_shape(app.bad_samples_shape->instance->shape,
    //           app.bad_samples, 0.003 * app.line_size);
    //   } else if (app.bad_samples.size() > 0 && it != app.added_points.end())
    //     clear_shape(app.bad_samples_shape->instance->shape);
  }
  continue_line(widget);
  draw_checkbox(widget, "Flipout", app.show_flipout_curve);
  if (draw_button(widget, "Reset")) {
    if (app.show_olr) {
      for (auto& spline : app.splines()) {
        if (spline.control_points.size() > 0) clear_spline(spline);
      }
      app.splines().clear();
    }

    if (app.show_wa_curve) {
      for (auto path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      app.added_paths.clear();
      for (auto points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_points.clear();
    }
    app.show_points = true;
  }
  // if (draw_button(widget, "export anchors")) {
  //   std::ofstream outfile0;
  //   std::ofstream outfile1;
  //   outfile0.open(
  //       "../WA_code_modified/anchors_" + path_basename(app.filename) +
  //       ".txt");
  //   outfile1.open("WA_anchors/anchors_" + path_basename(app.filename) +
  //   ".txt"); for (auto spline : app.splines()) {
  //     for (auto i = 0; i < spline.curves.size(); ++i) {
  //       auto polygon = get_control_polygon(spline, i);
  //       for (auto j = 0; j < 4; ++j) {
  //         //   auto point = force_point_on_vert(app.mesh, polygon[j]);
  //         outfile0 << polygon[j].face << " ";
  //         outfile0 << polygon[j].uv.x << " ";
  //         outfile0 << polygon[j].uv.y << "\n";
  //         outfile1 << polygon[j].face << " ";
  //         outfile1 << polygon[j].uv.x << " ";
  //         outfile1 << polygon[j].uv.y << "\n";
  //       }
  //     }
  //   }
  //   outfile0.close();
  //   outfile1.close();
  // }

  draw_separator(widget);
  // draw_checkbox(widget, "Weighted Averages", app.weighted_average);
  // draw_combobox(widget, "Mode", (int&)app.type_of_average, average_names);
  // if (draw_button(widget, "Compute Average")) {
  //   auto a    = mesh_point{3891, {0.14215076, 0.592992663}};
  //   auto b    = mesh_point{3628, {0.21376653, 0.781679987}};
  //   auto path = compute_geodesic_path(app.mesh, a, b);
  //   //   switch (app.type_of_average) {
  //   //     case 0: {
  //   //       auto p = weighted_average(
  //   //           app.mesh, app.control_points, app.input().parameter_t);
  //   //       auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
  //   //     } break;
  //   //     case 1: {
  //   //       vector<float> w;
  //   //       auto          p = zero3f;
  //   //       bernstein_polynomials(
  //   //           app.control_points.size(), app.input().parameter_t, w);
  //   //       for (auto i = 0; i < app.control_points.size(); ++i) {
  //   //         p += w[i] * eval_position(app.mesh, app.control_points[i]);
  //   //       }
  //   //       auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
  //   //     } break;
  //   //   }
  // }
  draw_separator(widget);
  if (draw_combobox(widget, "algorithm", (int&)app.bezier_params.algorithm,
          spline_algorithm_names)) {
    update_all_splines(app);
  }
  draw_slider(widget, "subdivisions", app.bezier_params.subdivisions, 1, 10);
  draw_slider(widget, "precision", app.bezier_params.precision, 0, 1.f);

  if (draw_slider(widget, "t parameter", app.input().parameter_t, 0.0f, 1.0f)) {
    if (app.added_paths.size() != 0) {
      for (auto& path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      app.added_paths.clear();
    }
    if (app.added_points.size() != 0) {
      for (auto& points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_points.clear();
    }
    if (!app.weighted_average) {
      auto curve_id = (app.input().selected_control_point / 3);
      curve_id      = clamp(curve_id, 0, (int)app.spline().curves.size());
      auto polygon  = get_control_polygon(app.spline(), curve_id);
      auto curve    = app.spline().curves[curve_id];
      if (app.bezier_params.algorithm ==
              spline_algorithm::de_casteljau_uniform ||
          app.bezier_params.algorithm ==
              spline_algorithm::de_casteljau_adaptive) {
        if (!app.show_construction)
          app.input().eval_point = eval_bezier_point(
              app.mesh, polygon, app.input().parameter_t);
        else {
          app.show_points      = false;
          auto [paths, points] = dc_construction(
              app.mesh, polygon, app.input().parameter_t);
          Added_Path*   paths_shape;
          Added_Points* points_shape;
          auto          rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
                       app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});

          for (auto i = 0; i < points.size(); ++i) {
            if (i < paths.size() - 1) {
              for (auto j = 0; j < 3; ++j) {
                paths_shape = add_path_shape(
                    app, paths[i][j], 0.0015 * app.line_size, rgb);
              }
              points_shape = add_points_shape(
                  app, points[i], 0.0030 * app.line_size, rgb);
            } else if (i == paths.size() - 1) {
              for (auto j = 0; j < 3; ++j) {
                paths_shape = add_path_shape(
                    app, paths[i][j], 0.0015 * app.line_size, {1, 0, 0});
              }
              points_shape = add_points_shape(
                  app, points[i], 0.0030 * app.line_size, {1, 0, 0});
            } else
              points_shape = add_points_shape(
                  app, points[i], 0.0030 * app.line_size, {0, 0, 0});
          }
        }

      } else if (app.bezier_params.algorithm ==
                     spline_algorithm::subdivision_uniform ||
                 app.bezier_params.algorithm ==
                     spline_algorithm::subdivision_adaptive) {
        if (!app.show_construction)
          app.input().eval_point = eval_spline_point(
              app.mesh, polygon, app.bezier_params, app.input().parameter_t);
      } else if (app.bezier_params.algorithm ==
                 spline_algorithm::de_casteljau_classic) {
        app.show_points      = false;
        auto [paths, points] = dc_classic_construction(
            app.mesh, polygon, app.input().parameter_t);
        Added_Path*   paths_shape;
        Added_Points* points_shape;
        auto          rgb = zero3f;
        for (auto i = 0; i < paths.size(); ++i) {
          if (i < 3)
            rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
                app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
          else if (i < 5)
            rgb = hsv_to_rgb({app.first_level_h / 360.f,
                app.first_level_s / 100.f, app.first_level_v / 100.f});
          else
            rgb = hsv_to_rgb({app.second_level_h / 360.f,
                app.second_level_s / 100.f, app.second_level_v / 100.f});

          paths_shape = add_path_shape(
              app, paths[i], 0.0015 * app.line_size, rgb);
        }
        for (auto i = 0; i < points.size(); ++i) {
          if (i == 0)
            rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
                app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
          else if (i == 1)
            rgb = hsv_to_rgb({app.first_level_h / 360.f,
                app.first_level_s / 100.f, app.first_level_v / 100.f});
          else if (i == 2)
            rgb = hsv_to_rgb({app.second_level_h / 360.f,
                app.second_level_s / 100.f, app.second_level_v / 100.f});
          else
            rgb = zero3f;

          points_shape = add_points_shape(
              app, points[i], 0.003 * app.line_size, rgb);
        }
      }
    } else if (app.show_gradient) {
      auto t       = app.input().parameter_t;
      auto polygon = get_control_polygon(app.spline(), 0);
      auto w       = vector<float>(polygon.size());
      bernstein_polynomials((int)polygon.size(), t, w);
      vector<int> vertices(polygon.size());
      for (int i = 0; i < 4; i++) {
        vertices[i] = vert_from_point(app.mesh.triangles, polygon[i]);
      }
      if (app.karcher_fields.size() == 0)
        std::tie(app.karcher_grads, app.karcher_fields) = compute_karcher_grad(
            *app.mesh.flipout, vertices, app.mesh.solver, app.mesh.triangles,
            app.mesh.positions, app.mesh.normals, app.mesh.Grad);

      app.vector_field = karcher_grad(
          app.karcher_grads, app.karcher_fields, w, app.mesh.normals, true);
      update_glvector_field(app, app.vector_field,
          app.vector_size_factor * app.vector_size, "vector_field");
    } else {
      app.show_points = false;
      switch (app.type_of_average) {
        case 0: {
          auto p = weighted_average(
              app.mesh, app.control_points, app.input().parameter_t);
          auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
        } break;
        case 1: {
          vector<float> w;
          auto          p = zero3f;
          bernstein_polynomials(
              app.control_points.size(), app.input().parameter_t, w);
          for (auto i = 0; i < app.control_points.size(); ++i) {
            p += w[i] * eval_position(app.mesh, app.control_points[i]);
          }
          auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
        } break;
        case 2: {
          auto p = weighted_average(
              app.mesh, app.control_points, app.input().parameter_t);
          auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
          auto gamma   = compute_geodesic_path(
                app.mesh, app.control_points[0], app.control_points[1]);
          add_path_shape(app, gamma, 0.0005, {1, 0, 0});
        } break;
        case 3: {
          vector<float> w;
          auto          p = zero3f;
          bernstein_polynomials(
              app.control_points.size(), app.input().parameter_t, w);
          for (auto i = 0; i < app.control_points.size(); ++i) {
            p += w[i] * eval_position(app.mesh, app.control_points[i]);
          }
          auto p_shape = add_points_shape(app, {p}, 0.003, {1, 0, 0});
          auto gamma   = compute_geodesic_path(
                app.mesh, app.control_points[0], app.control_points[1]);
          add_path_shape(app, gamma, 0.0005, {1, 0, 0});
        } break;
      }
      auto new_control_points = add_points_shape(
          app, app.control_points, 0.003, {0, 0, 1});
      // app.eval_points = eval_spline_point(
      //     app.mesh, polygon, app.bezier_params, 20);
    }
  }
  // draw_checkbox(widget, "Quadric Curve", app.quadric_curve);

  draw_checkbox(widget, "Show Construction", app.show_construction);
  if (draw_button(widget, "Show Algorithm")) {
    auto curve_id = (app.input().selected_control_point / 3);
    curve_id      = clamp(curve_id, 0, (int)app.spline().curves.size());
    auto polygon  = get_control_polygon(app.spline(), curve_id);

    if (app.added_paths.size() != 0) {
      for (auto& path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      app.added_paths.clear();
    }
    if (app.added_points.size() != 0) {
      for (auto& points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_points.clear();
    }
    if (app.bezier_params.algorithm == spline_algorithm::subdivision_uniform) {
      // control polygon
      auto [paths, points] = LR_algorithm(app.mesh, polygon);
      auto rgb             = hsv_to_rgb({app.control_polygon_h / 360.f,
                      app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
      for (auto i = 0; i < paths[0].size(); ++i) {
        add_path_shape(app, paths[0][i], 0.001 * app.line_size, rgb);
      }
      add_points_shape(app, points[0], 0.0020 * app.line_size, rgb);

      // intermediate points
      rgb = hsv_to_rgb({app.first_level_h / 360.f, app.first_level_s / 100.f,
          app.first_level_v / 100.f});

      auto q0      = eval_path_point(app.mesh, paths[0][1], 5 / 8.f);
      auto q1      = eval_path_point(app.mesh, paths[0][2], 0.25);
      auto gamma01 = compute_geodesic_path(app.mesh, q0, q1);
      add_points_shape(app, {q0, q1}, 0.002 * app.line_size, rgb);
      add_path_shape(app, gamma01, 0.001 * app.line_size, rgb);

      for (auto i = 2; i < paths[0].size() - 3; ++i) {
        q0      = eval_path_point(app.mesh, paths[0][i], 0.75);
        q1      = eval_path_point(app.mesh, paths[0][i + 1], 0.25);
        gamma01 = compute_geodesic_path(app.mesh, q0, q1);
        add_points_shape(app, {q0, q1}, 0.0020 * app.line_size, rgb);
        add_path_shape(app, gamma01, 0.001 * app.line_size, rgb);
      }
      q0 = eval_path_point(app.mesh, paths[0][paths[0].size() - 3], 0.75);
      q1 = eval_path_point(app.mesh, paths[0][paths[0].size() - 2], 3 / 8.f);
      gamma01 = compute_geodesic_path(app.mesh, q0, q1);
      add_points_shape(app, {q0, q1}, 0.0020 * app.line_size, rgb);
      add_path_shape(app, gamma01, 0.001 * app.line_size, rgb);
      // last level
      rgb = {1, 0, 0};  // hsv_to_rgb({app.second_level_h / 360.f,
                        // app.second_level_s / 100.f,
      // app.second_level_v / 100.f
      //});
      for (auto i = 0; i < paths[1].size(); ++i) {
        add_path_shape(app, paths[1][i], 0.0015 * app.line_size, rgb);
      }
      add_points_shape(app, points[1], 0.0030 * app.line_size, rgb);
    } else if (app.bezier_params.algorithm ==
               spline_algorithm::de_casteljau_uniform) {
      auto curve_id = (app.input().selected_control_point / 3);
      curve_id      = clamp(curve_id, 0, (int)app.spline().curves.size());
      auto polygon  = get_control_polygon(app.spline(), curve_id);
      auto [paths, points] = RDC_algorithm(app.mesh, polygon);
      if (app.added_paths.size() != 0) {
        for (auto& path : app.added_paths) {
          clear_shape(path->instance->shape);
        }
        app.added_paths.clear();
      }
      if (app.added_points.size() != 0) {
        for (auto& points : app.added_points) {
          clear_shape(points->instance->shape);
        }

        app.added_points.clear();
      }
      // control polygon
      auto rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
          app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
      for (auto i = 0; i < paths[0].size(); ++i) {
        add_path_shape(app, paths[0][i], 0.001 * app.line_size, rgb);
      }
      add_points_shape(app, points[0], 0.0020 * app.line_size, rgb);

      // intermediate level
      rgb = hsv_to_rgb({app.third_level_h / 360.f, app.third_level_s / 100.f,
          app.third_level_v / 100.f});
      auto p12   = eval_path_midpoint(app.mesh, paths[0][1]);
      auto gamma = compute_geodesic_path(app.mesh, points[1][1], p12);
      add_path_shape(app, gamma, 0.001 * app.line_size, rgb);
      gamma = compute_geodesic_path(app.mesh, p12, points[2][2]);
      add_path_shape(app, gamma, 0.001 * app.line_size, rgb);
      add_points_shape(app, {p12}, 0.0020 * app.line_size, rgb);

      // left
      rgb = hsv_to_rgb({app.first_level_h / 360.f, app.first_level_s / 100.f,
          app.first_level_v / 100.f});

      for (auto i = 0; i < paths[1].size(); ++i) {
        add_path_shape(app, paths[1][i], 0.0015 * app.line_size, rgb);
      }
      auto S = points[1].back();
      points[1].pop_back();
      add_points_shape(app, points[1], 0.0030 * app.line_size, rgb);

      // right
      rgb = hsv_to_rgb({app.second_level_h / 360.f, app.second_level_s / 100.f,
          app.second_level_v / 100.f});
      for (auto i = 0; i < paths[2].size(); ++i) {
        add_path_shape(app, paths[2][i], 0.0015 * app.line_size, rgb);
      }
      add_points_shape(app, {points[2][1], points[2][2], points[2][3]},
          0.0030 * app.line_size, rgb);

      add_points_shape(app, {S}, 0.0030 * app.line_size, {1, 0, 0});
    }
  }
  draw_separator(widget);

  {
    auto id = get_selected_anchor_point(app);
    if (id != -1) {
      auto s = (bool)app.spline().is_smooth[id];
      draw_checkbox(widget, "smooth", s);
      app.spline().is_smooth[id] = s;
    }
  }
  //    draw_checkbox(widget, "enforce tangents",
  //    app.input().enforce_tangents);
  draw_checkbox(widget, "Use Vector Heat", app.use_vector_heat);
  // draw_checkbox(widget, "parallel", app.bezier_params.parallel);
  draw_checkbox(widget, "show edges", app.show_edges);
  // app.shade_params.faceted = app.show_edges;
  if (draw_checkbox(widget, "show points", app.show_control_points)) {
    if (app.show_control_points) {
      if (app.input().selected_control_point != -1) {
        auto selected = app.input().selected_control_point;
        auto pos      = vector<vec3f>{};
        if (is_handle_control_point(selected)) {
          if (selected % 3 == 1) selected -= 1;
          if (selected % 3 == 2) selected += 1;
        }
        auto prev_id = (selected / 3) - 1;
        if (prev_id >= 0) {
          pos.push_back(eval_position(
              app.mesh, app.spline().control_points[selected - 1]));
        }
        auto next_id = (selected / 3);
        if (next_id < app.spline().curves.size()) {
          auto& next_curve = app.spline().curves[next_id];
          draw_shape(next_curve.tangents[0].shape);
          pos.push_back(eval_position(
              app.mesh, app.spline().control_points[selected + 1]));
        }
        pos.push_back(
            eval_position(app.mesh, app.spline().control_points[selected]));

        add_points_shape(app, pos, 0.0030 * app.line_size, vec3f{1, 0, 0});
      }
    } else {
      for (auto& point : app.added_points) {
        clear_shape(point->instance->shape);
      }
      app.added_points.clear();
    }
  }
  if (draw_slider(widget, "line size", app.line_size, 0, 3.0f)) {
    if (app.added_paths.size() != 0) {
      for (auto i = 0; i < app.added_paths.size(); ++i) {
        auto& curr_path = app.added_paths[i];
        if (curr_path != app.WA_curve_shape)
          update_path_shape(curr_path->instance->shape, app.mesh,
              curr_path->path, 0.0015 * app.line_size);
      }
    }
    if (app.added_points.size() != 0) {
      for (auto i = 0; i < app.added_points.size(); ++i) {
        auto& curr_points = app.added_points[i];
        if (curr_points != app.WA_samples_shape && !app.show_WA_samples)
          update_points_shape(curr_points->instance->shape, app.mesh,
              curr_points->points, 0.0030 * app.line_size);
        else if (curr_points == app.WA_samples_shape &&
                 app.WA_samples.size() > 0)
          update_points_shape(app.WA_samples_shape->instance->shape,
              app.WA_samples, 0.0030 * app.line_size);
      }
    }
  }

  if (draw_slider(widget, "curve size", app.curve_size, 0, 1.0f)) {
    for (auto& spline : app.splines()) {
      for (auto& curve : spline.curves) {
        //        set_polyline_shape(curve.shape->shape, curve.positions);
        if (!curve.shape) {
          curve.shape = add_shape(app.scene);
          add_instance(
              app.scene, identity3x4f, curve.shape, app.spline_material, false);
        }
        update_path_shape(
            curve.shape, app.mesh, curve.positions, 0.003f * app.curve_size);
        set_polyline_shape(
            curve.tangents[0].shape, curve.tangents[0].positions);
        set_polyline_shape(
            curve.tangents[1].shape, curve.tangents[1].positions);
      }
      spline.curves_to_update.clear();
    }
    // if (app.WA_curve.size() > 0)
    //   update_path_shape(app.WA_curve_shape->instance->shape, app.mesh,
    //       app.WA_curve, 0.005f * app.curve_size,
    //       app.avg_edge_length_factor * app.mesh.avg_edge_length);
  }
  if (draw_slider(widget, "scale factor", app.scale_factor, 0.001, 1.f)) {
    if (app.enable_jumps) {
      for (auto& spline : app.splines()) {
        for (auto& curve : spline.curves) {
          //        set_polyline_shape(curve.shape->shape,
          // curve.positions);
          if (!curve.shape) {
            curve.shape = add_shape(app.scene);
            add_instance(app.scene, identity3x4f, curve.shape,
                app.spline_material, false);
          }
          update_path_shape(curve, app.mesh, 0.003f * app.curve_size,
              app.scale_factor, app.cut_after_first_jump);
          set_polyline_shape(
              curve.tangents[0].shape, curve.tangents[0].positions);
          set_polyline_shape(
              curve.tangents[1].shape, curve.tangents[1].positions);
        }
        spline.curves_to_update.clear();
      }
    }
  }
  if (draw_checkbox(widget, "show jumps", app.enable_jumps)) {
    if (app.enable_jumps) {
      for (auto& spline : app.splines()) {
        for (auto& curve : spline.curves) {
          //        set_polyline_shape(curve.shape->shape,
          // curve.positions);
          if (!curve.shape) {
            curve.shape = add_shape(app.scene);
            add_instance(app.scene, identity3x4f, curve.shape,
                app.spline_material, false);
          }
          update_path_shape(curve, app.mesh, 0.003f * app.curve_size,
              app.scale_factor, app.cut_after_first_jump);
          set_polyline_shape(
              curve.tangents[0].shape, curve.tangents[0].positions);
          set_polyline_shape(
              curve.tangents[1].shape, curve.tangents[1].positions);
        }
        spline.curves_to_update.clear();
      }
    } else {
      for (auto& spline : app.splines()) {
        for (auto& curve : spline.curves) {
          //        set_polyline_shape(curve.shape->shape,
          // curve.positions);
          if (!curve.shape) {
            curve.shape = add_shape(app.scene);
            add_instance(app.scene, identity3x4f, curve.shape,
                app.spline_material, false);
          }
          update_path_shape(
              curve.shape, app.mesh, curve.positions, 0.003f * app.curve_size);
          set_polyline_shape(
              curve.tangents[0].shape, curve.tangents[0].positions);
          set_polyline_shape(
              curve.tangents[1].shape, curve.tangents[1].positions);
        }
        spline.curves_to_update.clear();
      }
    }
  }
  continue_line(widget);
  draw_checkbox(widget, "cut at jump", app.cut_after_first_jump);
  if (draw_checkbox(widget, "show control polygon", app.show_control_polygon)) {
    if (app.show_control_polygon) {
      app.show_points = false;
      if (app.added_paths.size() != 0) {
        for (auto& path : app.added_paths) {
          clear_shape(path->instance->shape);
        }
        for (auto& points : app.added_points) {
          clear_shape(points->instance->shape);
        }

        app.added_paths.clear();
        app.added_points.clear();
      }
      for (auto& spline : app.splines()) {
        for (auto i = 0; i < spline.curves.size(); ++i) {
          auto polygon = get_control_polygon(spline, i);
          // auto p0              = curve.tangents[0].path.start;
          // auto p1              = curve.tangents[0].path.end;
          // auto p2              = curve.tangents[1].path.end;
          // auto p3              = curve.tangents[1].path.start;
          auto rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
              app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
          if (!app.quadric_curve) {
            auto control_polygon = add_points_shape(app,
                {polygon[0], polygon[1], polygon[2], polygon[3]},
                0.003f * app.line_size, rgb);
            auto L0 = compute_geodesic_path(app.mesh, polygon[0], polygon[1]);
            auto L1 = compute_geodesic_path(app.mesh, polygon[1], polygon[2]);
            auto L2 = compute_geodesic_path(app.mesh, polygon[2], polygon[3]);
            if (!point_are_close(L0.start, L0.end))
              auto L0_curve = add_path_shape(
                  app, L0, 0.0015f * app.line_size, rgb);
            if (!point_are_close(L1.start, L1.end))
              auto L1_curve = add_path_shape(
                  app, L1, 0.0015f * app.line_size, rgb);
            if (!point_are_close(L2.start, L2.end))
              auto L2_curve = add_path_shape(
                  app, L2, 0.0015f * app.line_size, rgb);
          } else {
            auto control_polygon = add_points_shape(app,
                {polygon[0], polygon[1], polygon[3]}, 0.003f * app.line_size,
                rgb);
            auto L0 = compute_geodesic_path(app.mesh, polygon[0], polygon[1]);
            auto L1 = compute_geodesic_path(app.mesh, polygon[1], polygon[3]);
            if (!point_are_close(L0.start, L0.end))
              auto L0_curve = add_path_shape(
                  app, L0, 0.0015f * app.line_size, rgb);
            if (!point_are_close(L1.start, L1.end))
              auto L1_curve = add_path_shape(
                  app, L1, 0.0015f * app.line_size, rgb);
          }
        }
      }
    } else {
      app.show_points = true;
      for (auto& path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      for (auto& points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_paths.clear();
      app.added_points.clear();
    }
  }
  bool update_coords = false;
  if (app.input().selected_control_point != -1) {
    update_coords |= draw_checkbox(
        widget, "translate", app.input().translating);
    update_coords |= draw_checkbox(widget, "rotate", app.input().rotating);
    update_coords |= draw_checkbox(widget, "scale", app.input().scaling);
    if (update_coords) {
      update_local_coordinates(
          app.spline(), app.input().selected_control_point, app.mesh);
    }
  }
  if (app.temp_points.size())
    draw_slider(
        widget, "levels", app.temp_levels, 0, app.temp_points.size() - 1);

  if (is_anchor_control_point(app.input().selected_control_point) &&
      app.input().selected_control_point != -1) {
    if (draw_button(widget, "Delete point")) {
      auto  selected = app.input().selected_control_point;
      auto& points   = app.spline().control_points;
      points.erase(
          points.begin() + selected - 1, points.begin() + selected + 2);
      app.spline().curves.pop_back();
      app.spline().curves_to_update = make_set(app.spline().curves.size());
    }
  }
  // Load svg
  if (app.xxx_point.face != -1) {
    static string svg_filename = "test.svg";
    static string svg_dirname  = "data/";
    if (draw_filedialog_button(widget, "svg file", true, "svg file2",
            svg_filename, false, svg_dirname, "test", "*.svg")) {
      auto svg = load_svg(svg_filename);
      init_from_svg(app, svg);
    }
    // if (draw_filedialog_button(widget, "svg file (weights)", true, "svg
    // file2",
    //         svg_filename, false, svg_dirname, "test", "*.svg")) {
    //   auto svg = load_svg(svg_filename);
    //   init_from_svg_weights(app, svg);
    // }
  }
  draw_slider(widget, "OLR STEPS", app.OLR_steps, 1, 3);
  if (draw_button(widget, "Show OLR construction")) {
    auto curve_id = (app.input().selected_control_point / 3);
    curve_id      = clamp(curve_id, 0, (int)app.spline().curves.size());
    auto polygon  = get_control_polygon(app.spline(), curve_id);
    if (app.added_paths.size() != 0) {
      for (auto& path : app.added_paths) {
        clear_shape(path->instance->shape);
      }
      app.added_paths.clear();
    }
    if (app.added_points.size() != 0) {
      for (auto& points : app.added_points) {
        clear_shape(points->instance->shape);
      }

      app.added_points.clear();
    }
    auto prev = vector<mesh_point>{
        polygon[0], polygon[1], polygon[2], polygon[3]};
    auto curr = OLR(app.mesh, prev, 1);
    for (auto i = 2; i < app.OLR_steps; ++i) {
      prev = curr;
      curr = OLR(app.mesh, prev, i);
    }

    add_points_shape(app, prev, 0.0035 * app.line_size, vec3f{0, 0, 1});
    add_points_shape(app, curr, 0.0035 * app.line_size, vec3f{1, 0, 0});
    for (auto i = 0; i < prev.size() - 1; ++i) {
      auto path = compute_geodesic_path(app.mesh, prev[i], prev[i + 1]);
      add_path_shape(app, path, 0.003f * app.curve_size, vec3f{0, 0, 1});
    }

    for (auto i = 0; i < curr.size() - 1; ++i) {
      auto path = compute_geodesic_path(app.mesh, curr[i], curr[i + 1]);
      add_path_shape(app, path, 0.003f * app.curve_size, vec3f{1, 0, 0});
    }
  }
  // draw_slider(widget, "svg size", app.svg_size, 0, 1.0f);
  // draw_checkbox(widget, "EnvLight", app.envlight);
  // // if (draw_button(widget, "export scene")) {
  // //   export_yocto_scene(app, app.scene_name);
  // // }
  // continue_line(widget);

  // draw_separator(widget);
  // if (draw_combobox(widget, "DC Levels", (int&)app.level, levels_names)) {
  //   switch (app.level) {
  //     case 0: {
  //       app.hue   = app.control_polygon_h;
  //       app.sat   = app.control_polygon_s;
  //       app.value = app.control_polygon_v;
  //     } break;
  //     case 1: {
  //       app.hue   = app.first_level_h;
  //       app.sat   = app.first_level_s;
  //       app.value = app.first_level_v;
  //     } break;
  //     case 2: {
  //       app.hue   = app.second_level_h;
  //       app.sat   = app.second_level_s;
  //       app.value = app.second_level_v;
  //     } break;
  //   }
  // }

  // if (draw_slider(widget, "H", app.hue, 0, 360)) {
  //   if (app.added_paths.size() != 0) {
  //     if (app.level != 4) {
  //       auto& curr_point = app.added_points[app.level];
  //       auto  rgb        = hsv_to_rgb(
  //                   {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //       set_color(curr_point->instance->material, rgb);
  //     }

  //     switch (app.level) {
  //       case 0: {
  //         for (auto i = 0; i < 3; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.control_polygon_h = app.hue;
  //           set_color(curr_path->instance->material, rgb);
  //         }

  //       } break;
  //       case 1: {
  //         for (auto i = 3; i < 5; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.first_level_h = app.hue;
  //           set_color(curr_path->instance->material, rgb);
  //         }
  //       } break;

  //       case 2: {
  //         auto& curr_path = app.added_paths[5];
  //         auto  rgb       = hsv_to_rgb(
  //                    {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //         app.second_level_h = app.hue;
  //         set_color(curr_path->instance->material, rgb);

  //       } break;
  //         // case 3: {
  //         //   auto program = &app.shaders.at("polyline");
  //         //   bind_program(program);
  //         //   set_common_uniforms(app, program);
  //         //   set_uniform(program, "size", 3 * 0.0015f * app.line_size);
  //         //   for (auto& spline : app.splines()) {
  //         //     for (int i = 0; i < spline.curves.size(); i++) {
  //         //       auto& curve = spline.curves[i];
  //         //       set_uniform(
  //         //           program, "color", vec3f{app.red, app.green,
  //         app.blue});
  //         //       draw_shape(curve.shape->shape);
  //         //     }
  //         //   }
  //         // } break;
  //     }
  //   }
  // }
  // if (draw_slider(widget, "S", app.sat, 0, 100)) {
  //   if (app.added_paths.size() != 0) {
  //     if (app.level != 4) {
  //       auto& curr_point = app.added_points[app.level];
  //       auto  rgb        = hsv_to_rgb(
  //                   {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //       set_color(curr_point->instance->material, rgb);
  //     }
  //     switch (app.level) {
  //       case 0: {
  //         for (auto i = 0; i < 3; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.control_polygon_s = app.sat;
  //           set_color(curr_path->instance->material, rgb);
  //         }

  //       } break;
  //       case 1: {
  //         for (auto i = 3; i < 5; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.first_level_s = app.sat;
  //           set_color(curr_path->instance->material, rgb);
  //         }
  //       } break;

  //       case 2: {
  //         auto& curr_path = app.added_paths[5];
  //         auto  rgb       = hsv_to_rgb(
  //                    {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //         app.second_level_s = app.sat;
  //         set_color(curr_path->instance->material, rgb);

  //       } break;
  //         // case 3: {
  //         //   auto program = &app.shaders.at("polyline");
  //         //   bind_program(program);
  //         //   set_common_uniforms(app, program);
  //         //   set_uniform(program, "size", 0.0015f * app.line_size);
  //         //   for (auto& spline : app.splines()) {
  //         //     for (int i = 0; i < spline.curves.size(); i++) {
  //         //       auto& curve = spline.curves[i];
  //         //       set_uniform(
  //         //           program, "color", vec3f{app.red, app.green,
  //         app.blue});
  //         //       draw_shape(curve.shape->shape);
  //         //     }
  //         //   }
  //         // } break;
  //     }
  //   }
  // }
  // if (draw_slider(widget, "V", app.value, 0, 100)) {
  //   if (app.added_paths.size() != 0) {
  //     if (app.level != 4) {
  //       auto& curr_point = app.added_points[app.level];
  //       auto  rgb        = hsv_to_rgb(
  //                   {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //       set_color(curr_point->instance->material, rgb);
  //     }
  //     switch (app.level) {
  //       case 0: {
  //         for (auto i = 0; i < 3; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.control_polygon_v = app.value;
  //           set_color(curr_path->instance->material, rgb);
  //         }

  //       } break;
  //       case 1: {
  //         for (auto i = 3; i < 5; ++i) {
  //           auto& curr_path = app.added_paths[i];
  //           auto  rgb       = hsv_to_rgb(
  //                      {app.hue / 360.f, app.sat / 100.f, app.value /
  //                      100.f});
  //           app.first_level_v = app.value;
  //           set_color(curr_path->instance->material, rgb);
  //         }
  //       } break;

  //       case 2: {
  //         auto& curr_path = app.added_paths[5];
  //         auto  rgb       = hsv_to_rgb(
  //                    {app.hue / 360.f, app.sat / 100.f, app.value / 100.f});
  //         app.second_level_v = app.value;
  //         set_color(curr_path->instance->material, rgb);

  //       } break;
  //         // case 3: {
  //         //   auto program = &app.shaders.at("polyline");
  //         //   bind_program(program);
  //         //   set_common_uniforms(app, program);
  //         //   set_uniform(program, "size", 0.0015f * app.line_size);
  //         //   for (auto& spline : app.splines()) {
  //         //     for (int i = 0; i < spline.curves.size(); i++) {
  //         //       auto& curve = spline.curves[i];
  //         //       set_uniform(
  //         //           program, "color", vec3f{app.red, app.green,
  //         app.blue});
  //         //       draw_shape(curve.shape->shape);
  //         //     }
  //         //   }
  //         // } break;
  //     }
  //   }
  // }

  // draw_coloredit(widget, "mesh color", app.mesh_material->color);
  // draw_slider(widget, "roughness", app.mesh_material->roughness, 0, 1.0f);
  end_widget(widget);
};

int main(int num_args, const char* args[]) {
  auto app = App();

  bool   time       = true;
  string test       = "";
  bool   infolog    = true;
  bool   quiet      = false;
  bool   log_colors = true;
  string playback   = "";
  int    msaa       = 1;

  auto cli = make_cli("bezier", "interactive viewer for mesh processing");
  add_option(cli, "mesh", app.filename, "Model filenames", true);
  add_option(cli, "--test", app.testname, "Test filename", false);
  add_option(cli, "--time/--no-time", time, "Log times");
  add_option(cli, "--info/--no-info", infolog, "Log info");
  add_option(cli, "--envlight", app.envlight, "Enable environmental lighting");
  add_option(cli, "--quiet", quiet, "Disable logs");
  add_option(cli, "--colors/--no-colors", log_colors, "Colored logs");
  add_option(cli, "--msaa", msaa, "OpenGL multisample anti-aliasing");
  add_option(cli, "--playback", playback, "Playback recorded input session");
  parse_cli(cli, num_args, args);

  // Load model and init bvh for fast click-intersection.
  if (!load_mesh(app.filename, app.mesh, app.error)) print_fatal(app.error);
  app.mesh.flipout  = new flipout::flipout_mesh{};
  *app.mesh.flipout = make_flipout_mesh(app.mesh.triangles, app.mesh.positions);

  app.comparison = path_basename(app.filename);

  init_bvh(app);

  // Init window.
  auto win  = new gui_window();
  win->msaa = msaa;
  init_window(win, {1080, 720}, "mesh viewer", true);
  win->user_data = &app;

  init_gpu(app, app.envlight);

  init_widget(app.widget, win);
  load_editing(app, app.testname);

  app.show_points = false;
  if (app.added_paths.size() != 0) {
    for (auto& path : app.added_paths) {
      clear_shape(path->instance->shape);
    }
    for (auto& points : app.added_points) {
      clear_shape(points->instance->shape);
    }

    app.added_paths.clear();
    app.added_points.clear();
  }
  for (auto& spline : app.splines()) {
    for (auto i = 0; i < spline.curves.size(); ++i) {
      auto polygon = get_control_polygon(spline, i);
      // auto p0              = curve.tangents[0].path.start;
      // auto p1              = curve.tangents[0].path.end;
      // auto p2              = curve.tangents[1].path.end;
      // auto p3              = curve.tangents[1].path.start;
      auto rgb = hsv_to_rgb({app.control_polygon_h / 360.f,
          app.control_polygon_s / 100.f, app.control_polygon_v / 100.f});
      if (!app.quadric_curve) {
        auto control_polygon = add_points_shape(app,
            {polygon[0], polygon[1], polygon[2], polygon[3]},
            0.003f * app.line_size, rgb);
        auto L0 = compute_geodesic_path(app.mesh, polygon[0], polygon[1]);
        auto L1 = compute_geodesic_path(app.mesh, polygon[1], polygon[2]);
        auto L2 = compute_geodesic_path(app.mesh, polygon[2], polygon[3]);
        if (!point_are_close(L0.start, L0.end))
          auto L0_curve = add_path_shape(app, L0, 0.0015f * app.line_size, rgb);
        if (!point_are_close(L1.start, L1.end))
          auto L1_curve = add_path_shape(app, L1, 0.0015f * app.line_size, rgb);
        if (!point_are_close(L2.start, L2.end))
          auto L2_curve = add_path_shape(app, L2, 0.0015f * app.line_size, rgb);
      } else {
        auto control_polygon = add_points_shape(app,
            {polygon[0], polygon[1], polygon[3]}, 0.003f * app.line_size, rgb);
        auto L0 = compute_geodesic_path(app.mesh, polygon[0], polygon[1]);
        auto L1 = compute_geodesic_path(app.mesh, polygon[1], polygon[3]);
        if (!point_are_close(L0.start, L0.end))
          auto L0_curve = add_path_shape(app, L0, 0.0015f * app.line_size, rgb);
        if (!point_are_close(L1.start, L1.end))
          auto L1_curve = add_path_shape(app, L1, 0.0015f * app.line_size, rgb);
      }
    }
  }

  if (msaa > 1) set_ogl_msaa();

  if (playback.size()) {
    load_input_record(app.input_record, "apps/splinegui/test.recording");
    app.playback = true;
  }

  run_ui(win, draw);

  // TODO(giacomo): delete app
  clear_window(win);
}
