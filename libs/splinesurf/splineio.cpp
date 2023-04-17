#include "splineio.h"

#include <yocto/yocto_commonio.h>

#include "ext/json.hpp"

using json = nlohmann::json;

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------

inline bool save_json(const string& filename, const json& js, string& error) {
  return save_text(filename, js.dump(2), error);
}

inline bool load_json(const string& filename, json& js, string& error) {
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
    return true;
  } catch (std::exception& e) {
    return parse_error();
  }

  return true;
}

// support for json conversions
inline void to_json(json& js, const spline_algorithm& value) {
  js = spline_algorithm_names[(int)value];
}

inline void from_json(const json& js, spline_algorithm& value) {
  auto vs = js.get<string>();
  for (auto i = 0; i < spline_algorithm_names.size(); i++) {
    if (vs == spline_algorithm_names[i]) {
      value = (spline_algorithm)i;
      return;
    }
  }
  throw std::invalid_argument{"unknown spline_algorithm"};
}

bool load_bezier_params(const string& filename, vector<mesh_point>& points,
    bezier_params& params, string& error) {
  auto js = json{};
  if (!load_json(filename, js, error)) return false;
  try {
    points              = js["points"].get<vector<mesh_point>>();
    params.algorithm    = js["params"]["algorithm"];
    params.subdivisions = js["params"]["subdivisions"];
    params.precision    = js["params"]["precision"];
  } catch (std::exception& e) {
    error = "error parsing json: ("s + e.what() + ")"s;
    return false;
  }
  return true;
}

bool save_bezier_params(const string& filename,
    const vector<mesh_point>& control_points, const bezier_params& params,
    string& error) {
  auto js                      = json{};
  js["points"]                 = control_points;
  js["params"]["algorithm"]    = params.algorithm;
  js["params"]["subdivisions"] = params.subdivisions;
  js["params"]["precision"]    = params.precision;
  if (!save_json(filename, js, error)) return false;
  return true;
}

#define HEAVY 1
#if HEAVY
#include "karcher.h"
void init_mesh(bezier_mesh& mesh, bool with_opposite) {
  mesh.v2t = vertex_to_triangles(
      mesh.triangles, mesh.positions, mesh.adjacencies);
  mesh.solver = make_geodesic_solver(
      mesh.triangles, mesh.positions, mesh.adjacencies, mesh.v2t);
  mesh.angles = compute_angles(mesh.triangles, mesh.positions, mesh.adjacencies,
      mesh.v2t, mesh.total_angles, with_opposite);
  mesh.Grad   = init_riemannian_gradient_matrix(
        mesh.solver, mesh.angles, mesh.positions, mesh.normals);
  mesh.avg_edge_length = avg_edge_length(mesh.positions, mesh.triangles);
}
#endif

#include "stl_reader.h"

void load_mesh_stl(bezier_mesh& mesh, const string& filename) {
  std::vector<float>        coords, normals;
  std::vector<unsigned int> tris, solids;

  stl_reader::ReadStlFile(filename.c_str(), coords, normals, tris, solids);
  const size_t numTris = tris.size() / 3;
  for (size_t itri = 0; itri < numTris; ++itri) {
    for (size_t icorner = 0; icorner < 3; ++icorner) {
      float* c = &coords[3 * tris[3 * itri + icorner]];
    }

    float* n = &normals[3 * itri];
  }
  mesh.positions.assign(
      (vec3f*)coords.data(), (vec3f*)(&coords[0] + coords.size()));
  for (auto& p : mesh.positions) {
    std::swap(p.y, p.z);
    std::swap(p.x, p.z);
  }
  mesh.triangles.assign((vec3i*)tris.data(), (vec3i*)(&tris[0] + tris.size()));
  mesh.normals.assign(
      (vec3f*)normals.data(), (vec3f*)(&normals[0] + normals.size()));
}

bool load_mesh(const string& filename, bezier_mesh& mesh, string& error) {
  mesh = bezier_mesh{};

  vector<int>   points;
  vector<vec2i> lines;
  vector<vec3i> triangles;
  vector<vec4i> quads;
  vector<vec4i> quadspos;
  vector<vec4i> quadsnorm;
  vector<vec4i> quadstexcoord;
  vector<vec3f> positions;
  vector<vec3f> normals;
  vector<vec4f> colors;
  vector<float> radius;

#if 0
  mesh.positions = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  mesh.triangles = {{0, 1, 2}, {0, 2, 3}};
  mesh.normals   = {{0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
  mesh.texcoords = {{0, 1}, {1, 1}, {1, 0}, {0, 0}};
#else

  auto ext = path_extension(filename);
  if (ext == ".stl") {
    load_mesh_stl(mesh, filename);
  } else {
    if (!load_shape(filename, points, lines, mesh.triangles, quads, quadspos,
            quadsnorm, quadstexcoord, mesh.positions, normals, mesh.texcoords,
            colors, radius, error, false)) {
      return false;
    }
  }
  if (quads.size()) {
    mesh.triangles = quads_to_triangles(quads);
  }
  printf("%s: mesh has %ld  triangle\n", __FUNCTION__, mesh.triangles.size());
#endif

  // bumped_sphere(0.0001f, mesh.positions);
  // Normalize positions in the cube [-1, 1]^3
  auto bbox = invalidb3f;
  for (auto& p : mesh.positions) bbox = merge(bbox, p);
  auto center = (bbox.max + bbox.min) / 2;
  auto scale  = 1.0f / max(bbox.max - bbox.min);
  for (auto& p : mesh.positions) p = (p - center) * scale;

  mesh.normals     = compute_normals(mesh.triangles, mesh.positions);
  mesh.adjacencies = face_adjacencies(mesh.triangles);
  mesh.dual_solver = make_dual_geodesic_solver(
      mesh.triangles, mesh.positions, mesh.adjacencies);
#if HEAVY
  init_mesh(mesh, true);
#endif
  return true;
}

#define NANOSVG_ALL_COLOR_KEYWORDS
#define NANOSVG_IMPLEMENTATION
#include "../nanosvg/src/nanosvg.h"

Svg load_svg(const string& filename) {
  struct NSVGimage* image;
  image = nsvgParseFromFile(filename.c_str(), "px", 96);
  printf("size: %f x %f\n", image->width, image->height);
  auto size = vec2f{image->width, image->height};

  // Use...
  auto svg = Svg();
  for (auto shape = image->shapes; shape != NULL; shape = shape->next) {
    auto&        svg_shape = svg.emplace_back();
    unsigned int c;
    if (shape->fill.type == NSVG_PAINT_COLOR) {
      c = shape->fill.color;
    } else if (shape->fill.type >= NSVG_PAINT_LINEAR_GRADIENT) {
      c = shape->fill.gradient->stops[0].color;
    } else {
      c = 0;
    }
    float r         = ((c >> 16) & 0xFF) / 255.0;  // Extract the RR byte
    float g         = ((c >> 8) & 0xFF) / 255.0;   // Extract the GG byte
    float b         = ((c)&0xFF) / 255.0;
    svg_shape.color = yocto::pow(vec3f{b, g, r}, 2.2f);

    for (auto path = shape->paths; path != NULL; path = path->next) {
      auto& svg_path = svg_shape.paths.emplace_back();
      for (int i = 0; i < path->npts - 1; i += 3) {
        float* p     = &path->pts[i * 2];
        auto&  curve = svg_path.emplace_back();
        curve[0]     = vec2f{p[0], size.y - p[1]} / size.y;
        curve[1]     = vec2f{p[2], size.y - p[3]} / size.y;
        curve[2]     = vec2f{p[4], size.y - p[5]} / size.y;
        curve[3]     = vec2f{p[6], size.y - p[7]} / size.y;
        // printf("(%f %f) (%f %f) (%f %f) (%f %f)\n", curve[0].x, curve[0].y,
        //     curve[1].x, curve[1].y, curve[2].x, curve[2].y, curve[3].x,
        //     curve[3].y);
      }
    }
  }
  // Delete
  nsvgDelete(image);
  return svg;
}
