#pragma once

#include "spline.h"

bool load_mesh(const string& filename, bezier_mesh& mesh, string& error);

bool load_bezier_params(const string& filename, vector<mesh_point>& points,
    bezier_params& params, string& error);
bool save_bezier_params(const string& filename,
    const vector<mesh_point>& points, const bezier_params& params,
    string& error);

using Svg_Path = vector<array<vec2f, 4>>;
struct Svg_Shape {
  vec3f            color = {};
  vector<Svg_Path> paths = {};
};
using Svg = vector<Svg_Shape>;

Svg load_svg(const string& filename);

#include "ext/json.hpp"

namespace yocto {

using json = nlohmann::json;
using std::array;

// support for json conversions
inline void to_json(json& j, const vec2f& value) {
  nlohmann::to_json(j, (const array<float, 2>&)value);
}

inline void from_json(const json& j, vec2f& value) {
  nlohmann::from_json(j, (array<float, 2>&)value);
}

// support for json conversions
inline void to_json(json& js, const mesh_point& value) {
  js["face"] = value.face;
  js["uv"]   = value.uv;
}

inline void from_json(const json& js, mesh_point& value) {
  js.at("face").get_to(value.face);
  js.at("uv").get_to(value.uv);
}

}  // namespace yocto
