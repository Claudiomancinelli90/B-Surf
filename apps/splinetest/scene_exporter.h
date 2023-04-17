sceneio_shape* add_shape(sceneio_scene* scene, const string& name,
    const quads_shape& shape_data, int subdivisions = 0, float displacement = 0,
    sceneio_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->quads            = shape_data.quads;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
sceneio_shape* add_shape(sceneio_scene* scene, const string& name,
    const quads_fvshape& shape_data, int subdivisions = 0,
    float displacement = 0, sceneio_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->quadspos         = shape_data.quadspos;
  shape->quadsnorm        = shape_data.quadsnorm;
  shape->quadstexcoord    = shape_data.quadstexcoord;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
sceneio_shape* add_shape(sceneio_scene* scene, const string& name,
    const triangles_shape& shape_data, int subdivisions = 0,
    float displacement = 0, sceneio_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->triangles        = shape_data.triangles;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
sceneio_shape* add_shape(sceneio_scene* scene, const string& name,
    const lines_shape& shape_data, int subdivisions = 0, float displacement = 0,
    sceneio_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->lines            = shape_data.lines;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->radius           = shape_data.radius;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
sceneio_shape* add_shape(sceneio_scene* scene, const string& name,
    const points_shape& shape_data, int subdivisions = 0,
    float displacement = 0, sceneio_texture* displacement_tex = nullptr) {
  auto shape              = add_shape(scene, name);
  shape->points           = shape_data.points;
  shape->positions        = shape_data.positions;
  shape->normals          = shape_data.normals;
  shape->texcoords        = shape_data.texcoords;
  shape->radius           = shape_data.radius;
  shape->subdivisions     = subdivisions;
  shape->smooth           = subdivisions > 0 || displacement_tex;
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
  return shape;
}
sceneio_material* add_emission_material(sceneio_scene* scene,
    const string& name, const vec3f& emission, sceneio_texture* emission_tex) {
  auto material          = add_material(scene, name);
  material->emission     = emission;
  material->emission_tex = emission_tex;
  return material;
}

sceneio_instance* add_instance(sceneio_scene* scene, const string& name,
    const frame3f& frame, sceneio_shape* shape, sceneio_material* material) {
  auto instance      = add_instance(scene, name);
  instance->frame    = frame;
  instance->shape    = shape;
  instance->material = material;
  return instance;
}