if(YOCTO_OPENGL)
  add_definitions(-DYOCTO_OPENGL)
  find_package(OpenGL REQUIRED)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory(ext/glfw glfw)

  add_library(realtime
   gpu.h gpu.cpp
  window.h window.cpp
  gui.h gui.cpp
  # font.h font.cpp
  ext/imgui/imgui.cpp ext/imgui/imgui_draw.cpp ext/imgui/imgui_widgets.cpp
  ext/imgui/imgui_impl_glfw.cpp ext/imgui/imgui_impl_opengl3.cpp ext/glad/glad.c)

  target_include_directories(realtime PRIVATE .)
  target_include_directories(realtime PUBLIC ${CMAKE_SOURCE_DIR}/libs)
  target_compile_definitions(realtime PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
  set_target_properties(realtime PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  
  target_link_libraries(realtime glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(realtime dl)
  endif(UNIX AND NOT APPLE)
endif(YOCTO_OPENGL)

# warning flags
if(APPLE)
target_compile_options(realtime PRIVATE -Wall)
endif(APPLE)
if(MSVC)
target_compile_options(realtime PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
target_compile_options(realtime PUBLIC /EHsc)
endif(MSVC)

# if(APPLE)
#     target_compile_options(realtime PRIVATE -Wall)
#     include_directories(/usr/local/include /usr/local/include/freetype2)
#     link_directories(/usr/local/lib)
#     find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
#     set(GL_EXTRA_LIBRARIES)
# endif(APPLE)
