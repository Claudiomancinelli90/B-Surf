cmake_minimum_required (VERSION 3.9)
project (graphics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE “Release”)
endif()
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(graphics
    graphics/math.h
    graphics/common.h    graphics/commonio.h 
    graphics/modelio.h   graphics/modelio.cpp
    graphics/image.h     graphics/image.cpp
    graphics/bvh.h       graphics/bvh.cpp
    graphics/geometry.h  graphics/geometry.cpp
    # graphics/rendering.h graphics/rendering.cpp
)

target_include_directories(graphics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_definitions(graphics PUBLIC GRAPHICS_VERSION=1)

# 
#       How to link this library!!!
# 
# target_link_libraries(main opengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
