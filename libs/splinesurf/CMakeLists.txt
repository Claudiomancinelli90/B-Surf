add_library(splinesurf
  spline.h  spline.cpp
  splineio.h  splineio.cpp
  karcher.h karcher.cpp
  logging.h
)

target_include_directories(splinesurf PRIVATE .)
target_include_directories(splinesurf PUBLIC ${CMAKE_SOURCE_DIR}/libs)

set_target_properties(splinesurf PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_link_libraries(splinesurf yocto geometry-central)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(splinesurf Threads::Threads)
endif(UNIX AND NOT APPLE)

# warning flags
if(APPLE)
  target_compile_options(splinesurf PRIVATE -Wall)
endif(APPLE)
if(MSVC)
target_compile_options(splinesurf PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
target_compile_options(splinesurf PUBLIC /EHsc)
endif(MSVC)