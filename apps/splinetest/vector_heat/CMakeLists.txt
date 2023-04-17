add_library(vector_heat
  vector_heat.h  vector_heat.cpp
)



target_include_directories(vector_heat PRIVATE .)
target_include_directories(vector_heat PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_include_directories(vector_heat PUBLIC ${CMAKE_SOURCE_DIR}/libs/geometry-central/include)

set_target_properties(vector_heat PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_link_libraries(vector_heat yocto geometry-central splinesurf)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(vector_heat Threads::Threads)
endif(UNIX AND NOT APPLE)

# warning flags
if(APPLE)
  target_compile_options(vector_heat PRIVATE -Wall)
endif(APPLE)
if(MSVC)
target_compile_options(vector_heat PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
target_compile_options(vector_heat PUBLIC /EHsc)
endif(MSVC)