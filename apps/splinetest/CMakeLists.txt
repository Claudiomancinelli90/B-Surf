add_executable(splinetest 
vector_heat/vector_heat.h vector_heat/vector_heat.cpp
flipout/flipout.h flipout/flipout.cpp
splinetest.cpp app.h app.cpp)

set_target_properties(splinetest     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(splinetest    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_include_directories(splinetest PUBLIC ${CMAKE_SOURCE_DIR}/libs/geometry-central/include)

target_link_libraries(splinetest yocto splinesurf yocto_gui realtime)
