add_library(flipout
  flipout.h  flipout.cpp
)

# add_subdirectory(ext/flipout)

target_include_directories(flipout PRIVATE .)
target_include_directories(flipout PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_include_directories(flipout PUBLIC ${CMAKE_SOURCE_DIR}/libs/geometry-central/include)

set_target_properties(flipout PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_link_libraries(flipout yocto geometry-central)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(flipout Threads::Threads)
endif(UNIX AND NOT APPLE)

# warning flags
if(APPLE)
  target_compile_options(flipout PRIVATE -Wall)
endif(APPLE)
if(MSVC)
target_compile_options(flipout PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
target_compile_options(flipout PUBLIC /EHsc)
endif(MSVC)