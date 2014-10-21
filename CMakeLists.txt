cmake_minimum_required(VERSION 2.8)

project(HeightRidgeExtractor)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(GradientFilterTest gradient_filter_test.cc)

target_link_libraries(GradientFilterTest ${VTK_LIBRARIES})
