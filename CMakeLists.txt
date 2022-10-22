# CMakeList.txt: DemoForGeneticAlgorithm 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("DemoForGeneticAlgorithm")

# 将源代码添加到此项目的可执行文件。
file(GLOB_RECURSE source_files src/*.cpp)
add_executable (DemoForGeneticAlgorithm ${source_files} )

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET DemoForGeneticAlgorithm PROPERTY CXX_STANDARD 20)
endif()

include_directories(include)

# TODO: 如有需要，请添加测试并安装目标。
