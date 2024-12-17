cmake_minimum_required(VERSION 3.8)
project(nokolat2024_sim)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake_auto REQUIRED)
ament_auto_find_build_dependencies()



ament_auto_add_library(nokolat2024_sim src/airplane_eom.cpp include/nokolat2024_sim/airplane_eom.hpp)

ament_auto_add_executable(airplane_eom src/airplane_eom.cpp)



if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()
endif()

ament_auto_package()
