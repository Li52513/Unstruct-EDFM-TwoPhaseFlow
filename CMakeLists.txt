cmake_minimum_required(VERSION 3.16)
project(BoundaryDOFFramework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BNDF_ENABLE_WARNINGS "Enable extra warnings" ON)

if(BNDF_ENABLE_WARNINGS)
  if (MSVC)
    add_compile_options(/W4)
  else()
    add_compile_options(-Wall -Wextra -Wpedantic)
  endif()
endif()

add_library(bndf
  src/core/LinearSystem.cpp
  src/mesh/Mesh.cpp
  src/mesh/BoundaryDOF.cpp
  src/bc/TemperatureBC.cpp
  src/assembly/AssembleTemperature.cpp
  src/io/CSV.cpp
)

target_include_directories(bndf PUBLIC include)

add_executable(bndf_demo src/main.cpp)
target_link_libraries(bndf_demo PRIVATE bndf)

