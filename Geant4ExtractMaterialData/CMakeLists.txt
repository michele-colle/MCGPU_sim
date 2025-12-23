cmake_minimum_required(VERSION 3.16...3.27)
project(g4getdata)
find_package(Geant4 REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIR})
add_executable(Geant4ExtractMaterialData main.cpp)
target_link_libraries(Geant4ExtractMaterialData ${Geant4_LIBRARIES})