cmake_minimum_required(VERSION 3.16...3.27)
project(geant4_to_mcgpu)
find_package(Geant4 REQUIRED)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# DEVE essere presente il progetto CBCTSim nella cartella ../../CBCTSim
#
include_directories(${PROJECT_SOURCE_DIR}/../../CBCTSim/ICRPinclude
                    ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/../../CBCTSim/ICRPsrc/*.cc )
file(GLOB headers ${PROJECT_SOURCE_DIR}/../../CBCTSim/ICRPinclude/*.hh )


include(ExternalProject)
ExternalProject_Add(ICRPdata
  SOURCE_DIR ${PROJECT_BINARY_DIR}/ICRPdata
  URL https://cern.ch/geant4-data/datasets/examples/advanced/ICRP110Phantoms/ICRPdata.tar.gz
  URL_MD5 2aba5409d4e7b7a2f0328e7d3402f187
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

add_executable(main main.cpp ${sources} ${headers})
target_link_libraries(main ${Geant4_LIBRARIES} )
#add_dependencies(ICRP110phantoms ICRPdata)
