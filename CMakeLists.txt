cmake_minimum_required(VERSION 3.1.3)

add_custom_command ( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
                     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create-version.sh 
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_target(
  run ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

option(ENABLE_FLOW_GRAPH_TRACING "Enable flow graph tracing" OFF)
option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)
option(ENABLE_BENCHMARKS "Enable benchmarks" OFF)
option(ENABLE_FAST_OPTIMIZATONS "Enable more compiler optimizations" ON)
option(ENABLE_MPI "Enable algorithms that require MPI" OFF)
option(ENABLE_MPI_TIMING "Enable MPI debug logging (temporary)" OFF)
option(ENABLE_EPSILON_TIMING "Enable calculateEpsilonChange time logging (temporary)" OFF)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/Eigen/)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen/)

    # Silence some annoying warnings form Eigen with gcc >= 6.1
    SET(GCC_EIGEN_COMPILE_FLAGS "-Wno-ignored-attributes -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DEIGEN_UNROLLING_LIMIT=1000")
else()
    message(FATAL_ERROR "Cound not find Eigen submodule!")
endif()

find_package(Threads REQUIRED)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

set(TBB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/tbb)
include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS stdver=c++17 tbb_cpf=1)
find_package(TBB REQUIRED tbb_preview)

include(CheckCXXCompilerFlag)

# Check for -march=native support in the compiler
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

if(ENABLE_FAST_OPTIMIZATONS)
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -funroll-loops ")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_EIGEN_COMPILE_FLAGS}")

if(ENABLE_FLOW_GRAPH_TRACING)
    add_definitions(-DTBB_USE_THREADING_TOOLS)
    add_definitions(-DTBB_PREVIEW_FLOW_GRAPH_TRACE)
    add_definitions(-DTBB_PREVIEW_FLOW_GRAPH_FEATURES)
    add_definitions(-DTBB_PREVIEW_ALGORITHM_TRACE)
endif()

if (ENABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    add_definitions(-DMPI_ENABLED)

    if (ENABLE_MPI_TIMING)
        add_definitions(-DMPI_TIMING_ENABLED)
    endif()
endif()

if (ENABLE_EPSILON_TIMING)
    add_definitions(-DEPSILON_TIMING_ENABLED)
endif()

add_subdirectory(src)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(test)
endif()

if(ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
