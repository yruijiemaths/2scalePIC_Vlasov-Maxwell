cmake_minimum_required(VERSION 3.18)

project(PICfortran LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(PIC_METHOD "TSI-Sim" CACHE STRING "PIC time integrator")
set_property(CACHE PIC_METHOD PROPERTY STRINGS TSI-CPex TSI-CPim TSI-Sex TSI-Sim TSFV-CPimex TSFV-CP2imex)

set(PIC_GRID "auto" CACHE STRING "Spatial grid: auto, standard, or half-grid")
set_property(CACHE PIC_GRID PROPERTY STRINGS auto standard half-grid)

set(_method_names TSI-CPex TSI-CPim TSI-Sex TSI-Sim TSFV-CPimex TSFV-CP2imex)
list(FIND _method_names "${PIC_METHOD}" _method_index)
if(_method_index EQUAL -1)
  message(FATAL_ERROR "Unknown PIC_METHOD='${PIC_METHOD}'. Choose one of: ${_method_names}")
endif()

if(PIC_GRID STREQUAL "auto")
  if(PIC_METHOD MATCHES "^TSFV-CP(2)?imex$")
    set(_selected_grid "half-grid")
  else()
    set(_selected_grid "standard")
  endif()
else()
  set(_selected_grid "${PIC_GRID}")
endif()

if(PIC_METHOD MATCHES "^TSFV-CP(2)?imex$" AND NOT _selected_grid STREQUAL "half-grid")
  message(FATAL_ERROR "${PIC_METHOD} requires PIC_GRID=half-grid")
endif()
if(NOT PIC_METHOD MATCHES "^TSFV-CP(2)?imex$" AND NOT _selected_grid STREQUAL "standard")
  message(FATAL_ERROR "${PIC_METHOD} currently requires PIC_GRID=standard")
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)

set(FFTW_ROOT "" CACHE PATH "FFTW installation prefix")
find_path(FFTW_INCLUDE_DIR NAMES fftw3.f
  HINTS "${FFTW_ROOT}" ENV FFTW_ROOT ENV FFTW_HOME ENV FFTW_DIR
  PATH_SUFFIXES include)
find_library(FFTW_LIBRARY NAMES fftw3
  HINTS "${FFTW_ROOT}" ENV FFTW_ROOT ENV FFTW_HOME ENV FFTW_DIR
  PATH_SUFFIXES lib lib64)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW REQUIRED_VARS FFTW_INCLUDE_DIR FFTW_LIBRARY)

set(_method_source "${CMAKE_CURRENT_SOURCE_DIR}/src/methods/${PIC_METHOD}/main.F90")
set(_grid_dir "${CMAKE_CURRENT_SOURCE_DIR}/src/grids/${_selected_grid}")

add_executable(pic-vm
  "${CMAKE_CURRENT_SOURCE_DIR}/src/common/pic_config.F90"
  "${_grid_dir}/mesh.F90"
  "${_grid_dir}/interpolations.F90"
  "${_grid_dir}/particles.F90"
  "${_grid_dir}/poisson.F90"
  "${_grid_dir}/grid_output_adapter.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/common/pic_output.F90"
  "${_method_source}")

set_target_properties(pic-vm PROPERTIES
  OUTPUT_NAME "pic-vm.exe"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
target_include_directories(pic-vm PRIVATE "${CMAKE_BINARY_DIR}/modules" "${FFTW_INCLUDE_DIR}")
target_link_libraries(pic-vm PRIVATE MPI::MPI_Fortran "${FFTW_LIBRARY}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(pic-vm PRIVATE -std=legacy -ffree-line-length-none
    $<$<CONFIG:Debug>:-O0;-g;-fbacktrace;-fcheck=all>)
endif()

message(STATUS "PIC method: ${PIC_METHOD}")
message(STATUS "PIC grid:   ${_selected_grid}")
