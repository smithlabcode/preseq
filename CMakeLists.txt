# SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.30)
project(
  preseq
  VERSION 3.2.0
  DESCRIPTION
  "predict properties of genomic sequencing libraries"
  HOMEPAGE_URL https://github.com/smithlabcode/preseq
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED yes)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++23
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckCXXCompilerFlag)

include(GNUInstallDirs)

configure_file(data/config.h.in config.h)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
  -Werror
  -Wfatal-errors
)

file(READ data/LICENSE LICENSE_TEXT)
configure_file(data/LICENSES.inc.in license.h @ONLY)

if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

option(PRESEQ_TESTS "Enable testing (ctest or make test)" on)
if(PRESEQ_TESTS)
  include(cmake/tests.cmake)
endif()

add_subdirectory(src/libpreseq)
find_package(Threads REQUIRED)

find_package(ZLIB REQUIRED)
find_package(HTSLIB 1.20 REQUIRED)
add_subdirectory(src/bamxx)

# header only libs in subdirs
add_library(json INTERFACE)
target_include_directories(json INTERFACE src)
add_library(CLI11 INTERFACE)
target_include_directories(CLI11 INTERFACE src)

set(app_names
  lc_extrap
  gc_extrap
  c_curve
  bound_pop
  pop_size
)
set(obj_names
  ${app_names}
  load_data_for_complexity
  common
  Interval6
)
foreach(obj_name ${obj_names})
  add_library(${obj_name} OBJECT src/${obj_name}.cpp)
  target_include_directories(${obj_name} PUBLIC src)
  target_link_libraries(${obj_name} PUBLIC libpreseq::libpreseq_static)
endforeach()

add_executable(preseq src/preseq.cpp)
target_include_directories(preseq PUBLIC ${CMAKE_BINARY_DIR}) # for config.h
target_link_libraries(preseq PUBLIC
  common
  Interval6
  load_data_for_complexity
  ${app_names}
  libpreseq::libpreseq_shared
  json
  CLI11
  bamxx
  HTSLIB::HTSLIB
  ZLIB::ZLIB
)
