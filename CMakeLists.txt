# Copyright (C) 2026 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

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

if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

option(PRESEQ_TESTS "Enable testing (ctest or make test)" on)
if(PRESEQ_TESTS)
  include(cmake/tests.cmake)
endif()

find_package(Threads REQUIRED)

if(USE_HTSLIB)
  find_package(HTSLIB 1.20 REQUIRED)
  add_compile_definitions(HAVE_HTSLIB)
  add_subdirectory(src/bamxx)
endif()

add_library(common OBJECT src/common.cpp)
add_library(Interval6 OBJECT src/Interval6.cpp)
add_library(lc_extrap OBJECT src/lc_extrap.cpp)
add_library(gc_extrap OBJECT src/gc_extrap.cpp)
add_library(c_curve OBJECT src/c_curve.cpp)
add_library(bound_pop OBJECT src/bound_pop.cpp)
add_library(pop_size OBJECT src/pop_size.cpp)
add_library(continued_fraction OBJECT src/continued_fraction.cpp)
add_library(moment_sequence OBJECT src/moment_sequence.cpp)
add_library(load_data_for_complexity OBJECT src/load_data_for_complexity.cpp)

add_executable(preseq src/preseq.cpp)

# ADS: below, for config.h
target_include_directories(preseq PUBLIC ${CMAKE_BINARY_DIR})

target_link_libraries(preseq PRIVATE
  common
  continued_fraction
  moment_sequence
  Interval6
  load_data_for_complexity
  lc_extrap
  gc_extrap
  c_curve
  bound_pop
  pop_size
)

if(USE_HTSLIB)
  target_link_libraries(preseq PRIVATE
    bamxx
    HTSLIB::HTSLIB
  )
endif()
