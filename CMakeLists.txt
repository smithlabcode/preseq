# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
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
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED yes)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++17
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

add_subdirectory(src)
