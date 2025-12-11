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

if(USE_HTSLIB)
  find_package(HTSLIB REQUIRED)
endif()
find_package(Threads REQUIRED)

if(USE_HTSLIB)
  add_compile_definitions(HAVE_HTSLIB)
  add_subdirectory(bamxx)
endif()

add_executable(preseq
  preseq.cpp
  bound_pop.cpp
  c_curve.cpp
  common.cpp
  continued_fraction.cpp
  gc_extrap.cpp
  lc_extrap.cpp
  load_data_for_complexity.cpp
  moment_sequence.cpp
  pop_size.cpp
  Interval6.cpp
)
if(USE_HTSLIB)
  target_sources(preseq PRIVATE bam_record_utils.cpp)
endif()

# ADS: below, for config.h
target_include_directories(preseq PUBLIC ${CMAKE_BINARY_DIR})
if(USE_HTSLIB)
  target_link_libraries(preseq PUBLIC
    bamxx
    HTSLIB::HTSLIB
  )
endif()
