# SPDX-License-Identifier: GPL-3.0; Copyright 2026 Andrew D Smith

# Add test scripts and data

include(CTest)

## Ensure test data is available in the build dir
file(CREATE_LINK
  ${PROJECT_SOURCE_DIR}/data
  ${PROJECT_BINARY_DIR}/data
  SYMBOLIC
)

## Ensure test scripts are available in the build dir
file(CREATE_LINK
  ${PROJECT_SOURCE_DIR}/test
  ${PROJECT_BINARY_DIR}/test
  SYMBOLIC
)

find_program(BASH_PROGRAM bash)

## Add each test
add_test(NAME "lc_extrap histogram input" COMMAND bash test/lc_extrap_hist.sh)
