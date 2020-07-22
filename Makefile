# Copyright (C) 2011-2020 University of Southern California and
#                         Andrew D. Smith and Timothy Daley
#
# Authors: Timothy Daley and Andrew D. Smith
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

ifndef install_dir
install_dir := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
endif

all:
	@make -C src

install:
	@make -C src install_dir=$(install_dir) install

clean:
	@make -C src clean

distclean: clean
	@rm -rf $(install_dir)/bin

.PHONY: all distclean clean install
