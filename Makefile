#    This file is part of the methpipe system
#
#    Copyright (C) 2010 University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ROOT = $(shell pwd)


all:
	@make -C src ROOT=$(ROOT) OPT=1

install:
	@make -C src ROOT=$(ROOT) OPT=1 install

clean:
	@make -C src ROOT=$(ROOT) clean
.PHONY: clean

distclean: clean
	@rm -rf $(ROOT)/bin
	@rm -rf $(ROOT)/lib
	@rm -rf $(ROOT)/include
.PHONY: distclean
