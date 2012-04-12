This is the README file for the first release of Rational Function 
Complexity. Rational Function Complexity is a tool for estimating the 
yield of distinct reads from high-throughput sequencing technology.

CONTACT INFORMATION:
========================================================================
Timothy Daley
tdaley@usc.edu
http://smithlab.cmb.usc.edu/plone

SYSTEM REQUIREMENTS:
========================================================================
The Rational Function Complexity software will only run on UNIX-like 
operating systems, and was developed on Linux systems. The Rational
Function Complexity software requires a fairly recent C++ compiler 
(i.e. it must include tr1 headers).  Rational Function Complexity has 
been compiled and tested on Linux and OS X operating systems using GCC 
v4.1 or greater. Also, Rational Function Complexity will only run on 
64-bit machines.

INSTALLATION:
========================================================================
This should be easy: unpack the archive and change into the archive
directory. Then type 'make all'. The programs will be in the archive
directory. These can be moved around, and also do not depend on any
dynamic libraries, so they should simply work when executed. If the 
desired input is in .bam format, bamtools is required and is located
at '/bamtools_loc/.  Type 'make all BAMTOOLS_ROOT=/bamtools_loc/'
to make the programs.

USAGE EXAMPLES:
========================================================================
Each program included in this software package will print a list of
options if executed without any command line arguments. Many of the
programs use similar options (for example, output files are specified
with '-o'). For the most basic usage of lc_extrap to bound the number
of distinct reads and the expected yield, use the command:

  lc_extrap -o yield_estimates.txt -L size_estimates.txt input.bed

If the input file is in .bam format, use the command:

  lc_extrap -B -o yield_estimates.txt -L size_estimates.txt input.bed

The yield estimates will appear in yield_estimates.txt, and will 
be a column of future experiment sizes in total reads, a column of the 
corresponding expected distinct reads, followed by two columns giving
the corresponding confidence intervals.  The bounds on the number
of distinct reads will appear as size_estimates.txt, and will give
a column giving the upper and lower bounds, and corresponding confidence
intervals.

HISTORY
========================================================================
Rational Function Complexity was originally developed by Timothy Daley
and Andrew Smith at the University of Southern California.


LICENSE
========================================================================
The Rational Function Complexity software for estimating complexity
Copyright (C) 2012 Timothy Daley and Andrew D Smith and 
the University of Southern California

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.