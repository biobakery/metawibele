#!/usr/bin/env python

#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the SflE Scientific workFLow Environment for reproducible 
# research, authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Curtis Huttenhower, chuttenh@hsph.harvard.edu).
#
# If you use this environment, the included scripts, or any related code in your work,
# please let us know, sign up for the SflE user's group (sfle-users@googlegroups.com),
# pass along any issues or feedback, and we'll let you know as soon as a formal citation
# is available.
#######################################################################################

"""
Examples
~~~~~~~~

``data.pcl``::

	a	b
	c	d
	e	f

``Examples``::

	$ transpose.py < data.pcl
	a	c	e
	b	d	f

	$ echo "a	b	c" | transpose.py
	a
	b
	c

.. testsetup::

	from transpose import *
"""

import argparse
import csv
import sys

def transpose( aastrIn, ostm ):
	"""
	Outputs the matrix transpose of the input tab-delimited rows. 
	
	:param	aastrIn:	Split lines from which data are read.
	:type	aastrIn:	collection of string collections
	:param	ostm:		Output stream to which transposed rows are written.
	:type	ostm:		output stream

	>>> aastrIn = [list(s) for s in ("ab", "cd", "ef")]
	>>> transpose( aastrIn, sys.stdout ) #doctest: +NORMALIZE_WHITESPACE
	a	c	e
	b	d	f
	
	>>> transpose( [list("abc")], sys.stdout ) #doctest: +NORMALIZE_WHITESPACE
	a
	b
	c
	"""

	aastrLines = [a for a in aastrIn]
	csvw = csv.writer( ostm, csv.excel_tab )
	for iRow in range( len( aastrLines[0] ) ):
		csvw.writerow( [aastrLines[iCol][iRow] for iCol in range( len( aastrLines ) )] )

argp = argparse.ArgumentParser( prog = "transpose.py",
	description = """Transposes a tab-delimited text matrix.

The transposition process is robust to missing elements and rows of differing lengths.""" )
__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def main( ):
	args = argp.parse_args( )
	transpose( csv.reader( sys.stdin, csv.excel_tab ), sys.stdout )

if __name__ == "__main__":
	main( )
