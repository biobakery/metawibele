#!/usr/bin/env python

"""
MetaWIBELE: abundance_normalization module
Normalize abundance table within sample to get relative abundance

Copyright (c) 2019 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os
import os.path
import re
import argparse

try:
	from metawibele.common import utils
	#import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

description = """
Normalize abundance table
"""


# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------
def get_args():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', help='raw abundance table (tsv format)', required=True)
	parser.add_argument('-u', help='normalization scheme: copies per million [cpm], relative abundance [relab]; default=[cpm]',
	                    choices=["cpm", "relab"], default="cpm", required=True)
	parser.add_argument('-o', help='output normalized file', required=True)
	values = parser.parse_args()
	return values


# get_args


def normalize(table, cpm):
	divisor = 1e-6 if cpm else 1.0
	# compute totals by delim level
	totals_by_level = {}
	for i, row in enumerate(table.data):
		level = len(table.rowheads[i].split(utils.c_strat_delim))
		if level not in totals_by_level:
			totals_by_level[level] = [0 for k in range(len(table.colheads))]
		table.data[i] = [float(k) for k in row]
		totals_by_level[level] = [k1 + k2 for k1, k2 in zip(totals_by_level[level], table.data[i])]
	# check for sample / level combinations with zero sum
	for level in sorted(totals_by_level):
		totals = totals_by_level[level]
		for j, total in enumerate(totals):
			if total == 0:
				totals[j] = 1
				#print("WARNING: Column {} ({}) has zero sum at level {}".format(j + 1, table.colheads[j], level), file = sys.stderr)
				print("WARNING: Column {} ({}) has zero sum at level {}")
	# normalize
	for i, row in enumerate(table.data):
		level = len(table.rowheads[i].split(utils.c_strat_delim))
		# level=1 corresponds to the community level (no strata)
		#totals = totals_by_level[level] if levelwise else totals_by_level[1]
		totals = totals_by_level[1]
		table.data[i] = ["%.6g" % (row[j] / totals[j] / divisor) for j in range(len(totals))]


# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
	args = get_args()
	table = utils.Table(args.i)
	normalize(table, args.u)
	table.write(args.o)


if __name__ == "__main__":
	main()
