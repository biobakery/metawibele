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
	parser.add_argument("-m", 
						help="Normalize all levels by [community] total or [levelwise] total or [taxonwise] within taxon; default=[community]",
						choices=["community", "levelwise", "taxonwise"],
						default="community")
	parser.add_argument("-s",
						help="Include the special features UNMAPPED, UNINTEGRATED, and UNGROUPED; default=[y]",
						choices=["y", "n"],
						default="y")
	parser.add_argument('-o', help='output normalized file', required=True)
	values = parser.parse_args()
	return values


# get_args

# constants
c_special = [
	utils.c_unmapped, 
	utils.c_unintegrated, 
	utils.c_ungrouped,
	utils.c_msp_unknown,
]

def normalize(table, method, levelwise, special=True):
	if levelwise == 'levelwise':
		levelwise = True
		taxonwise = False
	if levelwise == 'community':
		levelwise = False
		taxonwise = False
	if levelwise == "taxonwise":
		taxonwise = True
		levelwise = False
	if method == "cpm":
		divisor = 1e-6
	else:
		divisor = 1.0
	#divisor = 1e-6 if cpm else 1.0
	
	# remove special features?
	if not special:
		test = [rowhead.split( util.c_strat_delim )[0] not in c_special for rowhead in table.rowheads]
		for flag, rowhead in zip( test, table.rowheads ):
			if not flag:
				print( "Excluding special feature:" + rowhead )
		table.rowheads = [rowhead for i, rowhead in enumerate( table.rowheads ) if test[i]]
		table.data = [row for i, row in enumerate( table.data ) if test[i]]
	
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
				#totals[j] = 1
				totals_by_level[level][j] = 1
				print("WARNING: Column {} ({}) has zero sum at level {}".format(j + 1, table.colheads[j], level))
	
	# compute totals by delim taxon
	totals_by_taxon = {}
	for i, row in enumerate(table.data):
		if not re.search(utils.c_strat_delim, table.rowheads[i]):
			continue
		taxon = table.rowheads[i].split(utils.c_strat_delim)[-1]
		if taxon not in totals_by_taxon:
			totals_by_taxon[taxon] = [0 for k in range(len(table.colheads))]
		table.data[i] = [float(k) for k in row]
		totals_by_taxon[taxon] = [k1 + k2 for k1, k2 in zip(totals_by_taxon[taxon], table.data[i])]
	# check for sample / taxon combinations with zero sum
	for taxon in totals_by_taxon:
		for j in range(len(totals_by_taxon[taxon])):
			if totals_by_taxon[taxon][j] == 0:
				totals_by_taxon[taxon][j] = 1

	# normalize
	for i, row in enumerate(table.data):
		taxon = table.rowheads[i].split(utils.c_strat_delim)[-1]
		if taxonwise:
			#print("Within-taxon level normalization")
			totals = totals_by_taxon 
			table.data[i] = ["%.6g" % (row[j] / totals[taxon][j] / divisor) for j in range(len(totals[taxon]))]
			continue
		
		level = len(table.rowheads[i].split(utils.c_strat_delim))
		totals = totals_by_level[level] if levelwise else totals_by_level[1]
		table.data[i] = ["%.6g" % (row[j] / totals[j] / divisor) for j in range(len(totals))]


# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
	args = get_args()
	table = utils.Table(args.i)
	normalize(table, args.u, args.m, args.s=="y")
	table.write(args.o)


if __name__ == "__main__":
	main()
