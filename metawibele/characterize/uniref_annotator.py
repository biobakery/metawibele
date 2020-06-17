#!/usr/bin/env python

"""
MetaWIBELE: uniref_annotator module
Annotate proteins to uniref DB and get UniRef90 and UniRef50 IDs

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


from __future__ import print_function
import os
import sys
import argparse
import csv
import re

from metawibele.characterize.utils import say, die, check_path, which, Hit, translate_fasta

"""
NOTES:
When we transition to DIAMOND 0.9+, change the search
to manually enforce query and subject coverage and then
only consider the top hit.
08/12/2018: modify function "reannotate" for UniRef90 and UniRef50 mapping info, including UniRef90_xxx-[1-9]
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_min_coverage    = 0.80
c_diamond_filters = ""
c_output_format   = "6 qseqid sseqid pident qlen qstart qend slen sstart send evalue"
#g_force_search    = False
g_force_search    = True

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

description = """
A script to annotate a fasta file of coding sequence against
HUMAnN2-formatted UniRef90/UniRef50 databases. Approximates the UniRef
clustering conventions in that query and target must exceed 80% mutual coverage
with the corresponding percent identity (e.g. 90 for UniRef90).
"""

def get_args( ):
    parser = argparse.ArgumentParser( description=description )
    parser.add_argument( "fasta", 
                         help="Sequences to annotate",
                         )
    parser.add_argument( "--seqtype",
                         choices=["cds", "nuc", "prot"],
                         metavar="<cds/nuc/prot>",
                         default="cds",
                         help="Sequence type [default: cds]",
                         )
    parser.add_argument( "--diamond",
                         metavar="<path>",
                         default="diamond",
                         help="Path to diamond binary [default: in PATH]",
                         )
    parser.add_argument( "--uniref90db",
                         metavar="<path>",
                         required=True,
                         help="Path to HUMAnN2-formatted UniRef90 database",
                         )
    parser.add_argument( "--temp",
                         metavar="<path>",
                         default=".",
                         help="Path for temp files [default: .]",
                         )
    parser.add_argument( "--out",
                         metavar="<path>",
                         default=None,
                         help="Path for output file [default: <fasta>.annotated]",
                         )
    parser.add_argument( "--force-search",
                         action="store_true",
                         help="Rerun searches, even if expected outputs exist",
                         )
    parser.add_argument( "--diamond-options",
                         metavar="<string>",
                         help="Additional options to pass to diamond, e.g. --threads",
                         )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utils 
# ---------------------------------------------------------------

def get_mode( path ):
    mode = None
    for test in "90 50".split( ):
        test = "uniref" + test
        if test in path.lower( ):
            mode = test
    if mode is None:
        die( "Could not infer mode from path: {}".format( path ) )
    return mode

def uniref_search( diamond=None, database=None, query=None, seqtype=None, temp=None, diamond_options=None ):
    if which( diamond ) is None:
        die( "<diamond> is not executable as: {}".format( diamond ) )
    for path in [database, query, temp]:
        check_path( path )
    binary = {"nuc":"blastx", "prot":"blastp"}[seqtype]
    mode = get_mode( database )
    results = os.path.split( query )[1]
    results = os.path.join( temp, results )
    results = ".".join( [results, mode, "hits"] )
    command = [
        diamond,
        binary,
        "--db", database,
        "--query", query,
        "--outfmt", c_output_format,
        "--tmpdir", temp,
        "--out", results,
        #"--id", get_mode( results ).replace( "uniref", "" ),
        c_diamond_filters,
        ]
    command = " ".join( [str( k ) for k in command] )
    command += (" " + diamond_options) if diamond_options is not None else ""
    if not os.path.exists( results ) or g_force_search:
        say( "Executing:\n ", command )
        os.system( command )
    else:
        say( "Using existing results file:\n ", results )
    return results

def parse_results( results ):
    say( "Parsing results file:\n ", results )
    check_path( results )
    mapping = {}
    mode = get_mode( results )
    min_pident = float( mode.replace( "uniref", "" ) )
    with open( results ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            h = Hit( row, config=c_output_format )
            if h.qseqid not in mapping:
                if float(h.pident) >= float(min_pident) and float(h.mcov) >= float(c_min_coverage):
                    uniref = h.sseqid.split( "|" )[0]
                    mapping[h.qseqid] = uniref
    return mapping

def trans_mapping( uniref90map, p_trans_map ):
    say( "Loading transitive mapping file:\n ", p_trans_map )
    check_path( p_trans_map )
    overrides = {}
    uniref90map_r = {}
    for header, uniref90 in uniref90map.items( ):
        # modify by yancong
        uniref90 = re.sub("-[0-9]+$", "", uniref90)
        uniref90map_r.setdefault( uniref90, set( ) ).add( header )
    with open( p_trans_map ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            uniref90, uniref50 = row
            headers = uniref90map_r.get( uniref90, set( ) )
            for h in headers:
                overrides[h] = uniref50
    return overrides

def reannotate( query=None, out=None, uniref90map=None, uniref50map=None, overrides=None ):
    say( "Writing new output file:\n ", out )
    oh = open( out, "w" )
    ntot, nmap90, ninf50, nmap50 = [0 for i in range( 4 )] 
    with open( query ) as fh:
        for line in fh:
            line = line.strip( )
            if line == "":
                continue
            elif line[0] != ">":
                print( line, file=oh )
            else:
                # diamond breaks the header on whitespace
                header = line[1:].split( )[0]
                ntot += 1
                uniref90code = "UniRef90_unknown"
                if header in uniref90map:
                    uniref90code = uniref90map[header]
                    nmap90 += 1
                uniref50code = "UniRef50_unknown"
                if header in overrides:
                    uniref50code = overrides[header]
                    ninf50 += 1
                elif header in uniref50map:
                    uniref50code = uniref50map[header]
                    nmap50 += 1
                print( "|".join( [line, uniref90code, uniref50code] ), file=oh )
    oh.close( )
    # report
    say( "Summary of annotations:" )
    say( "  Genes in input FASTA: {:,}".format( ntot ) )
    say( "  UniRef90 codes assigned: {:,} ({:.1f}%)".format(
            nmap90, 100 * nmap90 / float( ntot ) ) )
    say( "  UniRef50 codes assigned: {:,} ({:.1f}%)".format( 
            nmap50 + ninf50, 100 * (nmap50 + ninf50) / float( ntot ) ) )
    say( "  UniRef50 codes inferred from UniRef90 codes: {:,} ({:.1f}%)".format( 
            ninf50, 100 * ninf50 / float( ntot ) ) )
    # done
    return None

# ---------------------------------------------------------------
# main 
# ---------------------------------------------------------------
    
def main( ):
    args = get_args( )
    # global config
    global g_force_search
    if args.force_search:
        g_force_search = True
    # set defaults
    if args.out is None:
        args.out = args.fasta + ".annotated"
    # translate fasta?
    query = args.fasta
    if args.seqtype == "cds":
        query = os.path.split( query )[1]
        query = os.path.join( args.temp, query )
        query = query + ".translated"
        say( "Translating input fasta to:\n ", query )
        translate_fasta( args.fasta, query )
        args.seqtype = "prot"
    # perform uniref90 search
    uniref90hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref90db,
        query=query,
        seqtype=args.seqtype,
        temp=args.temp,
        diamond_options=args.diamond_options, )
    
    '''
	uniref90hits = query + ".uniref90.hits"
    uniref90map = parse_results( uniref90hits )
    # perform uniref50 search
    uniref50hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref50db,
        query=query,
        seqtype=args.seqtype,
        temp=args.temp,
        diamond_options=args.diamond_options, )
    uniref50hits = query + ".uniref50.hits"
    uniref50map = parse_results( uniref50hits )
    # override mappings?
    overrides = {}
    if args.transitive_map is not None:
        overrides = trans_mapping( uniref90map, args.transitive_map )
    # reannoate the fasta
    reannotate( 
        query=args.fasta, 
        out=args.out, 
        uniref90map=uniref90map, 
        uniref50map=uniref50map, 
        overrides=overrides, )
    # done
    '''
    say( "Finished successfully." )

if __name__ == "__main__":
    main( )
