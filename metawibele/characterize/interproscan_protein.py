#!/usr/bin/env python

"""
MetaWIBELE: interproscan_protein module
Extact functioanl annotation for each peptide from InterProScan results

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
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extact functional annotation for each peptide from InterProScan results
"""

def get_args():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
						default="interproscan.txt")
	parser.add_argument('-p', "--path",
	                    help='input the path of InterProScan folder',
	                    required=True)
	values = parser.parse_args()
	return values


#==============================================================
# collect InterProScan info
#==============================================================
def extract_interproscan_info (extension, interproscan_path):
	filelist = utilities.find_files(interproscan_path, extension, None)
	for myfile in filelist:
		#myfile = interproscan_path + "/" + samplelist + "/" + samplelist + ".interproscan.txt"
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			print ("OK!\t" + myfile)
			myout1 = re.sub(".interproscan.txt", ".signalp.signaling.tsv", myfile)
			myout2 = re.sub(".interproscan.txt", ".tmhmm.transmembrane.tsv", myfile)
			myout3 = re.sub(".interproscan.txt", ".phobius.signaling.tsv", myfile)
			myout4 = re.sub(".interproscan.txt", ".phobius.transmembrane.tsv", myfile)
			myout5 = re.sub(".interproscan.txt", ".interpro.PfamDomain.tsv", myfile)
			myout6 = re.sub(".interproscan.txt", ".interpro.SUPERFAMILY.tsv", myfile)
			myout7 = re.sub(".interproscan.txt", ".interpro.PROSITEPROFILES.tsv", myfile)
			myout8 = re.sub(".interproscan.txt", ".interpro.Gene3D.tsv", myfile)
			myout9 = re.sub(".interproscan.txt", ".interpro.PANTHER.tsv", myfile)
			myout10 = re.sub(".interproscan.txt", ".interpro.TIGRFAM.tsv", myfile)
			myout11 = re.sub(".interproscan.txt", ".interpro.SFLD.tsv", myfile)
			myout12 = re.sub(".interproscan.txt", ".interpro.ProDom.tsv", myfile)
			myout13 = re.sub(".interproscan.txt", ".interpro.Hamap.tsv", myfile)
			myout14 = re.sub(".interproscan.txt", ".interpro.SMART.tsv", myfile)
			myout15 = re.sub(".interproscan.txt", ".interpro.CDD.tsv", myfile)
			myout16 = re.sub(".interproscan.txt", ".interpro.PROSITEPATTERNS.tsv", myfile)
			myout17 = re.sub(".interproscan.txt", ".interpro.PRINTS.tsv", myfile)
			myout18 = re.sub(".interproscan.txt", ".interpro.PIRSF.tsv", myfile)
			myout19 = re.sub(".interproscan.txt", ".interpro.MobiDBLite.tsv", myfile)
			myout20 = re.sub(".interproscan.txt", ".interpro.Coils.tsv", myfile)
			open_out1 = open(myout1, "w")
			open_out2 = open(myout2, "w")
			open_out3 = open(myout3, "w")
			open_out4 = open(myout4, "w")
			open_out5 = open(myout5, "w")
			open_out6 = open(myout6, "w")
			open_out7 = open(myout7, "w")
			open_out8 = open(myout8, "w")
			open_out9 = open(myout9, "w")
			open_out10 = open(myout10, "w")
			open_out11 = open(myout11, "w")
			open_out12 = open(myout12, "w")
			open_out13 = open(myout13, "w")
			open_out14 = open(myout14, "w")
			open_out15 = open(myout15, "w")
			open_out16 = open(myout16, "w")
			open_out17 = open(myout17, "w")
			open_out18 = open(myout18, "w")
			open_out19 = open(myout19, "w")
			open_out20 = open(myout20, "w")

			open_file = open(myfile, "r")
			open_out1.write(utilities.PROTEIN_ID + "\tSP\tPrediction\tStart\tEnd\n")
			open_out2.write(utilities.PROTEIN_ID + "\tTM\tPrediction\tStart\tEnd\n")
			open_out3.write(utilities.PROTEIN_ID + "\tSP\tPrediction\tStart\tEnd\n")
			open_out4.write(utilities.PROTEIN_ID + "\tTM\tPrediction\tStart\tEnd\n")
			open_out5.write(utilities.PROTEIN_ID + "\tPfam\tDescription\tInterPro\tEvalue\n")
			open_out6.write(utilities.PROTEIN_ID + "\tSUPERFAMILY\tDescription\tInterPro\tEvalue\n")
			open_out7.write(utilities.PROTEIN_ID + "\tProSiteProfiles\tDescription\tInterPro\tEvalue\n")
			open_out8.write(utilities.PROTEIN_ID + "\tGene3D\tDescription\tInterPro\tEvalue\n")
			open_out9.write(utilities.PROTEIN_ID + "\tPANTHER\tDescription\tInterPro\tEvalue\n")
			open_out10.write(utilities.PROTEIN_ID + "\tTIGRFAM\tDescription\tInterPro\tEvalue\n")
			open_out11.write(utilities.PROTEIN_ID + "\tSFLD\tDescription\tInterPro\tEvalue\n")
			open_out12.write(utilities.PROTEIN_ID + "\tProDom\tDescription\tInterPro\tEvalue\n")
			open_out13.write(utilities.PROTEIN_ID + "\tHamap\tDescription\tInterPro\tEvalue\n")
			open_out14.write(utilities.PROTEIN_ID + "\tSMART\tDescription\tInterPro\tEvalue\n")
			open_out15.write(utilities.PROTEIN_ID + "\tCDD\tDescription\tInterPro\tEvalue\n")
			open_out16.write(utilities.PROTEIN_ID + "\tProSitePatterns\tDescription\tInterPro\tEvalue\n")
			open_out17.write(utilities.PROTEIN_ID + "\tPRINTS\tDescription\tInterPro\tEvalue\n")
			open_out18.write(utilities.PROTEIN_ID + "\tPIRSF\tDescription\tInterPro\tEvalue\n")
			open_out19.write(utilities.PROTEIN_ID + "\tMobiDBLite\tDescription\tInterPro\tEvalue\n")
			open_out20.write(utilities.PROTEIN_ID + "\tCoils\tDescription\tInterPro\tEvalue\n")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^#", line):
					continue
				info = line.split("\t")
				myid = info[0]
				mytype = info[3]
				myacc = info[4]
				mydec = info[5]
				start = info[6]
				end = info[7]
				myscore = info[8]
				mystatus = info[9]
				if len(info) < 12:
					interproacc = "NA"
					interprodec = "NA"
				else:
					interproacc = info[11]
					interprodec = info[12]
				if start == "":
					start = "NA"
				if end == "":
					end = "NA"
				if interproacc == "":
					interproacc = "NA"
				if interprodec == "":
					interprodec = "NA"
				if mydec == "":
					mydec = "NA"
				if interprodec != "NA":
					mydec = interprodec
				if mystatus != "T":	# not reliable prediction
					continue
				# SignalP
				if mytype == "SignalP_GRAM_NEGATIVE" or mytype == "SignalP_GRAM_POSITIVE":
					open_out1.write(myid + "\t" + mytype + "\t" + myacc + "\t" + start + "\t" + end + "\n")
				# TMHMM
				if mytype == "TMHMM":
					open_out2.write(myid + "\t" + mytype + "\t" + myacc + "\t" + start + "\t" + end + "\n")
				# Phobius
				if mytype == "Phobius":
					if re.search("SIGNAL_PEPTIDE", myacc): # signal peptide
						open_out3.write(myid + "\t" + myacc + "\t" + mydec + "\t" + start + "\t" + end + "\n")
					if re.search("TRANSMEMBRANE", myacc): # transmembrane
						open_out4.write(myid + "\t" + myacc + "\t" + mydec + "\t" + start + "\t" + end + "\n")
				# Pfam
				if mytype == "Pfam":
					open_out5.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# SUPERFAMILY
				if mytype == "SUPERFAMILY":
					open_out6.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# ProSiteProfiles
				if mytype == "ProSiteProfiles":
					open_out7.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# Gene3D
				if mytype == "Gene3D":
					open_out8.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# PANTHER
				if mytype == "PANTHER":
					open_out9.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# TIGRFAM
				if mytype == "TIGRFAM":
					open_out10.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# SFLD
				if mytype == "SFLD":
					open_out11.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# ProDom
				if mytype == "ProDom":
					open_out12.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# Hamap
				if mytype == "Hamap":
					open_out13.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# SMART
				if mytype == "SMART":
					open_out14.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# CDD
				if mytype == "CDD":
					open_out15.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# ProSitePatterns
				if mytype == "ProSitePatterns":
					open_out16.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# PRINTS
				if mytype == "PRINTS":
					open_out17.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# PIRSF
				if mytype == "PIRSF":
					open_out18.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# MobiDBLite
				if mytype == "MobiDBLite":
					open_out19.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
				# Coils
				if mytype == "Coils":
					open_out20.write(myid + "\t" + myacc + "\t" + mydec + "\t" + interproacc + "\t" + myscore + "\n") 
			# foreach line	
			open_out1.close()
			open_out2.close()
			open_out3.close()
			open_out4.close()
			open_out5.close()
			open_out6.close()
			open_out7.close()
			open_out8.close()
			open_out9.close()
			open_out10.close()
			open_out11.close()
			open_out12.close()
			open_out13.close()
			open_out14.close()
			open_out15.close()
			open_out16.close()
			open_out17.close()
			open_out18.close()
			open_out19.close()
			open_out20.close()
			open_file.close()
		# if file exist
	# foreach samplelist
# function extract_interproscan_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start interproscan_protein.py -p " + values.path + " ####\n")
	
	### collect Pfam info ###
	sys.stderr.write("Get InterProScan info ......starting\n")
	extract_interproscan_info (values.extension, values.path)
	sys.stderr.write("Get InterProScan info ......done\n")

	sys.stderr.write("### Finish interproscan_protein.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()
