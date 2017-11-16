#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171116"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"


import os
import sys
from Bio.Blast import NCBIXML

os.chdir("D:\\python_Code\\BlastXML")

infile = sys.argv[1]
outfile = sys.argv[2]

out_handle = open(outfile, "w")
infile_handle = open(infile, "r")

def retrieve_Mis_pos(HspM, label):
	p_index = -1

	posList = []
	while True:
		p_index = HspM.find(label, p_index + 1)
		if p_index == -1:
			break
		posList.append(p_index)
	if posList == []:
		return ["-"]

	return posList

out_handle.write("QuerySeq_ID\tquery_length\thits_alignmentID\thsp_matchlabel\tquery_start\tmismatch_Ofalignment\tquery_end\tMismatch_From3end\talignment_length\tright_mismatchNumber\tChromosomeID\tsubject_start\tsubject_end\tquery_seq\tsubject_seq\tHSP_bitScore\thsp_score\thsp_eValue\n")


blast_records = NCBIXML.parse(infile_handle)
for record in blast_records:
	## out_align.write(record.query + "### ### \n")
	for i , alignment in enumerate(record.alignments):
		try:
			## out_align.write(alignment.title + "--------------\n")
			for hsp in alignment.hsps:
				out_handle.write(record.query + "\t"+ str(record.query_length)+"\t" + alignment.title + "\t" )
				## print(dir(record))
				out_handle.write(hsp.match + "\t" )
				## out_align.write(str(alignment.length))
				
				if hsp.query_end < hsp.query_start:
					HSPQ_start = hsp.query_end
					HSPQ_end   = hsp.query_start
                      
				else:
					HSPQ_start = hsp.query_start
					HSPQ_end   = hsp.query_end
				
				if hsp.match.find(" ") == -1: ## MisPP is postions of mismatch in the hit query sequence
					misPL = "-"
					right_PL = "="
				else:
					if hsp.query_end < hsp.query_start:

						 revMatch = hsp.match[::-1]
						 misPs = retrieve_Mis_pos(revMatch, " ")                       
					else:

						 misPs = retrieve_Mis_pos(hsp.match, " ")
						 
					misPnum = [i+int(HSPQ_start) for i in misPs]
					misPL = ",".join([str(i+int(HSPQ_start)) for i in misPs])
					right_PL = ",".join([str(record.query_length +1 - i) for i in misPnum])
				##print(hsp.strand[1])
				
				if alignment.title.find("chromosome") != -1:
					z = alignment.title.find("chromosome")
					chrName = "chr" + alignment.title[z+11:].split()[0].strip(",")
				else:
					chrName = alignment.title
				
				if chrName.find("chr") != -1:
                    
					cIndex = chrName.find("chr")
					chrName = chrName[cIndex:]
                    
					out_handle.write(str(HSPQ_start) +"\t" + misPL + "\t" + str(HSPQ_end)+"\t" + right_PL + "\t" + str(hsp.align_length) + "\t" +str(record.query_length - HSPQ_end ) + "\t" + chrName +"\t" + str(hsp.sbjct_start) + "\t" + str(hsp.sbjct_end) + "\t"  + hsp.query + "\t" + hsp.sbjct + "\t")
				
					out_handle.write("{0}\t{1}\t{2}".format(hsp.bits, hsp.score, hsp.expect))
				out_handle.write("\n")
		except:
			out_handle.write("No Hits\tNone\n")

out_handle.close()  
infile_handle.close()
