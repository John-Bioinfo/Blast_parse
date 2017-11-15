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
from Bio.Blast import NCBIXML

os.chdir("D:\\python_Code\\BlastXML")
out_align = open("D:\\python_Code\\" + "out_a.xls", "w")
infile = open("test_Alignment.xml", "r")

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

out_align.write("QuerySeq_ID\tquery_length\thits_alignmentID\thsp_matchlabel\tquery_start\tmismatch_Ofalignment\tquery_end\tMismatch_From3end\talignment_length\tright_mismatchNumber\tChromosomeID\tsubject_start\tsubject_end\tquery_seq\tsubject_seq\tHSP_bitScore\thsp_score\thsp_eValue\n")


blast_records = NCBIXML.parse(infile)
for record in blast_records:
	## out_align.write(record.query + "### ### \n")
	for i , alignment in enumerate(record.alignments):
		try:
			## out_align.write(alignment.title + "--------------\n")
			for hsp in alignment.hsps:
				out_align.write(record.query + "\t"+ str(record.query_length)+"\t" + alignment.title + "\t" )
				## print(dir(record))
				out_align.write(hsp.match + "\t" )
				## out_align.write(str(alignment.length))
				misPs = retrieve_Mis_pos(hsp.match, " ")
				
				if misPs == ["-"]:          ## MisPP is postions of mismatch in the hit query sequence
					misPL = "-"
					right_PL = "="
				else:
					misPnum = [i+int(hsp.query_start) for i in misPs]
					misPL = ",".join([str(i+int(hsp.query_start)) for i in misPs])
					right_PL = ",".join([str(record.query_length +1 - i) for i in misPnum])
				##print(hsp.strand[1])
				
				if alignment.title.find("chromosome") != -1:
					z = alignment.title.find("chromosome")
                    
					chrName = "chr" + alignment.title[z+11:].split()[0].strip(",")
				else:
					chrName = alignment.title
				
				out_align.write(str(hsp.query_start) +"\t" + misPL + "\t" + str(hsp.query_end)+"\t" + right_PL + "\t" + str(hsp.align_length) + "\t" +str(record.query_length - hsp.query_end) + "\t" + chrName +"\t" + str(hsp.sbjct_start) + "\t" + str(hsp.sbjct_end) + "\t"  + hsp.query + "\t" + hsp.sbjct + "\t")
				
				out_align.write("{0}\t{1}\t{2}".format(hsp.bits, hsp.score, hsp.expect))
				out_align.write("\n")
		except:
			out_align.write("No Hits\tNone\n")

out_align.close()  
infile.close()
