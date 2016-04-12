"""
synonymRE.py

searches protein sequence for sites that can be synonymously mutated to create
restriction enzyme recognition sites.

Matt Rich, 04/2016
"""

from argparse import ArgumentParser
from fasta_iter import fasta_iter
from itertools import product

def main(seq, off, re, eight):
	
	#read RE sites file
	enzymes = {}	#dict of target:name
	for line in open(re, "r"):
		l = line.strip().split("\t")
		enzymes[l[1]] = l[0]

	#read off-target file
	other = "X"	#use X to delimit multiple sequences
	if off != None:
		for f in fasta_iter(off):
			other += f[1].upper()
			other += "X"
	
	#print header
	print "\t".join(["start", "enzyme", "wildtype", "synonym"])

	#remove enzymes with sites in off target sequences since these won't be
	#unique. This might not matter in practice (molecular biology-wise), 
	#but I'll put it in for now.
	re_lengths = set()
	for r in enzymes.keys():
		if r in other:
			del enzymes[r]
		elif len(r) > 6 and eight == False:
			del enzymes[r] 
		else:
			re_lengths.add(len(r))

	#enumerate all synonymous codons
	for l in re_lengths:
		for i in range(len(seq)-l):
			synonyms = enumerate_synonyms(seq[0:i], seq[i:i+l],
										  seq[i+l:])
	
			#check synonyms list for matches to restriction enzyme 
			for s in synonyms:
				for r in enzymes:
					if len(r) == l and r in s:
						print "\t".join([str(i), enzymes[r], seq[i:i+l], s])

#enumerate synonymous codons
def enumerate_synonyms(pre, seq, post):
	synonyms = []
	wt = translate_sequence(pre+seq+post)
	for s in product("ACTG", repeat=len(seq)):
		if translate_sequence(pre+"".join(s)+post) == wt:
			synonyms.append("".join(s))
	return synonyms

# lookup table for codon translation
def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	return lookup[codon.lower()]

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq
	
if __name__ == "__main__":
	
	parser = ArgumentParser()
	parser.add_argument('--seq', action = 'store', type = str, dest = 'sequence', 
		help = "sequence to be mutated")
	parser.add_argument('--offtarget', action = 'store', type = str, 
		dest = 'offtarget', help = "FASTA file containing off target sequenceto\
				search for non-unique cut sites. Contained target sequence will\
				be ignored.", default = None)
	parser.add_argument('--enzymes', action = 'store', type = str, 
		dest = 'enzymes', help = "tab-delimited file containing restriction \
				enzyme recognition sites (RE, site)")
	parser.add_argument('--8cut', action = 'store_true', dest = 'eightcut',
		help = 'use 8-cutters? (takes a while)', default=False)
	args = parser.parse_args()
	
	main(args.sequence.upper(), args.offtarget, args.enzymes, args.eightcut)	

