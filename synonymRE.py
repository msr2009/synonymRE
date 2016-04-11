"""
synonymRE.py

searches protein sequence for sites that can be synonymously mutated to create
restriction enzyme recognition sites.

Matt Rich, 04/2016
"""

def main(seq, off, re):
	
	#read RE sites file

	#read 
	

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
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--seq', action = 'store', type = str, dest = 'sequence', 
		help = "sequence to be mutated")
	parser.add_argument('--offtarget', action = 'store', type = str, dest = 'other', 
		help = "FASTA file containing off target sequence to search for 
				non-unique cut sites. Contained target sequence will be
				ignored.", default = None)
	parser.add_argument('--enzymes', action = 'store', type = str, 
		dest = 'enzymes', help = "tab-delimited file containing restriction 
				enzyme recognition sites ()")
	args = parser.parse_args()
	
	main(args.sequence, args.offtarget, args.enzymes)	

