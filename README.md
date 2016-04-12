#synonymRE

This script determines whether a protein coding sequence has possible synonymous mutations that create restriction enzyme cut sites.

For example, the sequence `ATGGCAGCTGCTAGCTAG` (MAAAS*) is two synonymous mutations away from `ATGgcggccgcTAGCTAG` (MAAAS*), which contains a NotI site (lowercase). 

`synonymRE.py` has three arguments:

  `--seq`:  The sequence to search for synonymous mutations to restriction enzyme recognition sites.
  
  `--offtarget`:  An optional sequence used for excluding non-unique restriction enzymes from the analysis. This step may not be necessary, depending on the downstream molecular biology required.
  
  `--enzymes`:  A tab-delimited file containing restriction enzymes and their recognition sites (in the format `NAME  SITE`, e.g., `EcoRI  GAATTC`). The file `re_list.txt` contained in this repository has a standard set of enzymes. 
  
The script outputs a list of possible synonymous mutations to recognition sequences, e.g.,

start	enzyme	wildtype	synonym
8	PsiI	CTATAA	TTATAA
29	PsiI	CTACAA	TTATAA
42	EcoRV 	GATATT	GATATC
45	ClaI 	ATTGAT	ATCGAT
50	PsiI	TTACAA	TTATAA
72	XbaI 	TCTAGG	TCTAGA
72	NruI 	TCTAGG	TCGCGA 
 
