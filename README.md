This script determines whether a protein coding sequence has possible synonymous mutations that create restriction enzyme cut sites.

For example, the sequence `ATGGCAGCTGCTAGCTAG` (MAAAS*) is two synonymous mutations away from `ATGgcggccgcTAGCTAG` (MAAAS*), which contains a NotI site (lowercase). 

`synonymRE.py` has four arguments:

  `--seq`:  The sequence to search for synonymous mutations to restriction enzyme recognition sites.
  
  `--offtarget`:  An optional sequence used for excluding non-unique restriction enzymes from the analysis. This step may not be necessary, depending on the downstream molecular biology required.
  
  `--enzymes`:  A tab-delimited file containing restriction enzymes and their recognition sites (in the format `NAME  SITE`, e.g., `EcoRI  GAATTC`). The file `re_list.txt` contained in this repository has a standard set of enzymes. 
  
  `--8cut`: Should restriction enzymes with 8-base recognition sequences be searched as well (default=False)? Because of the inefficient way I implemented the synonymous codon search (by enumerating all k-mers and filtering for synonymous sequences), 8-base cutters take significant amount of time to run.
  
  
