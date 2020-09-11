# sexPhase
A pipeline that was used to identify X-specific and Y-specific reads in Ficus hispida. The pipeline includes four scripts that are archived in the repo. 
- **Step 1:**
The first script `1.determineGeno.pl` was used to determine which genotype belong to X chromosomes and which belong to Y chromosome based on re-sequenced male and female samples.
- **Step 2:**
The second script `2.retrieveReads.py` identified pacbio reads that could cover phased genotypes from whatshap program and classified them into different haplotypes.
- **Step 3:**
The third one `3.identify_sex-phasedBLOCKs.pl` classified the SNP phased blocks into sex-specific block.
- **Step 4:**
The fourth one `4.filterBLOCK_reads.py` filtered and outputed phased pacbio reads that could be assigned into X and Y chromosomes. These phased pacbio reads were further used to de novo assemble sex chromosomes.
