# Zealexin_Dating
Estimation of duplication dates for zealexin genes in Ding et al. (2020).

# Citation
The study associated with the scripts and data in this repository is available as a preprint:

> Yezhang Ding, Philipp R. Weckwerth, Elly Poretsky, Katherine M. Murphy, James Sims, Evan Saldivar, Shawn A. Christensen, Si Nian Char, Bing Yang, Anh-dao Tong, Zhouxin Shen, Karl A. Kremling, Edward S. Buckler, Tom Kono, David R. Nelson, Jörg Bohlmann, Matthew G. Bakker, Martha M. Vaughan, Ahmed S. Khalil, Mariam Betsiashvili, Steven P. Briggs, Philipp Zerbe, Eric A. Schmelz, Alisa Huffaker bioRxiv 2020.03.04.977355; doi: https://doi.org/10.1101/2020.03.04.977355

# Contents
## Scripts
- `Translate.py`: Translate a nucleotide multi-FASTA file with the standard genetic code (NCBI table 1).
- `Translate_Align_Backtranslate.py`: Translate a nucleotide multi-FASTA file to amino acids with the standard genetic code, align them with clustal-omega, then back-translate the aligned sequences to nucleotides.

## Alignments
- `CYP71Z_Aln.fasta`: Back-translated alignment of coding sequences of the CYP71Z sequences from maize B73, maize PH207, and *Sorghum bicolor*.
- `Cyp81_Aln.fasta`: Back-translated alignment of coding sequences of the Cyp81 sequences from maize B73 and *Sorghum bicolor*.

## BEAST_XML
- `CYP71Z.xml`: XML control file for [BEAST2](https://www.beast2.org/) for the dating analysis of the CYP71Z gene cluster.
- `Cyp81.xml`: XML control file for BEAST2 for the dating analysis of the Cyp81 gene cluster.

## Trees
- `CYP71Z.trees.gz`: MCMC trees from BEAST2 (gzipped) for the CYP71Z gene cluster.
- `Cyp81.trees.gz`: MCMC trees from BEAST2 (gzipped) for the Cyp81 gene cluster.
- `CYP71Z_FullGrid.pdf`: Visualization of trees for the CYP71Z genes produced with [Densitree](https://www.cs.auckland.ac.nz/~remco/DensiTree/).
- `Cyp81_FullGrid.pdf`: Visualization of trees for the Cyp81 genes produced with Densitree.

