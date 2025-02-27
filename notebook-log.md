Name of article: A Genetic Analysis of Taoyuan Pig and Its Phylogenetic Relationship to Eurasian Pig Breeds

Link to article: https://pmc.ncbi.nlm.nih.gov/articles/PMC4341094/#sec7

The Taoyaun pig is a native Taiwan breed that was introduced around 1877.  It is believed that many of Taiwan’s pig breeds arrived with human migration from mainland China to Taiwan
In this study, mitochondrial DNA sequences of the D-loop region are used to evaluate the maternal lineage and genetic diversity between the Taoyuan pig and other Lower Changjiang River Basin and Central China Type pig breeds. Samples from 29 different pig breeds will be used.

The genotype data comes from NCBI genbank and used refernce mtDNA for the D-loop region. The links to the sequence regions were provided in the article. I collected all of the reads from each of the 29 breeds and entered them into a text file named Pigs_phylo_raw_data.fasfa 

This textfile then got entered into ClustalW and got aligned, and the aligned data can be found under the file pigs_phylo_raw_dataaligned.fasta

steps for alignment: 
in terminal: 
clustalw -ALIGN -INFILE=Pigs_phylo_raw_data.fasta -OUTFILE=Pigs_phylo_raw_dataaligned.fasta -OUTPUT=FASTA

alignment score: 2958718

Allignment method chosen: ClustalW
Description: a multiple alignment algorithm that aligns biological sequences (ex: DNA sequences) by following a stepwise process of: pairwise alignment, guide tree construction, and progressive alignment. 
Assumptions: sequences with higher similarity share a closer evolutionary relationship, evolutionary events such as substitutions, insertions, or deletions occur independently, gap penalties reflect biological realities
Limitations:early mistakes in pairwise alignment will follow through to the progressive alignment, gap pentalies are predefined and may not be evolutionarily accurate, ClustalW is less efficient than more modern MSA tools like MUSCLE for large data sets
Source: https://doi.org/10.1093/nar/22.22.4673

Parsimony-based Tree using R:

install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

library(ape)
library(adegenet)
library(phangorn)

> setwd("~/Desktop/pl-path/Finalproject")
> dna <- fasta2DNAbin(file="Pigs_phylo_raw_dataaligned.fasta")
 D <- dist.dna(dna, model="TN93")
> tre <- nj(D)
> tre <- ladderize(tre)
> plot(tre, cex=.6)
> title("A simple NJ tree")

Tree saved as "Simple NJ tree R.png"

Chosen algorithm: R distance-based methods
Description:NJ (neighbor joining) algorithm is a distance based method used to construct phylogenetic trees. NJ trees are based on a pairwise genetic distance matrix and is then processed using a hierachial clustering algorithm to infer phylogenetic relationships. 
Strengths: Fast, easy to use, integrates well with other R packages
Weaknesses: No built in model comparison, may be inaccurate depending on distance and clustering algorithm 
Assumptions: Assumes that the evolutionary distances can be accurately summarized in a pairwise matrix

