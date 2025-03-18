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

RAxML: (randomized axelerated Maxiumum Likelihood)

description: phylogenetic inference software that uses maxiumum likelihood methods to estimate evolutionary trees. 

Strengths: Best scoring tree for large data sets as well as a fast maximum likelihood algorithm

Weaknesses: Has shown to have a higher variance compared to other models such as IQ tree and does not support Bayesian inference or distance based methods. 

Assumptions: Mutation rate and base changes are constant across evolutionary space and time, sites evolve independently of one another, and all sites evolve at the same time. 


    Check for the corrected file:
    ./raxml-ng --check --msa pigs_phylo_raw_dataaligned.fasta --model GTR+G 
    Execution log saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs_phylo_raw_dataaligned.fasta.raxml.log
    
    Rerun with fixed file:
    ./raxml-ng --check --msa pigs_phylo_raw_dataaligned.fasta.raxml.reduced.phy --model GTR+G
    
    Infer the tree:
       ./raxml-ng --msa pigs_phylo_raw_dataaligned.fasta.raxml.reduced.phy --model GTR+G --prefix pigs --threads 2 --seed 2  
       
        Best ML tree with collapsed near-zero branches saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestTreeCollapsed
                Best ML tree saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestTree
                All ML trees saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.mlTrees
                Optimized model saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestModel

            pigs.raxml.bestTree
    
    in R:
        library(ape)
        > tre <-read.tree(text="(EF590178.1:0.000001,(((EF590165.1:0.002004,EF590169.1:0.000994):0.000001,EF590172.1:0.000001):0.000001,(((EF590161.1:0.000992,EF590146.1:0.004042):0.000001,EF590190.1:0.000001):0.000001,((((AF276925.1:0.000984,(AB015085.1:0.005062,(AB015094.1:0.013560,AB015087.1:0.000848):0.002132):0.002014):0.000001,AB015092.1:0.000001):0.000001,AF276932.1:0.000985):0.000986,((EF590149.1:0.000001,EF590141.1:0.000995):0.000994,EF590176.1:0.000001):0.000001):0.001977):0.000985):0.000985,AB015091.2:0.000987);")
            > tre2= root(tree, outgroup="AB015094.1")
            > plot(tre)
                Output:
                > 2025-03-13 14:11:16.652 R[50321:2286929] +[IMKClient subclass]: chose IMKClient_Modern
                2025-03-13 14:11:16.652 R[50321:2286929] +[IMKInputSession subclass]: chose IMKInputSession_Modern
                
            Other commands for plotting tree:
                edgelabels()
                ?root
                nodelabels()
                tre2= root(tree, outgroup="find from paper"
                nameing:  pdf("mytree.pdf")
                > dev.off()

