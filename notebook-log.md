Name of article: A Genetic Analysis of Taoyuan Pig and Its Phylogenetic Relationship to Eurasian Pig Breeds

Link to article: https://pmc.ncbi.nlm.nih.gov/articles/PMC4341094/#sec7

The Taoyaun pig is a native Taiwan breed that was introduced around 1877.  It is believed that many of Taiwan’s pig breeds arrived with human migration from mainland China to Taiwan
In this study, mitochondrial DNA sequences of the D-loop region are used to evaluate the maternal lineage and genetic diversity between the Taoyuan pig and other Lower Changjiang River Basin and Central China Type pig breeds. Samples from 29 different pig breeds will be used.

The genotype data comes from NCBI genbank and used refernce mtDNA for the D-loop region. The links to the sequence regions were provided in the article. I collected all of the reads from each of the 29 breeds and entered them into a text file named Pigs_phylo_raw_data.fasfa 


# Clustal W:

Allignment method chosen: ClustalW
Description: a multiple alignment algorithm that aligns biological sequences (ex: DNA sequences) by following a stepwise process of: pairwise alignment, guide tree construction, and progressive alignment. 
Strengths: Good sensitivity for both closely and distantly related sequences, position-specific gap penalties imrpove biological probability, sequence weighting reduces alighnement bias from overrepresented or cloosely related sequences, and it is widely adopted, well-documented, and available for multiple platforms. 
Limitations:early mistakes in pairwise alignment will follow through to the progressive alignment, gap pentalies are predefined and may not be evolutionarily accurate, ClustalW is less efficient than more modern MSA tools like MUSCLE for large data sets
Assumptions: sequences with higher similarity share a closer evolutionary relationship, evolutionary events such as substitutions, insertions, or deletions occur independently, gap penalties reflect biological realities, assumes progresive alignment is appropriate 
User choices: scoring matrices for nucleotides, gap opening and extension penalties for different datasets such as highly conserved or variable regions, can turn sequence weighting on and off
Source: https://doi.org/10.1093/nar/22.22.4673

dowloaded using tutorial on: http://www.clustal.org/clustal2/ -> clustalw-2.1-macosx.dmg
ClustalW 2.1 was downloaded

This textfile then got entered into ClustalW and got aligned, and the aligned data can be found under the file 
    pigs_phylo_raw_dataaligned.fasta

steps for alignment in terminal: 
```
clustalw -ALIGN -INFILE=Pigs_phylo_raw_data.fasta -OUTFILE=Pigs_phylo_raw_dataaligned.fasta -OUTPUT=FASTA
alignment score: 2958718
```


# Parsimony-based Tree using R:

Chosen algorithm: R distance-based methods
Description:NJ (neighbor joining) algorithm is a distance based method used to construct phylogenetic trees. NJ trees are based on a pairwise genetic distance matrix and is then processed using a hierachial clustering algorithm to infer phylogenetic relationships. 
Strengths: Fast, easy to use, integrates well with other R packages, commonly used for initial tree building, can be paired with bootstrap resampling to assess reliability
Weaknesses: No built in model comparison, may be inaccurate depending on distance and clustering algorithm, produces unrooted trees that require additional assumptions
Assumptions: Assumes that the evolutionary distances can be accurately summarized in a pairwise matrix, branch lengths are additive 
User choices: distance model (choose how to compute genetic distances), the tree inference agorithm, users can root the tree by giving an outgroup

Downloaded: Downloaded R for mac at this link: https://cran.rstudio.com/
    dowloaded R-4.1.1-arm64.pkg
    current version: version 4.4.2
In R:
download the ape package:
'''
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
'''

'''
library(ape)
library(adegenet)
library(phangorn)
'''

'''
setwd("~/Desktop/pl-path/Finalproject")
dna <- fasta2DNAbin(file="Pigs_phylo_raw_dataaligned.fasta")
 D <- dist.dna(dna, model="TN93")
tre <- nj(D)
tre <- ladderize(tre)
plot(tre, cex=.6)
title("A simple NJ tree")
'''

Tree saved as "Simple NJ tree R.png"


# RAxML: (randomized axelerated Maxiumum Likelihood)

description: phylogenetic inference software that uses maxiumum likelihood methods to estimate evolutionary trees. 
Strengths: Best scoring tree for large data sets as well as a fast maximum likelihood algorithm. Its main strength is its a fast maximum likelihood algroithrm that still produces good likelihood score. The program is constatnly being maintained and extended for growing datasets. 
Weaknesses: Has shown to have a higher variance compared to other models such as IQ tree and does not support Bayesian inference or distance based methods. 
Assumptions: Mutation rate and base changes are constant across evolutionary space and time, sites evolve independently of one another, and all sites evolve at the same time. 
User choices: Choice of substitution model, input data choices (FASTA or PHYLIP), starting tree choice (random or user-specified),number of boostrap replicates

    dowload: followed tutorial on github: https://github.com/amkozlov/raxml-ng
        dowloaded file for mac: raxml-ng_v1.0.2_macos_x86_64 and put it in my pl-path folder 
        version: 1.2.2

    Check for the corrected file:
    '''
    ./raxml-ng --check --msa pigs_phylo_raw_dataaligned.fasta --model GTR+G 
    Execution log saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs_phylo_raw_dataaligned.fasta.raxml.log
    '''
    
    Rerun with fixed file:
    '''
    ./raxml-ng --check --msa pigs_phylo_raw_dataaligned.fasta.raxml.reduced.phy --model GTR+G
    '''
    
    Infer the tree:
    '''
       ./raxml-ng --msa pigs_phylo_raw_dataaligned.fasta.raxml.reduced.phy --model GTR+G --prefix pigs --threads 2 --seed 2  
    '''
       
        Best ML tree with collapsed near-zero branches saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestTreeCollapsed
                Best ML tree saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestTree
                All ML trees saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.mlTrees
                Optimized model saved to: /Users/lindsaypropst/Software/raxml-ng_v1.2.2_macos (1)/pigs.raxml.bestModel

            pigs.raxml.bestTree
    
    in R:
    '''
        library(ape)
        > tre <-read.tree(text="(EF590178.1:0.000001,(((EF590165.1:0.002004,EF590169.1:0.000994):0.000001,EF590172.1:0.000001):0.000001,(((EF590161.1:0.000992,EF590146.1:0.004042):0.000001,EF590190.1:0.000001):0.000001,((((AF276925.1:0.000984,(AB015085.1:0.005062,(AB015094.1:0.013560,AB015087.1:0.000848):0.002132):0.002014):0.000001,AB015092.1:0.000001):0.000001,AF276932.1:0.000985):0.000986,((EF590149.1:0.000001,EF590141.1:0.000995):0.000994,EF590176.1:0.000001):0.000001):0.001977):0.000985):0.000985,AB015091.2:0.000987);")
            > tre2= root(tree, outgroup="AB015094.1")
            > plot(tre)
    '''
                Output:
                > 2025-03-13 14:11:16.652 R[50321:2286929] +[IMKClient subclass]: chose IMKClient_Modern
                2025-03-13 14:11:16.652 R[50321:2286929] +[IMKInputSession subclass]: chose IMKInputSession_Modern
                
                tree saved as "RAxML tree"


# MrBayes

Description: Program for Bayesian inference of a wide variety of phylogenic and evolutionary models using the MCMC (markov chain Monte Carlo). Basian inference estimates the probability under prior evidence 
Strengths: Accomindates data heterogeneity, increases chance of escaping local optima, speeds up convergence especially for large date, supports mixed models, incorporates model undertainty (estimate paramaters such as branch length and substitution rates)
Weaknesses: Choice of priors can change outcomes significantly => requires in-depth background of the program, bootstrap support values for ML or parsimony tend to be lower than posterior probabilities for Bayesian approaches 
assumptions: Assumes equal rates of evolution across sites, assumes sites evolve independently 
Users choices: set the model structure, link/unlink topologies, branch length, substitution rates, proportion of invariable sites, Gamma distributino 

downloaded using the following command:

'''
     brew tap brewsci/bio
     brew install mrbayes open-mpi
'''
     
     convert data to .nex file to run in MrBayes - aligned fasta file in R named Pigs_phylo_raw_dataaligned_copy_2.nexus

'''
    begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=1000000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
 outgroup AB015094.1;
 mcmc;
 sumt;
end;
'''
    
running Mrbayes in terminal:
    open finalproject in terminal
    
    '''
    which mb
    mb Pigs_phylo_raw_dataaligned_copy_2.nexus 
    '''

    Output saved in final project folder as "MrBayes.md"
    
    used the program icytree to turn the output into a tree, enter con.tre file 
    
    
#BUCKy

Description: BUCKy (bayesian untangling of concordance knots) uses the bayesian method to esitmade the primary concordance tree using posterior distribution from multiple loci. 
Strengths: Accounts for gene heterogeneity, concordance factors are meaningful and interpretable, it incorporates uncertnity, can incorporate many loci, and dirichlet process flexibility
Weaknesses: It is computationally intensive and requires significant computational effort, relies on accurate gene tree posteriors, only models ILS-driven discordance, sensitive to alpha parameter choice, and assumes independent loci
Assumptions: loci are unlinked and independently evolving, discordance is primarily caused by incomplete lineage sorting, dirichlet process approriately models gene tree clustering, and assumes convergence of MCMC chains
User choices: alpha paramater for the dirichlet process, MCMC iteraitions, burn-in length, gene tree posterior distributions

    note- my project does not involve multiple genes so I am using the "yeast" toy data by the program and this is not included in my final project. 
    
    program was downloaded at this link: https://pages.stat.wisc.edu/~ane/bucky/downloads.html and is version 1.4.4 from tarball. Path= /Users/lindsaypropst/Documents/bucky-1.4.4/data/yeast
    
    commands were follwed using the slides at this link: https://pages.stat.wisc.edu/~larget/AustinWorkshop/tutorial.pdf

inside yeast folder:

running mbsum
'''
mbsum -n 501 -o y000.in y000.run*.t
'''

running bucky
'''
./bucky -a 1 -k 4 -n 1000000 -c 4 -s1 23546 -s2 4564 -o yeast y*/*.in
'''
