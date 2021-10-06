# Simulations

This file contains the script used to simulate continuous epistatic traits from Grealey et al (in preparation) 2021. 

The intake for the file requires a file (either CSV or hdf5 format) which has SNPs in the rows and samples in the columns. 

an example command is run like so:

python3 simulate_cont_traits_fixed_aug_19.py --epistasis_level=0.1 --heritability=0.4 --epi_model= --NSNPs=10 --NEPIs=10 --outdir=test2_h40_epi10  --loadstr=data/data.h5 --seed 10

where:

--epistasis level fixes the proportion of epistatic heritbaility

--heritability is the heritability

--NSNPs is the number of SNPs that have a main effect only

--NEPIs is the number of SNPs that have both a main and an epistatic effect - must be a multiple of 2 otherwise they can't be placed into  epistatic pairs 

--seed fixes the random seed 

--epi_model is the epistatic model of interactions present within the phenotype for each pair from NEPIs number of SNPs.

--loadstr is the file (with extension) containing SNPs on the rows and samples on the columns. Currently the file should be in hdf5 format - https://pandas.pydata.org/docs/reference/api/pandas.read_hdf.html - however, there are functions in the code that can be used to read CSV.

--outdir is the directory where the files will be saved - it is okay if the directory is in existence, however, any files will be overwritten if they have the same savestring (Seed, number of SNPs etc). Best to keep these directories separate.

