#PHS code

Code for the PHS statistic from Toomajian et al. (2006) PLoS B, written by Shohei Takuno. This version is preliminary (unpublished), so caveat emptor! Also read the comments on the scripts. If a focal SNP is monomorphic, this script calculates only the second term of equation (1) in Toomajian et al. (2006) PLoS B. This may not be a good way.

As an example, I used the SNP data of chromosome 9 and 10 in the Andean population.

## Required data:

### phased_SNPs

You need the phased data without missing data. The rows indicates individual haplotpes and the columns indicate SNPs. The file name should be chrX.out, where X is 1 - 10.

### cM  
You need physical and genetic distance of all SNPs in phased-SNPs. In this file, the first column is SNP names, the second column is physical distance and the third column is genetic distance. The file name should be chrX.out, where X is 1 - 10. The number of SNPs in phased-SNPs and cM should be the same for every chromosome!!

The original cM files weren't sufficient (only had snp name in column 1, genetic dist in column 2).
I guessed at physical position by extracting it from their snp id (which wasn't always right, since some SNP names weren't position-based, to which I assigned a nearby value), and had to delete some of the candidate snps (file described in the next section) because they didn't have a matching physical position. -MCS

### candidate_SNPs  

This file contains the SNPs that you'd like to test the selection. The first column is chromosome number and the second column is physical position. The file name should be chrX.out, where X is 1 - 10.

## To run
Just run: 
	
	perl 1_ave1.pl
	
	perl 2_null1.pl
	
	perl 3_candidate1.pl
	
	perl 4_pval1.pl

Then, you find PHS_Pval.out that contains P-values of the PHS test. This file contains the statistics and P-values for both alleles at a SNP.
