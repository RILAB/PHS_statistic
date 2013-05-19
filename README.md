{\rtf1\ansi\ansicpg932\cocoartf949\cocoasubrtf540
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww15500\viewh13520\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\ql\qnatural\pardirnatural

\f0\fs28 \cf0 PHS test\
\
NOTE:  This version is preliminary (unpublished)!!  Please be careful!!\
NOTE:  Also read the comments on the scripts!!\
NOTE: If a focal SNP is monomorphic, this script calculates only the \
		second term of equation (1) in Toomajan et al. (2006) PLoS B\
		This may not be a good way.\
\
As an example, I used the SNP data of chromosome 9 and 10\
in the Andean population.\
\
Required data:\
	phased_SNPs: You need the phased data without missing data.\
				   The rows indicates individual haplotpes and the columns\
				   indicate SNPs.  	\
				  The file name should be chrX.out, where X is 1 - 10.\
	cM: You need physical and genetic distance of all SNPs in pahsed-SNPs.\
		In this file, the first column is SNP names, the second column is physical \
		distance and the third column is genetic distance.\
		The file name should be chrX.out, where X is 1 - 10.\
		NOTE:  The number of SNPs in phased-SNPs and cM should be the same for every\
		chromosome!!\
\
	candidate_SNPs: This file contains the SNPs that you'd like to test the selection.\
					The first column is chromosome number and the second column is\
					physical position.\
					The file name should be chrX.out, where X is 1 - 10.\
\
Just run\
perl 1_ave1.pl\
perl 2_null1.pl\
perl 3_candidate1.pl\
perl 4_pval1.pl\
\
Then, you find PHS_Pval.out that contains P-values of the PHS test.\
This file contains the statistics and P-values for both alleles at a SNP.\
\
\
shohei\
\
\
\
\
}