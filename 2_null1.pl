#! /usr/bin/perl
use strict ;
#use lib "/Users/stakuno/perl/module" ;
#use lib "/home/stakuno/module" ;
use lib "./module" ;
use stat ;
use format ;
use BeginPerlBioinfo ;
use evolve ;
srand ( time|$$ ) ;

# pairwide haplotype-sharing score
# calculate PHS in equation (1) in Toomajian et al. (2006) PLoSB
# for all SNPs as a null distribution given allele frequency

# NOTE: this program is preliminery 
# Null distribution should include the same number of PHS in
# the case that distance between alleles with the same frequency 
# is very close (especially in monomorphic sites)
# 3_null2.pl excludes such PHS values, but in the case of monomorphc SNPs
# auto-correlation between monomorphic sites may be observed 

# NOTE: if a SNP is monomorphic, this program calcualates 
# only the second term of equation (1) in Toomajian et al. (2006) PLoSB
# Not sure it is a correct way

my (@fl) = glob ("./phased_SNPs/*") ;
my ($n_fl) = scalar @fl ;
print "# of input files = $n_fl\n\n" ;

# open file for d_ij and sigma_ij
my ($openf) ;
my (%av) = () ;
my (%sd) = () ;
open (DAT, "d_ij.out")|| die "d_ij.out not found\n" ;
while ( $openf = <DAT> ){
	chomp $openf ;
	my (@tmp) = split ("\t", $openf) ;
	my ($k) = "$tmp[0]\t$tmp[1]" ;
	$av{$k} = $tmp[2] ;
	$sd{$k} = $tmp[3] ;

}
close (DAT) ;


my (@dt1) = () ;

my ($n_sam) = scalar &fileopen ("$fl[0]") ;
print "# of sample = $n_sam\n\n" ;

open (OUT, ">PHS_null.out")|| die "can't open output file, PHS_null.out\n" ;


# for each chromosome
foreach my $fl (@fl){

	my ($chr) = $fl ;
	$chr =~ s/\.\/phased\_SNPs\/// ;

	# open phased SNP file
	@dt1 = () ;
	open (DAT, "$fl")|| die "$fl not found\n" ;
	my ($i) = 0 ;
	while ( $openf = <DAT> ){
		chomp $openf ;
		@{$dt1[$i]} = split ("", $openf) ;
		$i++ ;
	}
	close (DAT) ;
	my ($n_dt1) = scalar @dt1 ;
	my ($c_dt1) = scalar @{$dt1[0]} ;
	print "  $fl, # of samples = $n_dt1, # of SNPs = $c_dt1\n" ;
	&TwoDarrayTest2 (\@dt1) ;
	if ( $n_sam != $n_dt1 ){print "Error\n# of samples are different among chromosomes $n_sam $n_dt1\n" ; exit ; }

	# open Genetic distance file
	my ($of_cm) = $fl ;
	$of_cm =~ s/phased\_SNPs/cM/ ;
	my (@hoge) = &fileopen ("$of_cm") ;
	my (@snp)  = &columnTab (\@hoge, \0) ;
	my ($n_snp) = scalar @snp ;
	my (@dist) = &columnTab (\@hoge, \1) ;
	my ($n_dist) = scalar @dist ;
	print "  # of SNPs in genetic distance = $n_dist\n" ;
	if ( $n_dist != $c_dt1 ){print "Error\n# of SNPs are different in SNP and genetic distance files, $n_dist $c_dt1\n" ; exit ; }
	if ( $n_snp != $n_dist ){print "Error\n#s of data are different between 1st and 2nd columns in a cM file." ; exit ; }
	#print "  $dist[0] $dist[10] $dist[20] $dist[30] $dist[50]\n" ;
	
	# calculate PHS at each site
	for ( my ($h)=0; $h < $n_dist; $h++ ){
		# extract a SNP in $h th column 
		my (@nt) = &column (\@dt1, \$h) ;
		my ($n_nt) = scalar @nt ;
		if ( $n_nt != $n_dt1 ){print "Error\n# of SNPs and sample size are different $n_nt $n_dt1\n" ; exit ; }
		my (@nt2) = grep (!/N/, @nt) ;
		my ($n_nt2) = scalar @nt2 ;
		if ( $n_nt2 <= 1 ){next ; }
		
		# check number of allele
		my (@al) = &union (@nt2) ;
		my ($n_al) = scalar @al ;
		if ( $n_al == 0 ){print "Error\n# of allele is zero\n" ; exit ; }
		if ( $n_al > 2 ){print "Error\n# of allele is >2\n" ; exit ; }
		#print "" ;
		
		# allele frequency
		my ($a1) = scalar grep (/^$al[0]$/, @nt2) ;
		my ($a2) = scalar grep (/^$al[1]$/, @nt2) ;
		if ( $n_nt2 != $a1+$a2 ){print "Error in allele frequency, $n_nt2 $a1 $a2\n" ; exit ; }
		
		my (@ttl) = () ;
		my (@hap1) = () ;
		my (@hap2) = () ;
		
		for ( my ($i)=0; $i < $n_sam-1; $i++ ){
		
			if ( $nt[$i] eq "N" ){next ; }
		
			for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
			
				if ( $nt[$j] eq "N" ){next ; }
			
				# get j_ij and sigma_ij
				my ($ky) = "$i\t$j" ;
				unless (exists $av{$ky} ){print "Error\nd_ij not found\n" ; exit ; }
				unless (exists $sd{$ky} ){print "Error\nsigma_ij not found\n" ; exit ; }

				my ($sts) = $dist[$h] ;
				# forward	
				for ( my ($p1)=($h+1); $p1 < $c_dt1; $p1++ ){
					if ( $dt1[$i][$p1] eq "N" || $dt1[$j][$p1] eq "N" ){next ; }
				
					if ( $dt1[$i][$p1] ne $dt1[$j][$p1] || $p1 == $c_dt1-1 ){					
						$sts = $dist[$p1-1] ; last ;
					}
				}
				# backword
				my ($end) = $dist[$h] ;
				for ( my ($p2)=($h-1); $p2 >= 0; $p2-- ){
					if ( $dt1[$i][$p2] eq "N" || $dt1[$j][$p2] eq "N" ){next ; }
				
					if ( $dt1[$i][$p2] ne $dt1[$j][$p2] || $p2 == 0 ){
						$end = $dist[$p2+1] ; last ;
					}
				}
				
				my ($hap_dist) = ($sts-$end-$av{$ky})/$sd{$ky} ;
				push (@ttl, $hap_dist) ;
				
				# allele 1
				if ( $n_al == 2 && $nt[$i] eq $al[0] && $nt[$j] eq $al[0] ){
					push (@hap1, $hap_dist) ;
				}
				# allele 2
				if ( $n_al == 2 && $nt[$i] eq $al[1] && $nt[$j] eq $al[1] ){
					push (@hap2, $hap_dist) ;
				}
				
			} # $j
		} # $i
		
		my ($whole) = &mean (@ttl) ;
		if ( scalar @ttl != $n_nt2*($n_nt2-1)/2 ){print "Error\n# of pairwise comparison 1\n" ; exit ; }
		
		if ( $n_al == 1 ){
			print OUT "$chr\t$snp[$h]\t$al[0]\t$a1\t$a1\t$whole\n" ;
		}
		elsif ( $n_al == 2 ){
			
			# allele 1
			unless ( $a1 == 1 ){
				if ( scalar @hap1 != $a1*($a1-1)/2 ){print "Error\n# of pairwise comparison 2\n" ; exit ; }
				
				my ($allele1) = &mean (@hap1) ;
				$allele1 -= $whole ;
				print OUT "$chr\t$snp[$h]\t$al[0]\t$n_nt2\t$a1\t$allele1\n" ;
			}
			# allele 2
			unless ( $a2 == 1 ){
				if ( scalar @hap2 != $a2*($a2-1)/2 ){print "Error\n# of pairwise comparison 3\n" ; exit ; }
				
				my ($allele2) = &mean (@hap2) ;
				$allele2 -= $whole ;
				print OUT "$chr\t$snp[$h]\t$al[1]\t$n_nt2\t$a2\t$allele2\n" ;
			}
		}
		
	}
	
	print "\n" ;
}


close (OUT) ;

		
exit ;

