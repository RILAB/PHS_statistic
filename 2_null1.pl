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
# calculate PHS in equation (1) in Toomajan (2006) PLoSB
# for all SNPs as a null distribution given allele frequency

# NOTE: this program is preliminery 
# Null distribution should include the same number of PHS in
# the case that distance between alleles with the same frequency 
# is very close (especially in monomorphic sites)

# NOTE: if a SNP is monomorphic, this program calcualates 
# only the second term of equation (1) in Toomajan (2006) PLoSB
# I'm not sure it is a correct way

my (@fl) = glob ("./phased_SNPs/*") ;
my ($n_fl) = scalar @fl ;
print "# of input files = $n_fl\n\n" ;

#my (@at) = qw (1 2 3 4 5 6 7 8 9 10) ;
#my (@at) = qw (9 10) ;

#my ($foc1) = $ARGV[0] ;
#my ($foc1) = "Andean" ;
#my ($foc1) = "Mex_Lowland" ;
#my ($foc1) = "Mex_Highland" ;
#my ($foc1) = "SA_Lowland" ;

#print "$foc1\n" ;

# open file for d_ij and sigma_ij
my ($openf) ;
my (%av) = () ;
my (%sd) = () ;
open (DAT, "d_ij.out")|| die "ppppppp\n" ;
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

open (OUT, ">PHS_null.out")|| die "ppp\n" ;


# for each chromosome
foreach my $fl (@fl){

	# open phased SNP file
	@dt1 = () ;
	open (DAT, "$fl")|| die "pppp\n" ;
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
	my (@dist) = &fileopen ("$of_cm") ;
	@dist = &columnTab (\@dist, \2) ;
	my ($n_dist) = scalar @dist ;
	print "  # of SNPs in genetic distance = $n_dist\n" ;
	if ( $n_dist != $c_dt1 ){print "Error\n# of SNPs are different in SNP and genetic distance files, $n_dist $c_dt1\n" ; exit ; }
	#print "  $dist[0] $dist[10] $dist[20] $dist[30] $dist[50]\n" ;
	
	# calculate PHS
	# cut first and last 5 SNPs
	for ( my ($h)=4; $h < $n_dist-5; $h++ ){
		# extract a SNP in $h th column 
		my (@nt) = &column (\@dt1, \$h) ;
		# check number of allele
		my (@al) = &union (@nt) ;
		my ($n_al) = scalar @al ;
		#print "" ;
	
		# biallelic
		if ( $n_al == 2 ){
			
			# allele frequency
			my ($a1) = scalar grep (/^$al[0]$/, @nt) ;
			my ($a2) = scalar grep (/^$al[1]$/, @nt) ;
			my (@ttl) = () ;
			my (@hap1) = () ;
			my (@hap2) = () ;

			for ( my ($i)=0; $i < $n_sam-1; $i++ ){
				for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
						# get j_ij and sigma_ij
						my ($ky) = "$i\t$j" ;
						unless (exists $av{$ky} ){print "Error\nd_ij not found\n" ; exit ; }
						unless (exists $sd{$ky} ){print "Error\nsigma_ij not found\n" ; exit ; }

						my ($sts) = $dist[$h] ;
						# forward	
						for ( my ($p1)=($h+1); $p1 < $c_dt1; $p1++ ){
							if ( $dt1[$i][$p1] ne $dt1[$j][$p1] || $p1 == $c_dt1-1 ){					
								$sts = $dist[$p1-1] ; last ;
							}
						}

						# backword
						my ($end) = $dist[$h] ;
						for ( my ($p2)=($h-1); $p2 >= 0; $p2-- ){
							if ( $dt1[$i][$p2] ne $dt1[$j][$p2] || $p2 == 0 ){
								$end = $dist[$p2+1] ; last ;
							}
						}
						my ($hap_dist) = ($sts-$end-$av{$ky})/$sd{$ky} ;
						push (@ttl, $hap_dist) ;	
						# allele 1
						if ( $nt[$i] eq $al[0] && $nt[$j] eq $al[0] ){
							push (@hap1, $hap_dist) ;
						}
						# allele 2
						if ( $nt[$i] eq $al[1] && $nt[$j] eq $al[1] ){
							push (@hap2, $hap_dist) ;
						}

				} # $j
			} # $i

			my ($whole) = &mean (@ttl) ;
			# allele 1
			unless ( $a1 == 1 ){
				my ($allele1) = &mean (@hap1) ;
				my ($stat1) = $allele1-$whole ;
				print OUT "$n_sam\t$a1\t$stat1\n" ;
			}

			# allele 2
			if ( $a1 != $a2 && $a2 > 1 ){
				my ($allele2) = &mean (@hap2) ;
				my ($stat2) = $allele2-$whole ;
				print OUT "$n_sam\t$a2\t$stat2\n" ;
			}

		} # if sample size 2

		elsif ( $n_al == 1 ) { # fixed polymorphisms
			my (@ttl) = () ;

			for ( my ($i)=0; $i < $n_sam-1; $i++ ){
				for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
					# ave and sd
					my ($ky) = "$i\t$j" ;
					unless (exists $av{$ky} ){print "Error\nd_ij not found\n" ; exit ; }
					unless (exists $sd{$ky} ){print "Error\nsigma_ij not found\n" ; exit ; }

					my ($sts) = $dist[$h] ;
					# forward       
					for ( my ($p1)=($h+1); $p1 < $c_dt1; $p1++ ){
						if ( $dt1[$i][$p1] ne $dt1[$j][$p1] || $p1 == $c_dt1-1 ){
							$sts = $dist[$p1-1] ; last ;
						}
					}

					# backword
					my ($end) = $dist[$h] ;
					for ( my ($p2)=($h-1); $p2 >= 0; $p2-- ){
						if ( $dt1[$i][$p2] ne $dt1[$j][$p2] || $p2 == 0 ){
							$end = $dist[$p2+1] ; last ;
						}
					}
					my ($hap_dist) = ($sts-$end-$av{$ky})/$sd{$ky} ;
					push (@ttl, $hap_dist) ;
				} # $j
			} # $i
	
			my ($whole) = &mean (@ttl) ;
			# allele 1
			print OUT "$n_sam\t$n_sam\t$whole\n" ;
		} # fixed
		else {
			print "Error\nMore than 2 alleles ($n_al alleles)\n" ; exit ; 
		}

	}
	
	print "\n" ;
}


close (OUT) ;

		
exit ;

