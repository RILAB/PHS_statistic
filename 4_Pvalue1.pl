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
# P-value pf PHS is calcuated given allele frequency (see Read Me).

# bin size of allele frequency
my ($bin) = 0.05 ;

# output file
open (OUT, ">PHS_Pvalue.out")|| die "Can't open output file, PHS_Pvalue.out\n" ;

# get P-values for PHSs given allele freqs
open (DAT, "PHS_null.out")|| die "PHS_null.out not found" ;
my ($openf) ;

while ( $openf = <DAT> ){
	chomp $openf ;
	my (@tmp) = split ("\t", $openf) ;
	my ($af) = $tmp[4]/$tmp[3] ;
	my ($phs) = $tmp[5] ;
	
	# get file for null distribution of PHS
	my ($fl) = "hoge" ;
	if ( $af == 1 ){
		$fl = "./null/null_1.out" ;
	}
	else {
		for ( my ($i)=0; $i < 1; $i+=$bin ){
			my ($s) = $i ;
			my ($e) = $i+$bin ;
			if ( $s <= $af && $af < $e ){
				$fl = "./null/null_$s\_$e\.out" ; last ;
			}
		}
		
	}
	
	# calculate P-value
	my ($pf) ;
	my ($ttl) = 0 ;
	my ($pvl) = 0 ;
	
	open (RES, "$fl")|| die "$fl not found\n" ;
	while ( $pf = <RES> ){
		chomp $pf ;
		$ttl++ ;
		if ( $phs <= $pf ){
			$pvl++ ;
		}
	}
	close (RES) ;

	$pvl /= $ttl ;
	print OUT "$openf\t$pvl\n" ;

}
close (DAT) ;


close (OUT) ;
		
exit ;

