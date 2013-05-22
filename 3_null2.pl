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

# bin size of allele frequency
my ($bin) = 0.05 ;

# make folder 'null' if it doesn't exist
my ($dir) = 'null' ;
unless ( -d $dir ){system ("mkdir null") ; }

# get all PHS values with allele freqs
open (DAT, "PHS_null.out")|| die "PHS_null.out not found" ;
my (@dt1) = () ;
my ($openf) ;
my ($i) = 0 ;
while ( $openf = <DAT> ){
	chomp $openf ;
	my (@tmp) = split ("\t", $openf) ;
	$dt1[$i][0] = $tmp[4]/$tmp[3] ;
	$dt1[$i][1] = $tmp[5] ;
	$i++ ;
}
close (DAT) ;
my ($n_dt1) = scalar @dt1 ;
my ($c_dt1) = scalar @{$dt1[0]} ;
print "# of dt1 = $n_dt1 $c_dt1\n" ;
&TwoDarrayTest2 (\@dt1) ;

# output null of PHS given allele frequency
for ( my ($i)=0; $i < 1; $i+=$bin ){
	my ($s) = $i ;
	my ($e) = $i+$bin ;
	print "$s $e\n" ;

	my (@phs) = () ;
	my ($tmp) = -999 ;
	for ( my ($j)=0; $j < $n_dt1; $j++ ){
		if ( $s <= $dt1[$j][0] && $dt1[$j][0] < $e ){
			if ( $tmp == $dt1[$j][1] ){next ; }
			push (@phs, $dt1[$j][1]) ;
			$tmp = $dt1[$j][1] ;
		}
	}
	#@phs = sort {$b <=> $a} @phs ;
	
	open (OUT, ">./null/null_$s\_$e\.out")|| die "Can't open output file, ./null/null_$s\_$e.out\n" ;
	foreach my $phs (@phs){
		print OUT "$phs\n" ;
	}
	close (OUT) ;
}

# output null for monomorphic SNPs
my (@phs) = () ;
my ($tmp) = -999 ;
for ( my ($j)=0; $j < $n_dt1; $j++ ){
	if ( $dt1[$j][0] == 1 ){
		if ( $tmp == $dt1[$j][1] ){next ; }
		push (@phs, $dt1[$j][1]) ;
		$tmp = $dt1[$j][1] ;
	}
}
#@phs = sort {$b <=> $a} @phs ;
	
open (OUT, ">./null/null_1.out")|| die "Can't open output file, ./null/null_1.out\n" ;
foreach my $phs (@phs){
	print OUT "$phs\n" ;
}
close (OUT) ;
		
exit ;

