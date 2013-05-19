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
# calculate d_ij and sigma_ij in equation (2) in Toomajian (2006) PLoSB

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

my ($openf) ;
my (@dt1) = () ;

my ($n_sam) = scalar &fileopen ("$fl[0]") ;
print "# of samples = $n_sam\n\n" ;

my (@res) = () ;
my ($h) = 0 ;
for ( my ($i)=0; $i < $n_sam-1; $i++ ){
	for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
		$res[$h][0] = $i ;
		$res[$h][1] = $j ;
		$h++ ;
	}
}

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
	
	
	# calculate d_ij
	my ($h) = 0 ;	
	for ( my ($i)=0; $i < $n_sam-1; $i++ ){
		for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
			if ( $res[$h][0] != $i ){print "Error\nhoge1\n" ; exit ; }
			if ( $res[$h][1] != $j ){print "Error\nhoge2\n" ; exit ; }

			my (@d) = () ;
			my (@resT) = () ;
			for ( my ($k)=0; $k < $c_dt1; $k++ ){
				if ( $dt1[$i][$k] eq $dt1[$j][$k] ){
					push (@d, $dist[$k]) ;
				}
				else {
					my ($n_d) = scalar @d ;
					if ( $n_d <= 1 ){}
					else {
						my ($distance) = &max (@d) - &min (@d) ;	
						push (@resT, $distance) ;
					}
					@d = () ;
				}
			}
			my ($n_d) = scalar @d ;
			if ( $n_d <= 1 ){}                                 
			else {
				my ($distance) = &max (@d) - &min (@d) ;
				push (@resT, $distance) ;
			}


			#print "@resT\n" ; exit ;
			@{$res[$h]} = (@{$res[$h]}, @resT) ;

			$h++ ;
		}
	}


	print "\n" ;
	#last ; 
} # @fl

open (OUT, ">d_ij.out")|| die "ooooo\n" ;

my ($h) = 0 ;
for ( my ($i)=0; $i < $n_sam-1; $i++ ){
	for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
		#my ($n_res) = scalar @{$res[$h]} ;
		#print "$n_res\n" ;
		if ( $res[$h][0] != $i ){print "Error\nfoo1\n" ; exit ; }
		if ( $res[$h][1] != $j ){print "Error\nfoo2\n" ; exit ; }
		my (@tmp) = @{$res[$h]} ;
		shift @tmp ; shift @tmp ;
		my ($mean) = &mean (@tmp) ;
		my ($sd) = &sd (@tmp) ;
		print OUT "$i\t$j\t$mean\t$sd\n" ;

		$h++ ;
	}
}

close (OUT) ;


		
exit ;


