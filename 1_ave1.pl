 #!/usr/bin/env perl
use strict ;
#use lib "/Users/stakuno/perl/module" ;
#use lib "/home/stakuno/module" ;
use lib "./module" ;
use stat ;
use format ;
use BeginPerlBioinfo ;
use evolve ;
srand ( time|$$ ) ;

# pairwide haplotype-sharing score (PHS)
# calculate d_ij and sigma_ij in equation (2) in Toomajian et al. (2006) PLoSB

# the order of individuals should be the same among chromosomes

my (@fl) = glob ("./phased_SNPs/*") ;
my ($n_fl) = scalar @fl ;
print "# of input files = $n_fl\n\n" ;

my ($openf) ;
my (@dt1) = () ;

my ($n_sam) = scalar &fileopen ("$fl[0]") ;
print "# of samples = $n_sam\n\n" ;

# set output array
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
	# A, T, G, C and N (missing data) are allowed
	@dt1 = () ;
	open (DAT, "$fl")|| die "$fl not found\n" ;
	my ($i) = 0 ;
	while ( $openf = <DAT> ){
		chomp $openf ;
		@{$dt1[$i]} = split ("", $openf) ;
		$i++ ;
	}
	close (DAT) ;
	my ($n_dt1) = scalar @dt1 ; # Number of samples
	my ($c_dt1) = scalar @{$dt1[0]} ; # Number of SNPs
	print "  $fl, # of samples = $n_dt1, # of SNPs = $c_dt1\n" ;
	&TwoDarrayTest2 (\@dt1) ;
	if ( $n_sam != $n_dt1 ){print "Error\n# of samples are different among chromosomes $n_sam $n_dt1\n" ; exit ; }
	
	# open Genetic distance file
	# 1st col: marker name
	# 2nd col: cM 
	my ($of_cm) = $fl ;
	$of_cm =~ s/phased\_SNPs/cM/ ;
	my (@dist) = &fileopen ("$of_cm") ;
	@dist = &columnTab (\@dist, \1) ;
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

				unless ( $dt1[$i][$k] =~ /[ATGCN]/ && $dt1[$j][$k] =~ /[ATGCN]/ ){
					print "Error\nBad nucleotide: $dt1[$i][$k] $dt1[$j][$k].\nThe data should contain only ATGC or N\n" ; exit ;
				}

				if ( $dt1[$i][$k] eq "N" || $dt1[$j][$k] eq "N" ){
					next ; 
				}
				elsif ( $dt1[$i][$k] eq $dt1[$j][$k] ){
					push (@d, $dist[$k]) ;
				}
				elsif ( $dt1[$i][$k] ne $dt1[$j][$k] || $k == ($c_dt1-1) ){
					if ( $k == ($c_dt1-1) && $dt1[$i][$k] eq $dt1[$j][$k] ){
						push (@d, $dist[$k]) ;
					}
				
					my ($n_d) = scalar @d ;
					if ( $n_d <= 1 ){}
					else {
						my ($distance) = &max (@d) - &min (@d) ;	
						push (@resT, $distance) ;
					}
					@d = () ;
				}
			}
			
			#print "@resT\n" ; exit ;
			@{$res[$h]} = (@{$res[$h]}, @resT) ;

			$h++ ;
		}
	}


	print "\n" ;
	#last ; 
} # @fl



open (OUT, ">d_ij.out")|| die "can't open output file, d_ij.out\n" ;

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







