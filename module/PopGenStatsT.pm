#
# $Id: PopGenStats.pm,v 1.1.1.1 2009/11/27 02:13:59 kteshima Exp $
#
#package PopGenStatsT ;

use strict ;
use warnings ;

##########  theta_and_D ##########
#
# theta	: pi within a population
# 	: Watterson's theta (1975) 
# 	: Fay and Wu's theta (2000) 
# 	: Zeng, Fu, Shi and Wu's theta_L (2006) 
# Test of neutrality 
#	: Tajima's D (1989) 
# 	: Fu and Li's D (1993) (using singleton, without outgroup & with outgroup)
# 	 	see also Simonsen, Churchill and Aquadro (1995) 
# 	: Fay and Wu's H (2000) (un-normalized)
# 	: Fay and Wu's H (2000) ( normalized )
# 		variance of H is obtained by Zeng, Fu, Shi and Wu (2006)
# 
# arguments are: 
# 	(1) array of seg
# 	(2) Number of seg sites
# 	(3) sample size
# 	(4) Length of simulated region
# 
##########  theta_and_D ##########

sub	theta_and_D{
	my ( $ss, $ll, $nn, $length ) = @_ ;
	# @$ss: array of seq
	# $$ll: # of seg sites
	# $$nn: # of sample
	# $$length: length of simulated region
	
	##############
	# theta_pi 
	# theta_w: Watterson
	# theta_H: Fay and Wu
	# theta_L: Zeng, Fu, Shi and Wu's
	#############
	my ($pi) = 0 ;
	my ($nsegl) = 0 ; # # of seg sites 
	my ($thetaH) = 0 ;
	my ($thetaL) = 0 ;
	my ($etas) = 0 ; # singleton
	my ($etas2) = 0 ; # singleton, derived
	my ($combination) = $$nn*($$nn-1) / 2 ;
	#print "n = $$nn\n" ;
	#print "c = $combination\n" ;
	
	# for each segregating site
	my ($nsnp) = 0 ; # # of derived alleles
	my (@ntD) = () ;
	for ( my ($itD)=0; $itD < $$ll; $itD++ ){
		
		# count the number of derived alleles
		$nsnp = 0 ;
		@ntD = map {substr ($_, $itD, 1)} @$ss ;
		for ( my ($jtD)=0; $jtD < $$nn; $jtD++ ){
			if ( $ntD[$jtD] == 1 ){
				$nsnp++ ;
			}
			elsif ( $ntD[$jtD] == 0 ){; }
			else {print "Error 1 in theta_and_D, $ntD[$jtD]\n" ; exit ; }
		}
		
		unless ( $nsnp == 0 || $nsnp == $$nn ){ 
			$pi += ($$nn-$nsnp)*$nsnp ;# pi
			$nsegl++ ; 		# W
			$thetaH += $nsnp*$nsnp ;	# H
			$thetaL += $nsnp ;	# L
		}
		if ( $nsnp == 1 || $nsnp == ($$nn-1) ){$etas++ ; } # eta
		if ( $nsnp == 1 ){$etas2++ ; } # eta, derived
	}
	$pi /= $combination ;
	$thetaH /= $combination ;
	$thetaL /= ($$nn-1) ;
	#print "Eta = $etas\n" ;
	
	# Tajima's D, Fu and Li's D and Fay and Wu's H (normalized)
	# define "a1" and "a2"
	my ($a1, $a2) = (0, 0) ;
	foreach my $i ( 1..($$nn-1) ){
		$a1 += 1/$i ;
		$a2 += 1/($i*$i) ;
	}
	#print "a1 = $a1\n" ; print "a2 = $a2\n" ;
	
	# define "b1" and "b2"
	my ($b1) = ($$nn+1)/(3*($$nn-1)) ;
	my ($b2) = 2*($$nn*$$nn+$$nn+3)/(9*$$nn*($$nn-1)) ;
	#print "b1 = $b1\n" ; print "b2 = $b2\n" ;
	
	# define "c1" and "c2"
	my ($c1) = $b1 - 1/$a1 ;
	my ($c2) = $b2 - ($$nn+2)/($a1*$$nn) + $a2/($a1*$a1) ;
	#print "c1 = $c1\n" ; print "c2 = $c2\n" ;
	
	# define "e1" and "e2"
	my ($e1) = $c1/$a1 ;
	my ($e2) = $c2/($a1*$a1+$a2) ;
	#print "e1 = $e1\n" ; print "e2 = $e2\n" ;
	
	# define "ud" and "vd" for Fu & Li without outgroup
	my ($vd) = $a2/$a1/$a1 - 2*(1+1/$a1-$a1+$a1/$$nn)/$$nn - 1/$$nn/$$nn ;
	$vd /= ($a1*$a1+$a2) ;
	my ($ud) = (($$nn-1)/$$nn - 1/$a1)/$a1 - $vd ;

	# define "udp" and "vdp" and $cn for Fu & Li with outgroup
	my ($cn) = ( $$nn == 2 ) ? 1: 2*($$nn*$a1-2*($$nn-1))/($$nn-1)/($$nn-2) ;
	my ($vdp) = 1 + $a1*$a1/($a2+$a1*$a1)*($cn-($$nn+1)/($$nn-1)) ;
	my ($udp) = $a1-1-$vdp ;
	
	# Tajima's D
	my ($tajimaD) = $pi - $nsegl/$a1 ;
	if ( $nsegl < 2 || ($e1*$nsegl + $e2*$nsegl*($nsegl-1)) == 0 ){ $tajimaD = -9999; }
	else { $tajimaD /= sqrt( $e1*$nsegl + $e2*$nsegl*($nsegl-1) ) ; }
	#my ($uuuu) = $nsegl/$a1 ;
	#print "$pi, $nsegl, $uuuu, $tajimaD\n" ;
	
	# Fu and Li's D (using S, without outgroup)
	my ($FuLiD) = $nsegl/$a1 - $etas*($$nn-1)/$$nn ;
	if (  $nsegl >= 2 && $ud*$nsegl + $vd*$nsegl*$nsegl > 0 ){ 
		$FuLiD /= sqrt ( $ud*$nsegl + $vd*$nsegl*$nsegl ) ; }
	else { $FuLiD = -9999 ; }
	
	# Fu and Li's D (using S, with outgroup)
	my ($FuLiD2) = $nsegl - $a1*$etas2 ;
	if (  $nsegl >= 2 && $udp*$nsegl + $vdp*$nsegl*$nsegl > 0 ){ 
		$FuLiD2 /= sqrt ( $udp*$nsegl + $vdp*$nsegl*$nsegl ) ; }
	else { $FuLiD2 = -9999 ; }
	
	# Fay and Wu's H (un-normalized) 
	my ($faywuH) = $pi - $thetaH ;
	
	# Fay and Wu's H (normalized)
	my ($normH) = $pi - $thetaL ;
	my ($theta1) = $nsegl/$a1 ;
	my ($theta2) = $nsegl*($nsegl-1)/($a1*$a1+$a2) ;
	# define "bn"
	my ($bn) = 0 ;
	foreach my $i (1..$$nn){
		$bn += 1/($i*$i) ;
	}
	my $varH = 
		#($$nn-2)*$pi/6/($$nn-1)  # is that right? 
		($$nn-2)*$theta1/6/($$nn-1)
		+ ( (18*$$nn*$$nn*(3*$$nn+2)*$bn - (88*($$nn**3)+9*$$nn*$$nn-13*$$nn+6))*$theta2 )
			/9/$$nn/($$nn-1)/($$nn-1) ;
	if ( $nsegl < 2 || $varH == 0 ){ $normH = -9999 ; }
	else             { $normH /= sqrt($varH) ; }
		
	my (%res) = () ;
	$res{"pi"} = $pi/$$length ;		       # theta_pi: per site
	$res{"thetaW"} = $nsegl/$$length/$a1 ; # theta_w: per site
	$res{"thetaH"} = $thetaH/$$length ;    # theta_H: per site
	$res{"thetaL"} = $thetaL/$$length ;    # theta_L: per site
	$res{"S"} = $nsegl ;                   # # of seg sites
	$res{"tajimaD"} = $tajimaD ;           # Tajima's D
	$res{"fuliD"} = $FuLiD ;               # Fu and Li's D without outgroup
	$res{"fuliD2"} = $FuLiD2 ;             # Fu and Li's D with outgroup
	$res{"faywuH"} = $faywuH ;             # Fay and Wu's H
	$res{"normalH"}= $normH;	           # normalized Fay and Wu's H
	
	return ( %res ) ;
}

##########  pi_B ##########
#
# Theta between allelic class
# 
# arguments are: 
# 	(1) array of seg (allelic class 1)
# 	(2) array of seg (allelic class 2)
# 	(3) Number of seg sites
# 	(4) sample size 1
# 	(5) sample size 2
# 	(6) Length of simulated region
# 
##########  theta_and_D ##########

sub	theta_B{
	my ( $ss1, $ss2, $ll, $nn1, $nn2, $length ) = @_ ;
	# @$ss1: array of seq1
	# @$ss2: array of seq2
	# $$ll: # of seg sites
	# $$nn1: # of sample1
	# $$nn2: # of sample2
	# $$length: length of simulated region
	
	# combi
	my ($combiB) = $$nn1*$$nn2 ;
	my ($combiW1) = ($$nn1*($$nn1-1))/2 ;
	my ($combiW2) = ($$nn2*($$nn2-1))/2 ;
#print "$combiB, $combiW1, $combiW2\n" ;
	
	my ($piB) = 0 ;
	my ($piW1) = 0 ;
	my ($piW2) = 0 ;
	
	my ($ns1, $ns2) ;
	my (@ntB) ;
	
	for ( my ($itB)=0; $itB < $$ll; $itB++ ){
		# count the number of derived alleles 1
		$ns1 = 0 ;
		@ntB = map {substr ($_, $itB, 1)} @$ss1 ;
		#print "@ntB\n" ;
		for ( my ($jtB)=0; $jtB < $$nn1; $jtB++ ){
			if ( $ntB[$jtB] == 1 ){
				$ns1++ ;
			}
			elsif ( $ntB[$jtB] == 0 ){; }
			else {print "Error 1 in theta_B, $ntB[$jtB]\n" ; exit ; }
		}
		
		# count the number of derived alleles 2
		$ns2 = 0 ;
		@ntB = map {substr ($_, $itB, 1)} @$ss2 ;
		#print "@ntB\n" ;
		for ( my ($jtB)=0; $jtB < $$nn2; $jtB++ ){
			if ( $ntB[$jtB] == 1 ){
				$ns2++ ;
			}
			elsif ( $ntB[$jtB] == 0 ){; }
			else {print "Error 2 in theta_B, $ntB[$jtB]\n" ; exit ; }
		}
		#print "$ns1, $ns2\n" ;
		
		$piB += ($$nn1-$ns1)*$ns2 + $ns1*($$nn2-$ns2) ;
		$piW1 += ($$nn1-$ns1)*$ns1 ;
		$piW2 += ($$nn2-$ns2)*$ns2 ;
		#last ;
	}

#print "$piB, $piW1, $piW2\n" ;
	$piB  /= $combiB ;
	$piW1 /= $combiW1 ;
	$piW2 /= $combiW2 ;
#print "$piB, $piW1, $piW2\n" ;
	$piB  /= $$length ;
	$piW1 /= $$length ;
	$piW2 /= $$length ;
#print "$piB, $piW1, $piW2\n" ;	

	# calculate Fstp	
	my ($Hw) = ($piW1+$piW2)/2 ;
	my ($Fstp) = 0 ;
	if ( (0.5*$Hw+0.5*$piB) > 0 ){
		$Fstp = 1 - $Hw/(0.5*$Hw+0.5*$piB) ;
	}
	else{
		$Fstp = -9999 ;
	}
#print "$Fstp\n" ;
	#exit ;
	
	# Fst based on Haplotype
	my ($HHw1) = 0 ;
	my ($HHw2) = 0 ;
	for ( my ($ii1)=0; $ii1 < $$nn1-1; $ii1++ ){
		for ( my ($ii2)=($ii1+1); $ii2 < $$nn1; $ii2++ ){
			if ( $$ss1[$ii1] ne $$ss1[$ii2] ){
				$HHw1++ ;
			}
		}
	}
$HHw1 /= $combiW1 ;
#print "$HHw1\n" ;
	for ( my ($ii1)=0; $ii1 < $$nn2-1; $ii1++ ){
		for ( my ($ii2)=($ii1+1); $ii2 < $$nn2; $ii2++ ){
			if ( $$ss2[$ii1] ne $$ss2[$ii2] ){
				$HHw2++ ;
			}
		}
	}
$HHw2 /= $combiW2 ;
#print "$HHw2\n" ;

	my ($HHb) = 0 ;
	for ( my ($ii1)=0; $ii1 < $$nn1; $ii1++ ){
		for ( my ($ii2)=0; $ii2 < $$nn2; $ii2++ ){
			if ( $$ss1[$ii1] ne $$ss2[$ii2] ){
				$HHb++ ;
			}
		}
	}
	$HHb /= $combiB ;
#print "$HHb\n" ;	
	
	my ($Fsth) = 0 ;
	my ($HHw) = ($HHw1+$HHw2)/2 ;
	if ( (0.5*$HHw+0.5*$HHb) > 0 ){
		$Fsth = 1 - $HHw/(0.5*$HHw+0.5*$HHb) ;
	}
	else{
		$Fsth = -9999 ;
	}
	
	my (%res) = () ;
	$res{"piW1"} = $piW1  ;
	$res{"piW2"} = $piW2  ;	
	$res{"piB"}  = $piB ;
	$res{"Fstp"} = $Fstp ;
	$res{"Fsth"} = $Fsth ;
	$res{"Deltapi"} = $piW2 - $piW1 ;
	
	return %res ;
	


}

=pod


######  between_pi #####
#
# calculate pi within a population
# arguments are: 
# 	(1) length of the region
# 	(2) string of allelic classes for each sample
# 	(3) Hash of sequences
# returns the value of "pi" per site

sub	between_pi{
	my ($len, $state, %data ) = @_;

	my @key = sort keys %data;
	my $totalss = scalar @key;	# total number of samples
	my @states = split //, $state;	# array of allelic class for each sample

	# split sequences into each sites
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	# count the number of samples in an ancestral allelic class
	my ($ssA, $ssD)=(0,0);
	foreach my $j ( 0..$#key ){
		if( $states[$j] eq 'a' ){ $ssA++; } 
	}
	# set the number of samples in the derived allelic class
	$ssD = $totalss-$ssA;

	# calculate between pi
	my $pi=0;
	if( $ssA*$ssD==0 ){
		$pi = -$len;
#print " ssA=$ssA or ssD=$ssD\n";
		return "NA";
	}else{
		# calculate pi for each segregating site
		foreach my $i ( 0..$#{$locusA{$key[0]}} ){	

			my ($ssA0, $ssD0)=(0,0);
			foreach my $j ( 0..$#key ){

				# count the number of sites with ancestral state
				# in each allelic class
				if( ${$locusA{$key[$j]}}[$i] eq '0' ){ 
					if( $states[$j] eq 'a' ){ 
						$ssA0++; 
					}else{
						$ssD0++; 
					}
				}
			}

			print " ssA=$ssA0 or ssD=$ssD0\n";

			# calculate pi
			$pi += ($ssA0*($ssD-$ssD0)+($ssA-$ssA0)*$ssD0)/($ssA*$ssD);
		}
		return ($pi/$len);
	}
}


######  between_dup_pi #####
#
# calculate pi between duplicates
# arguments are: 
# 	(1) length of the region
# 	(2) string of allelic classes for each sample
# 	(3) Hash of sequences
# returns the value of "pi" per site

sub	between_dup_pi{
	my ( $len, %data ) = @_;

	my @key = sort keys %data;
	my $ss = scalar @key;	# total number of samples

	# split sequences into each sites
	my %locusA = ();
	my %locusD = ();
	foreach my $i ( @key ){ 
		@{$locusA{$i}} = split //, ${$data{$i}}[0]; 
		@{$locusD{$i}} = split //, ${$data{$i}}[1]; 
	}
		
	# calculate pi between duplicates
	my $pi=0;

	foreach my $i ( 0..$#key ){
		foreach my $j ( 0..$#key ){

			# skip if genes are on the same chromosome
			if( $i eq $j ){ next; }

			my $diff=0;

			foreach my $k ( 0..$#{$locusA{$key[0]}} ){	

				# count the difference
				if( ${$locusA{$key[$i]}}[$k] ne ${$locusD{$key[$j]}}[$k] ){ 
					$diff++; 
				}
			}

#			print "i=$i\tj=$j\tnum of diff=$diff\n";

			$pi += $diff/$len;

		} # end  foreach my $j
	} # end  foreach my $i

	$pi /= ($ss*($ss-1));

	return ( $pi );
}

##########  Fst  ##########
#
# defined by Wrigt (1985) 
# calculated following Hudson, Slatkin and Maddison (1992)
#
# arguments are: 
# 	(1) length of the region
# 	(2) Hash of sequences
# returns the value of Fst 

sub	Fst{
	my ( $len, $state, %data ) = @_;

	my @key = sort keys %data;
	my $totalss = scalar @key;	# total number of samples
	my @states = split //, $state;	# array of allelic class for each sample
	my $ss = scalar @states;

	# make a list of sub-population
	my %stat;
	foreach my $st ( @states ){
		unless( exists $stat{$st} ){
			$stat{$st}=1;
		}
	}

#foreach my $st ( @states ){
#	print "$st ";
#}
#print "\n";
#
#foreach my $k ( keys %stat ){
#	print "$k ";
#}
#print "\n";

	# calculate divergence within populations
	my $Hw=0;
	my $combination=0;
	# for each population
	foreach my $k ( keys %stat ){

		# redefine sample size and data set
		my %subdata;
		my $sss=0;
		foreach my $i ( 0..$#states ){
			if( $states[$i] eq $k ){
				$subdata{$i} = $data{$i};
				$sss++;
			}
		}
		# calculate total diversity within a population
		$Hw += &theta_pi($len, %subdata)*$len*$sss*($sss-1)/2;
		$combination += $sss*($sss-1)/2;
	}
	$Hw /= $combination;

#print "Hw = $Hw\n";

#foreach my $k ( keys %stat ){
#	print "$k ";
#}
#print "\n";

	# calculate diversity between populations
	my $Hb=0;
	$combination=0;
	my @subpops = sort keys %stat;
	# for every combination of subpopulations
	foreach my $i ( 0..((scalar keys %stat)-2) ){
		foreach my $j ( $i+1..((scalar keys %stat)-1) ){

			# redefine data set
			my %subdata;
			my ($ss1, $ss2)=(0,0);
			foreach my $k ( 0..$#states ){
				if( $states[$k] eq $subpops[$i] || $states[$k] eq $subpops[$j] ){
					$subdata{$k} = $data{$k};
					if ( $states[$k] eq $subpops[$i] ){ $ss1++; }
					else{ $ss2++; }
				}
			}
#print "ss1=$ss1, ss2=$ss2\n";
			# calculate total number of difference between two subpopulations
			$Hb += &between_pi($len, $state, %subdata)*$len*$ss1*$ss2;
			$combination += $ss1*$ss2;
		}
	}
	# take average
#print "Hb=$Hb\n";
	if ( $combination<1 ){ return "NA"; } # if ss1==0 or ss2==0, then skip.
	$Hb /= $combination;
#print "Hb=$Hb\n";

	# calculate Fst
	my $Fst=0 ;
	if ( (0.5*$Hw+0.5*$Hb)>0 ){
		$Fst = 1-$Hw/(0.5*$Hw+0.5*$Hb) ;
	}
	else{
		$Fst = -999 ;
	}
 
	return $Fst ;
}


##########  Gst  ##########
#
# defined by Nei (1973) 
# calculated following Hudson, Boos and Kaplan (1992)
#



##########  Fu and Li's D* (with outgroup)  ##########
#
# calculate Fu and Li's D (1993) 
# see also Simonsen, Churchill and Aquadro (1995) 
# because the original paper contains typo in the equation.
#
# arguments are: 
# 	(1) length of the region
# 	(2) Hash of sequences; '0' is ancestral and '1' is derived
# returns the value of Fu and Li's D

sub	FuLiD_out{
	my ( $len, %data ) = @_;

	my @key = sort keys %data;
	my $ss = scalar @key;

	if ( $ss<2 ){ return "NA"; } # added, Oct 14, 2010

	# split strings of sequence into arrays
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	# count the number of segregating sites (S) per region
	my $S = 0;
	foreach my $i ( 0..$#{$locusA{$key[0]}} ){	

		# count the number of derived SNPs
		my $nsnp=0;
		foreach my $j ( 0..$#key ){	
			if( ${$locusA{$key[$j]}}[$i] eq '1' ){
				$nsnp++;
			}
		}
		# increment the number of segregating sites 
		# if the site is polymorphic
		unless( $nsnp==0 || $nsnp==$ss ){ $S++; }
	}
	if( $S==0 ){ return "NA"; }

	# count the number of singleton (eta)
	my $etas = 0;
	foreach my $i ( 0..$#{$locusA{$key[0]}} ){	

		# count the number of derived SNPs
		my $nsnp=0;
		foreach my $j ( 0..$#key ){	
			if( ${$locusA{$key[$j]}}[$i] eq '1' ){
				$nsnp++;
			}
		}
		# increment the number of segregating sites 
		# if the site is polymorphic
		if( $nsnp == 1 ){ $etas++; }
	}

	# calculate "an", "bn" and "cn"
	my ($an, $bn, $cn) = (0,0,0);
	foreach my $i (1..($ss-1)){
		$an += 1/$i;
		$bn += 1/($i*$i);
	}
	$cn = ($ss==2) ? 1: 2*($ss*$an-2*($ss-1))/($ss-1)/($ss-2);

	# define "ud" and "vd"
	my ($ud, $vd) = (0,0);
	$vd = 1 + $an*$an/($bn+$an*$an)*($cn-($ss+1)/($ss-1));
	$ud = $an-1-$vd;

	# calculate Fu and Li's D*
#	my $FuLiD = $S - $an*($S-$etas);
	my $FuLiD = $S - $an*$etas;
	$FuLiD /= sqrt( $ud*$S + $vd*$S*$S );
 
	return $FuLiD;
}



##########  Rm  ##########
#
# calculate minimum number of recombination 
# by Hudson and Kaplan (1985) 
#
# arguments are: 
# 	(1) length of the region
# 	(2) Hash of sequences
# returns the value of Tajima's D

sub	Rm{
	my ( $len, %data ) = @_;

	my @key = sort keys %data;
	my $ss = scalar @key;

	if ( $ss<2 ){ return "NA"; } # added Oct 14, 2010

	# split strings of sequence into arrays
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	my %fourGamMat;
	# Four gamete test
	foreach my $i ( 0..($#{$locusA{$key[0]}}-1) ){	
		foreach my $j ( ($i+1)..$#{$locusA{$key[0]}} ){	

			# count the number of each haplotype
			my ($f00, $f01, $f10, $f11) = (0, 0, 0, 0);
#			my ($nsnpA, $nsnpB) = (0, 0);
			foreach my $k ( 0..$#key ){	
				if( ${$locusA{$key[$k]}}[$i] eq '0' ){
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f00++;
					}else{
						$f01++;
					}
				}else{
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f10++;
					}else{
						$f11++;
					}
				}
			}
#print "($i,$j)\t($f00,$f01,$f10,$f11)\t";
#print ($f00*$f10*$f01*$f11);
#print "\n";
			if( $f00*$f10*$f01*$f11==0 ){
				${$fourGamMat{$i}}{$j}=0;
			}else{
				${$fourGamMat{$i}}{$j}=1;
			}
		}
	}

	# count the number of four-gamete-site-pairs
	my $Rm=0;
	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
			if( ${$fourGamMat{$i}}{$j}==1 ){
				$Rm++;
#				print "(", $i, ",", $j, ") ";
			}
		}
	}
#	print "\n";
#	print "number of site pairs with four gamete = ", $Rm, "\n";

#	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
#		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
#			print "(", $i, ",", $j, ") ";
#		}
#		print "\n";
#	}
#	print "\n";
#
#	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
#		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
#			print ${$fourGamMat{$i}}{$j}, " ";
#		}
#		print "\n";
#	}
#	print "\n";


	# remove four-gamete-site-pairs which completely include other four-gamete-site-pair
	foreach my $m ( sort {$a<=>$b} keys %fourGamMat ){
		foreach my $n ( sort {$a<=>$b} keys %{$fourGamMat{$m}} ){
			next unless( ${$fourGamMat{$m}}{$n}==1 );

			my $inclded=0;
			foreach my $i ( $m..($n-1) ){
				foreach my $j ( ($i+1)..$n ){

					next if ( $m==$i && $n==$j );
					if( ${$fourGamMat{$i}}{$j}==1 ){
						$inclded=1;
						last;
					}
				}
				last if( $inclded==1 );
			}

			if( $inclded==1 ){
				${$fourGamMat{$m}}{$n}=0;
			}
		}
	}

	# count the number of four-gamete-site-pairs
	$Rm=0;
	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
			if( ${$fourGamMat{$i}}{$j}==1 ){ 
				$Rm++; 
#				print "(", $i, ",", $j, ") ";
			}
		}
	}
#	print "\n";
#	print "number of site pairs with four gamete = ", $Rm, "\n";

	# remove four-gamete-pairs whose start point is included in the other four-gamete-pair
	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
			next unless( ${$fourGamMat{$i}}{$j}==1 );

			foreach my $m ( ($i+1)..($j-1) ){
				foreach my $n ( sort {$a<=>$b} keys %{$fourGamMat{$m}} ){
					next if( $n<$j );
					
					if( ${$fourGamMat{$m}}{$n}==1 ){
						${$fourGamMat{$m}}{$n}=0;
					}
				}
			}
		}
	}

	# count the number of four-gamete-site-pairs
	$Rm=0;
	foreach my $i ( sort {$a<=>$b} keys %fourGamMat ){
		foreach my $j ( sort {$a<=>$b} keys %{$fourGamMat{$i}} ){
			if( ${$fourGamMat{$i}}{$j}==1 ){ 
				$Rm++; 
#				print "(", $i, ",", $j, ") ";
			}
		}
	}
#	print "\n";
#	print "number of site pairs with four gamete = ", $Rm, "\n";
 
	return $Rm;
}

##########  Lewontin_D  ##########
#
# calculate Lewontin's D' (1964) 
# arguments are: 
# 	(1) length of the region
# 	(2) a string of SNP positions: eg. 0_1_2_5_...
# 	(3) Hash of sequences
# returns the value of D'

sub	Lewontin_D{
	my ( $len, $posit_ref, $hashdata_ref ) = @_;

	# split a string of positions into an array
	my @posit = @$posit_ref;
	#foreach my $i ( 0..(@$posit_ref-1) ){
	#	push @posit, @$posit_ref[$i];
	#}
	#split /_/, @$positions;

	my %data = %$hashdata_ref;
	my @key = sort keys %data;
	my $ss = scalar @key;

	# split strings of sequence into arrays
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	# count the number of segregating sites (S) per region
	my %listLD;
	foreach my $i ( 0..($#{$locusA{$key[0]}}-1) ){	
		foreach my $j ( ($i+1)..$#{$locusA{$key[0]}} ){	

			# calculate freqs of four gametes
			my ($f00, $f01, $f10, $f11) = (0, 0, 0, 0);
			my ($nsnpA, $nsnpB) = (0, 0);
			foreach my $k ( 0..$#key ){	
				if( ${$locusA{$key[$k]}}[$i] eq '0' ){
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f00++;
					}else{
						$f01++;
					}
				}else{
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f10++;
					}else{
						$f11++;
					}
				}
			}

			# count number -> frequency
			$f00 /= $ss;
			$f01 /= $ss;
			$f10 /= $ss;
			$f11 /= $ss;

			# calc freq of alleles in each locus
			my ($fa0, $fb0) = ($f00+$f01, $f00+$f10);
			
			# calc LD
			my $LD = $f00*$f11-$f01*$f10;

			my $LDmax;
			if($LD < 0){
				$LDmax = ($fa0*$fb0 < (1-$fa0)*(1-$fb0)) ? $fa0*$fb0 : (1-$fa0)*(1-$fb0);
			}else{
				$LDmax = ($fa0*(1-$fb0) < (1-$fa0)*$fb0) ? $fa0*(1-$fb0) : (1-$fa0)*$fb0;
			}

			# if one of two site is monomorphic, LDmax=0
			next if( $LDmax==0 );

			$LD /= $LDmax;

			# take absolute value of LD
			if($LD < 0){
				$LD*=-1;
			}

			# save result in a hash
			my $dist = $posit[$j]-$posit[$i];
			if( exists $listLD{$dist} ){
				push @{$listLD{$dist}}, $LD;
			}else{
				@{$listLD{$dist}}=();
				push @{$listLD{$dist}}, $LD;
			}
		}
	}
	return \%listLD;
}

##########  r^2  ##########
#
# calculate Hill and Robertson's (1968) 
# arguments are: 
# 	(1) length of the region
# 	(2) a string of SNP positions: eg. 0_1_2_5_...
# 	(3) Hash of sequences
# returns the value of D'

sub	r2{
	my ( $len, $posit_ref, $hashdata_ref ) = @_;

	# an array of positions
	my @posit = @$posit_ref;

	# a hash of sequence data
	my %data = %$hashdata_ref;
	my @key = sort keys %data;
	my $ss = scalar @key;
	# split strings of sequence into arrays
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	# count the number of segregating sites (S) per region
	my %listLD;
	foreach my $i ( 0..($#{$locusA{$key[0]}}-1) ){	
		foreach my $j ( ($i+1)..$#{$locusA{$key[0]}} ){	

			# calculate freqs of four gametes
			my ($f00, $f01, $f10, $f11) = (0, 0, 0, 0);
			my ($nsnpA, $nsnpB) = (0, 0);
			foreach my $k ( 0..$#key ){	
				if( ${$locusA{$key[$k]}}[$i] eq '0' ){
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f00++;
					}else{
						$f01++;
					}
				}else{
					if( ${$locusA{$key[$k]}}[$j] eq '0' ){
						$f10++;
					}else{
						$f11++;
					}
				}
			}

			# count number -> frequency
			$f00 /= $ss;
			$f01 /= $ss;
			$f10 /= $ss;
			$f11 /= $ss;

			# calc freq of alleles in each locus
			my ($fa0, $fb0) = ($f00+$f01, $f00+$f10);

			next if( $fa0==1 || $fa0==0 );
			next if( $fb0==1 || $fb0==0 );
			
			# calc LD
			my $LD = ($f00*$f11-$f01*$f10)*($f00*$f11-$f01*$f10);

			$LD /= ($fa0*(1-$fa0)*$fb0*(1-$fb0));

			# save result in a hash
			my $dist = $posit[$j]-$posit[$i];
			if( exists $listLD{$dist} ){
				push @{$listLD{$dist}}, $LD;
			}else{
				@{$listLD{$dist}}=();
				push @{$listLD{$dist}}, $LD;
			}
		}
	}
	return \%listLD;
}


##########  Dsum  ##########
#
# calculate Dsum defined by Innan (2003) 
# arguments are: 
# 	(1) reference of Hash of sequences
# returns the value of Dsum

sub	Dsum{
	my ( $hashdata_ref ) = @_;

	# a hash of sequence data
	my %data = %$hashdata_ref;
	my @key = sort keys %data;
	my $ss = scalar @key;

	# split strings of sequence into arrays
	my %locA = ();
	my %locB = ();
	foreach my $i ( @key ){
		@{$locA{$i}} = split //, ${$data{$i}}[0]; 
		@{$locB{$i}} = split //, ${$data{$i}}[1]; 
	}

	my $Dsum=0;
	# for each segregating site
	foreach my $i ( 0..$#{$locA{$key[0]}} ){

		my ($f00, $f01, $f10, $f11) = (0, 0, 0, 0);
		# for each chromosome
		foreach my $j ( 0..$#key ){	

			if( ${$locA{$key[$j]}}[$i] eq '0' ){
				if( ${$locB{$key[$j]}}[$i] eq '0' ){
					$f00++;
				}else{
					$f01++;
				}
			}else{
				if( ${$locB{$key[$j]}}[$i] eq '0' ){
					$f10++;
				}else{
					$f11++;
				}
			}
		}
		my $Dm = ($f00*$f11-$f01*$f10)/($ss*($ss-1));
#		print $Dm, "\n";
		$Dsum += $Dm;
	}
	return $Dsum;
}


##########  SFS of SNPs  ##########
#
# calculate SFS
# arguments are: 
# 	(1) reference of Hash of sequences
# returns the reference to an array of SFS 

sub	SFS{
	my ( $hashdata_ref ) = @_;

	# a hash of sequence data
	my %data = %$hashdata_ref;
	my @key = sort keys %data;
	my $ss = scalar @key;

	# split strings of sequence into arrays
	my %locA = ();
	foreach my $i ( @key ){
		@{$locA{$i}} = split //, $data{$i}; 
	}

	# initialize variable;
	my @sfs;
	#foreach my $i ( 0..$ss ){ $sfs[$i] = 0; }
	foreach my $i ( 0..100 ){ $sfs[$i] = 0; }

	# for each segregating site
	foreach my $i ( 0..$#{$locA{$key[0]}} ){

		my $der=0;
		# for each chromosome
		foreach my $j ( 0..$#key ){	

			# count the number of derived allele at this SNP site
			if( ${$locA{$key[$j]}}[$i] eq '1' ){
				$der++;
			}
		}
	
		# if there is no derived allele, skip.
		if( $der==0 ){ next; }
		# if erived allele is fixed, skip.
		if( $der==$ss ){ next; }

		# change count into %
		$sfs[int($der*100/$ss)]++;
	}
	return \@sfs;
}

##########  SFS of SNPs in duplicates  ##########
#
# calculate SFS
# arguments are: 
# 	(1) reference of Hash of sequences
# returns the array of reference to arrays 

sub	SFS_dup{
	my ( $hashdata_ref ) = @_;

	# a hash of sequence data
	my %data = %$hashdata_ref;
	my @key = sort keys %data;
	my $ss = scalar @key;

	# split strings of sequence into arrays
	my %locA = ();
	my %locB = ();
	foreach my $i ( @key ){
		@{$locA{$i}} = split //, ${$data{$i}}[0]; 
		@{$locB{$i}} = split //, ${$data{$i}}[1]; 
	}

	# initialize variable;
	my $fixed=0;
	my @specificA;
	my @specificB;
	my @sharedA;
	my @sharedB;
	my @sharedAll;
	foreach my $i ( 0..$ss ){
		$specificA[$i] = 0;
		$specificB[$i] = 0;
		$sharedA[$i] = 0;
		$sharedB[$i] = 0;
	}
	foreach my $i ( 0..2*$ss ){
		$sharedAll[$i] = 0;
	}

	# for each segregating site
	foreach my $i ( 0..$#{$locA{$key[0]}} ){

		my ($derA, $derB)=(0,0);
		# for each chromosome
		foreach my $j ( 0..$#key ){	

			if( ${$locA{$key[$j]}}[$i] eq '1' ){
				$derA++;
			}
			if( ${$locB{$key[$j]}}[$i] eq '1' ){
				$derB++;
			}
		}

print "$derA\t$derB\n";

		if( ($derA==$ss && $derB==0) || ($derA==0 && $derB==$ss) ){
			$fixed++;

		}elsif( ($derB==0 && 0<$derA ) ){
			$specificA[$derA]++;

		}elsif( ($derA==0 && 0<$derB ) ){
			$specificB[$derB]++;

		}elsif( (0<$derA && 0<$derB ) ){
			$sharedA[$derA]++;
			$sharedB[$derB]++;
			$sharedAll[$derA+$derB]++;
		}
	}
	return ( $fixed, \@specificA, \@specificB, \@sharedA, \@sharedB, \@sharedAll );
}


##########  nHaplo ##########
#
# calculate the number of haplotypes within a population
# arguments are: 
# 	(1) length of the region
# 	(2) Hash of sequences
# returns the number of haplotypes

sub	nHaplo{
	my ($len, %data ) = @_;

	# sample size 
	my @key = sort keys %data;
	my $ss = scalar @key;
	if ( $ss<2 ){ return ("NA", "NA"); }

	# create hash of haplotypes
	my %haplotype;
	my %freq;
	my @hkey=();
	my $nhaplo = scalar @hkey; 

	# classify haplotypes
	foreach my $id ( @key ){ # all sequences

		if( $nhaplo==0 ){		# first sample
			$haplotype{$id} = $data{$id};
			$freq{$id}=1;
			push @hkey, $id;
			$nhaplo++;

		}else{				# the next and after
			foreach my $hid ( @hkey ){	

				# haplotype already exists
				if( $data{$id} eq $haplotype{$hid} ){
					$freq{$hid}++;

				# new haplotype
				}else{
					# add new haplotype 
					$haplotype{$id} = $data{$id};
					$freq{$id}=1;

					# update the number and list of the types
					push @hkey, $id;
					$nhaplo++;

					last;
				}
			}

		}
	}
#print $nhaplo, "\n";

	#
	# calculate heterozygosity
	#

	# homozygosity
	my $homo=0;
	foreach my $hid ( @hkey ){
		$homo += ($freq{$hid}/$ss)*($freq{$hid}/$ss);
	}

#print $homo, "\t";
	# unbiased heterozygosity
	my $hetero=0;
	$hetero = (1.-$homo)*$ss/($ss-1);		
#print $hetero, "\n";

	return ( $nhaplo, $hetero);
}

=cut

1 ;
#####################
## end of the file ##
#####################
