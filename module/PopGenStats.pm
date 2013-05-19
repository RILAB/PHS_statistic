#
# $Id: PopGenStats.pm,v 1.1.1.1 2009/11/27 02:13:59 kteshima Exp $
#
package PopGenStats;

use strict;
use warnings;

##########  theta_and_D ##########
#
# theta	: pi within a population
# 	: Watterson's theta (1975) 
# 	: Fay and Wu's theta (2000) 
# 	: Zeng, Fu, Shi and Wu's theta_L (2006) 
# Test of neutrality 
#	: Tajima's D (1989) 
# 	: Fu and Li's D (1993) (without outgroup)
# 	 	see also Simonsen, Churchill and Aquadro (1995) 
# 	: Fay and Wu's H (2000) (un-normalized)
# 	: Fay and Wu's H (2000) ( normalized )
# 		variance of H is obtained by Zeng, Fu, Shi and Wu (2006)
# 
# arguments are: 
# 	(1) length of the region
# 	(2) Hash of sequences
# 
##########  theta_and_D ##########

sub	theta_and_D{
	my ( $len, %data ) = @_;

	my @key = sort keys %data;
	my $ss = scalar @key;

	#if ( $ss<2 ){ return "NA"; }

	# split sequences into sites
	my %locusA = ();
	foreach my $i ( @key ){
		@{$locusA{$i}} = split //, $data{$i}; 
	}

	##############
	# theta_pi 
	# theta_w: Watterson
	# theta_H: Fay and Wu
	# theta_L: Zeng, Fu, Shi and Wu's
	#############
	my $pi=0;
	my $nsegl=0;
	my $thetaH=0;
	my $thetaL=0;
	my $etas=0;
	my $combination = $ss*($ss-1)/2;
	# for each segregating site
	foreach my $i ( 0..$#{$locusA{$key[0]}} ){	

		my $nsnp=0;
		# count the number of derived SNP 
		foreach my $j ( 0..$#key ){	# sample

			if( ${$locusA{$key[$j]}}[$i] eq '1' ){
				$nsnp++;
			}
		}
		unless( $nsnp==0 || $nsnp==$ss ){ 
			$pi += ($ss-$nsnp)*$nsnp;# pi
			$nsegl++; 		# W
			$thetaH += $nsnp*$nsnp;	# H
			$thetaL += $nsnp;	# L
		}
		if( $nsnp==1 || $nsnp==($ss-1) ){ $etas++; }
	}
	$pi /= $combination;
	$thetaH /= $combination;
	$thetaL /= $combination;

	# calculate "a1" and "a2"
	my ($a1, $a2) = (0,0);
	foreach my $i (1..($ss-1)){
		$a1 += 1/$i;
		$a2 += 1/($i*$i);
	}
	# define "b1" and "b2"
	my ($b1, $b2) = (0,0);
	$b1 = ($ss+1)/(3*($ss-1));
	$b2 = 2*($ss*$ss+$ss+3)/(9*$ss*($ss-1));
	# define "c1" and "c2"
	my ($c1, $c2) = (0,0);
	$c1 = $b1-1/$a1;
	$c2 = $b2-($ss+2)/($a1*$ss)+$a2/($a1*$a1);
	# define "e1" and "e2"
	my ($e1, $e2) = (0,0);
	$e1 = $c1/$a1;
	$e2 = $c2/($a1*$a1+$a2);
	# define "ud" and "vd"
	my ($ud, $vd) = (0,0);
	$vd = $a2/$a1/$a1 - 2*(1+1/$a1-$a1+$a1/$ss)/$ss - 1/$ss/$ss;
	$vd /= ($a1*$a1+$a2);
	$ud = (($ss-1)/$ss - 1/$a1)/$a1 - $vd;
	# define "bn"
	my $bn=0;
	foreach my $i (1..$ss){
		$bn += 1/($i*$i);
	}

#print "sample size is $ss, n-segl is $nsegl.\n";
#print "a1,a2 is ($a1,$a2)\n";
#print "b1,b2 is ($b1,$b2)\n";
#print "c1,c2 is ($c1,$c2)\n";
#print "e1,e2 is ($e1,$e2)\n";
#print "vd,ud is ($vd,$ud)\n";

	# Tajima's D
	my $tajimaD = $pi - $nsegl/$a1;
	if( $nsegl<2 || ($e1*$nsegl + $e2*$nsegl*($nsegl-1))==0 ){ $tajimaD = "NA";}
	else{ $tajimaD /= sqrt( $e1*$nsegl + $e2*$nsegl*($nsegl-1) );}

	# Fu and Li's D
	my $FuLiD = $nsegl/$a1 - $etas*($ss-1)/$ss;
	if( $ud*$nsegl + $vd*$nsegl*$nsegl > 0 ){ 
		$FuLiD /= sqrt( $ud*$nsegl + $vd*$nsegl*$nsegl ); }
	else{ $FuLiD = "NA"; }

	# Fay and Wu's H (un-normalized) 
	my $faywuH = $pi - $thetaH;
	# normalized H
	my $theta2=$nsegl*($nsegl-1)/($a1*$a1+$a2);
	my $normH = $pi - $thetaL;
	my $varH = 
		($ss-2)*$pi/6/($ss-1) 
		+ ((18*$ss*$ss*(3*$ss+2)*$bn-(88*$ss*$ss*$ss+9*$ss*$ss-13*$ss+6))*$theta2)
			/9/$ss/($ss-1)/($ss-1);
	if( $varH==0 ){ $normH = "NA"; }
	else{ $normH /= sqrt($varH); }

#	return ($pi/$len);		# theta_pi: per site
#	return ($nsegl/$len/$a1);	# theta_w: per site
#	return ($thetaH/$len);		# theta_H: per site
#	return ($thetaL/$len);		# theta_L: per site
#	$tajimaD;
#	$FuLiD;
#	$faywuH;
#	$normH;

#	my $res = ($pi/$len);
	my %res=();
	$res{"pi"}= $pi/$len;		# theta_pi: per site
	$res{"thetaW"}= $nsegl/$len/$a1;	# theta_w: per site
	$res{"thetaH"}= $thetaH/$len;		# theta_H: per site
	$res{"thetaL"}= $thetaL/$len;		# theta_L: per site
	$res{"tajimaD"}= $tajimaD;	# Tajima's D
	$res{"fuliD"}= $FuLiD;		# Fu and Li's D
	$res{"faywuH"}= $faywuH;	# Fay and Wu's H
	$res{"normalH"}= $normH;	# normalized Fay and Wu's H

#	foreach my $k ( keys %res ){
#		print "$k\t$res{$k}\n";
#	}
	return( %res );
}




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
	my $Fst=0;
	if( (0.5*$Hw+0.5*$Hb)>0 ){
		$Fst = 1-$Hw/(0.5*$Hw+0.5*$Hb);
	}else{
		$Fst = -999;
	}
 
	return $Fst;
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

1;
#####################
## end of the file ##
#####################
