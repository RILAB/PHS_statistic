###################################
# Module for computing statistics #
###################################

# List of sub
# sum
# mean
# sd
# cov
# cc
# max
# min
# permutation
# columnTab
# union
# column
# columnSeq
# position
# log10


# calculating sum of sample
sub sum {
	use strict;
	my (@hogeSUM) = @_;	# argument of sample
	my ($resSUM) = 0;
	foreach my $tmpSUM (@hogeSUM){
		$resSUM += $tmpSUM;
	}
	
	return $resSUM;
}

# calculating sample means.
sub mean {
	use strict;
	my (@fooME) = @_;           # argument of sample
	my ($numME) = scalar @fooME;  # number of sample
	my ($sigmaME) = 0;          # sum of sample
	
	foreach my $tmpME (@fooME){
		$sigmaME += $tmpME;
	}
	
	return $sigmaME / $numME;
	
}# End "mean"


#Sub-routin for calculating unbiased standard deviation.
sub sd {
	use strict;
	my (@fooSD) = @_;            # argument of sample
	my ($numSD) = scalar @fooSD;   # number of sample
	my ($meanSD) = &mean (@fooSD); # mean of sample
	my ($sd) = 0;              #standard deviation of @foo
	
	for (my ($iSD)=0; $iSD < $numSD; $iSD++){
		$sd += ($fooSD[$iSD]-$meanSD) ** 2;
	}
	$sd = sqrt ($sd/($numSD-1));
	return $sd;
	
} # End "sd"


# Sub-routin for calculating Co-Variance between two samples.
sub cov {
	use strict;	
	my ($fooCV1, $fooCV2) = @_; # argument of two sample
	my ($numCV1) = scalar @$fooCV1; # number of each sample
	my ($numCV2) = scalar @$fooCV2; 
	# check whether sample of each sample are different or not.
	unless ($numCV1 == $numCV2){
		print "Error in sub cov\nArguments have different number of arrays\n"; exit;	
	}
	
	#calculating Co-Variance between @$foo1 & @$foo2.
	my ($cov) = 0; #sum of square between @$foo1 & @$foo2
	my ($meanCV1) = &mean (@$fooCV1);
	my ($meanCV2) = &mean (@$fooCV2);
	
	for (my ($iCV)=0; $iCV < $numCV1; $iCV++){
		$cov += ($$fooCV1[$iCV]-$meanCV1) * ($$fooCV2[$iCV]-$meanCV2);
	}
	
	return $cov / ($numCV1-1);
} # End "cov"


# subroutin for estimation correlation coefficient.
sub cc {
	use strict;
	my($fooCC1,$fooCC2) = @_;
	my ($numCC1) = scalar @$fooCC1; # number of each sample
	my ($numCC2) = scalar @$fooCC2; 
	# check whether sample of each sample are different or not.
	unless ($numCC1 == $numCC2){
		print "Error in sub cc\nArguments have different number of arrays\n"; exit;	
	}
	
	my ($covCC) = &cov (\@$fooCC1, \@$fooCC2);
	my ($sdCC1) = &sd (@$fooCC1);
	my ($sdCC2) = &sd (@$fooCC2);
	
	#calcuration and return of correlation coefficient

	return $covCC / ($sdCC1 * $sdCC2);
}

# determin a max value in array
sub max {
	use strict;
	my (@hogeMA) = @_;
	my ($maxMA) = $hogeMA[0];
	foreach (@hogeMA){
		if ($maxMA < $_){
			$maxMA = $_;
		}
	}
	return $maxMA;
}

# determin a min value in array
sub min{
	use strict;
	my (@hogeMI) = @_;
	my ($minMI) = $hogeMI[0];
	foreach (@hogeMI){
		if ($minMI > $_){
			$minMI = $_;
		}
	}
	return $minMI;
}

# Permutation an array
#sub permutation_old {
#	use strict;
#	my (@hogePM) = @_;
#	my ($n_hogePM) = scalar @hogePM;
#	my (@outPM) = ();
	
	#srand (time|$$);
	
#	for (my ($iPM)=0; $iPM < $n_hogePM; $iPM++){
#		my ($nPM) = scalar @hogePM;
#		my ($ranPM) = int (rand ($nPM));
#		my ($tmpPM) = $hogePM[$ranPM];
#		$hogePM[$ranPM] = $hogePM[0];
#		shift @hogePM;
#		push (@outPM, $tmpPM);
#	}
	
#	return @outPM;
#}

# Permutation an array
sub permutation {
	use strict ;
	my (@hogePM) = @_ ;
	my ($n_hogePM) = scalar @hogePM ;
	my (@outPM) = () ;
		
	for ( my ($iPM)=0; $iPM < $n_hogePM; $iPM++ ){
		my ($ranPM) = int (rand ($#hogePM+1)) ;
		push (@outPM, $hogePM[$ranPM]) ;
		splice @hogePM, $ranPM, 1 ;
	}
	
	return @outPM ;
}

# extract Column from an tab array.
# Two argues are needed
# When using, @res = &columnTab (\@tmp, \2);
sub columnTab {
	use strict;
	my ($hogeCT, $fooCT) = @_;
	my ($n_hogeCT) = scalar @$hogeCT;
	#print "$n_hoge\n";
	
	my (@outCT) = ();
	for (my ($il)=0; $il < $n_hogeCT; $il++){
		my (@tmpCT) = split ("\t", $$hogeCT[$il]);
		push (@outCT, $tmpCT[$$fooCT]);
	}

	return @outCT;
}

# Union an array
sub union {
	use strict;
	my (@hogeUN) = @_;
	my (%seen);
	@hogeUN = grep (!$seen{$_}++, @hogeUN);

	return @hogeUN;
}

# extract Column from 2-dimension array.
# Two argues are needed
# When using, @res = &column (\@tmp, \2);
sub column {
	use strict;
	my ($hogeCC, $fooCC) = @_;
	my ($n_hogeCC) = scalar @$hogeCC;
	my ($c_hogeCC) = scalar @{@$hogeCC[0]};
	#print "$n_hoge, $c_hoge\n";
	#print "$n_hoge, $$hoge[10][0], $$hoge[10][8]\n";
	
	if ($$fooCC >= $c_hogeCC){print "Argument is bigger than the number of columns in sub column.\n"; exit; }
	
	my (@outCC) = ();
	for (my ($il)=0; $il < $n_hogeCC; $il++){
		push (@outCC, $$hogeCC[$il][$$fooCC]);
	}

	return @outCC;
}

# extract Column from Sequences.
# Two argues are needed
# When using, @res = &columnTab (\@tmp, \2);
sub columnSeq {
	use strict;
	my ($hoge, $foo) = @_;
	#print "$n_hoge\n";
	
	my (@outS) = map {substr ($_, $$foo, 1)} @$hoge;
	
	return @outS;
}


# Search a given element in an array.
# Two argues are needed
# Note that this sub could not be used for numbers.
# When using, @res = &position (\@tmp, \$rrr);
sub position {
	use strict;
	my ($hogePO, $fooPO) = @_;
	my ($tmpPO);
	my ($ccPO) = 0;
	foreach $tmpPO (@$hogePO){
		if ($tmpPO eq $$fooPO){
			return $ccPO;
		}
		$ccPO++;
	}
	
	# not fount
	return -1;
	#print "Error, $$fooPO is not found\n";
	#exit;
}


# calculate common logarithm
sub log10{
    my ($xLG) = @_ ;
    return log ( $xLG ) / log ( 10 ) ;
}

sub log2{
	my ($x) = @_ ;
	return log ( $x ) / log ( 2 ) ;
}

sub median{
	return unless @_ ;
	return $_[0] unless @_ > 1 ;
	@_= sort {$a<=>$b} @_ ;
	return $_[$#_/2] if @_&1 ;
	my $mid= @_/2 ;
	return ($_[$mid-1]+$_[$mid])/2 ;
}


1;