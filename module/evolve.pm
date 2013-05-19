####################################
# Module for evolutionary analysis #
####################################

# List of sub
# cutgap
# codonSyn
# codonSub
# codonNG1
# codonNG2
# codonNG3
# kaks1
# jc69
# lensynsites

sub get_align_pep {
	use strict ;
	my (@pepAC) = @_ ;
	my ($n_pepAC) = scalar @pepAC ;
	#print "# of cds = $n_cdsAC\n" ;
	
	# align pep
	open (TMPAC, ">pep.fas")|| die "error tmp.fas\n" ;
	for ( my ($iAC)=0; $iAC < $n_pepAC; $iAC++ ){
		print TMPAC ">$iAC\n$pepAC[$iAC]\n" ;
	}
	close (TMPAC) ;
	system ("clustalw -align -INFILE=pep.fas -TYPE=PROTEIN -OUTFILE=pep.fas -OUTPUT=GDE -OUTORDER=INPUT -CASE=UPPER > a.out") ;
	
	my (@tmpAC) = () ;
	&fileopen2 (\@tmpAC, \ "pep.fas") ;
	my ($n_tmpAC) = scalar @tmpAC ;
	@pepAC = ("") x $n_pepAC ;
	my ($fAC1) = 0 ;
	for ( my ($iAC)=1; $iAC < $n_tmpAC; $iAC++ ){
		if ( $tmpAC[$iAC] =~ /\%/ ){
			$fAC1++ ;
		}
		else {
			$pepAC[$fAC1] = $pepAC[$fAC1].$tmpAC[$iAC] ;
		}
	}
	
	system ("rm a.out") ;
	system ("rm pep.fas") ;
	system ("rm pep.dnd") ;
	
	return @pepAC ;
}


sub get_align_DNA {
	my ($foo1, $foo2) = @_ ;
	$foo1 = uc $foo1 ;
	$foo2 = uc $foo2 ;
	open (TMP, ">tmp2.fas")|| die "iijjnn\n" ;
	print TMP ">a\n$foo1\n" ;
	print TMP ">b\n$foo2\n" ;
	close (TMP) ;
	system ("clustalw -align -INFILE=tmp2.fas -TYPE=DNA -OUTFILE=tmp.fas -OUTPUT=GDE -OUTORDER=INPUT -CASE=UPPER > a.out") ;
	
	my (@tmp) = &fileopen ("tmp.fas") ;
	my ($n_tmp) = scalar @tmp ;
	my (@res) = ("", "") ;
	my ($f1) = 0 ;
	foreach (@tmp){
		if ( $_ eq "#a" ){
			$f1 = 1 ;
		}
		elsif ( $_ eq "#b" ){
			$f1 = 2 ;
		}
		elsif ( $f1 == 1 ){
			$res[0] = $res[0].$_ ;
		}
		elsif ( $f1 == 2 ){
			$res[1] = $res[1].$_ ;			
		}
		else {
			print "Error in align_dna\n" ; exit ;
		}
	}
	if ( length $res[0] != length $res[1] ){
		print "Error2 in align_dna\n" ; exit; 
	}
	
	system ("rm a.out") ;
	system ("rm tmp.fas") ;
	system ("rm tmp2.fas") ;
	system ("rm tmp2.dnd") ;
	
	
	return @res ;
}


sub get_align_cds {
	use strict ;
	my (@cdsAC) = @_ ;
	my ($n_cdsAC) = scalar @cdsAC ;
	#print "# of cds = $n_cdsAC\n" ;
	
	# check length %3
	foreach (@cdsAC){
		if ( (length $_)%3 != 0 ){print "Error in sub get_align_cds, %3\n" ; }
	}
	
	# keep stop codon (if avaiable)
	my (@stopAC) = () ;
	foreach (@cdsAC){
		my ($lenAC) = length $_ ;
		my ($tmpAC) = substr ($_, $lenAC-3, 3) ;
		if ( &check_stop_codon (\$tmpAC) ){
			$_ = substr ($_, 0, $lenAC-3) ;
		}
		else {
			$tmpAC = "---" ;
		}
		push (@stopAC, $tmpAC) ;
	}
	#print "@stopAC\n" ;
	
	# translate
	my (@pepAC) = () ;
	foreach (@cdsAC){
		my ($tmpAC) = &mycodon ($_) ;
		push (@pepAC, $tmpAC) ;
	}
	
	# align pep
	open (TMPAC, ">pep.fas")|| die "error tmp.fas\n" ;
	for ( my ($iAC)=0; $iAC < $n_cdsAC; $iAC++ ){
		print TMPAC ">$iAC\n$pepAC[$iAC]\n" ;
	}
	close (TMPAC) ;
	system ("clustalw -align -INFILE=pep.fas -TYPE=PROTEIN -OUTFILE=pep.fas -OUTPUT=GDE -OUTORDER=INPUT -CASE=UPPER > a.out") ;
	
	my (@tmpAC) = () ;
	&fileopen2 (\@tmpAC, \"pep.fas") ;
	my ($n_tmpAC) = scalar @tmpAC ;
	@pepAC = ("") x $n_cdsAC ;
	my ($fAC1) = 0 ;
	for ( my ($iAC)=1; $iAC < $n_tmpAC; $iAC++ ){
		if ( $tmpAC[$iAC] =~ /\%/ ){
			$fAC1++ ;
		}
		else {
			$pepAC[$fAC1] = $pepAC[$fAC1].$tmpAC[$iAC] ;
		}
	}
	
	# align cds
	my (@resAC) = ("") x $n_cdsAC ;
	for ( my ($kAC)=0; $kAC < $n_cdsAC; $kAC++){
		# Divided for array
		my (@tmppAC) = split ("", $pepAC[$kAC]) ;
		my (@tmpcAC) = split ("", $cdsAC[$kAC]) ;
		
		# insert gap -> cds
		my (@outcAC) = () ; my ($jAC) = -1 ;
		foreach (@tmppAC){
			if ( $_ eq "-" ){
				push (@outcAC, "-"); push (@outcAC, "-"); push (@outcAC, "-");
			}
			else {
				push (@outcAC, $tmpcAC[++$jAC]); push (@outcAC, $tmpcAC[++$jAC]); push (@outcAC, $tmpcAC[++$jAC]);
			}
		}
		$resAC[$kAC] = join ("", @outcAC) ;
	}
	
	unless ( scalar grep (/\-\-\-/, @stopAC) == $n_cdsAC ){
		for ( my ($iAC)=0; $iAC < $n_cdsAC; $iAC++ ){
			$resAC[$iAC] = $resAC[$iAC].$stopAC[$iAC] ;
		}
	}
	
	#open (TMP, ">tmp.fas")|| die "pp\n" ;
	#print TMP "\>0\n$resAC[0]\n\>1\n$resAC[1]\n" ;
	#close (TMP) ;
	
	system ("rm a.out") ;
	system ("rm pep.fas") ;
	system ("rm pep.dnd") ;
	
	return @resAC ;
}

# check stop codon if stop, 1 return
sub check_stop_codon {
	use strict ;
	my ($hogeCSC) = @_ ;
	if ( length $$hogeCSC != 3 ){print "Error in check_stop_codon, $$hogeCSC\n"; }
	if ( $$hogeCSC eq "TAA" || $$hogeCSC eq "TAG" || $$hogeCSC eq "TGA" ){
		return 1 ;
	}
	else {
		return 0 ;
	}
}

# cut codons with gap and missing data 
sub cutNcod {
	use strict ;
	my (@hogeCNC) = @_ ;
	my ($n_hogeCNC) = scalar @hogeCNC ;
	
	my (@resCNC) = ("") x $n_hogeCNC ;
	
	my ($lenCNC) = length $hogeCNC[0] ;
	if ( $lenCNC%3 != 0 ){print "Error len3 in sub cutNcod, $lenCNC\n"; exit; }
	
	for ( my ($iCN)=0; $iCN < $lenCNC; $iCN+=3 ){
		my (@codCNC) = map {substr ($_, $iCN, 3)} @hogeCNC ;
		my (@tmpCNC) = grep (/^[ATGC][ATGC][ATGC]$/, @codCNC) ;
		my ($n_tmpCNC) = scalar @tmpCNC ;
		if ( $n_tmpCNC == $n_hogeCNC ){
			for ( my ($jCN)=0; $jCN < $n_hogeCNC; $jCN++ ){
				$resCNC[$jCN] = $resCNC[$jCN].$codCNC[$jCN] ;
			}
		}
	}
	
	return @resCNC ;

}

# cut gapped site from alignment of seq
# Note that argue can contain only seq data
sub cutgap {
	use strict;
	my (@hogeC) = @_;
	my ($n_hogeC) = scalar @hogeC;
	my (@hogeC2) = ();
	
	for (my ($iC)=0; $iC < $n_hogeC; $iC++){
		@{$hogeC2[$iC]} = split ("", $hogeC[$iC]);
	}
	my ($n_hogeC2) = scalar @hogeC2;
	my ($c_hogeC2) = scalar @{$hogeC2[0]};
	&TwoDarrayTest (@hogeC2);
	
	for (my ($iC)=0; $iC < $c_hogeC2; $iC++){
		my (@tmpC) = &column (\@hogeC2, \$iC);
		my (@ggC) = grep (/-/, @tmpC);
		my ($n_ggC) = scalar @ggC;
		if ($n_ggC > 0){
			for (my ($jC)=0; $jC < $n_hogeC2; $jC++){
				$hogeC2[$jC][$iC] = "-";
			}
		}
	}
	
	my (@resC) = ();
	my ($ttC);
	for (my ($iC)=0; $iC < $n_hogeC2; $iC++){
		$ttC = join ("", @{$hogeC2[$iC]});
		$ttC =~ s/\-//g;
		push (@resC, $ttC);
	}
	
	my ($lenC) = length $resC[0];
	for (my ($iC)=0; $iC < $n_hogeC2; $iC++){
		my ($lenTC) = length $resC[$iC];
		if ($lenC != $lenTC){print "Error Length in cutgap\n"; exit; }
	}
	
	return @resC;
	
}

# cut codon containing missing data
# Note that argue can contain only seq data
sub cutmiss {
	use strict;
	my (@hogeM) = @_;
	my ($n_hogeM) = scalar @hogeM;
	my (@hogeM2) = ();
	
	for (my ($iM)=0; $iM < $n_hogeM; $iM++){
		@{$hogeM2[$iM]} = split ("", $hogeM[$iM]);
	}
	my ($n_hogeM2) = scalar @hogeM2;
	my ($c_hogeM2) = scalar @{$hogeM2[0]};
	&TwoDarrayTest (@hogeM2);
	if ($c_hogeM2%3 != 0){print "Error in sub cutmiss, length of seq is not mupliplied by 3!!\n"; exit; }
	
	for (my ($iM)=0; $iM < $c_hogeM2; $iM+=3){
		my ($tmpM1) = scalar grep (/[RSNYWMBKVHrsnywmbkvh]/, &column (\@hogeM2, \$iM));
		my ($uuM) = $iM+1;
		my ($tmpM2) = scalar grep (/[RSNYWMBKVHrsnywmbkvh]/, &column (\@hogeM2, \$uuM));
		my ($uuM) = $iM+2;
		my ($tmpM3) = scalar grep (/[RSNYWMBKVHrsnywmbkvh]/, &column (\@hogeM2, \$uuM));
		
		if ($tmpM1+$tmpM2+$tmpM3 > 0){
			for (my ($jM)=0; $jM < $n_hogeM2; $jM++){
				$hogeM2[$jM][$iM]   = "N";
				$hogeM2[$jM][$iM+1] = "N";
				$hogeM2[$jM][$iM+2] = "N";
			}
		}
	}
	
	my (@resM) = ();
	my ($ttM);
	for (my ($iM)=0; $iM < $n_hogeM2; $iM++){
		$ttM = join ("", @{$hogeM2[$iM]});
		$ttM =~ s/[Nn]//g;
		push (@resM, $ttM);
	}
	
	my ($lenM) = length $resM[0];
	for (my ($iM)=0; $iM < $n_hogeM2; $iM++){
		my ($lenTM) = length $resM[$iM];
		if ($lenM != $lenTM){print "Error Length in cutmiss\n"; exit; }
	}
	
	return @resM;
	
}

# the number of potentially synonymous changes for codon
sub codonSyn {
	use strict;
	my($codonS) = @_;
	$codonS = uc $codonS;
	
	my(%genetic_codeS) = (
		'TCA' => 1.0,	# Serine
		'TCC' => 1.0,	# Serine
		'TCG' => 1.0,	# Serine
		'TCT' => 1.0,	# Serine
		'TTC' => 1/3,	# Phenylalanine
		'TTT' => 1/3,	# Phenylalanine
		'TTA' => 2/3,	# Leucine
		'TTG' => 2/3,	# Leucine
		'TAC' => 1.0,	# Tyrosine
		'TAT' => 1.0,	# Tyrosine
		#'TAA' => '*',    # Stop
		#'TAG' => '*',    # Stop
		'TGC' => 0.5,	# Cysteine
		'TGT' => 0.5,	# Cysteine
		#'TGA' => '*',    # Stop
		'TGG' => 0.0,	# Tryptophan
		'CTA' => 4/3,	# Leucine
		'CTC' => 1.0,	# Leucine
		'CTG' => 4/3,	# Leucine
		'CTT' => 1.0,	# Leucine
		'CCA' => 1.0,	# Proline
		'CCC' => 1.0,	# Proline
		'CCG' => 1.0,	# Proline
		'CCT' => 1.0,	# Proline
		'CAC' => 1/3,	# Histidine
		'CAT' => 1/3,	# Histidine
		'CAA' => 1/3,	# Glutamine
		'CAG' => 1/3,	# Glutamine
		'CGA' => 1.5,	# Arginine
		'CGC' => 1.0,	# Arginine
		'CGG' => 4/3,	# Arginine
		'CGT' => 1.0,	# Arginine
		'ATA' => 2/3,	# Isoleucine
		'ATC' => 2/3,	# Isoleucine
		'ATT' => 2/3,	# Isoleucine
		'ATG' => 0.0,	# Methionine
		'ACA' => 1.0,	# Threonine
		'ACC' => 1.0,	# Threonine
		'ACG' => 1.0,	# Threonine
		'ACT' => 1.0,	# Threonine
		'AAC' => 1/3,	# Asparagine
		'AAT' => 1/3,	# Asparagine
		'AAA' => 1/3,	# Lysine
		'AAG' => 1/3,	# Lysine
		'AGC' => 1/3,	# Serine
		'AGT' => 1/3,	# Serine
		'AGA' => 5/6,	# Arginine
		'AGG' => 2/3,	# Arginine
		'GTA' => 1.0,	# Valine
		'GTC' => 1.0,	# Valine
		'GTG' => 1.0,	# Valine
		'GTT' => 1.0,	# Valine
		'GCA' => 1.0,	# Alanine
		'GCC' => 1.0,	# Alanine
		'GCG' => 1.0,	# Alanine
		'GCT' => 1.0,	# Alanine
		'GAC' => 1/3,	# Aspartic Acid
		'GAT' => 1/3,	# Aspartic Acid
		'GAA' => 1/3,	# Glutamic Acid
		'GAG' => 1/3,	# Glutamic Acid
		'GGA' => 1.0,	# Glycine
		'GGC' => 1.0,	# Glycine
		'GGG' => 1.0,	# Glycine
		'GGT' => 1.0,	# Glycine
	);

	if (exists $genetic_codeS{$codonS}){
		return $genetic_codeS{$codonS};
	}
	else {
		print "Error in codonSyn, unknown or stop codon, $codonS\n"; exit;
	}
}

# number of total substituion in codon
sub codonSub {
	use strict;
	my (@hogeSU) = @_;
	my ($n_hogeSU) = scalar @hogeSU;
	if ($n_hogeSU != 2){print "Error number of array in codonSub, @hogeSU\n"; exit; }
	
	my ($lenSU1) = length $hogeSU[0];
	my ($lenSU2) = length $hogeSU[1];
	if ($lenSU1 != $lenSU2 || $lenSU1 != 3){print "Error length of codon in codonSub, @hogeSU\n"; exit; }
	
	my ($difSU) = 0;
	for (my ($qw)=0; $qw < $lenSU1; $qw++){
		if ( substr ($hogeSU[0], $qw, 1) ne substr ($hogeSU[1], $qw, 1) ){
			$difSU++;
		}
	}
	
	if ($difSU > 3){print "Error # of substitution > 3 in codonSub, @hogeSU\n"; exit; }
	
	return $difSU;
	
}

sub codonNG1 {
	use strict;
	my (@hogeNG1) = @_;
	my ($n_hogeNG1) = scalar @hogeNG1;
	if ($n_hogeNG1 != 2){print "Error number of array in codonNG1, @hogeNG1\n"; exit; }
	
	my ($pepNG1) = &codon2aa ($hogeNG1[0]);
	my ($pepNG2) = &codon2aa ($hogeNG1[1]);
	
	if ($pepNG1 eq "*" || $pepNG2 eq "*"){print "Stop codon are present in codonNG1, @hogeNG1\n"; exit;}
	if ($pepNG1 ne $pepNG2){
		return (0, 1);
	}
	elsif ($pepNG1 eq $pepNG2){
		return (1, 0);
	}
	else {
		print "Error in codonNG1, $pepNG1, $pepNG2\n"; exit;
	}
	#return $difSU;
	
}

sub codonNG2 {
	use strict;
	my (@hogeNG2) = @_;
	my ($n_hogeNG2) = scalar @hogeNG2;
	if ($n_hogeNG2 != 2){print "Error number of array in codonNG2, @hogeNG2\n"; exit; }
	
	# position of 2 mutaion
	my (@posiNG2) = (-1, -1);
	my ($jjNG2) = 0;
	if (substr ($hogeNG2[0], 0, 1) ne substr ($hogeNG2[1], 0, 1)){$posiNG2[$jjNG2] = 0; $jjNG2++;}
	if (substr ($hogeNG2[0], 1, 1) ne substr ($hogeNG2[1], 1, 1)){$posiNG2[$jjNG2] = 1; $jjNG2++;}
	if (substr ($hogeNG2[0], 2, 1) ne substr ($hogeNG2[1], 2, 1)){$posiNG2[$jjNG2] = 2; $jjNG2++;}
	if ($jjNG2 != 2|| $posiNG2[0] == -1 || $posiNG2[1] == -1){print "No two mutations, @hogeNG2\n"; exit;}
	#print "@posiNG2\n";
	
	# create intermediate codon
	# 1
	my (@tmpNG21) = split ("", $hogeNG2[0]);
	my (@tmpNG22) = @tmpNG21;
	$tmpNG21[$posiNG2[0]] = substr ($hogeNG2[1], $posiNG2[0], 1);
	my ($mmCod1) = join ("", @tmpNG21);
	# 2
	$tmpNG22[$posiNG2[1]] = substr ($hogeNG2[1], $posiNG2[1], 1);
	my ($mmCod2) = join ("", @tmpNG22);
	#print "$mmCod1, $mmCod2\n";
	
	# check stop codon
	my ($stopNG21) = 0;
	my ($stopNG22) = 0;
	if ($mmCod1 eq "TAA" || $mmCod1 eq "TAG" || $mmCod1 eq "TGA"){$stopNG21 = 1; }
	if ($mmCod2 eq "TAA" || $mmCod2 eq "TAG" || $mmCod2 eq "TGA"){$stopNG22 = 1; }
	#print "$stopNG21, $stopNG22\n";
	if ($stopNG21 + $stopNG22 == 2){print "Both pathway include stop codon in codonNG2, @hogeNG2\n"; exit; }
	
	# count syn and rep
	my ($pepNG21) =  &codon2aa ($hogeNG2[0]);
	my ($pepNG23) =  &codon2aa ($hogeNG2[1]);
	my (@difNG2) = (0, 0);
	my ($numPathNG2) = 0;
	# path1
	if ($stopNG21  == 0){
		my ($pepNG22) =  &codon2aa ($mmCod1);
		#print "$pepNG21, $pepNG22, $pepNG23\n";
		if ($pepNG21 eq $pepNG22){$difNG2[0]++; }
		elsif ($pepNG21 ne $pepNG22){$difNG2[1]++; }
		if ($pepNG22 eq $pepNG23){$difNG2[0]++; }
		elsif ($pepNG22 ne $pepNG23){$difNG2[1]++; }
		#print "@difNG2\n";
		$numPathNG2++;
	}
	# path2
	if ($stopNG22  == 0){
		my ($pepNG24) =  &codon2aa ($mmCod2);
		#print "$pepNG21, $pepNG24, $pepNG23\n";
		   if ($pepNG21 eq $pepNG24){$difNG2[0]++; }
		elsif ($pepNG21 ne $pepNG24){$difNG2[1]++; }
		   if ($pepNG24 eq $pepNG23){$difNG2[0]++; }
		elsif ($pepNG24 ne $pepNG23){$difNG2[1]++; }
		#print "@difNG2\n";
		$numPathNG2++;
	}
	@difNG2 = map {$_/$numPathNG2} @difNG2;
	#print "@difNG2\n";
	return @difNG2;

}

sub codonNG3 {
	use strict;
	my (@hogeNG3) = @_;
	my ($n_hogeNG3) = scalar @hogeNG3;
	if ($n_hogeNG3 != 2){print "Error number of array in codonNG3, @hogeNG3\n"; exit; }
	
	# count syn and rep
	my (@tmpNG31) = split ("", $hogeNG3[0]);
	my (@tmpNG32) = split ("", $hogeNG3[1]);
	my ($pepNG31) =  &codon2aa ($hogeNG3[0]);
	my ($pepNG32) =  &codon2aa ($hogeNG3[1]);
	my (@difNG3) = (0, 0);
	my ($numPathNG3) = 0;
	my ($mmCod31);
	my ($mmCod32);
	
	# path1
	$mmCod31 = "$tmpNG32[0]$tmpNG31[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG32[1]$tmpNG31[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path2
	$mmCod31 = "$tmpNG32[0]$tmpNG31[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG31[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path3
	$mmCod31 = "$tmpNG31[0]$tmpNG32[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG32[1]$tmpNG31[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path4
	$mmCod31 = "$tmpNG31[0]$tmpNG32[1]$tmpNG31[2]";
	$mmCod32 = "$tmpNG31[0]$tmpNG32[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path5
	$mmCod31 = "$tmpNG31[0]$tmpNG31[1]$tmpNG32[2]";
	$mmCod32 = "$tmpNG32[0]$tmpNG31[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	# path5
	$mmCod31 = "$tmpNG31[0]$tmpNG31[1]$tmpNG32[2]";
	$mmCod32 = "$tmpNG31[0]$tmpNG32[1]$tmpNG32[2]";
	#print "$mmCod31, $mmCod32\n";
	if ($mmCod31 ne "TAA" && $mmCod31 ne "TAG" && $mmCod31 ne "TGA" && $mmCod32 ne "TAA" && $mmCod32 ne "TAG" && $mmCod32 ne "TGA"){
		my ($mmPep31) = &codon2aa ($mmCod31);
		my ($mmPep32) = &codon2aa ($mmCod32);
		#print "$mmPep31, $mmPep32\n";
		   if ($pepNG31 eq $mmPep31){$difNG3[0]++; }
		elsif ($pepNG31 ne $mmPep31){$difNG3[1]++; }
		   if ($mmPep31 eq $mmPep32){$difNG3[0]++; }
		elsif ($mmPep31 ne $mmPep32){$difNG3[1]++; }
		   if ($mmPep32 eq $pepNG32){$difNG3[0]++; }
		elsif ($mmPep32 ne $pepNG32){$difNG3[1]++; }
		$numPathNG3++;
	}
	#print "@difNG3, $numPathNG3\n";
	@difNG3 = map {$_/$numPathNG3} @difNG3;
	#print "@difNG3, $numPathNG3\n";
	if ($numPathNG3 == 0){print "No pathway for stop codon in codonNG3, @hogeNG3\n"; exit; }
	return @difNG3;
}

# Nei-Gojobiri method
sub kaks1 {
	use strict;
	my (@seqN) = @_;
	@seqN = map uc, @seqN;
	#print "$seqN[0]\n\n$seqN[1]\n";

	# cut gap after checking frameshift mutation
	my ($ltmp1) = length $seqN[0];
	my ($ltmp2) = length $seqN[1];
	if ($ltmp1 != $ltmp2){print "Different length of seq in sub kaks1, $ltmp1, $ltmp2\n"; exit; }
=pod
	if ($ltmp1%3 != 0 || $ltmp2%3 != 0){print "Length of seq is not 3 times in kaks1, $ltmp1, $ltmp2\n"; exit; }
	for (my ($iN)=0; $iN < $ltmp1; $iN+=3){
		my ($tF1) = substr ($seqN[0], $iN, 3);
		my ($tF2) = substr ($seqN[1], $iN, 3);
		if ($tF1 =~ /-/ && $tF1 =~ /[ATGC]/){print "Frameshift of seq1 in sub kaks1, $iN\n"; exit; }
		if ($tF2 =~ /-/ && $tF2 =~ /[ATGC]/){print "Frameshift of seq2 in sub kaks1, $iN\n"; exit; }
	}

	#@seqN = &cutgap (@seqN);
	#print "$seqN[0]\n\n$seqN[1]\n";
=cut

	# cut start
	$ltmp1 = length $seqN[0];
	$ltmp2 = length $seqN[1];
	if ($ltmp1 != $ltmp2){print "Different length of seq in sub kaks1, $ltmp1, $ltmp2\n"; exit; }
	if ($ltmp1%3 != 0 || $ltmp2%3 != 0){print "Length of seq is not 3 times in kaks1, $ltmp1, $ltmp2\n"; exit; }
	my ($sN1) = substr ($seqN[0], 0, 3);
	my ($sN2) = substr ($seqN[1], 0, 3);
	if ($sN1 eq "ATG" || $sN2 eq "ATG"){
		#print "$sN1, $sN2\n";
		$seqN[0] = substr ($seqN[0], 3, $ltmp1-3);
		$seqN[1] = substr ($seqN[1], 3, $ltmp1-3);
	}
	$ltmp1 = length $seqN[0];
	$ltmp2 = length $seqN[1];
	if ($ltmp1 != $ltmp2){print "Different length of seq in sub kaks1 after excluding start codon, $ltmp1, $ltmp2\n"; exit; }
	if ($ltmp1%3 != 0 || $ltmp2%3 != 0){print "Length of seq is not 3 times in kaks1, $ltmp1, $ltmp2\n"; exit; }
	
	# cut last stop codon and check premature stop codon
	# premature

	for (my ($iN)=0; $iN < $ltmp1-3; $iN+=3){
		my ($tN1) = substr ($seqN[0], $iN, 3);
		my ($tN2) = substr ($seqN[1], $iN, 3);
		if ($tN1 eq "TAG" || $tN1 eq "TGA" || $tN1 eq "TAA"){
			print "Premature stop codon of seq1 in kaks1, $iN\n"; exit;
		}
		if ($tN2 eq "TAG" || $tN2 eq "TGA" || $tN2 eq "TAA"){
			print "Premature stop codon of seq2 in kaks1, $iN\n"; exit;
		}
	}


	# last stop codon
	my ($tN1) = substr ($seqN[0], $ltmp1-3, 3);
	my ($tN2) = substr ($seqN[1], $ltmp1-3, 3);
	if ($tN1 eq "TAG" || $tN1 eq "TGA" || $tN1 eq "TAA" || $tN2 eq "TAG" || $tN2 eq "TGA" || $tN2 eq "TAA"){
		$seqN[0] = substr ($seqN[0], 0, $ltmp1-3);
		$seqN[1] = substr ($seqN[1], 0, $ltmp1-3);
	}
	#print "$tN1, $tN2\n";
	#print "$seqN[0]\n\n$seqN[1]\n";

	# check length
	$ltmp1 = length $seqN[0];
	$ltmp2 = length $seqN[1];
	if ($ltmp1 != $ltmp2){print "Different length of seq after excluding start codon in kaks1, $ltmp1, $ltmp2\n"; exit; }
	if ($ltmp1%3 != 0 || $ltmp2%3 != 0){print "Length of seq is not 3 times in kaks1, $ltmp1, $ltmp2\n"; exit; }
	my ($LL) = $ltmp1;
	#print "$ltmp1, $ltmp2\n";
	
	# the number of potentially synonymous changes for codon
	my ($Lsyn1) = 0;
	my ($Lsyn2) = 0;
	for (my ($iN)=0; $iN < $LL; $iN+=3){
		my ($tN1) = substr ($seqN[0], $iN, 3);
		my ($tN2) = substr ($seqN[1], $iN, 3);
		my ($Lstmp) = &codonSyn ("$tN1");
		$Lsyn1 = $Lsyn1 + $Lstmp;
		$Lstmp = &codonSyn ("$tN2");
		$Lsyn2 = $Lsyn2 + $Lstmp;
		#print "$Lsyn1, $Lsyn2\n";
	}
	my ($Lrep1) = $LL - $Lsyn1;
	my ($Lrep2) = $LL - $Lsyn2;
	
	#print "$Lsyn1, $Lrep1, $Lsyn2, $Lrep2\n";
	
	my ($SS) = &mean (($Lsyn1, $Lsyn2));
	my ($NN) = &mean (($Lrep1, $Lrep2));
	#print "$SS, $NN\n";
	
	# count synonymous and nonsynonymous changes
	my ($dSS) = 0;
	my ($dNN) = 0;
	for (my ($iN)=0; $iN < $LL; $iN+=3){
		my ($tN1) = substr ($seqN[0], $iN, 3);
		my ($tN2) = substr ($seqN[1], $iN, 3);
		if ($tN1 ne $tN2){
			my ($ccNN) = &codonSub (($tN1, $tN2));
			my (@difSN);
			# number of substitution = 1
			if ($ccNN == 1){
			   @difSN = &codonNG1 (($tN1, $tN2));
			}
			# 2
			elsif ($ccNN == 2){
				#print "$tN1, $tN2\n";
				@difSN = &codonNG2 (($tN1, $tN2));
			}
			# 3
			elsif ($ccNN == 3){
				#print "$tN1, $tN2\n";
				@difSN = &codonNG3 (($tN1, $tN2));
			}
			else {
				print "Error in kaks1\n"; exit;
			}
			$dSS = $dSS + $difSN[0];
			$dNN = $dNN + $difSN[1];
		}
	}
	
	my ($pSS);
	my ($pNN);
	   if ($SS != 0){$pSS = $dSS/$SS; }
	elsif ($SS == 0){$pSS = 0; }
	   if ($NN != 0){$pNN = $dNN/$NN; }
	elsif ($NN == 0){$pNN = 0; }
	
	#print "$SS, $dSS, $pSS\n$NN, $dNN, $pNN\n";
	# 0 -> # of syn site
	# 1 -> # of syn change
	# 2 -> # of syn dif per site
	# 3 -> # of rep site
	# 4 -> # of rep change
	# 5 -> # of rep dif per site
	my ($resNG) = "$SS\t$dSS\t$pSS\t$NN\t$dNN\t$pNN";
	return $resNG;

}

# Jukes & Cantor (1969) correlation
sub jc69 {
	use strict;
	my ($hogeJC) = @_;
	if ($hogeJC == 0){
		return 0;
	}
	else {
		return -3/4 * log (1-4/3*$hogeJC);
	}
}

sub lensynsites {
	use strict;
	my ($seqLN) = @_;
	$seqLN = uc $seqLN;
	#print "$seqN[0]\n\n$seqN[1]\n";

	# cut gap after checking frameshift mutation
	my ($ltmp11) = length $seqLN;
	if ($ltmp11%3 != 0){print "Length of seq is not 3 times in lensynsites, $ltmp11\n"; exit; }
=pod
	for (my ($iL)=0; $iL < $ltmp11; $iL+=3){
		my ($tF1) = substr ($seqLN, $iL, 3);
		if ($tF1 =~ /-/ && $tF1 =~ /[ATGC]/){print "Frameshift of seq1 in sub lensynsites, $iL\n"; exit; }
	}
	#@seqN = &cutgap (@seqN);
	#print "$seqN[0]\n\n$seqN[1]\n";
=cut
	# cut start
	$ltmp11 = length $seqLN;
	if ($ltmp11%3 != 0){print "Length of seq is not 3 times in lensynsites, $ltmp11\n"; exit; }
	my ($sN1) = substr ($seqLN, 0, 3);
	if ($sN1 eq "ATG"){
		#print "$sN1\n";
		$seqLN = substr ($seqLN, 3, $ltmp11-3);
	}
	$ltmp11 = length $seqLN;
	if ($ltmp11%3 != 0){print "Length of seq is not 3 times in lensynsites, $ltmp11\n"; exit; }
	
	# cut last stop codon and check premature stop codon
=pod
	# premature
	for (my ($iL)=0; $iL < $ltmp11-3; $iL+=3){
		my ($tN1) = substr ($seqLN, $iL, 3);
		if ($tN1 eq "TAG" || $tN1 eq "TGA" || $tN1 eq "TAA"){
			print "Premature stop codon of seq1 in lensynsites, $iL, $tN1\n"; exit;
		}
	}
=cut
	# last stop codon
	my ($tN1) = substr ($seqLN, $ltmp11-3, 3);
	if ($tN1 eq "TAG" || $tN1 eq "TGA" || $tN1 eq "TAA"){
		$seqLN = substr ($seqLN, 0, $ltmp11-3);
	}
	#print "$tN1, $tN2\n";
	#print "$seqN[0]\n\n$seqN[1]\n";

	# check length
	my ($ltmp11) = length $seqLN;
	if ($ltmp11%3 != 0){print "Length of seq is not 3 times in lensynsites, $ltmp11\n"; exit; }
	my ($LL) = $ltmp11;
	#print "$ltmp1, $ltmp2\n";
	
	# the number of potentially synonymous changes for codon
	my ($LsynL1) = 0;
	for (my ($iL)=0; $iL < $LL; $iL+=3){
		my ($tN1) = substr ($seqLN, $iL, 3);
		my ($Lstmp) = &codonSyn ("$tN1");
		$LsynL1 = $LsynL1 + $Lstmp;
		#print "$Lsyn1, $Lsyn2\n";
	}
	#print "$Lsyn1, $Lrep1, $Lsyn2, $Lrep2\n";

	return $LsynL1;
}









1;