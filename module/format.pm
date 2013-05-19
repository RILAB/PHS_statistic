###############################################
# Module for changing FORMAT to UNIX and Perl #
###############################################

# linend
# fasta
# fileopen
# fileopenTab
# TwoDarrayTest
# fopen_hash_fas
# fopen_hash_two
# fopen_hash_multi

use encoding "utf8";

# Module to change line end format to UNIX(\n) from Windows(\r\n) and Macintosh(\r). 

sub linend{ 

	my ($directory) = @_;                # Catch the directory and file type (e.g. Users/stakumo/perl/*pl) 

	my (@file_name) = glob($directory);  # get file names
	my ($dir)='';
	my ($tmp) = '';

	foreach $dir(@file_name){

		open (OUT, "$dir")||die "error in sub-routin linend\n";
		my (@contents) = <OUT>;
		close (OUT);
		
		open (OUT, ">$dir")||die "error in sub-routin linend\n";
		
		
		foreach $tmp(@contents){
			$tmp =~ s/\r\n/\n/;  #change line end format, Windows to UNIX
			$tmp =~ s/\r/\n/g;   #change line end format, Macintosh to UNIX
			print OUT $tmp;
		}
		close (OUT);
	
	}

}


# Module for editing fasta format to use in perl

sub fasta{ 
	use strict;
	my ($directory) = @_;                # Catch the directory and file type (e.g. Users/stakumo/perl/*.fasta) 

	my (@file_name) = glob($directory);  # get file names
	my (@precode) = ();
	my ($dir) = '';
	my ($tmp) = '';
	my ($temp) = '';
	my ($a) = -1;
	my ($b) = 0;

	# Edit fasta format
	
	foreach $dir(@file_name){
	#print "$dir\n";
		open (OUT, "$dir")|| die "error in sub-routin fasta\n";
		my (@contents) = <OUT>;
		close (OUT);

		$a=-1;
		$b=0;
		@precode = ();
		$precode[0]='';
		
		foreach $tmp(@contents){
	
			if($b==0){
				$a=$a;
			}
			elsif($b==1){
				$a=$a+1;
			}
	
			if($tmp =~ />/){
				$a=$a+1;
				chomp $tmp;
				$precode[$a]=$tmp;
				$b=1;
			}
			else{
				chomp $tmp;
				$tmp =~ s/\s//;
				$precode[$a] = $precode[$a].$tmp;
				$b=0;
			}
		}

		open (OUT, ">$dir")||die "error in sub-routin fasta\n";
		foreach $temp(@precode){
			print OUT $temp,"\n";
		}
		close (OUT);
	}
}


# open file
sub fileopen {
	use strict;
	my ($hogeF) = @_;
	open (HOG, "$hogeF")|| die "error sub fileopen, $hogeF.\n";
	my (@fooF) = <HOG>;
	close (HOG);
	chomp @fooF;
	
	return @fooF;

}


# open file ver. 2 
sub fileopen2 {
	use strict;
	my ($hogeAR, $hogeF) = @_ ;
	open (HOG, "$$hogeF")|| die "Error sub fileopen2, $$hogeF.\n" ;
	@$hogeAR = <HOG> ;
	close (HOG) ;
	chomp @$hogeAR ;

}

# open file ver. 3
sub fileopen3 {
	use strict;
	my ($hogeAR, $hogeNN, $hogeFF) = @_ ;
	open (HOG, "$$hogeFF")|| die "Error sub fileopen3, $$hogeFF.\n" ;
	@$hogeAR = <HOG> ;
	close (HOG) ;
	chomp @$hogeAR ;
	$$hogeNN = scalar @$hogeAR ;

}

# open file as 2-dimensions array
sub fileopenTab {
	use strict;
	my ($hogeFT) = @_;
	open (HOG, "$hogeFT")|| die "error sub fileopenTab, $hogeFT.\n";
	my (@fooFT) = <HOG>;
	close (HOG);
	chomp @fooFT;
	my ($n_fooFT) = scalar @fooFT;
	
	# convert 2-dimensions array
	my (@fooFT2) = ();
	for (my ($ij)=0; $ij < $n_fooFT; $ij++){
		@{$fooFT2[$ij]} = split ("\t", $fooFT[$ij]);
	}
	
	return @fooFT2;

}

# open file as 2-dimensions array, ver. 2
sub fileopenTab2 {
	use strict;
	my ($hogeTB, $adTB) = @_;
	open (HOG, "$$adTB")|| die "error sub fileopenTab2, $$adTB.\n";
	
	my ($openfTB) ;
	my ($dTB) = 0 ;
	while ( $openfTB = <HOG> ){
		chomp $openfTB ;
	
		@{$$hogeTB[$dTB]} = split ("\t", $openfTB) ;

		$dTB++ ;
	} # while
	
	close (HOG) ;
	
}

# open file as 2-dimensions array, ver. 3
sub fileopenTab3 {
	use strict;
	my ($hogeTB, $nn_dt, $cc_dt, $adTB) = @_;
	open (HOG, "$$adTB")|| die "error sub fileopenTab2, $$adTB.\n";
	
	my ($openfTB) ;
	my ($dTB) = 0 ;
	while ( $openfTB = <HOG> ){
		chomp $openfTB ;
	
		@{$$hogeTB[$dTB]} = split ("\t", $openfTB) ;

		$dTB++ ;
	} # while
	
	close (HOG) ;
	
	$$nn_dt = scalar @$hogeTB ;
	$$cc_dt = scalar @{$$hogeTB[0]} ;
	
}

# test column numbers in 2-dimension array 
sub TwoDarrayTest{
	use strict;
	my (@hoge2d) = @_;
	my ($n_hoge2d) = scalar @hoge2d;
	my ($c_hoge2d) = scalar @{$hoge2d[0]};
	#print "$n_hoge, $c_hoge\n";
	for (my ($ik)=1; $ik < $n_hoge2d; $ik++){
		my ($ccc2d) = scalar @{$hoge2d[$ik]};
		if ($c_hoge2d != $ccc2d){
			print "Error in $ik th row.\n";
			exit;
		}
	}
	
}

# test column numbers in 2-dimension array ver.2
sub TwoDarrayTest2{
	use strict;
	my ($hoge2d) = @_;
	my ($n_hoge2d) = scalar @$hoge2d;
	my ($c_hoge2d) = scalar @{$$hoge2d[0]};
	#print "$n_hoge, $c_hoge\n";
	for (my ($ik)=1; $ik < $n_hoge2d; $ik++){
		my ($ccc2d) = scalar @{$$hoge2d[$ik]};
		if ($c_hoge2d != $ccc2d){
			print "Error in sub TwoDarrayTest2\n" ;
			print "Error in $ik th row.\n";
			exit;
		}
	}
	
}

# open fasta seq and return hash
sub fopen_hash_fas {
	use strict ;
	my ($hogeFHF) = @_ ;
	my (%hogeFHF) = () ;
	
	open (DAT, "$hogeFHF")|| die "Error in sub fopen_hash_fas.\nFile not found ($hogeFHF).\n" ;
	
	my ($openFHF) ;
	my ($dFHF) = 0 ;
	my ($idFHF) ;

	while ( $openFHF = <DAT> ){
		chomp $openFHF ;
		
		if ( $dFHF == 0 ){
			$idFHF = $openFHF ;
			$idFHF =~ s/^\>// ;
			$dFHF = 1 ;
		}
		elsif ( $dFHF == 1 ){
			$hogeFHF{$idFHF} = $openFHF ;
			$dFHF = 0 ;
		}
		else {
			print "Error in sub fopen_hash_fas.\n" ;
			print "Unexpected number, $dFHF\n" ;
			exit ;
		}
	}
	
	return ( %hogeFHF) ;
	
}

# open fasta seq and return hash
sub fopen_hash_fas2 {
	use strict ;
	my ($hogeFHF, $addD) = @_ ;
	
	open (DAT, "$$addD")|| die "Error in sub fopen_hash_fas2.\nFile not found ($hogeFHF).\n" ;
	
	my ($openFHF) ;
	my ($dFHF) = 0 ;
	my ($idFHF) ;

	while ( $openFHF = <DAT> ){
		chomp $openFHF ;
		
		if ( $dFHF == 0 ){
			$idFHF = $openFHF ;
			$idFHF =~ s/^\>// ;
			$dFHF = 1 ;
		}
		elsif ( $dFHF == 1 ){
			$$hogeFHF{$idFHF} = $openFHF ;
			$dFHF = 0 ;
		}
		else {
			print "Error in sub fopen_hash_fas.\n" ;
			print "Unexpected number, $dFHF\n" ;
			exit ;
		}
	}
	
	#return ( %hogeFHF) ;
	
}

# open Two tab-column file and return hash
sub fopen_hash_two {
	use strict ;
	my ($hogeFHT) = @_ ;
	my (%hogeFHT) = () ;
	
	open (DAT, "$hogeFHT")|| die "Error in sub fopen_hash_two.\nFile not found ($hogeFHT).\n" ;
	
	my ($openFHT) ;

	while ( $openFHT = <DAT> ){
		chomp $openFHT ;
		my (@tmpFHT) = split ("\t", $openFHT) ;
		my ($n_tmpFHT) = scalar @tmpFHT ;
		if ( $n_tmpFHT != 2 ){
			print "Error in sub fopen_hash_two.\n" ;
			print "Number of column is not 2, $n_tmpFHT\n" ;
			exit ;
		}
		
		$hogeFHT{$tmpFHT[0]} = $tmpFHT[1] ;
	}
	
	return ( %hogeFHT) ;
	
}


# open Multi tab-column file and return hash
sub fopen_hash_multi {
	use strict ;
	my ($hogeFHM) = @_ ;
	my (%hogeFHM) = () ;
	
	open (DAT, "$hogeFHM")|| die "Error in sub fopen_hash_multi.\nFile not found ($hogeFHM).\n" ;
	
	my ($openFHM) ;

	while ( $openFHM = <DAT> ){
		chomp $openFHM ;
		my (@tmpFHM) = split ("\t", $openFHM) ;
		my ($n_tmpFHM) = scalar @tmpFHM ;
		
		my ($keyFHM) = $tmpFHM[0] ;
		shift @tmpFHM ;
		my ($valFHM) = join ("\t", @tmpFHM) ;
		
		$hogeFHM{$keyFHM} = $valFHM ;
	}
	
	return ( %hogeFHM) ;
	
}

# open Multi tab-column file and return hash
sub fopen_hash_multi2 {
	use strict ;
	my ($hogeFHM, $added) = @_ ;
	#my (%hogeFHM) = () ;
	
	open (DAT, "$$added")|| die "Error in sub fopen_hash_multi2.\nFile not found ($$added).\n" ;
	
	my ($openFHM) ;

	while ( $openFHM = <DAT> ){
		chomp $openFHM ;
		my (@tmpFHM) = split ("\t", $openFHM) ;
		my ($n_tmpFHM) = scalar @tmpFHM ;
		
		my ($keyFHM) = $tmpFHM[0] ;
		shift @tmpFHM ;
		my ($valFHM) = join ("\t", @tmpFHM) ;
		
		$$hogeFHM{$keyFHM} = $valFHM ;
	}
	
	#return ( %hogeFHM) ;
	
}

# convert array into 2-D array
sub ConvertTwoDarray {
	my ($moto, $kekka) = @_ ;
	my ($N_moto) = scalar @$moto ;
	
	for ( my ($iCTD)=0; $iCTD < $N_moto; $iCTD++ ){
		@{$$kekka[$iCTD]} = split ("\t", $$moto[$iCTD]) ;
	
	}

}

1;