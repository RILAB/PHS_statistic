#!/usr/bin/env perl
use strict ;
use lib "./module" ;
#use lib "/home/shohei/module" ;
use stat ;
use format ;
use BeginPerlBioinfo ;
use evolve ;
use PopGenStatsT ;

# calulate P-values of PHS at candidate SNPs
# This program calcuates P-values given the frequency of 
# focused SNPs.  If the size of null distribution is small,
# please modify this script.



my (@nl) = () ;
my ($openf) ;
open (DAT, "PHS_null.out")|| die "ppp\n" ;
while ( $openf = <DAT> ){
	chomp $openf ;
	my (@tm) = split ("\t", $openf) ;
	my ($tm) = "$tm[1]\t$tm[2]" ;
	push (@nl, $tm) ;
}
close (DAT) ;

my (@dt1) = () ;
&fileopenTab2 (\@dt1, \ "PHS_candidate.out") ;
my ($n_dt1) = scalar @dt1 ;

open (OUT, ">PHS_Pval.out")|| die "ppp\n" ;

for ( my ($i)=0; $i < $n_dt1; $i++ ){
	my (@tmp) = grep (/^$dt1[$i][4]\t/, @nl) ;
	@tmp = &columnTab (\@tmp, \1) ;
	my ($n_tmp) = scalar @tmp ;
	my ($num) = scalar grep {$_ >= $dt1[$i][5]} @tmp ;
	$num /= $n_tmp ;
	
	my ($res) = join ("\t", @{$dt1[$i]}) ;
	print OUT "$res\t$num\n" ;

}

close (OUT) ;

exit ;

sub get_par {
	my ($comm, $n_nn, $n_aa, $n_dd, $rateu, $rater, $simreg, $n_rpp) = @_ ;
	
	my (@com) = split (/\s/, $$comm) ;
	my ($n_com) = scalar @com ;
	#print "# of com = $n_com\n" ;
	#print "$$comm\n" ;
	
	for ( my ($gg)=0; $gg < $n_com; $gg++ ){
		if  ( $gg == 0 && $com[$gg] ne "./ms" ){print "Error 1 in get par, $com[$gg]\n" ; exit ; }
		elsif ( $gg == 1 ){$$n_nn = $com[$gg] ; $$n_rpp = $com[++$gg] } # # of samples
		if ( $gg > 2 ){
			if    ( $com[$gg] eq "-t" ){$$rateu  = $com[++$gg] ; }
			elsif ( $com[$gg] eq "-r" ){$$rater  = $com[++$gg] ; $$simreg = $com[++$gg] ; }
			elsif ( $com[$gg] eq "-I" ){$gg++ ; $$n_aa = $com[++$gg] ; $$n_dd = $com[++$gg] ; }
		}
	}
	
}


sub round {
	my ($hoge) = @_;
	my (@foo) = split ("\\.", $hoge);
	my ($n_foo) = scalar @foo;
	if ($n_foo == 1){
		return $hoge;
	}
	elsif ($n_foo == 2){
		my ($bar) = substr ($foo[1], 0, 1);
		if ($bar <= 4){
			return $foo[0];
		}
		elsif ($bar >= 5){
			return $foo[0]+1;
		}
		else {
			print "Error1 in round, $hoge, @foo\n"; exit;
		}
	}
	else {
		print "Error2 in round, $hoge, @foo\n"; exit;
	}
	
}

sub round2 {
	my ($hoge) = @_;
	my (@foo) = split ("\\.", $hoge);
	my ($n_foo) = scalar @foo;
	if ($n_foo == 1){
		#print "Caution @foo, $n_foo\n";
		return $hoge;
	}
	elsif ($n_foo == 2){
		return $foo[0];
	}
	else {
		print "Error2 in round, $hoge, @foo\n"; exit;
	}
	
}

