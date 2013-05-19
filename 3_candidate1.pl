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
# calculate PHS for candidate SNPs


my (@fl) = glob ("./phased_SNPs/*") ;
my ($n_fl) = scalar @fl ;
print "# of input files = $n_fl\n\n" ;


#my (@tar) = () ;
#&fileopen2 (\@tar, \ "checkMX.out") ;
#my ($n_tar) = scalar @tar ;
#print "# of tar = $n_tar\n" ;

#my ($foc1) = $ARGV[0] ;
#print "$foc1\n" ;

my ($openf) ;

open (OUT, ">PHS_candidate.out")|| die "ppp\n" ;

#my (@at) = qw (1 2 3 4 5 6 7 8 9 10) ;
#my (@at) = qw (9 10) ;

my ($n_sam) = scalar &fileopen ("$fl[0]") ;


foreach my $fl (@fl){

	# SNP
	my (@dt1) = () ;
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
	print " $fl # of individuals = $n_dt1 # of SNPs = $c_dt1\n" ;
	&TwoDarrayTest2 (\@dt1) ;
	if ( $n_sam != $n_dt1 ){print "E2 sam $n_sam $n_dt1\n" ; exit ; }

	my ($of_tar) = $fl ;
	$of_tar =~ s/phased/candidate/ ;
	my (@tar2) = &fileopen ("$of_tar") ;
	my ($n_tar2) = scalar @tar2 ;
	print "  $n_tar2\n" ;

	my ($openf) ;
	my (%av) = () ;
	my (%sd) = () ;
	open (DAT, "d_ij.out")|| die "ppppppp\n" ;
	while ( $openf = <DAT> ){
		chomp $openf ;
		my (@tmp) = split ("\t", $openf) ;
		my ($k) = "$tmp[0]\t$tmp[1]" ;
		$av{$k} = $tmp[2] ;
		$sd{$k} = $tmp[3] ;
	}
	close (DAT) ;

	# cM and physical distance
	my ($of_cm) = $fl ;
	$of_cm =~ s/phased\_SNPs/cM/ ;
	my (@w1) = &fileopen ("$of_cm") ;
	my (@dist) = &columnTab (\@w1, \2) ;
	my ($n_dist) = scalar @dist ;
	my (@mega) = &columnTab (\@w1, \1) ;
	my ($n_mega) = scalar @mega ;
	print " $n_dist $n_mega\n" ;
	if ( $n_dist != $c_dt1 || $n_dist != $n_mega ){print "E1 num $n_dist $c_dt1\n" ; exit ; }

	# calculate PHS
	foreach my $tar (@tar2){
		my (@tar3) = split ("\t", $tar) ;

		# get position
		my ($ppp) = &position (\@mega, \$tar3[1]) ;
		if ( $ppp == -1 ){print "E posi\n" ; exit ; }
		if ( $tar3[1] != $mega[$ppp] ){print "EEE\n" ; exit ; }
		my ($h) = $ppp ;

		my (@nt) = &column (\@dt1, \$h) ;
                my (@al) = &union (@nt) ;
                my ($n_al) = scalar @al ;

		# biallelic
                if ( $n_al == 2 ){
			my ($a1) = scalar grep (/^$al[0]$/, @nt) ;
	                my ($a2) = scalar grep (/^$al[1]$/, @nt) ;
	                my (@ttl) = () ;
	                my (@hap1) = () ;
	                my (@hap2) = () ;

			for ( my ($i)=0; $i < $n_sam-1; $i++ ){
                        for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
                                # ave and sd
                                my ($ky) = "$i\t$j" ;
                                unless (exists $av{$ky} ){print "E1\n" ; exit ; }
                                unless (exists $sd{$ky} ){print "E2\n" ; exit ; }

                                my ($sts) = $dist[$h] ;
                                # forward       
                                for ( my ($p1)=($h+1); $p1 < $c_dt1; $p1++ ){
                                        if ( $dt1[$i][$p1] ne $dt1[$j][$p1] || $p1 == $c_dt1-1 ){                                       
                                                $sts = $dist[$p1-1] ; last ;
                                        }

                                }

                                # backword
                                my ($end) = $dist[$h] ;
                                for ( my ($p2)=($h-1); $p2 >= 0; $p2-- ){
                                        if ( $dt1[$i][$p2] ne $dt1[$j][$p2] || $p2 == 0 ){
                                                $end = $dist[$p2+1] ; last ;
                                        }
                                }
                                my ($hap_dist) = ($sts-$end-$av{$ky})/$sd{$ky} ;
                                push (@ttl, $hap_dist) ;        
                                if ( $nt[$i] eq $al[0] && $nt[$j] eq $al[0] ){
                                        push (@hap1, $hap_dist) ;
                                }
                                if ( $nt[$i] eq $al[1] && $nt[$j] eq $al[1] ){
                                        push (@hap2, $hap_dist) ;
                                }

                        } # $j
                	} # $i
			
			my ($whole) = &mean (@ttl) ;
			# allele 1
	                unless ( $a1 == 1 ){
		                my ($allele1) = &mean (@hap1) ;
        		        my ($stat1) = $allele1-$whole ;
        		        print OUT "$tar\t$al[0]\t$n_sam\t$a1\t$stat1\n" ;
                	}

                	# allele 2
                	if ( $a2 > 1 ){
                        	my ($allele2) = &mean (@hap2) ;
                        	my ($stat2) = $allele2-$whole ;
                        	print OUT "$tar\t$al[1]\t$n_sam\t$a2\t$stat2\n" ;
                	}



		} # biallic
		# fix
		elsif ( $n_al == 1 ) {
			my (@ttl) = () ;
                	for ( my ($i)=0; $i < $n_sam-1; $i++ ){
                        for ( my ($j)=($i+1); $j < $n_sam; $j++ ){
                                # ave and sd
                                my ($ky) = "$i\t$j" ;
                                unless (exists $av{$ky} ){print "E1\n" ; exit ; }
                                unless (exists $sd{$ky} ){print "E2\n" ; exit ; }

                                my ($sts) = $dist[$h] ;
                                # forward       
                                for ( my ($p1)=($h+1); $p1 < $c_dt1; $p1++ ){
                                        if ( $dt1[$i][$p1] ne $dt1[$j][$p1] || $p1 == $c_dt1-1 ){
                                                $sts = $dist[$p1-1] ; last ;
                                        }

                                }

                                # backword
                                my ($end) = $dist[$h] ;
                                for ( my ($p2)=($h-1); $p2 >= 0; $p2-- ){
                                        if ( $dt1[$i][$p2] ne $dt1[$j][$p2] || $p2 == 0 ){
                                                $end = $dist[$p2+1] ; last ;
                                        }
                                }
                                my ($hap_dist) = ($sts-$end-$av{$ky})/$sd{$ky} ;
                                push (@ttl, $hap_dist) ;
                        } # $j
                	} # $i
        
               		 my ($whole) = &mean (@ttl) ;
                	# allele 1
                	print OUT "$tar\t$al[0]\t$n_sam\t$n_sam\t$whole\n" ;

 
		} # fix
		else {
			print "hogehgoe $n_al\n"; exit ; 
		}


	}


} # @at
#print "$cc\n" ;

close (OUT) ;


		
exit ;
