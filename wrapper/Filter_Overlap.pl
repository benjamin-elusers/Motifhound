#! /usr/bin/perl
#------------------------------------------------------------+
# Date    : 30/07/2012  (Last Time modified : 11/03/2013)    |
# Authors : B. Dubreuil                                      |
# Contact : dubreuil.benjamin@hotmail.com                    |
#------------------------------------------------------------+

use strict;
use warnings;
BEGIN
{
	use File::Basename;
	use Cwd qw(getcwd chdir abs_path);
	use FindBin qw($Bin $Script);
	$ENV{'DIR'}=abs_path("$Bin");
	use if (exists $ENV{'SCRIPTSDIR'}), lib => $ENV{'SCRIPTSDIR'};
	use if (not exists $ENV{'SCRIPTSDIR'} and exists $ENV{'DIR'}), lib => $ENV{'DIR'};
	if(not exists $ENV{'DIR'} and not exists $ENV{'SCRIPTSDIR'}){
		die "No Environment Variable for locating Scripts Directory. Please refer to documentation.\n";
	}
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File_Utils;
use String_Print_Utils;

sub load_all_motifs_in_Hash {

	my ($Motifs_file,$Smin,$Smax)=@_;

	my %All_MOTIFS=();
	my $nb_tot_motifs=0; my $s=3;
	my $nb_motifs_cutoff=100;
	my $nb_motifs=0;
	my ($fname,$dir,$ext)  = fileparse($Motifs_file,qr{\..*});
	print "\n    o Loads all motifs (from length $Smin to $Smax)...";
	for($s=$Smin; $s <= $Smax ; $s++){
		$nb_motifs=0;
		my ($Ref_MOTIFS) = load_motifs("$dir/$fname".".L$s");
		my %MOTIFS = %{$Ref_MOTIFS}; my @Mots = sort { $MOTIFS{$b}->{"Pval"} <=> $MOTIFS{$a}->{"Pval"} } (keys (%MOTIFS));
		foreach my $motif ( @Mots ){
			$nb_tot_motifs++;
			if( $nb_motifs < $nb_motifs_cutoff){
				$nb_motifs++;
				$All_MOTIFS{$motif}->{"Domain_Homolog"}=0;
				$All_MOTIFS{$motif}->{"S_set"}= $MOTIFS{$motif}->{"S_set"};
				$All_MOTIFS{$motif}->{"S_pop"}=$MOTIFS{$motif}->{"S_pop"};
				$All_MOTIFS{$motif}->{"N_set"}= $MOTIFS{$motif}->{"N_set"};
				$All_MOTIFS{$motif}->{"N_pop"}= $MOTIFS{$motif}->{"N_pop"};
				$All_MOTIFS{$motif}->{"Pval"}=$MOTIFS{$motif}->{"Pval"};
			}
		}
	}
	my @final_Motifs = keys(%All_MOTIFS);
	my $nb_top_motifs = ($#final_Motifs+1);
	print "$nb_tot_motifs motifs from Length $Smin to $Smax. \n";
	print "      $nb_top_motifs retained motifs (top 100 motifs for each Length from $Smin to $Smax). \n";
	return(\%All_MOTIFS,$nb_top_motifs);
}

sub log10 {
# Computes the log10 of a given number $n
	my $n = shift;
	return (log($n)/log(10));
}

sub intersect {
#     returns an integer which corresponds to the relative positions between 2 intervals (I1 & I2) :
#  1. rc= -2  Length(Overlap) = L(I1 inter I2)          |     Overlap between I2 and I1 ( I2 starts in the middle of I1 )                 
#  2. rc= -1  Length(Overlap) = 1                       |     I1 & I2 Adjacents ( the last pos of I1 = the first pos of I2 )              
#  3. rc= 0   Length(Overlap) = 0                       |     No overlap ( No positions shared between I1 & I2 but I1 before I2)          
#  4. rc= 1   Length(Overlap) = L(I2)                   |     I2 included in I1 ( I1 before I2 but both ends at the same position )       
#  5. rc= 2   Length(Overlap) = L(I2)                   |     I2 inside I1 ( I1 contains I2 )                                             
#  6. rc= 3   Length(Overlap) = L(I1)                   |     I1 included in I2 ( I1 & I2 begins at the same position but L(I2) > L(I1) ) 
#  7. rc= 4   Length(Overlap) = L(I1) = L(I2)           |     I1 equal I2 ( perfect Overlap )                                             
#  8. rc= 5   Length(Overlap) = L(I1 inter I2)          |     I2 included in I1 (I1 & I2 begins at the same position but L(I1) > L(I2) )  
#  9. rc= 6   Length(Overlap) = L(I1)                   |     I1 inside I2 ( I2 contains I1 )                                             
# 10. rc= 7   Length(Overlap) = L(I1)                   |     I1 included in I2 ( I2 before I1 but ends both at the same position )       
# 11. rc= 8   Length(Overlap) = L(I1 inter I2)          |     Overlap between I2 and I1 ( I1 before I2 )                                  
# 12. rc= 9   Length(Overlap) = 1                       |     I1 & I2 Adjacents ( the last pos of I2 = the first position of I1 )         
# 13. rc= 10  Length(Overlap) = 0                       |     No overlap ( No positions shared between I1 & I2 but I2 before I1)          

	my($s1,$e1,$s2,$e2)=@_;

	my $v1 = ($s1 > $s2) + ($s1 >= $s2);
	my $v2 = ($e1 > $e2) + ($e1 >= $e2);
#	print "m1 $s1-$e1 m2 $s2-$e2\t\t v1=$v1 v2=$v2 ";

	my $v3=0;
	if($v1*3+$v2 == 0){
		$v3 += ($e1 > $s2) + ($e1 >= $s2);
#		print "v3=$v3";
	}elsif($v1*3+$v2 == 8){
		$v3 -= ($s1 > $e2) + ($s1 >= $e2);
#		print "v3=$v3";
	}
#	print "\n";
	return ($v1*3+$v2-$v3);
}

sub Compute_Overlap {
	
	my ($m1,$m2,$Ref_Seq)=@_;
	
	my %Seq=%{$Ref_Seq}; my $i=0; my $i_orf=1;
	my ($S1,$S2,$shared)=(0,0,0);
	my $L1 = length($m1); my $L2 = length($m2);
	my @ORFS = keys (%Seq);
#	print "$m1 et $m2 \n";

	foreach my $orf (@ORFS){
#		if($i_orf % 100 ){ print "ORF #$i_orf/".($#ORFS+1)."\r"; }
		my $seq = $Seq{$orf}->{"Seq"};
		if($seq =~ /$m1/ and $seq !~ /$m2/){
			$S1+=$L1;
		}elsif($seq !~ /$m1/ and $seq =~ /$m2/){
			$S2+=$L2;
		}elsif($seq =~ /$m1/ and $seq =~ /$m2/){
			my ($s1,$e1)=match_positions($m1,$seq); my ($s2,$e2)=match_positions($m2,$seq);
			my $intersect = intersect($s1,$e1,$s2,$e2);
			$i++;
			if( $intersect == 0 or $intersect == 10){								# No overlap
				$S1+=$L1; $S2+=$L2;
			}elsif ( $intersect == 4 ){												# Equals
				my $L = $L1 = $L2;
				$S1+=$L; $S2+=$L; $shared+=$L;
			}elsif ( $intersect == -1 or $intersect == 9 ){							# Adjacents
				$S1+=$L1; $S2+=$L2; $shared+=1;
			}elsif ( $intersect == -2 ){											# Overlap
				$S1+=$L1; $S2+=$L2; $shared+=abs($e1-$s2+1);
			}elsif ( $intersect == 8 ){												# Overlap
				$S1+=$L1; $S2+=$L2; $shared+=abs($s1-$e2+1);
			}elsif ( $intersect == 3 or $intersect == 6 or $intersect == 7 ){		# Overlap (included)
				$S1+=$L1; $S2+=$L2; $shared+=$L1;
			}elsif ( $intersect == 1 or $intersect == 2 or $intersect == 5 ){		# Overlap (included)
				$S1+=$L1; $S2+=$L2; $shared+=$L2;
			}
		}
		$i_orf++;
	}
#	print STDERR "\n";

	return($S1,$S2,$shared);
}

my ($script,$dir) = fileparse($0);
# Silent options
my $help = ""; my $quiet= "";
# File path
my $resfile = "No_overlap.res"; my $Seq_file = ""; my $Motifs_file = ""; 
my $Threshold = 50;
# Directory path
my $logdir = "./";
my @Size;

# Processing command line options
GetOptions(
	"help!"         => \$help,
	"nov|quiet!"    => \$quiet,
	"logdir=s"      => \$logdir,
	"resfile=s"     => \$resfile,
	"Motifs_file=s" => \$Motifs_file,
	"Seq=s"         => \$Seq_file,
	"Threshold:i"   => \$Threshold,
	"Size=i{2}"     => \@Size
	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
---------------------------------------------------------------------------------------------------------------------------------------
     Options         |  Arguments                          |  Description                                                              
---------------------------------------------------------------------------------------------------------------------------------------
 *  --Seq            |  <File directory>                   |  Sequence File's Location                                                 
                     |                                     |                                                                           
 *  --resfile        |  <File directory>                   |  Location of Results Files                                                
                     |                                     |                                                                           
 *  --Motifs_file    |  <File directory>                   |  Motif Enumeration File's Location                                        
                     |                                     |                                                                           
    --logdir         |  <File or Path directory>           |  Location or Name of the results files                                    
                     |                                     |                                                                           
    --Size           |  <Integer value>   <Integer value>  |  Minimum length and Maximum length of the Motifs                          
                     |                                     |                                                                           
    --Threshold      |  <Integer value>                    |  Min. percentage to consider two motifs as overlapping motifs             
                     |                                     |                                                                           
    --nov            |                                     |  Avoid printing message (Non verbose mode)                                
                     |                                     |                                                                           
    --help           |                                     |  Prints the help description and exit                                     
                     |                                     |                                                                           
---------------------------------------------------------------------------------------------------------------------------------------
(*) Mandatory options                                                                                                                  \n";

# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#########################################################################################################################################\n";
	print "	#                                    HYPERMOTIF - A de novo Method For Linear Motif Discovery -                                         #\n";
	print "	#------------------------------------------------------------------------------------------------------------------------------- -------#\n";
	print "	# This program takes at most 8 parameters :                                                                                             #\n";
	print "	#    1)  --help                                    + Description of Program Usage                                    +  (Optional)      #\n";
	print "	#    2)  --quiet,--nov                             + Non verbose Mode                                                +  (Optional)      #\n";
	print "	#    3)  --Motifs_file                             + Motifs Enumeration file                                         +  (Mandatory)     #\n";
	print "	#    4)  --Seq                                     + File containing the sequences where occurs those motifs         +  (Mandatory)     #\n";
	print "	#    5)  --resfile                                 + Name of the results file or result directory                    +  (Mandatory)     #\n";
	print "	#    6)  --logdir                                  + Name of the Log directory                                       +  (Optional)      #\n";
	print "	#    7)  --Threshold   (50% by default)            + Threshold (%) for removing overlapping motifs                   +  (Optional)      #\n";
	print "	#    8)  --Size        (From 3 to 10 by default)   + Min and Max Motifs Sizes (2 values required between 3 to 12 AA) +  (Optional)      #\n";
	print "	#########################################################################################################################################\n";
	exit(0);
}


if($Seq_file ne "" and $Motifs_file ne ""){

	# Non verbose mode 
	if ($quiet){
		# Redirection of printed messages to a log file ('Overlap.log')
		open(STDOUT,">$logdir/Overlap.log") or die ("Unable to redirect STDOUT to a file...");
	}

	# Existence of the Sequence file
	print "- Sequence file             :    \"$Seq_file\" => "; my ($rc_seq,$TEXT) = File_verification($Seq_file); print $TEXT;
	if($rc_seq != 0){ die "# Something wrong happened when verifying the location of the Sequence file #\n"; }
	# Existence of the Motif Enumeration file
	print "- Motif Enumeration file    :    \"$Motifs_file\" => "; my ($rc_motif,$TEXT2) = File_verification($Motifs_file); print $TEXT2;
	if($rc_motif != 0){ die "# Something wrong happened when verifying the location of the Motif Enumeration file #\n"; }

	my ($Ref_MOTIFS,$nb_motifs) = load_motifs($Motifs_file);
	my %MOTIFS = %{$Ref_MOTIFS};
	my $SEQUENCES = Load_FASTA($Seq_file);

	open(TEST,">test.txt");
	my %new_MOTIFS=(); my $i=1; my $j=1; my @motifs=keys(%MOTIFS); my $chosen_motif;
	foreach my $motif (sort @motifs){
		my @Overlap; $j=1;
		foreach my $motif2 (sort @motifs){
			$Overlap[$i]=0;
			if( $i % 10 == 0 and $j == $nb_motifs ){ print "/=====OVERLAP  FILTER=====/   Applying Overlap Filter...motif $i against motif $j (".($#motifs+1)." motifs)          \r"; }
			if($motif ne $motif2){
				my ($s1,$s2,$shared) = Compute_Overlap($motif,$motif2,$SEQUENCES);
				if($Threshold <= (100*($shared/($s1))) and $Threshold <= (100*($shared/($s2)))){
					print TEST "$Threshold $motif $motif2 ".(100*($shared/($s1)))." ".(100*($shared/($s2)))."\n";
					$Overlap[$i]=1;
					if($MOTIFS{$motif}->{"Pval"} > $MOTIFS{$motif2}->{"Pval"}){
						$chosen_motif = $motif2;
					}elsif($MOTIFS{$motif}->{"Pval"} < $MOTIFS{$motif2}->{"Pval"}){
						$chosen_motif = $motif;
					}elsif($MOTIFS{$motif}->{"Pval"} == $MOTIFS{$motif2}->{"Pval"}){
						if($MOTIFS{$motif}->{"Cdef"} > $MOTIFS{$motif2}->{"Cdef"}){
							$chosen_motif=$motif2;
						}elsif($MOTIFS{$motif}->{"Cdef"} < $MOTIFS{$motif2}->{"Cdef"}){
							$chosen_motif=$motif;
						}elsif($MOTIFS{$motif}->{"Cdef"} == $MOTIFS{$motif2}->{"Cdef"}){
							if($MOTIFS{$motif}->{"S_set"} > $MOTIFS{$motif2}->{"S_set"}){
								$chosen_motif=$motif;
							}elsif($MOTIFS{$motif}->{"S_set"} < $MOTIFS{$motif2}->{"S_set"}){
								$chosen_motif=$motif2;
							}
						}
					}
					$new_MOTIFS{$chosen_motif}->{"Len"}	= $MOTIFS{$motif2}->{"Len"};	$new_MOTIFS{$chosen_motif}->{"Cdef"}	= $MOTIFS{$motif2}->{"Cdef"};
					$new_MOTIFS{$chosen_motif}->{"S_set"}	= $MOTIFS{$motif2}->{"S_set"};	$new_MOTIFS{$chosen_motif}->{"S_pop"}	= $MOTIFS{$motif2}->{"S_pop"};
					$new_MOTIFS{$chosen_motif}->{"N_set"}	= $MOTIFS{$motif2}->{"N_set"};	$new_MOTIFS{$chosen_motif}->{"N_pop"}	= $MOTIFS{$motif2}->{"N_pop"};
					$new_MOTIFS{$chosen_motif}->{"Pval"}	= $MOTIFS{$motif2}->{"Pval"};
					}
				}
			$j++;
		}
		$MOTIFS{$motif}->{"OL"}=\@Overlap;
		$i++;
	}
	
	my @new_MOTIFS = keys(%new_MOTIFS);
	print TEST ($#new_MOTIFS+1);
	open(NO_OVERLAP,">>$resfile"); my $i_mot=1;
	foreach my $motif (sort { $new_MOTIFS{$a}->{"Pval"} <=> $new_MOTIFS{$b}->{"Pval"} } (@new_MOTIFS)){
		if( $i_mot % 1 == 0 ){ print TEST "    o Writing the results...Motifs($i_mot/".($#new_MOTIFS+1).")\n"; }
		my @l=@{$MOTIFS{$motif}->{"OL"}};
		if( sum(@l) < 1){
			my ($k,$m,$n,$N) = ($MOTIFS{$motif}->{"S_set"},$MOTIFS{$motif}->{"S_pop"},$MOTIFS{$motif}->{"N_set"},$MOTIFS{$motif}->{"N_pop"});
			my $Pval = $MOTIFS{$motif}->{"Pval"};
			print NO_OVERLAP $motif."\t".$k."\t".$m."\t".$n."\t".$N."\t".$Pval."\n";
		}else{
			print TEST "$motif\n";
		}
		$i_mot++;
	}
	close(NO_OVERLAP);
	close(TEST);
	print "\n";
}else{
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, consult the help to launch the program :"," ","#",">");
	print "\n".format_string(70,"#        $0       --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
