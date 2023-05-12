#! /usr/bin/perl
#------------------------------------------------------------+
# Date    : 30/07/2012                                       |
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
	if( exists $ENV{'SCRIPTSDIR'}){ $ENV{'DATADIR'}=abs_path($ENV{'SCRIPTSDIR'}."/../Data"); }else{ $ENV{'DATADIR'}="./"; }
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File_Utils;
use String_Print_Utils;

# Silent options
my $help = ""; my $quiet= "";
# Shared options
my $Seq_file = ""; 
# Homology options
my $Homology=""; my $percId=40; my $evalue=0.000001; my $NB_THREADS = 1;
my $Blast_file=abs_path($ENV{'DATADIR'}."Blast/Blast_Yeast_Proteome_sorted.txt");
# Disorder options
my $Disorder=""; my $minCF=5;
my $Diso_file=abs_path($ENV{'DATADIR'}."Disorder/YEAST.dat");

my $resdir = "./Seqs_Filtered.faa"; my $logdir="./";

# Processing command line options
GetOptions(
	"help!"        => \$help,
	"nov|quiet!"   => \$quiet,
	"out:s"        => \$resdir,
	"logdir:s"     => \$logdir,
	"Seq=s"        => \$Seq_file,
	"Diso:s"       => \$Diso_file,
	"Blast:s"      => \$Blast_file,
	"p|percId:i"   => \$percId,
	"e|evalue:f"   => \$evalue,
	"thr:i"        => \$NB_THREADS,
	"minCF:i"      => \$minCF,
	"Homology|H!"  => \$Homology,
	"Disorder|D!"  => \$Disorder
	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
--------------------------------/*       SHARED OPTIONS       */--------------------------------
			 *   --Seq           <File Directory>
			     --out           <Path or File Directory>
			     --logdir        <Path or File Directory>
			     --help OR --nov

--------------------------------/* HOMOLOGY FILTERING OPTIONS */--------------------------------
			 1   --H
			     --Blast         <File Directory>
			     --percID        <Integer Value>
			     --evalue        <Floating Value>
			     --thr           <Integer Value>

--------------------------------/* DISORDER FILTERING OPTIONS */--------------------------------
			 2   --D
			 2   --Diso_file     <File Directory>
			     --minCF         <Integer Value>

			(1) Required options for Homology Filtering
			(2) Required options for Disorder Filtering
			(*) Mandatory options\n";

# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	##################################################################################################################################################\n";
	print "	#                                 Sequence Filtering - Replacing Amino Acids by X in Homologuous or Ordered Regions -                            #\n";
	print "	#------------------------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 13 parameters :                                                                                                     #\n";
	print "	#   1)  --help                                            + Description of Program Usage                                      +  (Optional)      #\n";
	print "	#   2)  --quiet,--nov                                     + Non verbous Mode                                                  +  (Optional)      #\n";
	print "	#   3)  --Seq                                             + Fasta Formatted File of Proteins Sequences                        +  (Mandatory)     #\n";
	print "	#   5)  --Diso          ( Disorder )                      + File that contains Disorder Predictions                           +  (Mandatory)     #\n";
	print "	#   4)  --Blast         ( Homology )                      + File that contains Blast Results                                  +  (Optional)      #\n";
	print "	#   6)  --out                                             + Location of the Output Directory                                  +  (Optional)      #\n";
	print "	#   7)  --logdir                                          + Location of the Log Directory                                     +  (Optional)      #\n";
	print "	#   8)  --percId        (40    by default - Homology )    + Minimum Percentage of Identity to get rid of Homologuous regions  +  (Optional)      #\n";
	print "	#   9)  --evalue        (10e-5 by default - Homology )    + Minimum Evalue to consider the alignment of homologuous regions   +  (Optional)      #\n";
	print "	#  10)  --thr           (4     by default - Homology )    + Number of Threads (number of CPUs) to use for the Blast           +  (Optional)      #\n";
	print "	#  11)  --minCF         (5     by default - Disorder )    + Minimal Confidence Value for Disordered Regions (between 0 and 9) +  (Optional)      #\n";
	print "	#  12)  --H             ( Homology )                      + Activate the Homology Filter which masks the homologuous regions  +  (Optional)      #\n";
	print "	#  13)  --D             ( Disorder )                      + Activate the Disorder Filter which masks the ordered regions      +  (Optional)      #\n";
	print "	##################################################################################################################################################\n";
	exit(0);
}


if($Seq_file ne "" and ($Homology or $Disorder) ){
	$Seq_file=abs_path($Seq_file); $resdir = abs_path($resdir);
	my ($Seqname,$Seqdir,$Seqext)=fileparse($Seq_file,qr{\..*});
	my ($outname,$outdir)=fileparse($resdir);

	my ($fh_HOMO,$fh_DISO)=(\*STDERR,\*STDERR);
	# Quiet mode - Homology.log -
	if ($quiet and $Homology){
		$fh_HOMO=();
		# Redirect displayed output to a log file ('Homology.log') to the user defined directory
		open ($fh_HOMO,">".$logdir."/Homology.log") or die ("Unable to write in file $logdir/Homology.log...");
	}
	# Quiet mode - Disorder.log -
	if ($quiet and $Disorder){
		$fh_DISO=();
		# Redirect displayed output to a log file ('Disorder.log') to the user defined directory
		open ($fh_DISO,">".$logdir."/Disorder.log") or die ("Unable to write in file $logdir/Disorder.log...");
	}
	
	# File verifications #
	# Existence of the Sequence file
	my ($rc_seq,$TEXT);
	for my $fh ($fh_HOMO,$fh_DISO){ print $fh "- Sequence file       :    \"".$Seq_file."\" => "; ($rc_seq,$TEXT) = File_verification($Seq_file); print $fh $TEXT; }
	if($rc_seq != 0){ die "# Something wrong happened when verifying the location of the Sequence file #\n"; }
	my $Ref_SEQUENCES = Load_FASTA($Seq_file);
	my %SEQUENCES = %{$Ref_SEQUENCES};

	my @IDS_SEQ = (keys %SEQUENCES);
	if($Homology){
		# Existence of the BLAST results file
		print $fh_HOMO "- BLAST file          :    \"$Blast_file\" => "; my ($rc_blast,$TEXT2) = File_verification($Blast_file); print $fh_HOMO $TEXT2;
		# No Blast results file - Launching Blast File 1 against itself #
		if($rc_blast != 0){ 
			if(!(-e "$Seq_file.psd" and -e "$Seq_file.psi" and -e "$Seq_file.phr" and -e "$Seq_file.pin" and -e "$Seq_file.psq")){
				print $fh_HOMO "Creation of the database...";
				`formatdb -i $Seq_file -p T -o T -l $logdir/formatdb.log`;
				print $fh_HOMO "OK\n";
			}
			if(!(-f -s -e -T "$Blast_file")){
				print $fh_HOMO "Blastall \"$Seqname.faa\" against \"$Seqname.faa\"...";
				`blastall -p blastp -i $Seq_file -d $Seq_file -a $NB_THREADS -m 8 -o $Blast_file -e $evalue 2>> "$logdir"/Homology.log`;
				print $fh_HOMO "OK\n";
			}

			my ($blastname,$blastdir,$blastext)=fileparse($Blast_file,qr{\..*});
			`sort $Blast_file | awk -v percID=$percId -F '\t' '{ if(\$3 > percID && \$1 != \$2){ print \$0 } }' > $blastdir"$blastname".dat`;
			unlink("$Seq_file.pin","$Seq_file.psd","$Seq_file.phr","$Seq_file.psi","$Seq_file.psq",$Blast_file);
			$Blast_file=$blastdir.$blastname.".dat";
			print $fh_HOMO "- New BLAST file      :    \"$Blast_file\" => "; my ($rc_blast2,$TEXT3) = File_verification($Blast_file); print $fh_HOMO $TEXT3;
		}
		my %data_Blast = %{Load_data_BLAST($Blast_file)}; 
#		unlink($Blast_file);
		
		my @ID1 = sort { $data_Blast{$a} cmp $data_Blast{$b} } (keys (%data_Blast));
		print $fh_HOMO "/=====HOMOLOGY FILTER=====/   Applying Homology Filter...\r";

		# Filter the homologuous regions based on the percentage of identity between two sequences
		my %ID_Filtered=(); my $cpt_IDS=0;
#		open(OUT,">$Blast_file");
		foreach my $ID1 (@ID1){
			if(($cpt_IDS+1) % 100 == 0){
				print $fh_HOMO "/=====HOMOLOGY FILTER=====/   Applying Homology Filter...(".($cpt_IDS+1)."/".($#ID1+1).")                            \r";
			}
			my @ID2 = sort { $data_Blast{$ID1}->{$a} cmp $data_Blast{$ID1}->{$b} } (keys (%{$data_Blast{$ID1}}));
			foreach my $ID2 (@ID2){
				foreach my $cpt (sort {$data_Blast{$ID1}->{$ID2}->{$b}->{"percId"} <=> $data_Blast{$ID1}->{$ID2}->{$a}->{"percId"}} (keys (%{$data_Blast{$ID1}->{$ID2}}))){
#					print OUT "$ID1\t$ID2\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"percId"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"alnLength"}."\t";
#					print OUT $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_start"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_end"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID2_start"}."\t";
#					print OUT $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID2_end"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"E-value"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"Score"};
					my $start = $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_start"};
					my $end = $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_end"};
					my $replacement=my_print("X",$end-$start);
	#				print STDERR "$ID1\t$ID2\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"percId"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"alnLength"}."\t";
	#				print STDERR $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_start"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID1_end"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID2_start"}."\t";
	#				print STDERR $data_Blast{$ID1}->{$ID2}->{$cpt}->{"ID2_end"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"E-value"}."\t".$data_Blast{$ID1}->{$ID2}->{$cpt}->{"Score"};
	#				print STDERR ">$ID1\n".$SEQUENCES{$ID1}->{"Seq"}."\n";
	#				print STDERR "Replacement: $replacement Pos:$start-$end\n";
			
					substr($SEQUENCES{$ID1}->{"Seq"},$start-1,$end-$start,$replacement);
					delete($data_Blast{$ID1}->{$ID2}->{$cpt}); if( exists $data_Blast{$ID2}->{$ID1}->{$cpt} ){ delete($data_Blast{$ID2}->{$ID1}->{$cpt}); }
	#				if( exists $ID_Filtered{$ID1}->{$cpt})){
	#					$ID_Filtered{$ID1}->{$cpt}=1;
	#				}
				}
			}
			$cpt_IDS++;
		}


		print $fh_HOMO "/=====HOMOLOGY FILTER=====/   Applying Homology Filter...(".($cpt_IDS)."/".($#ID1+1).") IDS Filtered. => Finished\r";
		print $fh_HOMO "/=====HOMOLOGY FILTER=====/   Applying Homology Filter...OK                                                           \n";
#		close(OUT);
	}

	if($Disorder){
		# Existence of the DISO results file
		print $fh_DISO "- Disorder file    :    \"$Diso_file\" => "; my ($rc_diso,$TEXT2) = File_verification($Diso_file); print $fh_DISO $TEXT2;
		if($rc_diso != 0){ die "# Something wrong happened when verifying the location of the Disorder Predictions result file #\n"; }

		if($minCF > 9){ $minCF=9; }elsif($minCF < 0){ $minCF=0; }
		my $Ref_data_Diso = Load_data_DISO($Diso_file);
		my %data_Diso = %{$Ref_data_Diso};

		print $fh_DISO "/=====DISORDER FILTER=====/   Applying Disorder Filter...\r";

		my $i=0; my @IDS=keys(%data_Diso); my $cpt_IDS=0;
		foreach my $ID (keys(%data_Diso)){
			if(($cpt_IDS+1) % 100 == 0){
				print $fh_DISO "/=====DISORDER FILTER=====/   Applying Disorder Filter...".($cpt_IDS+1)."/".($#IDS_SEQ+1)." IDS Filtered.              \r";
			}
			if( $ID =~ /\|([A-Z0-9]+)\|/ ){ $ID = $1; };
			if(exists $SEQUENCES{$ID}){
	#			print STDERR "$ID\n";
#				if($ID eq "O00141"){
#					print STDERR ">$k | $ID LenDis ".length($data_Diso{$k}->{"Seq"})." LenSeq ".length($SEQUENCES{$ID}->{"Seq"})."                   \n";
#					print STDERR $data_Diso{$k}->{"Seq"}."\n";
#					print STDERR $SEQUENCES{$ID}->{"Seq"}."\n";
#					print STDERR $data_Diso{$k}->{"Dis"}."\n";
#					<STDIN>;
#				}

		#		my @dis=split($data_Diso{$k}->{"Dis"},"");
				my @seq=split("",$SEQUENCES{$ID}->{"Seq"});
		#		my @CF=split("",$data_Diso{$k}->{"CF_Dis"});
				my @disoseq=split("",$data_Diso{$ID}->{"Dis"});
				my @tabdisonew=();
				for($i=0;$i<($#disoseq+1);$i++){ $tabdisonew[$i]=int($disoseq[$i]); }
				for($i=5;$i<($#disoseq-4);$i++){
					my @before = @disoseq[$i-5..$i-1]; my @after = @disoseq[$i+1..$i+5];
					my $Window = sum(@before) + sum(@after) + $disoseq[$i];
					$tabdisonew[$i]=$Window;
		#			print "$k $seq[$i] $i Value Window = ".$Window." (@before+@after)\n";
		#			print "$k WINDOW     : ".$seq[$i-1].$seq[$i-2].$seq[$i-3].$seq[$i-4].$seq[$i-5]."|".$seq[$i]."|".$seq[$i+1].$seq[$i+2].$seq[$i+3].$seq[$i+4].$seq[$i+5]."\n";
		#			print "$k WINDOW     : ".$disoseq[$i-1].$disoseq[$i-2].$disoseq[$i-3].$disoseq[$i-4].$disoseq[$i-5]."|".$disoseq[$i]."|".$disoseq[$i+1].$disoseq[$i+2].$disoseq[$i+3].$disoseq[$i+4].$disoseq[$i+5]."\n";
		#			print "\n";
		#			if($disoseq[$i] < 5){ 
		#				if($CF[$i] < $minCF){
		#					$seq[$i]="X";
		#				}
		#			}else{
		#				$seq[$i]="X";
		#			}
				}

				my $cpt=0;

				for($i=0;$i<($#tabdisonew+1);$i++){ 
		#			print "$i $tabdisonew[$i] $seq[$i]\n"; 
					if($i < 5 and $tabdisonew[$i] == 0 ){ 
						$seq[$i] = "X"; 
					}
					if($i >= 5 and $i < ($#disoseq-4) and $tabdisonew[$i] < 6 ){ 
						$seq[$i] = "X"; 
					}
					$cpt++;
					if($i >= ($#disoseq-4) and $tabdisonew[$i] == 0 ){ 
						$seq[$i] = "X"; 
					}
				}
				$cpt_IDS++;
				$SEQUENCES{$ID}->{"Seq"}=join("",@seq) unless not exists($SEQUENCES{$ID});

#				if($ID eq "O00141"){
#					print STDERR ">$k | $ID LenDis ".length($data_Diso{$k}->{"Seq"})." LenSeq ".length($SEQUENCES{$ID}->{"Seq"})."  LenDisNew    ".($#tabdisonew+1)."                \n";
#					print STDERR $data_Diso{$k}->{"Seq"}."\n";
#					print STDERR $SEQUENCES{$ID}->{"Seq"}."\n";
#					print STDERR "@tabdisonew"."\n";
#					<STDIN>;
#				}
			}
		}
		print $fh_DISO "/=====DISORDER FILTER=====/   Applying Disorder Filter...".($cpt_IDS)."/".($#IDS_SEQ+1)." IDS Filtered => Finished\r";
		print $fh_DISO "/=====DISORDER FILTER=====/   Applying Disorder Filter...OK                                                      \n";
	}

	# Write the sequences without the homologuous regions (replaced by X)
	open(OUT2,">".$outdir.$outname);
	foreach my $ID ( sort { $a cmp $b } (keys (%SEQUENCES))){
		print OUT2 ">$ID\n".fasta_format($SEQUENCES{$ID}->{"Seq"});
	}
	close(OUT2);
}elsif($Disorder eq "" and $Homology eq "" or $Seq_file eq ""){
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, read the help :"," ","#",">");
	print "\n".format_string(70,"#        $0       --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
