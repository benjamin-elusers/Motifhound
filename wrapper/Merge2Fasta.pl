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
	use if (exists $ENV{'SCRIPTSDIR'}), lib => $ENV{'SCRIPTSDIR'};
	if (not exists $ENV{'SCRIPTSDIR'}) {	die "No Environment Variable for locating Scripts Directory ( \$SCRIPTSDIR ). Please refer to documentation.\n"; }
}
use File::Basename;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use File_Utils;
use String_Print_Utils;

sub List_shared_IDS {

	my($file1,$file2)=@_;
	
	my $Ref_Seq1 = Load_FASTA($file1); my $Ref_Seq2 = Load_FASTA($file2);

	my %Seq1 = %{$Ref_Seq1}; my %Seq2 = %{$Ref_Seq2};
	my @IDS_Seq1= sort{ $a cmp $b } (keys %Seq1 );
	my @IDS_Seq2= sort{ $a cmp $b } (keys %Seq2 );
	my @IDS_shared=();

	foreach my $ID ( @IDS_Seq1 ){	push(@IDS_shared,$ID); }
	foreach my $ID ( @IDS_Seq2 ){	push(@IDS_shared,$ID); }

	my %tmp=();
	foreach (@IDS_shared){
		if((exists($Seq1{$_}->{"Seq"})) and (exists($Seq2{$_}->{"Seq"}))){ 
			$tmp{$_}++ 
		}
	}
	@IDS_shared = sort{ $a cmp $b } (keys %tmp);

	my $nb_max_seq = max($#IDS_Seq1,$#IDS_Seq2);

	print STDERR "==========================================> IDS shared = ".($#IDS_shared+1)." <==========================================\n";
#	print STDERR " Max nb of sequences = ".($nb_max_seq+1)." F1 = ".($#IDS_Seq1+1)." F2 = ".($#IDS_Seq2+1)."\n";
	my @IDS_Set=(); my @IDS_Prot=(); my $Ref_Set=(); my $Ref_Prot=();
	if    ( $nb_max_seq == ($#IDS_Seq1) ){
		@IDS_Prot = @IDS_Seq1; @IDS_Set = @IDS_Seq2; $Ref_Prot = $Ref_Seq1; $Ref_Set = $Ref_Seq2;
		print STDERR " Proteome (".($#IDS_Seq1+1)." IDS)  : $file1 \n Set      (".($#IDS_Seq2+1)." IDS)   : $file2 \n";
	}elsif( $nb_max_seq == ($#IDS_Seq2) ){
		@IDS_Prot = @IDS_Seq2; @IDS_Set = @IDS_Seq1; $Ref_Prot = $Ref_Seq2; $Ref_Set = $Ref_Seq1;
		print STDERR " Proteome (".($#IDS_Seq2+1)." IDS)  : $file2 \n Set      (".($#IDS_Seq1+1)." IDS)   : $file1 \n";
	}
	print STDERR "---------------------------------------------------------------------------------------------------------\n";

	return(\@IDS_shared,\@IDS_Set,\@IDS_Prot,$Ref_Set,$Ref_Prot);
}

my $help = ""; my $quiet= ""; 
my $oneline="";
my $all=""; my $shared= ""; my $not_shared= ""; my $proteome= "";
my $Seq_file1 = ""; my $Seq_file2 = "";
my $resfile = "New_Merged.faa"; my $logdir = ".";
my %Merged=();
my $i=0; my $cpt_IDS=0; my $k="";

# Processing command line options
GetOptions(
	"all"                =>         \$all,
	"help!"              =>        \$help,
	"nov|quiet!"         =>       \$quiet,
	"shared!"            =>      \$shared,
	"logdir:s"           =>      \$logdir,
	"resfile:s"          =>     \$resfile,
	"1Ls|one-line-seq!"  =>     \$oneline,
	"proteome!"          =>    \$proteome,
	"Seq1=s"             =>   \$Seq_file1,
	"Seq2=s"             =>   \$Seq_file2,
	"not-shared!"        =>  \$not_shared,

	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]

			   --help                                        -     Display help    -
			   --nov                                         -   Non-verbose Mode  -
			 * --Seq1                <File Directory>        -     FASTA File 1    -
			 * --Seq2                <File Directory>        -     FASTA File 2    -
			   --one-line-seq                                -   Writing options   -
			   --resfile             <Path Directory>        -  FASTA File output  -
			   --logdir              <Path Directory>        -    Log directory    -
			   --shared                                      - IDS saving options -
			   --not-shared                                  - IDS saving options -
			   --proteome                                    - IDS saving options -
			   --all                                         - IDS saving options -
			(*) Mandatory options\n";

# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#########################################################################################################################################\n";
	print "	#                                       Merge2Fasta - Merge two FASTA formatted files in one -                                          #\n";
	print "	#---------------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 11 parameters :                                                                                            #\n";
	print "	#    1)  --help                              + Description of Program Usage                                             +  (Optional)   #\n";
	print "	#    2)  --quiet,--nov                       + Non verbous Mode                                                         +  (Optional)   #\n";
	print "	#    3)  --Seq1                              + 1st FASTA Formatted File of Proteins Sequences                           +  (Mandatory)  #\n";
	print "	#    4)  --Seq2                              + 2nd FASTA Formatted File of Proteins Sequences                           +  (Mandatory)  #\n";
	print "	#    5)  --one-line-seq                      + Write each sequence in a single line with only one carriage return (\\n) +  (Optional)   #\n";
	print "	#    6)  --shared                            + Write only shared IDS from two FASTA formatted file                     +  (Optional)   #\n";
	print "	#    7)  --not-shared                        + Write only not shared IDS from two FASTA formatted files                +  (Optional)   #\n";
	print "	#    8)  --all                               + Write all IDS from two FASTA formatted files                            +  (Optional)   #\n";
	print "	#    9)  --proteome                          + Write IDS from the smallest file if IDS are shared between the two     +  (Optional)   #\n";
	print "	#                                            |  FASTA formatted files, and IDS from the largest file if IDS are not   |               #\n";
	print "	#                                            |  shared between the two FASTA formatted files                            |               #\n";
	print "	#   10)  --resfile                           + Name of the new Merged FASTA formatted file                              +  (Optional)   #\n";
	print "	#   11)  --logdir                            + Directory for saving log files                                           +  (Optional)   #\n";
	print "	#########################################################################################################################################\n";
	exit(0);
}

if($Seq_file1 ne "" and $Seq_file2 ne ""){

	# quiet mode
	if ($quiet){
		# Redirection of printed messages to a log file ('Disorder_output.log') to the user defined directory
		open(STDERR,">$logdir/Merge2Fasta.log") or die ("Unable to redirect STDERR to a file...");
	}

	# Existence of the 1st Sequence file
	print STDERR "- 1st Sequence file    : \"$Seq_file1\" => "; my ($rc_seq1,$TEXT1) = File_verification($Seq_file1); #print STDERR $TEXT1;
	if($rc_seq1 != 0){ die "# Something wrong happened when verifying the location of the 1st Sequence file #\n"; }
	# Existence of the 2nd Sequence file
	print STDERR "- 2nd Sequence file    : \"$Seq_file2\" => "; my ($rc_seq2,$TEXT2) = File_verification($Seq_file2); #print STDERR $TEXT2;
	if($rc_seq2 != 0){ die "# Something wrong happened when verifying the location of the 2nd Sequence file #\n"; }

	my ($Ref_shared_IDs,$Ref_IDS_Set,$Ref_IDS_Prot,$Ref_Set,$Ref_Prot) = List_shared_IDS($Seq_file1,$Seq_file2);

	my %Set = %{$Ref_Set}; my %Prot = %{$Ref_Prot};
	my @IDS_shared = @{$Ref_shared_IDs} ;
	my @IDS_Set    = @{$Ref_IDS_Set} ;
	my @IDS_Prot   = @{$Ref_IDS_Prot} ;
	
	if($all){
		print STDERR "  == Concatenate the two fasta files and replace the IDS from the largest file with IDS from the smallest file ==  \n\n";
		$cpt_IDS=0; $k="";
		my %All_IDS=();
		foreach $k (@IDS_Prot){ $All_IDS{$k}=1; }
		foreach $k (@IDS_Set){  $All_IDS{$k}=1; }
		my @all_IDs = keys %All_IDS;
		foreach $k (keys %All_IDS){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Saving ID #$k# ".($cpt_IDS+1)."/".($#all_IDs+1)."       \r";
			}
			
			if(exists($Prot{$k})){
				$Merged{$k}->{"Seq"} = $Prot{$k}->{"Seq"};
			}
			
			if(exists($Set{$k})){
				$Merged{$k}->{"Seq"} = $Set{$k}->{"Seq"};
			}
			$cpt_IDS++;
		}
		print STDERR "\n-> Saving ID : Finished\n";
		# Turn off the other saving options
		$proteome="";
		$shared="";
		$not_shared="";
	}

	if($proteome){
		print STDERR "  == Replace IDS from the largest sequence file with IDS from the smallest if they are shared ==  \n\n";
		$cpt_IDS=0; $k="";
		my %All_IDS=();
		foreach $k (@IDS_Prot){ $All_IDS{$k}=1; }
		foreach $k (@IDS_Set){  $All_IDS{$k}=1; }
		my @all_IDs = keys %All_IDS;
		foreach $k (keys %All_IDS){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Saving ID #$k# ".($cpt_IDS+1)."/".($#all_IDs+1)."       \r";
			}
			if(exists($Prot{$k})){
				$Merged{$k}->{"Seq"} = $Prot{$k}->{"Seq"};
				if(exists($Set{$k})){
					$Merged{$k}->{"Seq"} = $Set{$k}->{"Seq"};
				}
			}
			$cpt_IDS++;
		}
		print STDERR "\n-> Saving ID : Finished\n";
		# Turn off the other saving options
		$shared="";
		$not_shared="";
	}

	if($shared){
		print STDERR "  == Retains IDS only if they are shared between the two FASTA formatted files ==  \n\n";
		$cpt_IDS=0; $k="";
		foreach $k (@IDS_shared){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Saving ID #$k# ".($cpt_IDS+1)."/".($#IDS_shared+1)."       \r";
			}
			$Merged{$k}->{"Seq"} = $Set{$k}->{"Seq"};
			$cpt_IDS++;
		}
		print STDERR "\n-> Saving ID : Finished\n";
		# Turn off the other saving options
		$not_shared="";
		$cpt_IDS=0; $k="";
	}

	if($not_shared){
		print STDERR "  == Retains IDS only if they are not shared between the two FASTA formatted files ==  \n\n";
		$cpt_IDS=0; $k="";
		my %All_IDS=();
		foreach $k (@IDS_Prot){ if(not exists($Set{$k})){ $All_IDS{$k}=1; } }
		foreach $k (@IDS_Set){ if(not exists($Prot{$k})){ $All_IDS{$k}=1; } }
		my @IDS_not_shared = keys %All_IDS;
		foreach $k (keys %All_IDS){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Saving ID #$k# ".($cpt_IDS+1)."/".($#IDS_not_shared+1)."       \r";
			}
			$Merged{$k}->{"Seq"} = $Prot{$k}->{"Seq"} unless exists($Prot{$k});
			$Merged{$k}->{"Seq"} = $Set{$k}->{"Seq"}  unless exists($Set{$k});
			$cpt_IDS++;
		}
		print STDERR "\n-> Saving ID : Finished\n";
	}

	my $k2=""; $cpt_IDS=0;
	open(OUT,">".$resfile);
	my @IDS_Merged = sort { $a cmp $b } (keys %Merged);
	if($oneline){
		foreach $k2 (@IDS_Merged){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Writing ID #$k2# ".($cpt_IDS+1)."/".($#IDS_Merged+1)."       \r";
			}
			print OUT ">$k2\n".join("\n",$Merged{$k2}->{"Seq"})."\n";
			$cpt_IDS++;
		}
	}else{
		foreach $k2 ( sort { $a cmp $b } (keys %Merged) ){
			if( ($cpt_IDS+1) % 100 == 0){
				print STDERR "Writing ID #$k2# ".($cpt_IDS+1)."/".($#IDS_Merged+1)."       \r";
			}
			print OUT ">$k2\n".fasta_format($Merged{$k2}->{"Seq"});
			$cpt_IDS++;
		}
	}
	close(OUT);
	print STDERR "\n-> Writing IDS : Finished\n";
}else{
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, read the help :"," ","#",">");
	print "\n".format_string(70,"#        $0       --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
