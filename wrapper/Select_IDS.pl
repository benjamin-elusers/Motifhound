#!/usr/bin/perl
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
use autodie qw(open close); # open/close succeed or die
use File_Utils;
use String_Print_Utils;


sub Select_Sequence {

	my ($Ref_Seq,$IDS_file,$Mapfile) = @_;
	my @IDS_selected=@{load_IDS($IDS_file)};
	my %Sequences = %{$Ref_Seq};

	my %Selected; my %Conversion;my %Rejected;
	my @keep;
	my $i=1; my $cpt = 0;
	if($Mapfile ne ""){
		print  "- ID Mapping file      :    \"".$Mapfile."\" => "; my ($rc_map,$TEXT) = File_verification($Mapfile); print $TEXT;
		if($rc_map == 0){ 
			open(CONV,"$Mapfile") or die "Unable to open file $Mapfile for : $!\n";
			while(<CONV>){ if($_ !~ /[\*]/){ push(@keep,$_); } }# print STDERR $_; 
			close(CONV);
			foreach (@keep){
				chomp($_);
				my @c=split("\t",$_); my @Uni=split("_",$c[1]);
				$Conversion{$c[0]}->{"UniProtAC"}=$Uni[1];
			}
		}
	}

	print "\n".format_string(30,"#","#","#",">")."\n";
	foreach my $ID (@IDS_selected){
		if(not exists $Rejected{$ID} ){
			printf("%3d : #%s# ->",$i,$ID);
			if(exists $Sequences{$ID}){
				if(not exists($Selected{$ID})){ $cpt++;	$Selected{$ID}->{"Seq"} = $Sequences{$ID}->{"Seq"}; print " OK\n"; }
				elsif ( exists $Selected{$ID} ){ print " already Added (ID duplicated?)\n"; }
			}elsif(exists $Conversion{$ID} and not exists $Sequences{$ID}){
				my $ID_homolog=$Conversion{$ID}->{"UniProtAC"};
				if(not exists($Selected{$ID_homolog})){ $cpt++;	$Selected{$ID_homolog}->{"Seq"} = $Sequences{$ID_homolog}->{"Seq"}; print " OK\n"; }
				elsif ( exists $Selected{$ID_homolog} ){ print " already Added (ID duplicated?)\n"; }
			}elsif(not exists $Conversion{$ID} and not exists $Sequences{$ID}){
				$Rejected{$ID}->{"Seq"} = "";
				print " Rejected (Not found in FastaFile)\n";
			}
		}else{
			if ( exists $Rejected{$ID} ){ printf('%3d',$i); print ": #$ID# -> already Rejected\n"; }
		}
		$i++;
	}
	print "\n".format_string(30,"#","#","#",">")."\n";
	print " $cpt sequence(s) has/have been selected over the ".($#IDS_selected+1)." IDS, in the input list of IDS.\n";
	return(\%Selected,\%Rejected);
}

my $help = ""; my $quiet= "";
my $Fastafile = ""; my $resdir = "./Seqs.faa"; my $IDS_file = "";
my $logdir = "./"; my $order=1;
my $Mapfile="";
# Processing command line options
GetOptions(
	"help!"							=> \$help,
	"nov|quiet!"					=> \$quiet,
	"out:s"							=> \$resdir,
	"logdir:s"     					=> \$logdir,
	"Fastafile=s" 					=> \$Fastafile,
	"IDS=s"							=> \$IDS_file,
	"alphabetical|ord|alpha|a!"		=> \$order,
	"Map|m:s"						=> \$Mapfile
	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
			 * --Fastafile     <File Directory>
			 * --IDS           <File Directory>
			   --Map           <File Directory>
			   --out           <File or Path Directory>
			   --logdir        <Path Directory>
			   --alpha
			   --help OR --nov
			(*) Mandatory options\n";

# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#####################################################################################################################################\n";
	print "	#                           Selection of IDs - Pick the sequences corresponding to the IDs from a list -                            #\n";
	print "	#-----------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 8 parameters :                                                                                         #\n";
	print "	#    1)  --help                              + Description of Program Usage                                      +  (Optional)      #\n";
	print "	#    2)  --quiet,--nov                       + Non verbous Mode                                                  +  (Optional)      #\n";
	print "	#    3)  --Fastafile                         + Fasta Formatted File of Proteins Sequences                        +  (Mandatory)     #\n";
	print "	#    4)  --IDS                               + File that contains a list of IDS                                  +  (Optional)      #\n";
	print "	#    5)  --out                               + Location of the Output File                                       +  (Optional)      #\n";
	print "	#    6)  --logdir                            + Location of the Log Directory                                     +  (Optional)      #\n";
	print "	#    7)  --a            (Active by default)  + Write IDs in alphabetical order                                   +  (Optional)      #\n";
	print "	#    8)  --m                                 + Mapping Table (used only if an ID was not found in the FASTA)     +  (Optional)      #\n";
	print "	#####################################################################################################################################\n";
	exit(0);
}

# Si le mode non verbeux est active
if ($quiet){
	$logdir=abs_path($logdir);
	# on redirige les affichages a l'ecran dans un fichier log_seq_output.txt dans le repertoire courant
	open (STDOUT,">$logdir"."/Selection_IDs.log") or die ("Unable to redirect STDOUT to a file...");
}

if ($Fastafile ne ""){
	$Fastafile=abs_path($Fastafile); $IDS_file = abs_path($IDS_file);
	print  "- Sequence file       :    \"".$Fastafile."\" => "; my ($rc_seq,$Fastafile,$N) = guess($Fastafile);
	if($rc_seq  != 1){ print "# Something wrong happened when verifying the location of the Sequence file #\n"; exit($rc_seq);  }
	print  "- List IDs file       :    \"".$IDS_file."\" => "; my ($rc_list,$IDS_file,$n) = guess($IDS_file);
	if($rc_list != 2){ print "# Something wrong happened when verifying the location of the list IDS file #\n"; exit($rc_list); }

	my ($Ref_S) = Load_FASTA($Fastafile);
	my @IDS = keys(%{$Ref_S});

	my @ids = @{load_IDS($IDS_file)};
	print format_string(90,"=","=","\n",">");
	print "  Reading the list of $n IDs from the List IDs file.\n";
	print "  Selection of the corresponding sequences from the $N sequences in the Sequence File.\n";
	print format_string(90,"=","=","\n",">");

	print format_string(80,"-","-","\n",">");
	print " Creation of a FASTA file which contains all selected sequences (ID listed below).\n";
	print format_string(80,"-","-","\n",">");

	my ($Ref_Selected,$Ref_Rejected) = Select_Sequence($Ref_S,$IDS_file,$Mapfile);
	my ($fout,$fdir)=fileparse(abs_path($resdir));
	print " ".format_string(76+length($fout),"=","=","\n",">");
	print "  Location of the output file : $fdir$fout\n";
	print " ".format_string(76+length($fout),"=","=","\n",">");
	Write_Seq_to_FASTA($Ref_Selected,$fdir.$fout,$order);
	my @rejected_IDs= keys (%{$Ref_Rejected});
	if($#rejected_IDs+1 >= 1){
		Write_Seq_to_FASTA($Ref_Rejected,$fdir."/Rejected_IDs.txt",$order);
	}
}else{
	print "\n".format_string(30+length($0),"#","#","#",">");
	print "\n".format_string(30+length($0),"#    For more details, read the help :"," ","#",">");
	print "\n".format_string(30+length($0),"#        $0       --help"," ","#",">");
	print "\n".format_string(30+length($0),"#","#","#",">")."\n\n";
	exit(0);
}
