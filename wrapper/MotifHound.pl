#! /usr/bin/perl
#------------------------------------------------------------+
# Date    : 30/07/2012  (Last Time modified : 11/03/2013)    |
# Authors : B. Dubreuil                                      |
# Contact : dubreuil.benjamin@hotmail.com                    |
#------------------------------------------------------------+
use strict;
use warnings;
use diagnostics;
BEGIN
{
	use File::Basename;
	use Cwd qw(getcwd chdir abs_path);
	use FindBin qw($Bin $Script);
	$ENV{'MOTIFHOUNDDIR'}=abs_path("$Bin/..");
	$ENV{'SCRIPTSDIR'}=$ENV{'MOTIFHOUNDDIR'}."/Scripts";
	$ENV{'DATADIR'}=$ENV{'MOTIFHOUNDDIR'}."/Data";
	$ENV{'SRCDIR'}=$ENV{'MOTIFHOUNDDIR'}."/Sources";
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use autodie qw(open close); # open/close succeed or die
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Path qw(mkpath);
use Benchmark;
use lib $ENV{'SCRIPTSDIR'};
use File_Utils;
use String_Print_Utils;

# Silent options
my $help = ""; my $quiet= "";
# Minimum Occurence for the Motif Enumeration
my $MinOcc=3;
# Min and Max Number of non-wildcard positions (by default)
my $Dmax=10; my $Dmin=3;
my @Def=($Dmin,$Dmax);
# User defined Sizes (2 values: Size Min and Size Max)
# Size max and min (by default)
my $Smax=10; my $Smin=3;
my @Size=($Smin,$Smax);
# Number of Sequences in the Set and in the Proteome
my $nseq=0; my $Nseq=0;
# File or Path directories
## (Mandatory) Inputs
my $Setfile = ""; my $Proteome = ""; 
## (Optional) Output directory
my $Workdir = "./";
## Location of precomputed Datafiles
my $Map="";
my $Disofile=abs_path($ENV{'DATADIR'}."/Disorder/YEAST.dat");
my $Blastfile=abs_path($ENV{'DATADIR'}."/Blast/Blast_Yeast_Proteome_sorted.txt");
my $Pfam_annot=abs_path($ENV{'DATADIR'}."/Domains/Yeast_Pfam_annotations.txt");
my $Superfamily_annot=abs_path($ENV{'DATADIR'}."/Domains/Saccharomyces_cerevisiae_SUPERFAMILY_domains.txt");
my $Genefile=abs_path($ENV{'DATADIR'}."/Genes/YEAST.data");
# Homology, Overlap and Disorder Filters
my $Homology_Filter = ""; my $Overlap_Filter = ""; my $Disorder_Filter = "";
# Method to look for motif in Proteome and counting motifs occurences
my $Scan = 1; my $OnePerSeq = 1;
# Number of simultaneous threads to launch
my $THR=1;
# Output Activation 
my $HTML = ""; my $Graph = "";
# Processing command line options
my @tab_command = split("--",join(" ",@ARGV));
my $Command="Command: $0\nOptions: \n"; foreach(@tab_command){ if($_ ne ""){ $Command.=" --$_\n"; } }

GetOptions(
	"help!"              => \$help,
	"nov|quiet!"         => \$quiet,

	"Setfile=s"          => \$Setfile,
	"Proteome=s"         => \$Proteome,

	"Map:s"              => \$Map,
	"Disorder|D!"        => \$Disorder_Filter,
	"Disofile:s"         => \$Disofile,
	"Homology|H!"        => \$Homology_Filter,
	"Blastfile:s"        => \$Blastfile,
	"Overlap|O!"         => \$Overlap_Filter,

	"Threads:i"          => \$THR,

	"Size:i{2}"          => \@Size,
	"nonW|nw:i{1,2}"     => \@Def,
	"MinOcc:i"           => \$MinOcc,

	"Scan!"              => \$Scan,
	"Enumeration|Enum!"  => sub { $Scan = 0 },

	"one!"               => \$OnePerSeq,
	"all!"               => sub { $OnePerSeq = 0 },

	"Gene_annot:s"       => \$Genefile,
	"Pfam_annot:s"       => \$Pfam_annot,
	"Superfam_annot:s"   => \$Superfamily_annot,

	"Graphical|graphic!" => \$Graph,
	"HTML|html!"         => \$HTML,
	"Workdir|WD:s"       => \$Workdir)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
-------------------------------------------------------------------------------------------------------------------------------------------------------
     Options         |  Arguments                          |  Description                                                                              
-------------------------------------------------------------------------------------------------------------------------------------------------------
    --help           |                                     |  Display the help description and exit                                                    
                     |                                     |                                                                                           
 *  --Proteome       |  <File directory>                   |  Location of the file that contains the Proteome sequences which should contain the       
                     |                                     |   sequences of the Setfile (background model)                                             
                     |                                     |                                                                                           
 *  --Setfile        |  <File directory>                   |  Location of the Setfile that contains Protein sequences of interest                      
                     |                                     |                                                                                           
    --Map            |  <File directory>                   |  Location of the ID Mapping file                                                          
                     |                                     |                                                                                           
    --WD             |  <Path directory>                   |  Location of the Working directory that will contain all outputs                          
                     |                                     |                                                                                           
    --Size           |  <Integer value>   <Integer value>  |  Minimum and Maximum length for the motif Enumeration                                     
    --nonW           |  <Integer value>   <Integer value>  |  Minimum and Maximum number of non-wildcard positions for the motif Enumeration           
                     |                                     |                                                                                           
    --MinOcc         |  <Integer value>                    |  Minimum number of  occurences for a motif to be considered                               
                     |                                     |                                                                                           
    --H --O --D      |                                     |  Activate the filters that remove Overlapping motifs, (--O)                               
                     |                                     |                       that remove Homologous regions, (--H)                               
                     |                                     |                       that remove Ordered regions,    (--D)                               
                     |                                     |                                                                                           
                     |                                     |  Select a method to look for motifs in Proteome :                                         
    --Scan           |  ( active by default )              |   (--Scan) Scan the Proteome and only take into account the Motifs from the Set           
    --Enum           |                                     |   (--Enum) Enumerates both Set and Proteome files                                         
                     |                                     |                                                                                           
                     |                                     |  Select a counting method for motif occurrences :                                         
    --one            |  ( active by default )              |   (--one) Counting only one occurrence per Sequence                                       
    --all            |                                     |   (--all) Counting all occurrences in Sequences                                           
                     |                                     |                                                                                           
    --graphic        |                                     |  Make a graphical ouptut for seeing where motifs appear in the sequences                  
                     |                                     |                                                                                           
    --HTML           |                                     |  Make an HTML output                                                                      
                     |                                     |                                                                                           
    --nov            |                                     |  Only print Logo, progression and time usage                                              
                     |                                     |                                                                                           
                     |                                     |  By default, data files are available for S. Cerevisiae proteome :                        
    --Disofile       |  <File directory>                   |      Disorder predictions file obtained from Disopred 2.0                                 
    --Blastfile      |  <File directory>                   |      Output file obtained from BLAST listing homologuous regions                          
    --Pfam_annot     |  <File directory>                   |      Domains Annotation file obtained from PFAM Database                                  
    --Superfam_annot |  <File directory>                   |      Domains Annotation file obtained from SUPERFAMILY Database                           
    --Gene_annot     |  <File directory>                   |      File containing Names and Descriptions of genes                                      
                     |                                     |                                                                                           
    --THR            |  <Integer value>                    |  Number of simultaneous threads to launch ( do not exceed twice times the number of CPUs  
                     |                                     |                                                                                           
-------------------------------------------------------------------------------------------------------------------------------------------------------
(*) Mandatory options                                                                                                                                  \n";

# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#####################################################################################################################################################\n";
	print "	##                                           MOTIFHOUND - A de novo Method For Linear Motif Discovery -                                            ##\n";
	print "	##-------------------------------------------------------------------------------------------------------------------------------------------------##\n";
	print "	## This program takes at most 24 parameters :                                                                                                      ##\n";
	print "	##     1) --help                                     + Description of Program Usage                                             +  (Optional)      ##\n";
	print "	##     2) --quiet,--nov                              + Non verbous Mode                                                         +  (Optional)      ##\n";
	print "	##     3) --Setfile                                  + File of Proteins Sequences of Interest (FASTA or list of IDS)            +  (Mandatory)     ##\n";
	print "	##     4) --Proteome                                 + Fasta Formatted File containing Proteome Sequences                       +  (Mandatory)     ##\n";
	print "	##     5) --Map                                      + File of ID mapping                                                       +  (Optional)      ##\n";
	print "	##     6) --Disorder,D     (inactive by default)     + Activate the Disorder Filter                                             +  (Optional)      ##\n";
	print "	##     7) --Overlap,O      (inactive by default)     + Activate the Overlap  Filter                                             +  (Optional)      ##\n";
	print "	##     8) --Homology,H     (inactive by default)     + Activate the Homology Filter                                             +  (Optional)      ##\n";
	print "	##     9) --Size           (From 3 to 10 by default) + Min and Max Motifs Sizes (2 values required between 3 to 12 AA)          +  (Optional)      ##\n";
	print "	##    10) --nonW           (From 3 to 10 by default) + Min and Max non-wildcard positions allowed (1 or 2 values required)      +  (Optional)      ##\n";
	print "	##    11) --MinOcc         (3 by default)            + Minimum Motifs Occurence                                                 +  (Optional)      ##\n";
	print "	##    12) --Workdir,--WD                             + Location of the Working Directory (by default current location)          +  (Optional)      ##\n";
	print "	##    13) --one            (activate by default)     + Counting only one occurrence per Sequence for each motif                 +  (Optional)      ##\n";
	print "	##    14) --all            (inactive by default)     + Counting all occurrences in Sequences for each motif                     +  (Optional)      ##\n";
	print "	##    15) --Scan           (activate by default)     + Scan the Proteome with motifs from the set                               +  (Optional)      ##\n";
	print "	##    16) --Enum           (inactive by default)     + Enumerates Motifs from both Setfile and Proteome                         +  (Optional)      ##\n";
	print "	##    17) --Graphical                                + Graphical Output for the top 100 Motifs                                  +  (Optional)      ##\n";
	print "	##    18) --HTML                                     + HMTL Output for top 100 over-rerpresented Motifs                         +  (Optional)      ##\n";
	print "	##    19) --Disofile       (by default for yeast)    + Disorder predictions file obtained by DISOPRED                           +  (Optional)      ##\n";
	print "	##    20) --Blastfile      (by default for yeast)    + Blast output file listing homologuous regions in S. Cerevisiae           +  (Optional)      ##\n";
	print "	##    21) --Pfam_annot     (by default for yeast)    + Pfam domains annotations file for S. Cerevisiae                          +  (Optional)      ##\n";
	print "	##    22) --Superfam_annot (by default for yeast)    + Superfamily domains annotations file for S. Cerevisiae                   +  (Optional)      ##\n";
	print "	##    23) --Gene_annot     (by default for yeast)    + File containing Names and Descriptions of S. Cerevisiae genes            +  (Optional)      ##\n";
	print "	##    24) --THR            ( 2 by default )          + Number of Simultaneous threads to launch                                 +  (Optional)      ##\n";
	print "	#####################################################################################################################################################\n";
	exit(0);
}

sub Logo {
	my ($t0,$Command)=@_;
	print STDERR format_string(90," ****     ****            **   **   **** **      **                                **"," ","\n",">");
	print STDERR format_string(90,"/**/**   **/**           /**  //   /**/ /**     /**                               /**"," ","\n",">");
	print STDERR format_string(90,"/**//** ** /**  ******  ****** ** ******/**     /**  ******  **   ** *******      /**"," ","\n",">");
	print STDERR format_string(90,"/** //***  /** **////**///**/ /**///**/ /********** **////**/**  /**//**///**  ******"," ","\n",">");
	print STDERR format_string(90,"/**  //*   /**/**   /**  /**  /**  /**  /**//////**/**   /**/**  /** /**  /** **///**"," ","\n",">");
	print STDERR format_string(90,"/**   /    /**/**   /**  /**  /**  /**  /**     /**/**   /**/**  /** /**  /**/**  /**"," ","\n",">");
	print STDERR format_string(90,"/**        /**//******   //** /**  /**  /**     /**//****** //****** ***  /**//******"," ","\n",">");
	print STDERR format_string(90,"//         //  //////     //  //   //   //      //  //////   ////// ///   //  ////// "," ","\n",">");
	print STDERR format_string(90,"====================================================================================="," ","\n",">");
	print STDERR format_string(90,"-------------------| Author: Dubreuil Benjamin & Kelil Abdellali |-------------------"," ","\n",">");
	print STDERR format_string(90,"-------------------| Mailto: dubreuil.benjamin\@hotmail.com       |-------------------"," ","\n",">");
	print STDERR format_string(90,"-------------------|         Abdellali.Kelil\@umontreal.ca        |-------------------"," ","\n",">");
	print STDERR format_string(90,"-------------------| From MICHNIK LAB  -  Universite de Montreal |-------------------"," ","\n",">");
	print STDERR format_string(90,"-------------------| http://michnick.bcm.umontreal.ca/           |-------------------"," ","\n",">");
	print STDERR "\n";
	print "$Command\n";
	print "				   RUN START : ".localtime($t0)."\n";
}

# File and Directory organization
sub File_and_Directory_Managment {

	my ($Proteome,$Setfile,$Workdir)=@_;

	my $Scriptsdir=abs_path($ENV{'SCRIPTSDIR'});

	my $oldProteome=$Proteome;
	if ((-e -s -f -T "$oldProteome") and !(-d $oldProteome)){ $Proteome=$oldProteome; }else{ print "### \"$Proteome\" is not found. ###\n"; exit(2); }
	my $oldSetfile=$Setfile;
	if ((-e -s -f -T "$oldSetfile") and !(-d $oldSetfile)){ $Setfile=$oldSetfile; }else{ print "### \"$Setfile\" is not found. ###\n"; exit(2); }

	my ($Motifhound,$path_to_Motifhound)=fileparse($0,qr{\..*});
	if(!(-d abs_path($ENV{'DATADIR'}))){
		print STDERR "# \"Data\" directory not found at ".$ENV{'DATADIR'}.". #\n";
		print STDERR "# Have you moved MotifHound.pl from the original directory ? #\n";
		exit(1);
	}elsif(!(-d abs_path($ENV{'SRCDIR'}))){ 
		print STDERR "# \"Sources\" directory not found at ".$ENV{'SRCDIR'}.". #\n";
		print STDERR "# Have you moved MotifHound.pl from the original directory ? #\n";
		exit(1);
	}

	$Workdir=abs_path("$Workdir") unless !(-d $Workdir);
	# If the user has specified a Working Directory (unless it corresponds to the current directory, or the project directory)
	if(!(-d $Workdir)){
		mkpath($Workdir,0,0766) or die "# Cannot create directory \"$Workdir\"  : $! #\n";
		$Workdir=abs_path("$Workdir");
#		chdir($Workdir) or die "# Unable to move to directory \#$Workdir\" for : $! #\n";
	}elsif(-d $Workdir and abs_path($Workdir) ne abs_path($ENV{'PWD'}) ){
		$Workdir=abs_path("$Workdir");
#		chdir($Workdir) unless !(-d $Workdir);
	}elsif(abs_path($Workdir) eq abs_path($Scriptsdir) ){ 
		$Workdir=abs_path("$Workdir/../");
#		chdir("$Workdir/../");
#		$Workdir=getcwd;
	}else{ 
		$Workdir=abs_path("$Scriptsdir/../");
#		chdir("$Scriptsdir/../") or die "# Unable to change directory to \#$Scriptsdir\" for : $! #\n";
#		$Workdir=getcwd;
	}

	my ($setname,$setdir) = fileparse($Setfile); my $setext="";
	if($Setfile =~ /(\..*)$/){ ($setname,$setdir,$setext) = fileparse($Setfile,qr{\..*}); }
	my $Execdir=abs_path($Workdir."/Executables");
	my $Logdir=$Workdir."/Log/$setname";
	my $Enumdir=abs_path($Workdir."/Enumeration");
	my $Overdir=abs_path($Workdir."/Over_Representation");
	# Make new directories for the Results ouptut, statistics on Time usage and Memory usage, Log files, and Executables
	my @DIR_to_create=("$Enumdir","$Overdir","$Execdir","$Logdir");
	foreach(@DIR_to_create){
		mkpath($_,0,0766) or die "# Cannot create directory \"$_\"  : $! #\n" unless -d $_;
	}

	# Change all created directories to match the working directory
	unlink ($Logdir."/*.log");
	my %Final=( "Proteome"        => abs_path($Proteome)          , "Setfile"     => abs_path($Setfile)      , "Mapfile"     => abs_path($Map)      ,
			    "Workdir"         => $Workdir                     , "Logdir"      => $Logdir                 ,
			    "Enumdir"         => $Enumdir                     , "Overdir"     => $Overdir                ,
			    "Blastfile"       => abs_path($Blastfile)         , "Disofile"    => abs_path($Disofile)     ,
			    "Superfam_annot"  => abs_path($Superfamily_annot) , "Pfam_annot"  => abs_path($Pfam_annot)   ,
			    "Gene_file"       => abs_path($Genefile)          , "Execdir"     => $Execdir                );

	return(\%Final);
}

sub Print_Parameters {

	my($Ref_AllPaths,$Reference_Proteome,$Ref_Filters,$Ref_Size,$Ref_Def,$Smax,$Smin,$K,$Graph,$HTML,$Scan,$OnePerSeq)=@_;

	my %AllPaths=%{$Ref_AllPaths};
	my $Workdir=$AllPaths{"Workdir"}; my $Logdir=$AllPaths{"Logdir"};
	my $Setfile=$AllPaths{"Setfile"}; my $Proteome=$AllPaths{"Proteome"};
	my @Size=@{$Ref_Size}; my @Def=@{$Ref_Def}; 

	my $nseq=0; my $Nseq=0; 
	my $Path_rc_set=""; my $Path_rc_prot=""; my $Path_rc_dir=""; my $rc=""; my $Text="";
	print "************************************************** PARAMETERS *******************************************************\n\n";
	# Scan or Enumeration
		print " 			=================================================================			 \n";
	if($Scan){ 
		print " 			| Looking for motifs in Proteome by :   Scanning the Proteome   |			 \n";
	}else{
		print " 			| Looking for motifs in Proteome by :  Enumerating the Proteome |			 \n";
	}
		print " 			-----------------------------------------------------------------			 \n";
	if($OnePerSeq){ 
		print " 			|           Counting only one occurrence per Sequence           |			 \n";
	}else{
		print " 			|             Counting all occurrences in Sequences             |			 \n";
	}
		print " 			=================================================================			 \n";

	# Filters Status
	my ($D,$H,$O)=("","","");
	my $Homology_Filter=@{$Ref_Filters}[0]; my $Disorder_Filter=@{$Ref_Filters}[1]; my $Overlap_Filter=@{$Ref_Filters}[2];
	if($Homology_Filter){ $H=format_string(9,"Yes"," ","|",">"); @{$Ref_Filters}[0]="Yes";}else{ $H=format_string(9,"No"," ","|",">"); @{$Ref_Filters}[0]="No";}
	if($Disorder_Filter){ $D=format_string(9,"Yes"," ","|",">"); @{$Ref_Filters}[1]="Yes";}else{ $D=format_string(9,"No"," ","|",">"); @{$Ref_Filters}[1]="No";}
	if($Overlap_Filter) { $O=format_string(9,"Yes"," ","|",">"); @{$Ref_Filters}[2]="Yes";}else{ $O=format_string(9,"No"," ","|",">"); @{$Ref_Filters}[2]="No";}

	# Minimum number of occurences for any motif
	print " 			|                  Minimum Motif Occurence : $K                  |			\n";

	# Motif Size (From Smin to Smax)
	if(($#Size+1) == 2){
		$Smax=max($Size[0],$Size[1]); if ($Smax > 12){ $Smax=12;}
		$Smin=min($Size[0],$Size[1]); if ($Smin < 3){ $Smin=3;}
		if ($Smax < $Smin){ $Smax=$Smin;}
		if ($Smin > $Smax){ $Smin=$Smax;}
	}

	if(($#Def+1) == 2){
		$Dmax=max($Def[0],$Def[1]); if ($Dmax > 12){ $Dmax=12;} if($Dmax>$Smax){ $Dmax=$Smax; }
		$Dmin=min($Def[0],$Def[1]); if ($Dmin < 3){ $Dmin=3;} if($Dmin<$Smin){ $Dmin=$Smin; }
		if ($Dmax < $Dmin){ $Dmax=$Dmin;}
		if ($Dmin > $Dmax){ $Dmin=$Dmax;}
	}elsif(($#Def+1) == 1){
		$Dmax=$Smax;
		$Dmin=$Def[0]; if ($Dmin < 3){ $Dmin=3;} if($Dmin<$Smin){ $Dmin=$Smin; }
		if ($Dmax < $Dmin){ $Dmax=$Dmin;}
		if ($Dmin > $Dmax){ $Dmin=$Dmax;}
	}
	

	my $sizefromto=format_string(26,"   From ".$Smin." to ".$Smax." AAs"," ","|",">");
	my $deffromto=format_string(26,"   From ".$Dmin." to ".$Dmax." non-wild."," ","|",">");
	print " 			+-------------------------+          +--------------------------+			\n";
	print " 			| Filters activated :     |          | Motif Size :             |			\n";
	print " 			|    -Homology...$H".     "          |".$sizefromto            ."			\n";
	print " 			|    -Disorder...$D".     "          | Non-wildcard positions : |			\n";
	print " 			|    -Overlap ...$O".     "          |".$deffromto             ."			\n";
	print " 			+-------------------------+          +--------------------------+			\n";

	# Additional Outputs
	my($h,$g)=("","");
	if($HTML){  $h=format_string(6,"Yes"," ","|",">"); }else{ $h=format_string(6,"No"," ","|",">"); }
	if($Graph){ $g=format_string(5,"Yes"," ","|",">"); }else{ $g=format_string(5,"No"," ","|",">"); }
	print " 			|     HTML output : $h".  "          |  Graphical output : $g"."			\n";
	print " 			=================================================================			 \n";
	print "\n";

	# File verifications, Merging/Filtering Sequences from Proteome/Setfile
	print "- Current directory     :	\"".$ENV{'PWD'}."\"\n";
	print "- Working directory     :	\"".$Workdir ."\" => ";   ($Path_rc_dir,$AllPaths{"Workdir"})=guess($Workdir);
	print "- Proteome              :	\"".$Proteome."\" => ";   ($Path_rc_prot,$AllPaths{"Proteome"},$Nseq)=guess($Proteome);
	print "- Set of Sequences      :	\"".$Setfile ."\" => ";   ($Path_rc_set,$AllPaths{"Setfile"},$nseq)=guess($Setfile); 
	$AllPaths{"Proteome_old"}=$Proteome;
	# If Sequence filters are enabled
	$AllPaths{"Setfile_old"} = $Setfile;
	# If the Setfile is FASTA Formatted
	if($Path_rc_set == 1){
		$AllPaths{"Proteome"} = Merge_Setfile_Proteome($AllPaths{"Proteome"},$AllPaths{"Setfile"},$Logdir);
		my ($setname,$setdir,$setext) = fileparse($AllPaths{"Setfile_old"},qr{\..*});
		my $Setfile=$setdir.$setname.".IDS";
		write_IDS(Load_FASTA($AllPaths{"Setfile_old"}),$Setfile);
		($Path_rc_set,$AllPaths{"Setfile"},$nseq)=List_IDS_to_Fasta($Setfile,$AllPaths{"Proteome"},$AllPaths{"Mapfile"},$Logdir,$ENV{'SCRIPTSDIR'});
		if($Setfile !~ /([\.faa])$/){ unlink($Setfile); } #print "Removing old Proteome   : $oldProteome\n";   } 
	# If the Setfile is a list of IDS
	}elsif($Path_rc_set == 2){
		($Path_rc_set,$AllPaths{"Setfile"},$nseq)=List_IDS_to_Fasta($AllPaths{"Setfile"},$AllPaths{"Proteome"},$AllPaths{"Mapfile"},$Logdir,$ENV{'SCRIPTSDIR'});
	}

	if(@{$Ref_Filters}[0] eq "Yes" or @{$Ref_Filters}[1] eq "Yes"){
		my $Path_rc_filt=""; my $Path_rc_merg="";
		($Path_rc_filt,$AllPaths{"Setfile"},$nseq)=Sequence_Masking($AllPaths{"Setfile"},$Logdir,$Ref_Filters,$AllPaths{"Blastfile"},$AllPaths{"Disofile"});
		$AllPaths{"Proteome"} = Merge_Setfile_Proteome($AllPaths{"Proteome"},$AllPaths{"Setfile"},$Logdir);
		my ($ext_Diso,$ext_Homo)=("","");
		if(@{$Ref_Filters}[0] eq "Yes"){ $ext_Homo="_H"; }
		if(@{$Ref_Filters}[1] eq "Yes"){ $ext_Diso="_D"; }
		my ($Protname,$Protdir,$Protext) = fileparse($AllPaths{"Proteome"},qr{\..*});
		my $Proteome_filtered=$Protdir.$Protname.$ext_Homo.$ext_Diso.$Protext;
		fcopy($AllPaths{"Proteome"},$Proteome_filtered); if($AllPaths{"Proteome"} !~ /(\.faa|\.fasta|\.fa)$/){ unlink($AllPaths{"Proteome"}); } $AllPaths{"Proteome"}=$Proteome_filtered;
		print "- New Proteome (with filtered sequences from the set) :	\"".$AllPaths{"Proteome"}."\" => "; ($Path_rc_merg,$AllPaths{"Proteome"},$Nseq)=guess($AllPaths{"Proteome"}); 
	}
	
	print "\n";
	print " Number of sequences in the Set      : $nseq.\n";
	print " Number of sequences in the Proteome : $Nseq.\n";
	print "\n";

	print "*********************************************** END OF PARAMETERS ***************************************************\n";
	return($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,\%AllPaths,$Ref_Filters,$Scan);
}

sub Compilation_status {
# Compile C source code and create an executable (binary file).
# If compilation fails, it appends the compilation error to an error log file and stops the run.
	my ($code,$Ref_AllPaths,$fh) = @_;
	my %AllPaths=%{$Ref_AllPaths}; my $Execdir=$AllPaths{"Execdir"}; my $Logdir=$AllPaths{"Logdir"};
	# Name of the source code
	my $Name=$code; $Name=~ tr/_/ /;
	print $fh "__".$Name."__\n";
	my $execmode = 0711;
	# Remove executables and compile the source code
	# unlink("$Execdir/$code.exe");

	# Printing Compilation Command and the status of compilation
	if (!(-e -x "$Execdir/".$code.".exe")){ 
		print $fh "gcc -Wall -ansi -pedantic -g $code.c -o $code.exe -lm -lJudy -lpthread...";
		system("gcc -Wall -ansi -pedantic -g ".$ENV{'SRCDIR'}."/$code.c -o $Execdir/$code.exe -lm -lJudy -lpthread 2>> $Logdir/compilation_status.log");
		# Give the permission to execute the compiled source code
		chmod $execmode, "$Execdir/".$code.".exe";
		if (!(-e -x "$Execdir/".$code.".exe")){ 
			print $fh "Compilation Failed.\n";
			print STDERR "# Check error log file for compilation errors \"$Logdir/compilation_status.log\". #\n";
			print "# Check also permissions of the binary file created ! #\n";
			exit(1);
		}else{
			print $fh "OK\n";
		}
	}
}

sub Merge_Setfile_Proteome {

	my($Proteome,$Setfile,$Logdir)=@_;
	my $rand=""; my $newProteome="";
	my $Scriptsdir=$ENV{'SCRIPTSDIR'};
	print "[Create Background Sequence Model] (see Merge2Fasta.log)\r";
	my ($oldProteome,$olddir,$oldext)=fileparse($Proteome,qr{\..*});
	$newProteome=$olddir.$oldProteome.".".int(rand(1000000))."tmp";
#	while( not -f $newProteome){ $newProteome=$olddir.$oldProteome.".".int(rand(1000000))."tmp"; }
	`perl $Scriptsdir/Merge2Fasta.pl --Seq1 $Proteome --Seq2 $Setfile --all --resfile $newProteome --logdir $Logdir --nov`;
	if(!(-f ($newProteome))){ print "# Something wrong happened when merging the Setfile with the Proteome file #\n"; exit(1); }
	print "[Create Background Sequence Model] (see Merge2Fasta.log) => OK\n";
	return ($newProteome);
}

#sub Compute_AA_frequency {

#	my($Setfile,$Logdir)=@_;
#	my $Scriptsdir=$ENV{'SCRIPTSDIR'};
#	print "[Compute Amino Acid frequencies] (see AA_Freq.log)\r";
#	my ($Setname,$Setdir,$Setext)=fileparse($Setfile,qr{\..*});
#	my $SetFreq="$Setdir/$Setname.freq";
#	`perl $Scriptsdir/AA_Freq.pl --Fastafile $Setfile --out $SetFreq --logdir $Logdir --nov`;
#	print "[Compute Amino Acid frequencies] (see AA_Freq.log) => OK\n";
#	return ($SetFreq);
#}

sub Sequence_Masking{

		my($oldSequence,$Logdir,$Ref_Filters,$Blastfile,$Disofile)=@_;
		my ($oldname,$oldir,$oldext)=fileparse($oldSequence,qr{\..*});

		my $Scriptsdir=$ENV{'SCRIPTSDIR'};
		my ($ext_Homo,$opt_Homo,$ext_Diso,$opt_Diso)=("","","","");
		if(@{$Ref_Filters}[0] eq "Yes"){ $ext_Homo="_H"; $opt_Homo="--H"; }
		if(@{$Ref_Filters}[1] eq "Yes"){ $ext_Diso="_D"; $opt_Diso="--D"; }

		my $Path_rc=""; my $Nseq="";
		# Create a new Proteome with the ordered regions replaced by 'X'
		my $Sequence=$oldir.$oldname.$ext_Homo.$ext_Diso.".faa";
		if(!(-f $Sequence)){
			# If Proteome does not exists and is a file, it launchs the Script which removes the ordered regions and/or the homologuous regions from the Proteome
			print "[Filtering Sequences]\r";
			`perl $Scriptsdir/Sequence_Filtering.pl --Seq $oldSequence $opt_Homo $opt_Diso --Diso $Disofile --Blast $Blastfile --percId 40 --out $Sequence --thr $THR --logdir $Logdir --nov `;
		}
		print "[Filtering Sequences] => OK\n";
		print "- Filtered Sequence file          :	\"".$Sequence."\" => "; ($Path_rc,$Sequence,$Nseq)=guess($Sequence); 
		if($Sequence !~ /(\.faa|\.fasta|\.fa)$/){ unlink($Sequence); }
		return($Path_rc,$Sequence,$Nseq);
}

sub Filter_Overlap{

	my ($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,$MinOcc,$Ref_AllPaths,$Ref_Filters)=@_;

	my %AllPaths=%{$Ref_AllPaths};
	my $Scriptsdir = $ENV{'SCRIPTSDIR'};
	my $Proteome   = $AllPaths{"Proteome"};   my $Setfile = $AllPaths{"Setfile"}; 
	my $Resultfile = $AllPaths{"Resultfile"}; my $Logdir  = $AllPaths{"Logdir"}; 

	my ($resname,$resdir,$resext)=fileparse($Resultfile,qr{\..*});
	my $resfile = "$resdir/$resname"."_O.res";

	print "[Filtering Overlapping motifs]\r";
	header_over($resfile,$AllPaths{"Setfile_old"},$Setfile,$AllPaths{"Proteome_old"},$Proteome,$MinOcc,$Nseq,$nseq,$Ref_Filters,$Dmin,$Dmax);
	`perl $Scriptsdir/Filter_Overlap.pl --Seq $Setfile --Motifs_file $Resultfile --Size $Smin $Smax --resfile $resfile --logdir $Logdir --nov `;
	footer_over($resfile);
	print "[Filtering Overlapping motifs] => OK\n";

	return($Ref_AllPaths);
}

sub Enumeration {

	my ($Smin,$Smax,$Dmin,$Dmax,$nseq,$MinOcc,$Ref_AllPaths,$Seqfile,$old_Seqfile,$Ref_Filters,$fh,$Ref_times)= @_;

	my %AllPaths=%{$Ref_AllPaths}; my $Resdir=$AllPaths{"Enumdir"}; my $Execdir=$AllPaths{"Execdir"}; 
	my($t0,$t1,$td);


	# my $Setfile=$AllPaths{"Setfile"};
	my (@path_to_Seqfile)   = fileparse($Seqfile,qr{\..*});

	my $S_enu = $Resdir."/$path_to_Seqfile[0]";
	mkpath($S_enu,0,0766) or die "# Cannot create directory \"$S_enu\"  : $! #\n" unless -d $S_enu;

	my ($Homo,$Diso)=("",""); $Homo="_H" unless @{$Ref_Filters}[0] eq "No";	$Diso="_D" unless @{$Ref_Filters}[1] eq "No";

	my @DIR1=($S_enu."/Data_User_Input", $S_enu."/Statistics", $S_enu."/Log");
	foreach(@DIR1){
		mkpath($_,0,0766) or die "# Cannot create directory \"$_\"  : $! #\n" unless -d $_;
	}

	my $New_Seqfile=$S_enu."/Data_User_Input/$path_to_Seqfile[0].faa";
	print STDERR "Copying Input files to the Working directory...(\"".$S_enu."/Data_User_Input/\")\n";
	fcopy($Seqfile,$New_Seqfile) or die "Copy failed: $!";;
#	my $Freqfile=Compute_AA_frequency($New_Seqfile,$AllPaths{"Logdir"}); $AllPaths{"Freqfile"} = $Freqfile;
	# If $Seqfile has been renamed since input
	if($Seqfile !~ /([\.faa])$/){ unlink($Seqfile); }

	my $s=0; 
	my $str="SIZE :   ";
	for($s=$Smin;$s<=$Smax;$s++){
		$str.=" $s "; print STDERR format_string(40,$str," ",".",">")."\r";
		my $DIR=("$S_enu/Statistics/L$s"); mkpath($DIR,0,0766) or die "# Cannot create directory \"$DIR\"  : $! #\n" unless -d $DIR;
		
		my $dmin=$Dmin; if($Dmin > $s){ $dmin=3;  }
		my $dmax=$Dmax; if($Dmax > $s){ $dmax=$s; }

		open($fh,">","$S_enu/Statistics/L$s/Time_usage.txt");
		$t0 = Benchmark->new;
		my @times = times;
		my $resfile="$S_enu/result_K$MinOcc".$Homo.$Diso.".L$s";

#		my $Cmd="$Execdir/Enumeration.exe $New_Seqfile $resfile $nseq $s $MinOcc";
#		print STDERR "\n$Cmd\n";

		if(-e -s -f -T "$resfile"){
			my $footer = Extract_Last_Line($resfile);
			if( $footer ne "# END OF ENUMERATION #" ){
				# If the file exists, is not empty, is text file but the Enumeration is not finished
				# Overwrites the result file and Launchs again the enumeration on the Set
				header_enum($resfile,$old_Seqfile,$New_Seqfile,$MinOcc,$nseq,$Ref_Filters,$dmin,$dmax);
				`$Execdir/Enumeration.exe $New_Seqfile $resfile $nseq $s $MinOcc $dmin $dmax`; # $Freqfile
				footer_enum($resfile);
			}
		}else{
			# If the file doesn't exists, is empty, is not a text file, Launch the Enumeration on the Set
			header_enum($resfile,$old_Seqfile,$New_Seqfile,$MinOcc,$nseq,$Ref_Filters,$dmin,$dmax);
			`$Execdir/Enumeration.exe $New_Seqfile $resfile $nseq $s $MinOcc $dmin $dmax`; # $Freqfile
			footer_enum($resfile);
		}
		# Write time Usage in a file
		$t1 = Benchmark->new;
		$td = timediff($t1, $t0);
		print $fh "Enumeration took: ",timestr($td),"\n";
		#$Ref_times=write_Time_usage(\*$fh,\@times);
		close($fh);
	}
	print STDERR "\n-----------------------------------------\n";
	return($Ref_times,$New_Seqfile);
}

sub Scan_Proteome {

	my ($Smin,$Smax,$Dmin,$Dmax,$Nseq,$Ref_AllPaths,$Proteome,$old_Proteome,$Ref_Filters,$Ref_times,$OnePerSeq)= @_;
	my %AllPaths=%{$Ref_AllPaths}; my $Resdir=$AllPaths{"Enumdir"}; my $Execdir=$AllPaths{"Execdir"};
	my (@path_to_Proteome)  = fileparse($Proteome,qr{\..*});
	my (@path_to_Setfile)   = fileparse($AllPaths{"Setfile"},qr{\..*});
	my($t0,$t1,$td);
	
	my $S_enu = $Resdir."/$path_to_Setfile[0]";
	my $P_enu = $Resdir."/$path_to_Setfile[0]/$path_to_Proteome[0]";
	my $Scan_Log   = $P_enu."/Log";
	my ($Homo,$Diso)=("",""); $Homo="_H" unless @{$Ref_Filters}[0] eq "No"; $Diso="_D" unless @{$Ref_Filters}[1] eq "No";

	print STDERR "Copying Input files to the Working directory...(\"$S_enu/Data_User_Input/\")\n";
	my $New_Proteome="$S_enu/Data_User_Input/$path_to_Proteome[0].faa";
	fcopy($Proteome,$New_Proteome) or die "Copy failed: $!";
	# If $Proteome has been renamed since input
	if($Proteome !~ /(\.faa|\.fasta|\.fa)$/){ unlink($Proteome); }

	my $s=0; 
	my $str="SIZE :   ";
	for($s=$Smin;$s<=$Smax;$s++){
		$str.=" $s "; print STDERR format_string(40,$str," ",".",">")."\r";
		my @DIR=("$P_enu/Statistics/L$s",$Scan_Log); foreach(@DIR){ mkpath($_,0,0766) or die "# Cannot create directory \"$_\" : $! #\n" unless -d $_; }
		my $Set_resfile="$S_enu/result_K$MinOcc".$Homo.$Diso.".L$s";
		my $resfile="$P_enu/result_K$MinOcc".$Homo.$Diso.".L$s";
		open(PROTTIME,">","$P_enu/Statistics/L$s/Time_usage.txt");
		my @times = times;
		$t0 = Benchmark->new;
		my $dmin=$Dmin; if($Dmin > $s){ $dmin=3;  }
		my $dmax=$Dmax; if($Dmax > $s){ $dmax=$s; }

		# Scan the Proteome with the list of enumerated motifs in the Set and count their occurences
#		my $Cmd="$Execdir/Scan_Proteome.exe $New_Proteome $Set_resfile $Nseq $s $resfile";
#		print STDERR "\n$Cmd\n";
		if(-e -s -f -T "$resfile"){
			my $footer = Extract_Last_Line($resfile);
			if( $footer ne "# END OF ENUMERATION #" ){
				header_enum($resfile,$old_Proteome,$New_Proteome,$MinOcc,$Nseq,$Ref_Filters,$dmin,$dmax);
				`$Execdir/Scan_Proteome.exe $New_Proteome $Set_resfile $Nseq $s $resfile $OnePerSeq `;
				footer_enum($resfile);
			}
		}else{
			header_enum($resfile,$old_Proteome,$New_Proteome,$MinOcc,$Nseq,$Ref_Filters,$dmin,$dmax);
			`$Execdir/Scan_Proteome.exe $New_Proteome $Set_resfile $Nseq $s $resfile $OnePerSeq `;
			footer_enum($resfile);
		}
		$t1 = Benchmark->new;
		$td = timediff($t1, $t0);
		print PROTTIME "Proteome Scan took: ",timestr($td),"\n";
#		$Ref_times=write_Time_usage(\*PROTTIME,$Ref_times);
		close(PROTTIME);
	}
	print STDERR "\n-----------------------------------------\n";
	return($Ref_times,$New_Proteome);
}

sub Over_Representation {

	my ($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,$MinOcc,$Ref_AllPaths,$Ref_Filters,$Scan,$Ref_times,$OnePerSeq)= @_;
	my %AllPaths=%{$Ref_AllPaths}; my $Execdir=$AllPaths{"Execdir"}; my $Resdir=$AllPaths{"Overdir"}; my $Enumdir=$AllPaths{"Enumdir"};
	my $Setfile=$AllPaths{"Setfile"};   my (@path_to_Setfile)  = fileparse($Setfile,qr{\..*});  
	my $Proteome=$AllPaths{"Proteome"}; my (@path_to_Proteome) = fileparse($Proteome,qr{\..*}); 
	my($t0,$t1,$td);
	my ($Homo,$Diso)  = ("","");   $Homo="_H" unless @{$Ref_Filters}[0] eq "No";   $Diso="_D" unless @{$Ref_Filters}[1] eq "No";

	my $S_enu     = abs_path($Enumdir."/$path_to_Setfile[0]/");
	my $P_enu     = abs_path($Enumdir."/$path_to_Setfile[0]/$path_to_Proteome[0]/");
	my $OvR_Res   = $Resdir."/$path_to_Setfile[0]";
	my $OvR_Stats = $OvR_Res."/Statistics";
	my $OvR_Log   = $OvR_Res."/Log";

	my $s=0;
	print STDERR "-----------------------------------------\n";
	print STDERR "       STEP OF OVER-REPRESENTATION       \n";
	print STDERR "-----------------------------------------\n";
	my $str="SIZE :   ";
	for($s=$Smin;$s<=$Smax;$s++){
		$str.=" $s "; print STDERR format_string(40,$str," ",".",">")."\r";
		my @DIR = ("$OvR_Stats/L$s","$OvR_Log"); foreach(@DIR){ mkpath($_,0,0766) or die "# Cannot create directory \"$_\" : $! #\n" unless -d $_; }

		my $Set_resfile="$S_enu/result_K$MinOcc".$Homo.$Diso.".L$s";
		my $Proteome_resfile="$P_enu/result_K$MinOcc".$Homo.$Diso.".L$s";
		my $Over_resfile="$OvR_Res/".$path_to_Setfile[0].".L$s";
		my $sorted_Over_resfile="$OvR_Res/sorted_".$path_to_Setfile[0].".L$s";

		my $dmin=$Dmin; if($Dmin > $s){ $dmin=3;  }
		my $dmax=$Dmax; if($Dmax > $s){ $dmax=$s; }

		open(TIME,">","$OvR_Stats/L$s/Time_usage.txt");
		my @times = times;
		$t0 = Benchmark->new;

		if(-e -s -f -T "$sorted_Over_resfile"){
			if(-f $Set_resfile and -f $Proteome_resfile){
				my $footer = Extract_Last_Line($sorted_Over_resfile);
				if( $footer ne "# END OF OVER REPRESENTATION #" ){
					# If the file exists, is not empty, is text file but the Over Representation is not finished
					# Overwrites the result file and Launchs again the Over Representation from the Set and the Proteome
					`$Execdir/Over_Representation.exe $Proteome_resfile $Set_resfile $Over_resfile $s $nseq $Nseq $OnePerSeq`;
				}
			}else{ print "# Something wrong happened during the Enumeration step. #\n"; exit(1); }
		}else{
			# If the file doesn't exists, is empty, is not a text file, Launch the Over Representation from the Set and the Proteome
			`$Execdir/Over_Representation.exe $Proteome_resfile $Set_resfile $Over_resfile $s $nseq $Nseq $OnePerSeq`;
		}

		if(-f $Over_resfile ) {
			my ($Ref_Motifs,$nb_motifs)=load_motifs($Over_resfile);
			my %Motifs = %{$Ref_Motifs};
			header_over($sorted_Over_resfile,$AllPaths{"Setfile_old"},$Setfile,$AllPaths{"Proteome_old"},$Proteome,$MinOcc,$Nseq,$nseq,$Ref_Filters,$dmin,$dmax);
			open(OUT,">>$sorted_Over_resfile");
#			foreach (sort {$Motifs{$a}->{'Pval_adj'} <=> $Motifs{$b}->{'Pval_adj'} || $Motifs{$b}->{'Cdef'} <=> $Motifs{$a}->{'Cdef'} } keys %Motifs) {
#				print OUT "$_\t$Motifs{$_}->{'S_set'}\t$Motifs{$_}->{'S_pop'}\t$Motifs{$_}->{'N_set'}\t$Motifs{$_}->{'N_pop'}\t$Motifs{$_}->{'Pval'}\t$Motifs{$_}->{'nbdegen'}\t$Motifs{$_}->{'Pval_adj'}\t$Motifs{$_}->{'nwild'}\n";
			foreach (sort {$Motifs{$a}->{'Pval'} <=> $Motifs{$b}->{'Pval'} || $Motifs{$b}->{'Cdef'} <=> $Motifs{$a}->{'Cdef'} } keys %Motifs) {
				print OUT "$_\t$Motifs{$_}->{'S_set'}\t$Motifs{$_}->{'S_pop'}\t$Motifs{$_}->{'N_set'}\t$Motifs{$_}->{'N_pop'}\t$Motifs{$_}->{'Pval'}\n";
			}
			close(OUT);
			footer_over($sorted_Over_resfile);
			unlink("$Over_resfile");
		}

		$AllPaths{"Resultfile"}=$sorted_Over_resfile; $AllPaths{"Resultdir"}=$OvR_Res;
		$t1 = Benchmark->new;
		$td = timediff($t1, $t0);
		print TIME "Over Representation took: ",timestr($td),"\n";
#		$Ref_times=write_Time_usage(\*TIME,$Ref_times);	
		close(TIME);
	}
	print STDERR "\n-----------------------------------------\n";
	return(\%AllPaths);
}

if($Setfile ne "" and $Proteome ne ""){ 

	my ($Ref_AllPaths) = File_and_Directory_Managment($Proteome,$Setfile,$Workdir);
	my %AllPaths=%{$Ref_AllPaths};
	my $Scriptsdir = $ENV{'SCRIPTSDIR'};

	my $Reference_Proteome = $ENV{'DATADIR'}."/Yeast_Proteome/Yeast_proteome_5761_IDS.faa";
	my @Filters = ($Homology_Filter,$Disorder_Filter,$Overlap_Filter); my $Ref_Filters=\@Filters;
	# If the non verbous mode is active
	if ($quiet){
		# Redirection of printed messages to a log file ('MotifHound.log') in the Log directory Directory
		open (STDOUT,">".$AllPaths{"Logdir"}."/MotifHound.log") or die "# Unable to open MotifHound.log #\n";
	}

	# MotifHound in ASCII Art
	my $t0 = time;
	my @start=times;
	Logo($t0,$Command);

	# Chech enabled parameters and print a summary of the parameters status
	($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,$Ref_AllPaths,$Ref_Filters,$Scan)=Print_Parameters($Ref_AllPaths,$Reference_Proteome,\@Filters,\@Size,\@Def,$Smax,$Smin,$MinOcc,$Graph,$HTML,$Scan,$OnePerSeq);

	# Compilations of C/C++ Source code
	open(COMPIL,">".$AllPaths{"Logdir"}."/compilation_status.log");
	print COMPIL "- Sources Directory : \"".$ENV{'SRCDIR'}."\" \n";
	print COMPIL "---------------------------------- \n";
	print COMPIL "Compilations of the source codes : \n";
	print COMPIL "---------------------------------- \n";
	Compilation_status("Enumeration",$Ref_AllPaths,\*COMPIL);
	Compilation_status("Over_Representation",$Ref_AllPaths,\*COMPIL);
	if($Scan){ Compilation_status("Scan_Proteome",$Ref_AllPaths,\*COMPIL); }
	close(COMPIL);

	print STDERR "\n"; 
	my $Ref_times="";
	print STDERR "-----------------------------------------\n";
	print STDERR "        STEP OF MOTIF ENUMERATION        \n";
	print STDERR "-----------------------------------------\n";
	%AllPaths=%{$Ref_AllPaths};
	my $fh;
	print STDERR "=>  Enumeration (Setfile Enumeration)  <=\n";
	($Ref_times,$AllPaths{"Setfile"})      = Enumeration($Smin,$Smax,$Dmin,$Dmax,$nseq,$MinOcc,$Ref_AllPaths,$AllPaths{"Setfile"},$AllPaths{"Setfile_old"},$Ref_Filters,$fh,$Ref_times);
	if($Scan == 0){
		print STDERR "=> Enumeration (Proteome Enumeration)  <=\n";
		($Ref_times,$AllPaths{"Proteome"}) = Enumeration($Smin,$Smax,$Dmin,$Dmax,$Nseq,$MinOcc,$Ref_AllPaths,$AllPaths{"Proteome"},$AllPaths{"Proteome_old"},$Ref_Filters,$fh,$Ref_times);
	}elsif($Scan == 1){
		print STDERR "=>     Enumeration (Proteome Scan)     <=\n";
		($Ref_times,$AllPaths{"Proteome"}) = Scan_Proteome($Smin,$Smax,$Dmin,$Dmax,$Nseq,$Ref_AllPaths,$AllPaths{"Proteome"},$AllPaths{"Proteome_old"},$Ref_Filters,$Ref_times,$OnePerSeq);
	}
	print STDERR "\n";
	($Ref_AllPaths) = Over_Representation($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,$MinOcc,\%AllPaths,$Ref_Filters,$Scan,$Ref_times,$OnePerSeq);
	print STDERR "\n";

	if(@{$Ref_Filters}[2] eq "Yes"){
		($Ref_AllPaths)=Filter_Overlap($Smin,$Smax,$Dmin,$Dmax,$nseq,$Nseq,$MinOcc,$Ref_AllPaths,$Ref_Filters);
	}
	print STDERR "\n";

	%AllPaths=%{$Ref_AllPaths};
	my $Proteome   = $AllPaths{"Proteome"};   my $Setfile   = $AllPaths{"Setfile"};
	my $Resultfile = $AllPaths{"Resultfile"}; my $Resultdir = $AllPaths{"Resultdir"};
	my $Logdir     = $AllPaths{"Logdir"};     my $Genefile  = $AllPaths{"Gene_file"}; my $Pfam_annot = $AllPaths{"Pfam_annot"};

	print STDERR "Result Directory : ".$Resultdir."\n";
	my($fname,$dir,$ext)=fileparse($Resultfile,qr{\..*});

	if($HTML){
		print STDERR "Producing HTML Output...(".abs_path($Resultdir."/$fname").".html)\r";
		`$Scriptsdir/Generate_HTML.pl --Seq $Setfile --Proteome $Proteome --Motifs $Resultfile --Gene_annot $Genefile --Pfam_annot $Pfam_annot --Size $Smin $Smax --resfile $Resultdir/$fname.html --logdir $Resultdir/Log`;
		print STDERR "Producing HTML Output...OK                                                                                                                                                                              \n";
	}

	if($Graph){
		print STDERR "Producing Graphical Output...(".abs_path($Resultdir."/$fname").".ps)\r";
		`$Scriptsdir/Draw_motif_image.pl --Seq $Setfile --Proteome $Proteome --Motifs $Resultfile --Gene_annot $Genefile --Pfam_annot $Pfam_annot --resfile $Resultdir/$fname.ps --logdir $Logdir --nov`;
		print STDERR "Producing Graphical Output...OK                                                                                                                                                                         \n";
	}
	my $t1=time;
	print "				   RUN END : ".localtime($t1)."\n";
	print STDERR "===========================================================================================================================\n";
	print STDERR "RUN START    : ".localtime($t0)."\n";
	print STDERR "RUN END      : ".localtime($t1)."\n";
#	write_Time_usage(\*STDERR,\@start);
	print STDERR "The program ran for ", time() - $^T, " seconds.\n";
	print STDERR "===========================================================================================================================\n";
	print STDERR "\n";
}else{
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, read the help :"," ","#",">");
	print "\n".format_string(70,"#        $0                   --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
