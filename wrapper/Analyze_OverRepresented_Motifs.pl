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
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use List::MoreUtils qw(uniq);
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use Data::Dumper;

sub File_verification {
# Test wheter the path given corresponds to an existing file and if is not empty or not readable
# => returns 0 if the file is not empty, readable and if it's an ordinary file
# => returns -1 if the file is empty, not readable or if it's not an ordinary file

# File Test Operators
#-e  file exists
#-f  file is an ordinary file (not a directory or a device)
#-s  file's size greater to 0
#-r  file has the permission of read      (for the user who has run this command)
#-w  file has the permission of write     (for the user who has run this command)
#-x  file has the permission of execution (for the user who has run this command)

	my ($path) = @_;
	if(-e $path and !(-d $path)){
		if("-f -s -r $path"){
			return(0,"OK [Plain File]\n");
		}elsif(!("-f -s -r $path")){
			return(-1,"File is empty or is not readable.\n");
		}
	}else{
		return (-1,"Path does not exist and/or corresponds to a directory.\n");
	}
}

sub format_string {

	my ($len_to_max,$str,$spacer,$end,$side) = @_;
	
	my $j = 0; my $i = 0;

	my $fin =""; my $new_str = "";
	my $len = length($str);
	if( $end ne ""){ $fin = $end; }
		
	if($side eq ">"){
		for($i=$len; $i<$len_to_max;$i+=length($spacer)){
			$str.="$spacer";
		}
	}elsif($side eq "<"){
		for($j=0; $j<($len_to_max-$len);$j+=length($spacer)){
			$new_str.="$spacer";
		}
		$new_str.=$str;
		$str=$new_str;
	}
	$str.="$fin";
	return($str);
}

sub Count_lines($){

	my ($filename) = (@_);
	my $lines = 0; my $buffer;
	open my($fh), '<'.$filename or die "Can't open $filename: $!";
	while( $buffer = <$fh> ) {
		$lines += ( $buffer =~ /\n$/ and $buffer !~ /^#/ );
	}
	close($fh);
	
	return ($lines)
}

sub Count_cdef {
# Count fixed positions (non wildcards) in a string 
	my ($motif) = @_;
	my $flag_cdef=2;
	my @tab = split("",$motif);

	for(my $i=1;$i<$#tab;$i++){
		if ($tab[$i] ne "." and $tab[$i] ne "X"){
			$flag_cdef++;
		}
	}
	return ($flag_cdef);
}

sub Count_wild {
# Count fixed positions (non wildcards) in a string 
	my ($motif) = @_;
	my $flag_wild=0;
	my @tab = split("",$motif);

	for(my $i=1;$i<$#tab;$i++){
		if ($tab[$i] eq "." or $tab[$i] eq "X"){
			$flag_wild++;
		}
	}
	return ($flag_wild);
}

sub load_motifs {
# Load the content of over representation results file
	my ($file,$NR) = (@_);

	my $cpt=0;
	my $motif=""; my $line="";
	my @c;
	my %H=();
	open(IN, "<$file") or print("Could not open $file\n") ;
	while( $line = <IN>){
		if( $line =~ /^[^#]/ and $cpt <= int($NR/5)){
			chop($line);
			@c = split("\t",$line);
			my @pos=split("",$c[0]);
			$cpt++;

			$H{$c[0]}->{"Len"} = length($c[0]);
			$H{$c[0]}->{"Cdef"} = Count_cdef($c[0]);
			$H{$c[0]}->{"S_set"} = $c[1];
			$H{$c[0]}->{"S_pop"} = $c[2];
			$H{$c[0]}->{"N_set"} = $c[3];
			$H{$c[0]}->{"N_pop"} = $c[4];
			$H{$c[0]}->{"Pval"} = $c[5];
#			$H{$c[0]}->{"nbdegen"} = $c[6];
#			$H{$c[0]}->{"Pval_adj"} = $c[7];
			$H{$c[0]}->{"nwild"} = Count_wild($c[0]);
			$H{$c[0]}->{"Pos"}   = \@pos;
		}
	}
	close(IN);
	return (\%H,$cpt);
}

# Silent options
my $help = ""; my $quiet= "";
my $Motifs_file = "";
my $logdir   = "./";
my $Size=3;
my $resfile  = "./Best_Motifs.";

GetOptions(
	"help!"              => \$help,
	"nov|quiet!"         => \$quiet,
	"Motifs=s"           => \$Motifs_file,
	"resfile|out:s"      => \$resfile,
	"logdir:s"           => \$logdir,
	"Size=i"             => \$Size,
	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
-------------------------------------------------------------------------------------------------------------------------------------
     Options         |  Arguments                 |  Description                                                                     
-------------------------------------------------------------------------------------------------------------------------------------
                     |                            |                                                                                  
 *  --Motifs         |  <File directory>          |  Location of the file that contains all the over-represented Motifs              
                     |                            |                                                                                  
 *  --Size           |  <Integer value>           |  Length of the Motifs                                                            
                     |                            |                                                                                  
    --out            |  <File or Path directory>  |  Location or Name of the Results file                                            
                     |                            |                                                                                  
    --logdir         |  <Path directory>          |  Location of the Log Directory                                                   
                     |                            |                                                                                  
    --nov            |                            |  Only print Logo, progression and time usage                                     
                     |                            |                                                                                  
    --help           |                            |  Prints the help description and exit                                            
                     |                            |                                                                                  
-------------------------------------------------------------------------------------------------------------------------------------
(*) Mandatory options                                                                                                                \n";
# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#########################################################################################################################################\n";
	print "	#                         Analyze Over-Represented Motifs (effect of Wildcards on p-value and motif occurences)                         #\n";
	print "	#---------------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 6 parameters :                                                                                             #\n";
	print "	#   1)  --help                                     + Description of Program Usage                                    +  (Optional)      #\n";
	print "	#   2)  --quiet,--nov                              + Non verbous Mode                                                +  (Optional)      #\n";
	print "	#   3)  --Motifs                                   + Over Representation Result File                                 +  (Mandatory)     #\n";
	print "	#   4)  --resfile,--out                            + Location of the Result file                                     +  (Optional)      #\n";
	print "	#   5)  --logdir                                   + Location of the Log Directory                                   +  (Optional)      #\n";
	print "	#   6)  --Size                                     + Motif Length                                                    +  (Optional)      #\n";
	print "	#########################################################################################################################################\n";
	exit(0);
}

if($Motifs_file ne ""){

	# Non verbose mode 
	if ($quiet){
		# Redirection of printed messages to a log file ('Analyze_Output.log')
		open (STDOUT,">$logdir/filter_wildcards.log") or die ("Unable to redirect STDOUT to a file");
	}

	# Existence of the Motif Enumeration file
	my ($Motifsname,$Motifsdir,$Motifsext)  = fileparse(abs_path($Motifs_file),qr{\..*});
	print "- Over Represented Motifs file    :  \"$Motifs_file\" => "; my ($rc_motif,$TEXT) = File_verification($Motifs_file); print $TEXT;
	if($rc_motif != 0){ die "# Something wrong happened when verifying the location of the over-represented Motifs file #\n"; }

	my $nbLines=Count_lines($Motifs_file);
#	print " Number of Motifs : $nbLines.\n";
	my ($Ref_MOTIFS,$nb_motifs)=load_motifs($Motifs_file,$nbLines);
	my %MOTIFS = %{$Ref_MOTIFS};
	my $nb=0; my@tab; my @Motifs_tab;
	my @List_sorted = sort { $MOTIFS{$b}->{"Cdef"} <=> $MOTIFS{$a}->{"Cdef"} || $MOTIFS{$b}->{"S_set"} <=> $MOTIFS{$a}->{"S_set"} || $MOTIFS{$a}->{"S_pop"} <=> $MOTIFS{$b}->{"S_pop"} } (keys %MOTIFS);
#	print Dumper(\@List_sorted);
	foreach my $m (0..$#List_sorted){
		@{$Motifs_tab[$m]} = @{$MOTIFS{$List_sorted[$m]}->{'Pos'}};
#		printf("%3d: %s (%d)\t %d \t %d \t %.5e\n",$nb,$m,$MOTIFS{$m}->{'Cdef'},$MOTIFS{$m}->{'S_set'},$MOTIFS{$m}->{'S_pop'},$MOTIFS{$m}->{'Pval'});
		$nb++;
	}
	print " Number of Motifs : $nb_motifs.\n";

#	my @test=(1,2,3,4,5,6,7,8,9,10);

#	my $ind=0;
#	foreach (@test){
#		if($_ % 2 ==0){
#			splice(@test,$ind,1);  # splice ARRAY, OFFSET, LENGTH, LIST
#		}
#		print $_."\n";
#		$ind++;
#	}

#	print " Motifs tab ".Dumper(\@Motifs_tab);

	$nb=0; my $c=0; my $i=0; my $j=0; my @tmp; my @tmp2;
	my @Final_Motifs; my @Motifs_tab_tmp;
	while ( $#Motifs_tab > -1){
		@tmp=@{shift(@Motifs_tab)};
#	for($i=0; $i <= $#Motifs_tab ; $i++){
#		print "$nb".@{$_}[0]."\n";
		my $L = $#Motifs_tab+1;
		$j=0;# $cpt=1;
		@{$Final_Motifs[$#Final_Motifs+1]} = @tmp;
		for($j=0; $j < $L ; $j++){
			my @tmp2=@{shift(@Motifs_tab)};
#		for($j=0; $j <= $#Motifs_tab ; $j++){
			my $w1 = $MOTIFS{join('',@tmp)}->{"nwild"};
			my $w2 = $MOTIFS{join('',@tmp2)}->{"nwild"};
			if(join('',@tmp) ne join('',@tmp2) and ($w1+$w2 != 0) ){
				if($tmp[0] eq $tmp2[0] and $tmp[$Size-1] eq $tmp2[$Size-1] ){
					my $occ1 = $MOTIFS{join('',@tmp)}->{"S_set"};
					my $occ2 = $MOTIFS{join('',@tmp2)}->{"S_set"};
					my $occCdef1 = $occ1 - $w1;
					my $occCdef2 = $occ2 - $w2;
#					print " $i @tmp [$occ1 - $occCdef1] #*# $j @tmp2  [$occ2 - $occCdef2]\n";
					if( $occCdef1 <= $occCdef2 ){
						@{$Motifs_tab_tmp[$#Motifs_tab_tmp+1]} = @tmp2;
	#					splice(@Motifs_tab,$j,1);  # splice ARRAY, OFFSET, LENGTH, LIST
					}
				}else{
					@{$Motifs_tab_tmp[$#Motifs_tab_tmp+1]} = @tmp2;
				}
			}else{
					@{$Motifs_tab_tmp[$#Motifs_tab_tmp+1]} = @tmp2;
			}
		}
#		print Dumper(@Motifs_tab_tmp);
		@Motifs_tab=@Motifs_tab_tmp;
		@Motifs_tab_tmp=();
#		print Dumper(@Motifs_tab);
		$i++;
	}

	$nb=0; my $M; my %Final;
	foreach (uniq(@Final_Motifs)){
		$M=join('',@{$_});
		printf("%3d: %s (%d)\t %d \t %d \t %.5e\n",$nb+1,$M,$MOTIFS{$M}->{'Cdef'},$MOTIFS{$M}->{'S_set'},$MOTIFS{$M}->{'S_pop'},$MOTIFS{$M}->{'Pval'});
		$Final{$M}{'L'}    = $MOTIFS{$M}->{'Len'} ;
		$Final{$M}{'C'}    = $MOTIFS{$M}->{'Cdef'} ;
		$Final{$M}{'k'}    = $MOTIFS{$M}->{'S_set'} ;
		$Final{$M}{'m'}    = $MOTIFS{$M}->{'S_pop'} ;
		$Final{$M}{'n'}    = $MOTIFS{$M}->{'N_set'} ;
		$Final{$M}{'N'}    = $MOTIFS{$M}->{'N_pop'} ;
		$Final{$M}{'P'}    = $MOTIFS{$M}->{'Pval'} ;
#		$Final{$M}{'d'}    = $MOTIFS{$M}->{'nbdegen'} ;
#		$Final{$M}{'Padj'} = $MOTIFS{$M}->{'Pval_adj'} ;
#		$Final{$M}{'w'}    = $MOTIFS{$M}->{'nwild'} ;

		$nb++;
	}
	print "Nb of motifs : $nb_motifs (before) $nb (after)\n";
	
	my($resname,$resdir,$resext) = fileparse($resfile,qr{\..*});
#	print abs_path($resdir).$resname."L$Size\n";
	open(OUT,">".abs_path($resdir)."/$resname".".L$Size.tmp");

#	foreach (sort {$Final{$a}->{'Padj'} <=> $Final{$b}->{'Padj'} || $Final{$a}->{'P'} <=> $Final{$b}->{'P'} || $Final{$b}->{'C'} <=> $Final{$a}->{'C'} } keys %Final) {
	foreach (sort {$Final{$a}->{'P'} <=> $Final{$b}->{'P'} || $Final{$b}->{'C'} <=> $Final{$a}->{'C'} } keys %Final) {
#		print $Final{$_}->{'P'}." $_ ".($Final{$_}->{'C'})." ".$Final{$_}->{'k'}." ".$Final{$_}->{'m'}."\n";
		print OUT "$_\t$Final{$_}->{'k'}\t$Final{$_}->{'m'}\t$Final{$_}->{'n'}\t$Final{$_}->{'N'}\t$Final{$_}->{'P'}\n";
#		print OUT "$_\t$Final{$_}->{'k'}\t$Final{$_}->{'m'}\t$Final{$_}->{'n'}\t$Final{$_}->{'N'}\t$Final{$_}->{'P'}\t$Final{$_}->{'d'}\t$Final{$_}->{'Padj'}\t$Final{$_}->{'w'}\n";
	}
	close(OUT);
#			printf("%3d: %s (%d)\t %d \t %d \t %.5e\n",$nb,$m,$MOTIFS{$m}->{'Cdef'},$MOTIFS{$m}->{'S_set'},$MOTIFS{$m}->{'S_pop'},$MOTIFS{$m}->{'Pval'});
#			@Motifs_tab=sort { $MOTIFS{$b}->{"Cdef"} <=> $MOTIFS{$a}->{"Cdef"} } (keys %MOTIFS);
}else{
	my ($bin_name,$bin_dir) = fileparse($0);
	print format_string(70,"#","#","#\n",">");
	print format_string(70,"#  For more details, read the help :"," ","#\n",">");
	print format_string(70,"#    $bin_name    --help"," ","#\n",">");
	print format_string(70,"#","#","#\n\n\n",">");
	exit(0);
}
