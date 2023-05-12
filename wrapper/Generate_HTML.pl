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
	$ENV{'DATADIR'}=abs_path($ENV{'SCRIPTSDIR'}."/../Data");
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use List::MoreUtils qw(uniq);
use File_Utils;
use String_Print_Utils;
use Data::Dumper;
sub Font_Color_motifs($){

	my ($motif)=@_;
	
	my $font_motif="";
	my @splitted_motif = split('',$motif);
	
	foreach my $pos (@splitted_motif){
		if($pos eq '.'){
			$font_motif.="<font size=2 color=green>$pos</font>";
		}elsif($pos ne '.'){
			$font_motif.="<font size=2 color=red>$pos</font>";
		}
	}

	return($font_motif);
}

sub which_min {

	my ($Ref_tab) = @_;
	my @tab = @{$Ref_tab};
	my $ind = 0; my $min=0;

	foreach $ind (1..$#tab) {
		$min=$ind if $tab[$ind]<$tab[$min];
	}

	return ($min);
}

sub which {

	my ($Ref_tab,$val) = @_;
	my @tab = @{$Ref_tab}; my $ind = 0;
	my @list_ind;
	foreach $ind (0..$#tab) {
		push(@list_ind,$ind) if $tab[$ind]==$val;
	}
	return (@list_ind);
}

sub log10 {
	my $n = shift;
	return (log($n)/log(10));
}

sub Filter_Wildcards {
	my ($Ref_Motifs,$s,$scripts,$logdir) = @_;
	`$scripts/Analyze_OverRepresented_Motifs.pl --Motifs $Ref_Motifs --out $Ref_Motifs --logdir $logdir --Size $s --nov`;
	return("$Ref_Motifs.tmp");
}

sub load_all_motifs_in_Hash {

	my ($Motifs_file,$Smin,$Smax,$scripts,$logdir,$msg)=@_;

	my %All_MOTIFS=();
	my $nb_tot_motifs=0; my $s=3;
	my $nb_motifs_cutoff=100;
	my $nb_motifs=0; my $nb=0; my $perc=0;
	my ($fname,$dir,$ext)  = fileparse($Motifs_file,qr{\..*});
	print "\n    o Loads all motifs (from length $Smin to $Smax)...";
	for($s=$Smin; $s <= $Smax ; $s++){
		$nb_motifs=0;
		print STDERR "$msg [Loading Motif of Length $s - ".sprintf("%2d",$perc)."%]\r"; $perc=int((($s-$Smin+1)/($Smax-$Smin))*100);
		my $Ref_MOTIFS = Filter_Wildcards("$dir/$fname".".L$s",$s,$scripts,$logdir);
		my ($Ref_HMOTIFS,$nb) = load_motifs($Ref_MOTIFS); unlink("$dir/$fname".".L$s.tmp");
		my %MOTIFS = %{$Ref_HMOTIFS};
		my @Mots = sort { $MOTIFS{$a}->{"Pval"} <=> $MOTIFS{$b}->{"Pval"} || $MOTIFS{$b}->{"Cdef"} <=> $MOTIFS{$a}->{"Cdef"}} (keys (%MOTIFS));
		foreach my $motif ( @Mots ){
			if( $nb_motifs < $nb_motifs_cutoff){
				$All_MOTIFS{$motif}->{"Len"}= $MOTIFS{$motif}->{"Len"};
				$All_MOTIFS{$motif}->{"Cdef"}= $MOTIFS{$motif}->{"Cdef"};
				$All_MOTIFS{$motif}->{"S_set"}= $MOTIFS{$motif}->{"S_set"};
				$All_MOTIFS{$motif}->{"S_pop"}=$MOTIFS{$motif}->{"S_pop"};
				$All_MOTIFS{$motif}->{"N_set"}= $MOTIFS{$motif}->{"N_set"};
				$All_MOTIFS{$motif}->{"N_pop"}= $MOTIFS{$motif}->{"N_pop"};
				$All_MOTIFS{$motif}->{"Pval"}=$MOTIFS{$motif}->{"Pval"};
#				$All_MOTIFS{$motif}->{"nbdegen"}=$MOTIFS{$motif}->{"nbdegen"};
#				$All_MOTIFS{$motif}->{"Pval_adj"}=$MOTIFS{$motif}->{"Pval_adj"};
#				$All_MOTIFS{$motif}->{"nwild"}=$MOTIFS{$motif}->{"nwild"};
				$nb_tot_motifs++;	$nb_motifs++;
			}
		}
	}
	if( $nb_tot_motifs == 0 ){ print STDERR "# NO MOTIFS FOUND ! #\n-------------------\nPossible reasons :\n-------------------\n- The set might not contain enough sequences\n- The minimum number of occurence was set too high.\n- If sequence filtering was activated, the set might contain only homologuous and/or only ordered proteins.\n"; exit(-1);}
	return(\%All_MOTIFS,$nb_tot_motifs);
}

sub load_Evolution {

	my ($Evofile,$Seq_file)=@_;
	my %Evo=();
	my ($fname,$dir,$ext)  = fileparse($Evofile,qr{\..*});
#	print STDERR "EVO $Evofile"; #<STDIN>;
	unless (open(EVO,"$Evofile")){		print "Unable to open Evolution file ($Evofile)\n"; 	return(\%Evo);		}
	my $SEQUENCES = Load_FASTA($Seq_file);

	my $line; my $i=0;
	my @c;
	while($line=<EVO>){
		if($i % 1000000 == 0){ print "  o Reading Data...[Evolutionary Data]  $i    \r"; }
		chomp($line);
		@c=split(" ",$line);
		if( exists $SEQUENCES->{$c[0]} ){
			if( length($SEQUENCES->{$c[0]}->{'Seq'}) == $c[1] ) {
				$Evo{$c[0]}->{"nAA"}=$c[1];
				$Evo{$c[0]}->{"min"}=$c[3];
				$Evo{$c[0]}->{"max"}=$c[4];
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"AA"}=$c[6];
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"Score"}=$c[7]+abs($c[3]);
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"Conf"}=$c[8];
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"Std"}=$c[9];
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"AlignedAA"}=$c[10];
				$Evo{$c[0]}->{"pos"}->{$c[5]-1}->{"Species"}=$c[11];
			}else{
				print STDERR "ERROR EVOLUTION DATA SEQUENCE NOT IDENTICAL TO FASTA SEQUENCE\n";
				my @seqtmp=split('',$SEQUENCES->{$c[0]}->{'Seq'});
				print STDERR " LEN SEQ FASTA  ".length($SEQUENCES->{$c[0]}->{'Seq'})."  ".($#seqtmp+1)." LEN SEQ EVO ".$c[1]."\n";
				print STDERR ">$c[0]\n".$SEQUENCES->{$c[0]}->{'Seq'}."\n";
				exit(0);
			}
		}
		$i++;
	}
	close(EVO);
	return(\%Evo);
}


# Silent options
my $help = ""; my $quiet= "";
my $Seq_file = ""; my $Motifs_file = ""; my $Proteome = "";
my $resfile  = "./Best_Motifs.html";
my $logdir   = "./";
my $Pfam_annot=abs_path($ENV{'DATADIR'}."/Domains/Yeast_Pfam_annotations.txt");
my $Superfamily_annot=abs_path($ENV{'DATADIR'}."/Domains/Saccharomyces_cerevisiae_SUPERFAMILY_domains.txt");
my $Gene_annot=abs_path($ENV{'DATADIR'}."/Genes/YEAST.data");
my $Evofile="";
my @Size=(3,10);

GetOptions(
	"help!"              => \$help,
	"nov|quiet!"         => \$quiet,
	"Seq=s"              => \$Seq_file,
	"Motifs=s"           => \$Motifs_file,
	"Proteome=s"         => \$Proteome,
	"resfile|out:s"      => \$resfile,
	"logdir:s"           => \$logdir,
	"Size=i{2}"          => \@Size,
	"Gene_annot:s"       => \$Gene_annot,
	"Pfam_annot:s"       => \$Pfam_annot,
	"Superfam_annot:s"   => \$Superfamily_annot,
	"Evo_file:s"         => \$Evofile
	)
or die "=========================\n /!\\ Incorrect Usage /!\\\n=========================\n
USAGE: ".$0." [Options]
----------------------------------------------------------------------------------------------------------------------------------------------
     Options         |  Arguments                          |  Description                                                                     
----------------------------------------------------------------------------------------------------------------------------------------------
 *  --Seq            |  <File directory>                   |  Location of the Set file                                                        
                     |                                     |                                                                                  
 *  --Motifs         |  <File directory>                   |  Location of the file that contains all the over-represented Motifs              
                     |                                     |                                                                                  
 *  --Proteome       |  <File directory>                   |  Location of the file that contains the Proteome sequences which should contain  
                     |                                     |   the sequences of the Setfile                                                   
                     |                                     |                                                                                  
    --Size           |  <Integer value>   <Integer value>  |  Minimum length and Maximum length of the Motifs                                 
                     |                                     |                                                                                  
    --out            |  <File or Path directory>           |  Location or Name of the Results file                                            
                     |                                     |                                                                                  
    --logdir         |  <Path directory>                   |  Location of the Log Directory                                                   
                     |                                     |                                                                                  
    --Pfam_annot     |  <File directory>                   |  Domains Annotation file obtained from PFAM Database                             
    --Superfam_annot |  <File directory>                   |  Domains Annotation file obtained from SUPERFAMILY Database                      
    --Gene_annot     |  <File directory>                   |  File containing Names and Descriptions of genes                                 
    --Evo_file       |  <File directory>                   |  Evolution Rate per site for protein orthologuous to at least 12 mammals species 
                     |                                     |                                                                                  
    --nov            |                                     |  Only print Logo, progression and time usage                                     
                     |                                     |                                                                                  
    --help           |                                     |  Prints the help description and exit                                            
                     |                                     |                                                                                  
----------------------------------------------------------------------------------------------------------------------------------------------
(*) Mandatory options                                                                                                                         \n";
# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#########################################################################################################################################################\n";
	print "	#                                      Generate HTML - Generating HTML output for HYPERMOTIF -                                                          #\n";
	print "	#-------------------------------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 11 parameters :                                                                                                            #\n";
	print "	#   1)  --help                                     + Description of Program Usage                                                    +  (Optional)      #\n";
	print "	#   2)  --quiet,--nov                              + Non verbous Mode                                                                +  (Optional)      #\n";
	print "	#   3)  --Seq                                      + File of Proteins Sequences of Interest (FASTA)                                  +  (Mandatory)     #\n";
	print "	#   4)  --Proteome                                 + Fasta Formatted File containing Proteome Sequences                              +  (Mandatory)     #\n";
	print "	#   5)  --Motifs                                   + Over Representation Result File                                                 +  (Mandatory)     #\n";
	print "	#   6)  --resfile,--out                            + Location of the Result file                                                     +  (Optional)      #\n";
	print "	#   7)  --logdir                                   + Location of the Log Directory                                                   +  (Optional)      #\n";
	print "	#   8)  --Size        (From 3 to 10 by default)    + Min and Max Motifs Sizes (2 values required between 3 to 12 AA)                 +  (Optional)      #\n";
	print "	#   9)  --Pfam_annot     (by default for yeast)    + Pfam domains annotations file for S. Cerevisiae                                 +  (Optional)      #\n";
	print "	#  10)  --Superfam_annot (by default for yeast)    + Superfamily domains annotations file for S. Cerevisiae                          +  (Optional)      #\n";
	print "	#  11)  --Gene_annot     (by default for yeast)    + File containing Names and Descriptions of S. Cerevisiae genes                   +  (Optional)      #\n";
	print "	#  12)  --Evo_file       (by default for human)    + Evolution Rate per site for protein orthologuous to at least 12 mammals species +  (Optional)      #\n";
	print "	#########################################################################################################################################################\n";
	exit(0);
}

if($Seq_file ne "" and $Motifs_file ne "" and $Proteome ne ""){ 

	# Non verbose mode 
	if ($quiet){
		# Redirection of printed messages to a log file ('HTML_Output.log')
		open (STDOUT,">$logdir/HTML_Output.log") or die ("Unable to redirect STDOUT to a file");
	}

	# Existence of the Setfile
	my ($fseqname,$fseqdir,$fseqext)  = fileparse(abs_path($Seq_file),qr{\..*});
	print  "- Set file                        :  \"$Seq_file\" => "; my ($rc_seq,$TEXT) = File_verification($Seq_file); print  $TEXT;
	if($rc_seq != 0){ die "# Something wrong happened when verifying the location of the Setfile #\n"; }
	print "  o Reading Data...[Set Sequences]                                      \n";
	my $SEQUENCES = Load_FASTA($Seq_file);

	# Existence of the Proteome file
	my ($protname,$protdir,$protext)  = fileparse(abs_path($Proteome),qr{\..*});
	print  "- Proteome Sequence File          :  \"$Proteome\" => "; my ($rc_prot,$TEXT2) = File_verification($Proteome); print $TEXT2;
	if($rc_prot != 0){ die "# Something wrong happened when verifying the location of the Proteome file #\n"; }
	print "  o Reading Data...[Proteome Sequences]                                 \n";
	my $PROTEOME = Load_FASTA($Proteome);

	# Existence of the Motif Enumeration file
	my ($Motifsname,$Motifsdir,$Motifsext)  = fileparse(abs_path($Motifs_file),qr{\..*});
	print  "- Over Represented Motifs file    :  \"$Motifs_file\" => "; my ($rc_motif,$TEXT3) = File_verification($Motifs_file); print $TEXT3;
	if($rc_motif != 0){ die "# Something wrong happened when verifying the location of the over-represented Motifs file #\n"; }
	my ($fout,$dir,$ext)  = fileparse($resfile,qr{\..*});

	print "  o Reading Data...[Overrepresented Motifs]                             \n";
	my ($Ref_All_MOTIFS,$nb_motifs)=load_all_motifs_in_Hash($Motifs_file,min(@Size),max(@Size),$ENV{'SCRIPTSDIR'},$logdir,"Producing HTML Output...($dir"."$fout.html)");
	my %All_MOTIFS = %{$Ref_All_MOTIFS};


	open(HTML_all,">".$dir.$fout.".html");
	my $header = "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n<HTML xmlns=\"http://www.w3.org/1999/xhtml\" lang=\"fr\" xml:lang=\"fr\">\n<head>\n<title> Over Represented Motifs </title>\n<meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\" />\n<script src=\"".$ENV{'SCRIPTSDIR'}."/sorttable.js\"></script>\n<link rel=\"stylesheet\" type=\"text/css\" href=\"".$ENV{'SCRIPTSDIR'}."/menu.css\" />\n</head>\n\n<body bgcolor=#fffccc  margin:10px 20px 30px 10px;>\n<div class=\"titre\"><font size=5><b> Set File : <i> $fseqname </i></b></font><BR><font size=2>Location : $fseqdir</font></div>\n<div class=\"contenu\">\n";
	my $table_header = "<table style=\" position:absolute; left:5px; right:5px; border-width:1px; border-style:dashed; border-color:black;\" class=\"sortable\">\n<th><div class=\"cell-en-tete1\"><b>#</b></div></th>\n<th><div class=\"cell-en-tete2\"><b>MOTIF</b></div></th>\n<th><div class=\"cell-en-tete3\"><b># Set</b></div></th>\n<th><div class=\"cell-en-tete4\"><b>ID/Gene <BR>Set</b></div></th>\n<th><div class=\"cell-en-tete3\"><b># Proteome</b></div></th>\n<th><div class=\"cell-en-tete4\"><b>ID/Gene Proteome</b></div></th>\n<th><div class=\"cell-en-tete5\"><b>P-value</b></div>\n<th><div class=\"cell-en-tete6\"><b>Domains</b></div></th>\n<th><div class=\"cell-en-tete6\"><b>Motifs matches (Set)</b></div></th>\n<th><div class=\"cell-en-tete6\"><b>Evolution per Occurences</b></div></th>\n<th><div class=\"cell-en-tete6\"><b>Average Evolution</b></div></th>\n<th><div class=\"cell-en-tete6\"><b>Ratio Evolution</b></div></th>\n";
	print HTML_all $header.$table_header;


#	print "$ENV{DATADIR}/Genes/YEAST.data";
	print "  o Reading Data...                                                     \n";
	print "  o Reading Data...[Gene Annotation]                                    \n";
	my %Data_ID = %{load_Gene_ID_desc($Gene_annot)};
	print "  o Reading Data...[Pfam Domains]                                       \n";
	my %Data_domain = %{load_PFAM_annot($Pfam_annot)};
	print "  o Reading Data...[Evolutionary Data]                                  \r";
	my %Data_Evo  = %{load_Evolution($Evofile,$Seq_file)};
	print "  o Reading Data...DONE                                                 \n";
#	my $Ref_Sharing_domains = load_Sharing_PFAM_Domains($ENV{'DATADIR'}."/Domains/IDs_sharing_same_Pfam_domains.txt",$Pfam_annot,$ENV{'SCRIPTSDIR'});
#	my %Data_sharing_domain  = %{$Ref_Sharing_domains};

	my $i=1; my $cpt=1; my $NSeq=0;
	my $i_motif=1; my $left_motif=$nb_motifs; my $i_ID=1;
	my $lien_domain=""; my $lien_ID="";

	my $Uniprot_link = "href=\"http://www.uniprot.org/uniprot/";
	my $SGD_link = "href=\"http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=";
	my $PFAM_link="http://pfam.sanger.ac.uk/search/jump?entry=";

	print "    o Search motifs in Proteome Sequences...\n";
	$i_motif=1;
	my @IDS_PROT = keys %{$PROTEOME};
	foreach my $motif (sort {$All_MOTIFS{$b}->{"Pval"} <=> $All_MOTIFS{$a}->{"Pval"}} (keys (%All_MOTIFS))){
		$i_ID=1; my %ID_cl=();
		foreach my $ID (@IDS_PROT){
#			if( $i_motif % 50 == 0 and $i_ID == ($#IDS_PROT+1) ){ print "\r        Motif ($i_motif/$left_motif)"; }
			$i_ID++;
			if(exists $PROTEOME->{$ID}){
				my $seq = $PROTEOME->{$ID}->{"Seq"};
				if( $seq =~ /$motif/){
				my @start=@{match_all_start_positions($motif,$seq)}; my @end=@{match_all_end_positions($motif,$seq)};
				($ID_cl{$ID}->{"start"},$ID_cl{$ID}->{"end"})=(join(";",@start),join(";",@end))
				}
			}
		}
		$All_MOTIFS{$motif}->{"ID_Prot"}={%ID_cl};
		$i_motif++;
	}
	print "\n\n";
	
	$i_motif=1; my $Gene_name=""; my $perc=0;
	print "    o Search motifs in Set Sequences...\n";
	my @IDS = keys %{$SEQUENCES}; my @ALL_MOTIFS = keys (%All_MOTIFS);
	foreach my $motif (sort {$All_MOTIFS{$b}->{"Pval"} <=> $All_MOTIFS{$a}->{"Pval"}} (@ALL_MOTIFS)){
		$perc=int(($i_motif/($#ALL_MOTIFS+1))*100);
		if($perc % 10 == 0){ print STDERR "Producing HTML Output...($dir"."$fout.html) [Loading Data in a table - ".sprintf("%2d",$perc)."%]\r";}
		$i_ID=0; my $i_motif_ID=1; my %ID_set_cl=(); my %all_dom_cl=();
		foreach my $ID (@IDS){
			$i_ID++; my %dom_cl=();
#			if( $i_motif % 50 == 0 and $i_ID == ($#IDS+1) ){ print "\r        Motif ($i_motif/$left_motif)"; }
			my $seq = $SEQUENCES->{$ID}->{"Seq"};
#				print STDERR "\n".$motif."\n".$ID."\n".$seq."\n";
			if( $seq =~ /$motif/){
				my @SEQ=split('',$seq);
				my $aa_pos=0; my $aa_dis=0; my $sum=0;
				if(exists($Data_Evo{$ID})){
					foreach my $aa ( @SEQ ){
						if($aa ne "X"){
							if( $aa eq $Data_Evo{$ID}->{"pos"}->{$aa_pos}->{"AA"} ){
								$sum+=$Data_Evo{$ID}->{"pos"}->{$aa_pos}->{"Score"};
								$aa_dis++;
							}else{
								print STDERR  ($aa_pos)." ".$Data_Evo{$ID}->{"pos"}->{$aa_pos}->{'AA'}." ".$aa." | SEQUENCE DISORDER Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID\n\n";
							}
						}
						$aa_pos++;
					}
					$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"disorder"}->{"mean"}=$sum/$aa_dis;
					$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"disorder"}->{"nAA"}=$aa_dis;
				}
				my @start=@{match_all_start_positions($motif,$seq)}; my @end=@{match_all_end_positions($motif,$seq)};
				($ID_set_cl{$ID}->{"start"},$ID_set_cl{$ID}->{"end"})=(join(";",@start),join(";",@end));
				$Gene_name=$ID; if(exists $Data_ID{$ID}){ $Gene_name=$Data_ID{$ID}->{"Gene"}; }
				$ID_set_cl{$ID}->{"Gene_name"} = $Gene_name;

#				print STDERR  "        o Calculating evolution score...\r";
				my $position=0;
				for($position=0;$position<$#start+1;$position++){
					my $start=$start[$position]; my $end=$end[$position];
					my $len_C_before=10;
					my $C_before=$start-10; if($C_before < 0){ $len_C_before=$len_C_before-abs($C_before); $C_before=0;}

#					print STDERR " $ID $start $end $motif BEF:".substr($seq,$start-10,10)." MOT: ".substr($seq,$start,$end-$start)." AFT: ".substr($seq,$end,10)."\n";
					if(exists($Data_Evo{$ID})){
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"Evo"}="";

						my @pos_motif=split('',$motif);
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"motif"}->{"def"}=0; $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"motif"}->{"wild"}=0;
						my $ind_p=0;
	#					print STDERR "###$motif###\n";
						foreach my $p (@pos_motif){
							if( $p ne '.'){
								if( $p eq $Data_Evo{$ID}->{"pos"}->{$start+$ind_p}->{"AA"} ){
									$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"motif"}->{"def"}+=$Data_Evo{$ID}->{"pos"}->{$start+$ind_p}->{"Score"};
								}else{ 
									print STDERR  $start+$ind_p." ".$Data_Evo{$ID}->{"pos"}->{$start+$ind_p}->{'AA'}." ".$p." | MOTIF Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID pos_start $start pos_end $end\n";
								}
							}else{
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"motif"}->{"wild"}+=$Data_Evo{$ID}->{"pos"}->{$start+1+$ind_p}->{"Score"};
							}
						
	#						print STDERR $start+1+$ind_p." ".$p."  ".$Data_Evo{$ID}->{"pos"}->{$start+1+$ind_p}->{"Score"}." ".$Data_Evo{$ID}->{"pos"}->{$start+1+$ind_p}->{'AA'}." WILD ".$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"Evolution"}->{"motif_wild"};
	#						print STDERR " DEF ".$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"Evolution"}->{"motif_def"}."\n";
						
							$ind_p++;
	#						<STDIN>;
						}
					
						

						my @pos_contextA=split('',substr($seq,$end,10));
						my $Cont_before=$C_before; if($Cont_before < 0){ $Cont_before=0; }
						my @pos_contextB=split('',substr($seq,$C_before,$len_C_before));

						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"before"}->{"Score"}=0; $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"after"}->{"Score"}=0;
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"before"}->{"Score"}=0; $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"after"}->{"Score"}=0;

						$ind_p=0;
						foreach my $p (@pos_contextB){
							if( $p eq "X"){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"before"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$Cont_before+$ind_p}->{"Score"};
							}elsif( $p eq $Data_Evo{$ID}->{"pos"}->{$Cont_before+$ind_p}->{"AA"} ){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"before"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$Cont_before+$ind_p}->{"Score"};
							}else{
								print STDERR Dumper(@pos_contextB);
								print STDERR Dumper($All_MOTIFS{$motif});
								print STDERR Dumper(@pos_contextA);
								print STDERR Dumper($seq); 
								print STDERR  $Cont_before+$ind_p." ".$Data_Evo{$ID}->{"pos"}->{$Cont_before+$ind_p}->{'AA'}." ".$p." | CONTEXT BEFORE Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID pos_start $start pos_end $end\n";
								exit(0);
							}
							$ind_p++;
						}
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"before"}->{"nAA"}=$ind_p;


						$ind_p=0;
						foreach my $p (@pos_contextA){
							if( $p eq "X"){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"after"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"Score"};
							}elsif( $p eq $Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"AA"} ){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"after"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"Score"};
							}else{
								print STDERR  $end+$ind_p." ".$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{'AA'}." ".$p." | CONTEXT AFTER Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID pos_start $start pos_end $end\n";
							}
							$ind_p++;
						}
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"context"}->{"disorder"}->{"after"}->{"nAA"}=$ind_p;

						# FINDING THE DISORDER STRETCH BEFORE AND AFTER
						my $pos_stretch=0; my $stretchB="";  my $stretchA="";
						for($pos_stretch=$start; $pos_stretch > 0; $pos_stretch--){
							if( substr($seq,$pos_stretch-1,1) ne 'X' ){
#								print STDERR "STRETCH BEFORE pos $pos_stretch AA ".substr($seq,$pos_stretch-1,1)."\n";
								$stretchB.=substr($seq,$pos_stretch-1,1);
							}else{
								last; 
							}
						}

						for($pos_stretch=$end; $pos_stretch < length($seq)+1; $pos_stretch++){
							if( substr($seq,$pos_stretch,1) ne 'X' ){
								$stretchA.=substr($seq,$pos_stretch,1);
							}else{
								last; 
							}
						}

						my @pos_stretchB=split('',reverse($stretchB));
						my @pos_stretchA=split('',$stretchA);


						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"before"}->{"seq"}=reverse($stretchB); $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"after"}->{"seq"}=$stretchA;

						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"before"}->{"Score"}=0; $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"after"}->{"Score"}=0;
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"before"}->{"Score"}=0; $All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"after"}->{"Score"}=0;

						$ind_p=0;
						foreach my $p (@pos_stretchB){
							if( $p eq "X"){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"before"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"Score"};
							}elsif( $p eq $Data_Evo{$ID}->{"pos"}->{$start-length($stretchB)+$ind_p}->{"AA"} ){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"before"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$start-length($stretchB)+$ind_p}->{"Score"};
							}else{
								print STDERR  $start-length($stretchB)+$ind_p." ".$Data_Evo{$ID}->{"pos"}->{$start-length($stretchB)+$ind_p}->{'AA'}." ".$p." | STRETCH BEFORE Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID pos_start $start pos_end $end\n\n";
							}
							$ind_p++;
						}
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"before"}->{"nAA"}=$ind_p;

						$ind_p=0;
						foreach my $p (@pos_stretchA){
							if( $p eq "X"){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"after"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"Score"};
							}elsif( $p eq $Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"AA"} ){
								$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"after"}->{"Score"}+=$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{"Score"};
							}else{
								print STDERR  $end+$ind_p." ".$Data_Evo{$ID}->{"pos"}->{$end+$ind_p}->{'AA'}." ".$p." | STRETCH AFTER Not the same amino acid (SHIFT IN SEQUENCE POSITIONS ???) M $motif SET $fseqname SEQ $ID pos_start $start pos_end $end\n\n";
							}
							$ind_p++;
						}
						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"stretch"}->{"disorder"}->{"after"}->{"nAA"}=$ind_p;


						$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"occ"}->{$start}->{"Evo"} = "<font size=2 color=blue> >$ID ( #AA : ".$Data_Evo{$ID}->{'nAA'}." #AA(disorder) : ".$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"disorder"}->{"nAA"}.") mean(disorder) :".sprintf("%.3f",$All_MOTIFS{$motif}->{"Evo_ID"}->{$ID}->{"disorder"}->{"mean"})." <BR></font>\n"; #<BR> min(Evolution)-max(Evolution) : ".$Data_Evo{$ID}->{'min'}."-".$Data_Evo{$ID}->{'max'}."<BR></font>\n";
					}

					$ID_set_cl{$ID}->{"Seq"}->{$start} = ">$ID <font size=4 color=blue>".$ID_set_cl{$ID}->{"Gene_name"}." (".length($seq)." res.) </font>\n<BR><b>&nbsp;&nbsp;&nbsp;<font size=4 color=purple>".($start+1)."</font></b>...".substr($seq,$C_before,$len_C_before)."<i><b><font size=4 color=purple>".substr($seq,$start,$end-$start)."</font></b></i>".substr($seq,$end,10)."...<b><font size=4 color=purple>".($end+1)."</font></b><BR>\n";
					if(exists $Data_domain{$ID}){
						foreach my $id ( keys (%{$Data_domain{$ID}})){
							foreach my $cpt ( keys (%{$Data_domain{$ID}->{$id}})){
								my $dom_start = $Data_domain{$ID}->{$id}->{$cpt}->{"ali_start"}; my $dom_end = $Data_domain{$ID}->{$id}->{$cpt}->{"ali_end"};
								if(exists $all_dom_cl{$id}){
									$all_dom_cl{$id}->{"Count"}++;
									$all_dom_cl{$id}->{"ID"}.=" $ID";
								}else{
									$all_dom_cl{$id}->{"Desc"}=$Data_domain{$ID}->{$id}->{$cpt}->{"Annot"};
									$all_dom_cl{$id}->{"Count"}=1;
									$all_dom_cl{$id}->{"ID"}=" $ID";
								}
								if($start >= $dom_start and $end <= $dom_end){
									if(exists $dom_cl{$id}->{$start}){
										$dom_cl{$id}->{"Count"}++;
									}else{
										$dom_cl{$id}->{"Type"}=$Data_domain{$ID}->{$id}->{$cpt}->{"Type"};
										$dom_cl{$id}->{"Start"}=$Data_domain{$ID}->{$id}->{$cpt}->{"ali_start"};
										$dom_cl{$id}->{"End"}=$Data_domain{$ID}->{$id}->{$cpt}->{"ali_end"};
										$dom_cl{$id}->{"E-value"}=$Data_domain{$ID}->{$id}->{$cpt}->{"E-value"};
										$dom_cl{$id}->{"Count"}=1;
										$dom_cl{$id}->{"Desc"}=$Data_domain{$ID}->{$id}->{$cpt}->{"Annot"};
									}
								}
							}
						}
					}
					$ID_set_cl{$ID}->{"Domain"}= {%dom_cl};
				}
			}
		}
		$All_MOTIFS{$motif}->{"ID_Set"} = {%ID_set_cl};
		$All_MOTIFS{$motif}->{"Domains_list"}={%all_dom_cl};
		$i_motif++;
	}
	$perc=int($i_motif/$#ALL_MOTIFS+1)*100;
	print STDERR "Producing HTML Output...($dir"."$fout.html) [Loading Data in a table - ".sprintf("%2d",$perc)."%]\r"; 

#	print "\n\n";
#	print "    o Filters Homologous Domains...\n";
#	$i_motif=1; my $left_motif2=0; my $i_dom=1;
#	foreach my $k (sort {$All_MOTIFS{$b}->{"Pval"} <=> $All_MOTIFS{$a}->{"Pval"}} (keys (%All_MOTIFS))){
#		my @DOMS = keys (%{$All_MOTIFS{$k}->{"Domains_list"}});	$i_dom=1;
#		my %ID_ban;
#		$All_MOTIFS{$k}->{"S_set_Dom"} = $All_MOTIFS{$k}->{"S_set"};
#		$All_MOTIFS{$k}->{"S_pop_Dom"} = $All_MOTIFS{$k}->{"S_pop"};
#		$All_MOTIFS{$k}->{"N_set_Dom"} = $All_MOTIFS{$k}->{"N_set"};
#		$All_MOTIFS{$k}->{"N_pop_Dom"} = $All_MOTIFS{$k}->{"N_pop"};
#		foreach my $dom (@DOMS){
#			if( $i_motif % 50 == 0 ){ print "\r       Motif ($i_motif/$left_motif)"; }
#			my @IDS = sort(split(" ",$All_MOTIFS{$k}->{"Domains_list"}->{$dom}->{"ID"}));
##				print "\n$k $dom ".$All_MOTIFS{$k}->{"Domains_list"}->{$dom}->{"Desc"}." ".$All_MOTIFS{$k}->{"Domains_list"}->{$dom}->{"Count"}." ".($#IDS+1)."\n";
#			if( ($#IDS+1) > 1){
#				my @count_IDS_sharing_domains=();
#				foreach my $ID (@IDS){
#					if(exists $Data_sharing_domain{$ID}){
#						push(@count_IDS_sharing_domains,$Data_sharing_domain{$ID}->{"#ID"});
##							print "      $ID ".$Data_sharing_domain{$ID}->{"#ID"}." ".$Data_sharing_domain{$ID}->{"ID_sharing_same_dom"};
#					}else{
#						push(@count_IDS_sharing_domains,0);
##							print "      $ID 0\n";
#					}
#				}
#				my $min_sharing_IDS = min(@count_IDS_sharing_domains);
##					print "MINIMUM ".$IDS[which_min(\@count_IDS_sharing_domains)]." ===>  $min_sharing_IDS \n";
#				my @list=which(\@count_IDS_sharing_domains,$min_sharing_IDS);
##					print "LIST IDS MIN "; print @list; print "\n";
#				my $rand=rand($#list+1); my $sel=$list[$rand];
##					print "SELECTION ".$IDS[$sel]."\n";
#				foreach (my $j=0;$j<($#IDS+1);$j++){
#					if($j != $sel){
#						$ID_ban{$IDS[$j]}=1;
#						$All_MOTIFS{$k}->{"S_set_Dom"}--;	$All_MOTIFS{$k}->{"S_pop_Dom"}--;
#						$All_MOTIFS{$k}->{"N_set_Dom"}--;	$All_MOTIFS{$k}->{"N_pop_Dom"}--;
#					}
#				}
#			}
#			$i_dom++;
#		}
#		$All_MOTIFS{$k}->{"ID_Banned"}={%ID_ban};
#		$i_motif++;
#	}
	print "\n\n";

#	open(SELECTED,">Best_selected_motifs.txt");
#	open(HTML_all,">".$dir.$fout.".html");
	open(TOP,">$dir"."Top_motifs.txt");
#	open(FINAL,">$dir"."Final_motifs.txt");
	$i=1;@IDS = keys %{$SEQUENCES}; $perc=0; @ALL_MOTIFS=keys (%All_MOTIFS);
	foreach my $k (sort { $All_MOTIFS{$a}->{"Pval"} <=> $All_MOTIFS{$b}->{"Pval"} || $All_MOTIFS{$b}->{"Cdef"} <=> $All_MOTIFS{$a}->{"Cdef"} || $All_MOTIFS{$b}->{"Len"} <=> $All_MOTIFS{$a}->{"Len"} } (@ALL_MOTIFS)) {
		$perc=int(($i/($#ALL_MOTIFS+1))*100);
		if($perc % 10 ==0){ print STDERR "Producing HTML Output...($dir"."$fout.html) [Writing Results HTML table - ".sprintf("%2d",$perc)."%]\r"; }
#		print STDERR "$k ".$All_MOTIFS{$k}->{"Pval"}." ".$All_MOTIFS{$k}->{"Cdef"}." ".$All_MOTIFS{$k}->{"Len"}."\n";
#		if( $All_MOTIFS{$k}->{"Domain_Homolog"} == 0 ){
#		if( $i % 50 == 0 ){ print "    o Writing the results...Motifs($i/$left_motif)\r"; }
		print HTML_all "<TR valign=middle>\n";
		print HTML_all "<TD><div class=\"display-cell1\"> $i </div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell2\"> $k </div></TD>\n";
#		print HTML_all "<TD><div class=\"display-cell3\"> <a \" title=\" Occurences in non-homologuous domain : ".$All_MOTIFS{$k}->{"S_set_Dom"}."\">".$All_MOTIFS{$k}->{"S_set"}."</a></div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell3\"> ".$All_MOTIFS{$k}->{"S_set"}."</div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell4\"> ";
		my @ID_set = sort (keys (%{$All_MOTIFS{$k}->{"ID_Set"}}));
#				print "$k long: ".($#ID_set+1)."\n";
		foreach my $ID (@ID_set){
#			if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
				if(exists $Data_ID{$ID}){
					print HTML_all "<a ".$Data_ID{$ID}->{"Gene"}."\" title=\"".$Data_ID{$ID}->{"Desc"}."\">".$Data_ID{$ID}->{"Gene"}."</a><BR> ";
				}else{
					print HTML_all "<a href=".$Uniprot_link.$ID." \"#\">".$ID."</a><BR> ";
				}
#			}
		}
		print HTML_all " </div></TD>\n";
#		print HTML_all "<TD><div class=\"display-cell3\"> <a \" title=\" Occurences in non-homologuous domain : ".$All_MOTIFS{$k}->{"S_pop_Dom"}."\">".$All_MOTIFS{$k}->{"S_pop"}."</a></div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell3\"> ".$All_MOTIFS{$k}->{"S_pop"}."</div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell4\">";
		foreach my $ID (sort (keys (%{$All_MOTIFS{$k}->{"ID_Prot"}}))){
#			if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
				if(exists $Data_ID{$ID}){
					print HTML_all "<a ".$Data_ID{$ID}->{"Gene"}."\" title=\"".$Data_ID{$ID}->{"Desc"}."\">".$Data_ID{$ID}->{"Gene"}."</a><BR> ";
				}else{
					print HTML_all "<a href=".$Uniprot_link.$ID." \"#\">".$ID."</a><BR> ";
				}
#			}
		}
		print HTML_all " </div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell5\"> ".sprintf("%.3e",$All_MOTIFS{$k}->{"Pval"})." </div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell6\">\n";
		foreach my $ID (sort (keys (%{$All_MOTIFS{$k}->{"ID_Set"}}))){
			foreach my $id (sort (keys (%{$All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}}))){
#				if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
					my $Type = $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}->{$id}->{"Type"}; my $Annot = $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}->{$id}->{"Desc"}; 
					my $Start_dom = $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}->{$id}->{"Start"}; my $End_dom = $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}->{$id}->{"End"};
					my $Eval=$All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Domain"}->{$id}->{"E-value"};
					print HTML_all $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Gene_name"}."&nbsp;<a href=\"$PFAM_link".$id."&redirect=1\" title=\"$id $Type\"> $Annot </a> $Start_dom-$End_dom <font size=\"2\" color=red> $Eval</font><BR>\n";
#				}
			}
		}

		print HTML_all "\n</div></TD>\n";
		print HTML_all "<TD><div class=\"display-cell6\">";
		foreach my $ID (sort (keys (%{$All_MOTIFS{$k}->{"ID_Set"}}))){
			foreach my $st (sort (keys (%{$All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Seq"}}))){
#				if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
#					print STDERR " $st \n"; <STDIN>;
					if(exists $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Seq"}->{$st}){
						print HTML_all $All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"Seq"}->{$st};
					}
#				}
			}
		}
		
		

		print HTML_all "<TD><div class=\"display-cell6\">";
		foreach my $ID (sort (keys (%{$All_MOTIFS{$k}->{"Evo_ID"}}))){
#			if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
				foreach my $sta (sort (keys (%{$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}}))){
					if(exists $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"Evo"}){
						print HTML_all $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"Evo"};
					}
					if(exists $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"motif"}){
						if(exists $All_MOTIFS{$k}->{"ID_Set"}->{$ID}){
							my @M = split('',$k); my $color="";
	#						my $s=$All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"start"}+1; my $e=$All_MOTIFS{$k}->{"ID_Set"}->{$ID}->{"end"}+1; 
							my $ind_m=0;
							print HTML_all "<table class=\"motiftable\">\n";
							foreach (@M){
								if($_ eq '.'){ $color="green"; }else{ $color="red"; }
								print HTML_all "<TH><font size=2 color=$color> <b><i>$_</i></b></font></TH>\n";
							}
							print HTML_all "<TR>";
							foreach (@M){
								print HTML_all "<TD><font size=2 color=purple> <b><i>".sprintf("%.2f",$Data_Evo{$ID}->{"pos"}->{$ind_m+$sta}->{"Score"})."</i></b></font></TD>\n";
								$ind_m++;
							}
							print HTML_all "</table>\n";
						}
						my $mean_def = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"motif"}->{"def"} / $All_MOTIFS{$k}->{"Cdef"};
						my $mean_wild = 0;
						if( $All_MOTIFS{$k}->{"Len"} != $All_MOTIFS{$k}->{"Cdef"} ){ $mean_wild = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"motif"}->{"wild"} / ($All_MOTIFS{$k}->{"Len"}-$All_MOTIFS{$k}->{"Cdef"}) ; }
						print HTML_all "<font size=2 color=purple> M(Pos. Def.) = ".sprintf("%.2f",$mean_def)." M(Wildcards) = ".sprintf("%.2f",$mean_wild)."<BR></font>\n";
					}

	#				if(exists $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"context"}){
	#					print HTML_all "<font size=2 color=green> S(Context Before)= ".$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"context"}->{"before"}."</font>\n";
	#					print HTML_all "<font size=2 color=red> S(Context After)= ".$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"context"}->{"after"}."<BR></font>\n";
	#				}
				
						my $mean_contextB =0;
						if( $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"before"}->{"nAA"} != 0 ){
							$mean_contextB = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"before"}->{"Score"} / $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"before"}->{"nAA"};
						}
						my $mean_contextA =0;
						if( $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"after"}->{"nAA"} != 0 ){
							$mean_contextA = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"after"}->{"Score"} / $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"after"}->{"nAA"};
						}
						print HTML_all "<font size=2 color=green> M(Context Before)= ".sprintf("%.2f",$mean_contextB)."</font>\n<font size=2 color=red> M(Context After)= ".sprintf("%.2f",$mean_contextA)."<BR></font>\n";
						print HTML_all "<font size=2 color=black> M(Context After+Before)= ".sprintf("%.2f",(($mean_contextB+$mean_contextA)/2))."<BR></font>\n";
				}
#			}
		}
		print HTML_all " </div></TD>\n";

		my $occ=0;
		my $sum_score_motif_stretch_def=0; my $sum_score_motif_stretch_wild=0; my $sum_score_motif_stretch=0;
		my $sum_aa_motif_stretch_def=0; my $sum_aa_motif_stretch_wild=0; my $sum_aa_motif_stretch=0;
		my $mean_motif_stretch_def=0; my $mean_motif_stretch_wild=0; my $aa_motif_stretch=0;


		my $sum_score_motif_protein_def=0; my $sum_score_motif_protein_wild=0; my $sum_score_motif_protein=0;
		my $sum_aa_motif_protein_def=0; my $sum_aa_motif_protein_wild=0; my $sum_aa_motif_protein=0;
		my $mean_motif_protein_def=0; my $mean_motif_protein_wild=0; my $aa_motif_protein=0;


		my $mean_motif_def=0; my $mean_motif_wild=0; my $mean_motif=0;
		my $mean_context=0; my $mean_stretch=0;
		my $mean_prot=0;

		my $ratio_motif_stretch_def=0;
		my $ratio_motif_stretch_wild=0;
		my $ratio_motif_stretch=0;

		my $ratio_motif_protein_def=0;
		my $ratio_motif_protein_wild=0;
		my $ratio_motif_protein=0;


		my $len=$All_MOTIFS{$k}->{"Len"}; my $cdef=$All_MOTIFS{$k}->{"Cdef"}; my $wild = $len - $cdef; if($wild == 0){ $wild++; }
		print HTML_all "<TD><div class=\"display-cell6\">";
		foreach my $ID (sort (keys (%{$All_MOTIFS{$k}->{"Evo_ID"}}))){
#			if(!(exists $All_MOTIFS{$k}->{"ID_Banned"}->{$ID})){
			my $occ_ID=0; 
			my $sum_motif_def_ID=0; my $sum_motif_wild_ID=0; my $sum_motif_ID=0;
			my $sum_aa_context_ID=0; my $sum_score_context_ID=0; my $sum_score_stretch_ID=0; my $sum_aa_stretch_ID=0;
			$mean_prot+=$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"disorder"}->{"mean"};

			foreach my $sta (sort (keys (%{$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}}))){
#				print STDERR " ---|| ($occ) $ID $k $sta [$occ_ID] ||--- \n";
				$sum_motif_def_ID += $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"motif"}->{"def"};
				$sum_motif_wild_ID += $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"motif"}->{"wild"};
				$sum_motif_ID+= ( $sum_motif_def_ID + $sum_motif_wild_ID );

	
				my $Score_context_before = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"before"}->{"Score"};
				my $Score_context_after = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"after"}->{"Score"};
				
				my $aa_context_before = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"before"}->{"nAA"};
				my $aa_context_after = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"context"}->{"disorder"}->{"after"}->{"nAA"};

				my $Score_stretch_before = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"disorder"}->{"before"}->{"Score"};
				my $Score_stretch_after = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"disorder"}->{"after"}->{"Score"};
				
				my $aa_stretch_before = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"disorder"}->{"before"}->{"nAA"};
				my $aa_stretch_after = $All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"disorder"}->{"after"}->{"nAA"};

				$sum_score_context_ID += ($Score_context_before + $Score_context_after);
				$sum_aa_context_ID += ($aa_context_before + $aa_context_after);
				if( $sum_aa_context_ID == 0 ){$sum_aa_context_ID=1;}

				$sum_score_stretch_ID += ($Score_stretch_before + $Score_stretch_after);
				$sum_aa_stretch_ID += ($aa_stretch_before + $aa_stretch_after);
				if( $sum_aa_stretch_ID == 0 ){$sum_aa_stretch_ID=1;}
#				print STDERR "$k $sta $ID occID $occ_ID score bef $Score_stretch_before score after $Score_stretch_after aa bef $aa_stretch_before aa after $aa_stretch_after\n"."stretch before ".$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"before"}->{"seq"}."\n"."stretch after ".$All_MOTIFS{$k}->{"Evo_ID"}->{$ID}->{"occ"}->{$sta}->{"stretch"}->{"after"}->{"seq"}."\n"."sequence ".$SEQUENCES->{$ID}->{'Seq'}."\n";
#				<STDIN>;
				$occ_ID++;
			}
#			print STDERR " $ID SUM aa $sum_aa_context_ID SUM score $sum_score_context_ID\n";
#			<STDIN>;


			$mean_motif_def+=(($sum_motif_def_ID/$cdef)/$occ_ID); $mean_motif_wild+=(($sum_motif_wild_ID/$wild)/$occ_ID); $mean_motif+=(($sum_motif_ID/$len)/$occ_ID);
			$mean_context+=(($sum_score_context_ID/$sum_aa_context_ID)/$occ_ID);	$mean_stretch+=(($sum_score_stretch_ID/$sum_aa_stretch_ID)/$occ_ID);

#			$sum_score_motif_stretch_def+=$mean_motif_def+$mean_stretch; 
#			$sum_score_motif_stretch_wild+=$mean_motif_wild+$mean_stretch;
#			$sum_score_motif_stretch+=($sum_motif_def_ID/$occ_ID)+($sum_motif_wild_ID/$occ_ID)+($sum_score_stretch_ID/$occ_ID); 
#			$sum_aa_motif_stretch_def+=($sum_motif_def_ID/$occ_ID)+($sum_score_stretch_ID/$occ_ID); 
#			$sum_aa_motif_stretch_wild+=($sum_motif_wild_ID/$occ_ID)+($sum_score_stretch_ID/$occ_ID); 
#			$sum_aa_motif_stretch+=($sum_motif_def_ID/$occ_ID)+($sum_motif_wild_ID/$occ_ID)+($sum_score_stretch_ID/$occ_ID); 
#			$mean_motif_stretch_def += ; $mean_motif_stretch_wild += ; $mean_motif_stretch +=  ;
#			$mean_motif_protein_def += ; $mean_motif_protein_wild += ; $mean_motif_protein +=  ;
#				print STDERR "$ID $k $mean_context\n"; <STDIN>;
#			}
			$occ++;
		}

		if($occ !=0){
			$mean_motif_def=($mean_motif_def/$occ); $mean_motif_wild=($mean_motif_wild/$occ); $mean_motif=($mean_motif/$occ);
			$mean_stretch=($mean_stretch/$occ); 
			$mean_context=($mean_context/$occ); 
			$mean_prot=($mean_prot/$occ);
		
			print HTML_all "<b>Average Evolution for ".Font_Color_motifs($k)." :<b> (<font size=2 color=purple> All Pos = ".sprintf("%.2f",$mean_motif)."</font>)<BR>\n";
			print HTML_all "<font size=2 color=red   > Defined Pos.   = ".sprintf("%.2f",$mean_motif_def)."</font>";
			print HTML_all "<font size=2 color=green > Wildcards Pos. = ".sprintf("%.2f",$mean_motif_wild)." </font><BR>\n";
			print HTML_all "<font size=2 color=black > Context (AAs [-10,+10])    = ".sprintf("%.2f",$mean_context)."</font>";
			print HTML_all "<font size=2 color=blue  > Proteins (k=$occ)          = ".sprintf("%.2f",$mean_prot)."</font><BR>\n";
	
			$All_MOTIFS{$k}->{"Evolution"}->{"mean"}->{"motif"}->{"all"}=$mean_motif;
			$All_MOTIFS{$k}->{"Evolution"}->{"mean"}->{"motif"}->{"def"}=$mean_motif_def;
			$All_MOTIFS{$k}->{"Evolution"}->{"mean"}->{"motif"}->{"wild"}=$mean_motif_wild;
			$All_MOTIFS{$k}->{"Evolution"}->{"mean"}->{"context"}->{"disorder"}=$mean_context;
			$All_MOTIFS{$k}->{"Evolution"}->{"mean"}->{"protein"}->{"disorder"}=$mean_prot;
		

			$ratio_motif_stretch_def=($mean_motif_def/$mean_stretch);
			$ratio_motif_stretch_wild=($mean_motif_wild/$mean_stretch);
			$ratio_motif_stretch=($mean_motif/$mean_stretch);

			$ratio_motif_protein_def=($mean_motif_def/$mean_prot);
			$ratio_motif_protein_wild=($mean_motif_wild/$mean_prot);
			$ratio_motif_protein=($mean_motif/$mean_prot);
		}
		print HTML_all " </div></TD>\n";

		print HTML_all "<TD><div class=\"display-cell6\">";
		if($occ !=0){
			print HTML_all "<div class=\"EvoRate\">\n";
			print HTML_all "<font size=2 color=purple>S(motif/stretch)=".sprintf("%.2f",$ratio_motif_stretch)."</font><BR>\n";
			print HTML_all "<font size=2 color=red>S(def/stretch)=".sprintf("%.2f",$ratio_motif_stretch_def)."</font><BR>\n";
			print HTML_all "<font size=2 color=green>S(wild/stretch)=".sprintf("%.2f",$ratio_motif_stretch_wild)."</font><BR>\n";
			print HTML_all " </div>\n";
			print HTML_all "<div class=\"EvoRate\">\n";
			print HTML_all"<font size=2 color=purple>S(motif/prot)=".sprintf("%.2f",$ratio_motif_protein)."</font><BR>\n";
			print HTML_all"<font size=2 color=red>S(def/prot)=".sprintf("%.2f",$ratio_motif_protein_def)."</font><BR>\n";
			print HTML_all"<font size=2 color=green>S(wild/prot)=".sprintf("%.2f",$ratio_motif_protein_wild)."</font><BR>\n";
			print HTML_all " </div>\n";
		}
#		print HTML_all "<BR>\n";
#		print HTML_all "<font size=2 color=purple>S(motif/prot)=".$ratio_motif_protein."</font><BR>\n";
#		print HTML_all "<font size=2 color=red>S(def/prot)=".$ratio_motif_protein_def."</font><BR>\n";
#		print HTML_all "<font size=2 color=green>S(wild/prot)=".$ratio_motif_protein_wild."</font><BR>\n";
		print HTML_all " </div></TD>\n";

#			print SELECTED "$k ".$All_MOTIFS{$k}->{"S_set"}." ".$All_MOTIFS{$k}->{"S_pop"}." ".$All_MOTIFS{$k}->{"N_set"}." ".$All_MOTIFS{$k}->{"N_pop"}." ".$All_MOTIFS{$k}->{"Pval"}."\n";

#		if($All_MOTIFS{$k}->{"S_set_Dom"}>=3){
#			print FINAL "$k\t".$All_MOTIFS{$k}->{"S_set_Dom"}."\t".$All_MOTIFS{$k}->{"S_pop_Dom"}."\t".$All_MOTIFS{$k}->{"N_set_Dom"}."\t".$All_MOTIFS{$k}->{"N_pop_Dom"}."\n";
#		}

		# Write Header to the final selection of motifs
#		if($i==1){ print TOP "# Mot\tOcc.set\tOcc.back\tnseq.set\tnseq.back\tPval\tComb.AA\tnb.degen.patt\tPval.adj\tnwild\tL.Mot\tDef.pos\tSetname\tratio.evo.Mot(def/stretch)\tratio.evo.Mot(def/prot)\tratio.evo.Mot(wild/stretch)\tratio.evo.Mot(wild/prot)\tratio.evo(Mot/stretch)\tratio.evo(Mot/prot)\tmean.evo(def)\tmean.evo(wild)\tmean.evo(Mot)\tmean.evo(context)\tmean.evo(stretch)\tmean.evo(prot)\n"; }
		if($i==1){ print TOP "# Mot\tOcc.set\tOcc.back\tnseq.set\tnseq.back\tPval\tL.Mot\tDef.pos\tSetname\tratio.evo.Mot(def/stretch)\tratio.evo.Mot(def/prot)\tratio.evo.Mot(wild/stretch\tratio.evo.Mot(wild/prot)\tratio.evo(Mot/stretch)\tratio.evo(Mot/prot)\tmean.evo(def)\tmean.evo(wild)\tmean.evo(Mot)\tmean.evo(context)\tmean.evo(stretch)\tmean.evo(prot)\n"; }
		print TOP "$k\t$All_MOTIFS{$k}->{'S_set'}\t$All_MOTIFS{$k}->{'S_pop'}\t$All_MOTIFS{$k}->{'N_set'}\t$All_MOTIFS{$k}->{'N_pop'}\t$All_MOTIFS{$k}->{'Pval'}\t";
#		print TOP "$k\t$All_MOTIFS{$k}->{'S_set'}\t$All_MOTIFS{$k}->{'S_pop'}\t$All_MOTIFS{$k}->{'N_set'}\t$All_MOTIFS{$k}->{'N_pop'}\t$All_MOTIFS{$k}->{'Pval'}\t$All_MOTIFS{$k}->{'nbdegen'}\t";
		print TOP "$All_MOTIFS{$k}->{'Len'}\t$All_MOTIFS{$k}->{'Cdef'}\t$fseqname\t";
#		print TOP "$All_MOTIFS{$k}->{'Pval_adj'}\t$All_MOTIFS{$k}->{'nwild'}\t$All_MOTIFS{$k}->{'Len'}\t$All_MOTIFS{$k}->{'Cdef'}\t$fseqname\t";
		print TOP "$ratio_motif_stretch_def\t$ratio_motif_protein_def\t$ratio_motif_stretch_wild\t$ratio_motif_protein_wild\t$ratio_motif_stretch\t$ratio_motif_protein\t";
		print TOP "$mean_motif_def\t$mean_motif_wild\t$mean_motif\t$mean_context\t$mean_stretch\t$mean_prot\n";
		$i++;
	}
	$perc=int(($i/($#ALL_MOTIFS+1))*100);
	print STDERR "Producing HTML Output...($dir"."$fout.html) [Writing Results HTML table - ".sprintf("%2d",$perc)."%] => OK \n"; 
	
	print "\n\n";
	print HTML_all "\n</table>\n</div><div class=\"contenu2\"><font size=5 color=\"blue\"><b> ".($#IDS+1)." Sequences</b></font> </div>\n</body>\n   </html>";
	close(HTML_all); 
#	close(FINAL);
	close(TOP);
}else{
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, consult the help to launch the program :"," ","#",">");
	print "\n".format_string(70,"#        $0       --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
