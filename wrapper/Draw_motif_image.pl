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
	use if (exists $ENV{'SCRIPTSDIR'}), lib => $ENV{'SCRIPTSDIR'};
	if (not exists $ENV{'SCRIPTSDIR'}) {	die "No Environment Variable for locating Scripts Directory ( \$SCRIPTSDIR ). Please refer to documentation.\n"; }
	$ENV{'DATADIR'}=abs_path($ENV{'SCRIPTSDIR'}."/../Data");
}
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;
use Tk;
use Tk::TableMatrix;
use File_Utils;
use String_Print_Utils;

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

# Silent options
my $help = ""; my $quiet= "";

my $Seq_file = ""; my $Motifs_file = ""; my $Proteome = "";
my $resfile  = "./Best_Motifs.html";
my $logdir   = "./";
my $Pfam_annot=abs_path($ENV{'DATADIR'}."Domains/Yeast_Pfam_annotations.txt");
my $Superfamily_annot=abs_path($ENV{'DATADIR'}."Domains/Saccharomyces_cerevisiae_SUPERFAMILY_domains.txt");
my $Gene_annot=abs_path($ENV{'DATADIR'}."Genes/YEAST.data");

GetOptions(
	"help!"              => \$help,
	"nov|quiet!"         => \$quiet,
	"Seq=s"              => \$Seq_file,
	"Motifs=s"           => \$Motifs_file,
	"Proteome=s"         => \$Proteome,
	"resfile|out:s"      => \$resfile,
	"logdir:s"           => \$logdir,
	"Gene_annot:s"       => \$Gene_annot,
	"Pfam_annot:s"       => \$Pfam_annot,
	"Superfam_annot:s"   => \$Superfamily_annot
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
    --Pfam_annot     |  <File directory>                   |  Domains Annotation file obtained from PFAM Database                             
    --Superfam_annot |  <File directory>                   |  Domains Annotation file obtained from SUPERFAMILY Database                      
    --Gene_file      |  <File directory>                   |  File containing Names and Descriptions of genes                                 
                     |                                     |                                                                                  
    --out            |  <File or Path directory>           |  Location or Name of the Results file                                            
                     |                                     |                                                                                  
    --logdir         |  <Path directory>                   |  Location of the Log Directory                                                   
                     |                                     |                                                                                  
    --nov            |                                     |  Only print Logo, progression and time usage                                     
                     |                                     |                                                                                  
    --help           |                                     |  Prints the help description and exit                                            
                     |                                     |                                                                                  
----------------------------------------------------------------------------------------------------------------------------------------------
(*) Mandatory options                                                                                                                         \n";
# If the help parameter is active, the program does not run but the text above is printed to the screen.
if ($help) {
	print "	#########################################################################################################################################\n";
	print "	#                               VisualMotifs - Draw a graphical representation for a list of motifs -                                   #\n";
	print "	#---------------------------------------------------------------------------------------------------------------------------------------#\n";
	print "	# This program takes at most 10 parameters :                                                                                            #\n";
	print "	#   1)  --help                                     + Description of Program Usage                                    +  (Optional)      #\n";
	print "	#   2)  --quiet,--nov                              + Non verbous Mode                                                +  (Optional)      #\n";
	print "	#   3)  --Seq                                      + File of Proteins Sequences of Interest (FASTA)                  +  (Mandatory)     #\n";
	print "	#   4)  --Proteome                                 + Fasta Formatted File containing Proteome Sequences              +  (Mandatory)     #\n";
	print "	#   5)  --Motifs                                   + Over Representation Result File                                 +  (Mandatory)     #\n";
	print "	#   6)  --resfile,--out                            + Location of the Result file                                     +  (Optional)      #\n";
	print "	#   7)  --logdir                                   + Location of the Log Directory                                   +  (Optional)      #\n";
	print "	#   8)  --Pfam_annot     (by default for yeast)    + Pfam domains annotations file for S. Cerevisiae                 +  (Optional)      #\n";
	print "	#   9)  --Superfam_annot (by default for yeast)    + Superfamily domains annotations file for S. Cerevisiae          +  (Optional)      #\n";
	print "	#  10)  --Gene_annot     (by default for yeast)    + File containing Names and Descriptions of S. Cerevisiae genes   +  (Optional)      #\n";
	print "	#########################################################################################################################################\n";
	exit(0);
}

if($Seq_file ne "" and $Motifs_file ne "" and $Proteome ne "") {

	# Quiet mode
	if ($quiet){
		# Redirection of printed messages to a log file ('Visualization.log') to the user defined directory
		open (STDOUT,">$logdir/Motif_Visualization.log") or die ("Unable to redirect STDOUT to a file...");
	}

	# File verifications #
	my @path_to_Seq_file=fileparse($Seq_file,qr{\..*});

	# Existence of the Setfile
	print "- Set file                        :  \"$Seq_file\" => "; my ($rc_seq,$TEXT) = File_verification($Seq_file); print $TEXT;
	if($rc_seq != 0){ die "# Something wrong happened when verifying the location of the Setfile #\n"; }
	my $SEQUENCES = Load_FASTA($Seq_file);

	# Existence of the Proteome file
	print "- Proteome Sequence File          :  \"$Proteome\" => "; my ($rc_prot,$TEXT2) = File_verification($Proteome); print $TEXT2;
	if($rc_prot != 0){ die "# Something wrong happened when verifying the location of the Proteome file #\n"; }
	my $PROTEOME = Load_FASTA($Proteome);

	# Existence of the BLAST results file
	print "- Motifs over-representation file  :    \"$Motifs_file\" => "; my ($rc_motifs,$TEXT3) = File_verification($Motifs_file); print $TEXT3;
	if($rc_motifs != 0){ die "# Something wrong happened when verifying the location of the over-represented Motifs file #\n"; }
	my ($Ref_MOTIFS)=load_motifs($Motifs_file);
	my %MOTIFS = %{$Ref_MOTIFS};


	# LOADING THE MOTIFS, THE SEQUENCES FROM THE SET AND FROM THE PROTEOME, SOME INFORMATION ABOUT THE ORFS, THE PFAM ANNOTATIONS FOR DOMAINS, THE 

	my $Ref_Data_orf = load_data_ORF($Gene_annot);
	my %Data_orf = %{$Ref_Data_orf};
	my $Ref_Data_domain = load_PFAM_annot($Pfam_annot);
	my %Data_domain  = %{$Ref_Data_domain};
	my $Ref_Sharing_domains = load_Sharing_PFAM_Domains($ENV{'DATADIR'}."Domains/Yeast_ORFs_sharing_same_Pfam_domains.txt");
	my %Data_sharing_domain  = %{$Ref_Sharing_domains};

	my $nb_motifs=0;
	my $i=1; my $cpt=1; my $NSeq=0;
	my $i_motif=1; my $left_motif=0; my $i_orf=1;
	my ($seq, $ORF, $motif)=("","","");

	my $Gene_name="";
	$i_motif=1;
	# FIND THE MOTIF IN EVERY SEQUENCES IN THE PROTEOME AND SEARCH FOR THE DOMAINS INSIDE EACH ORF AND FOR THOSE WHICH CONTAIN THE MOTIF
	my @ORFS = keys %{$PROTEOME};
	foreach $motif (sort {$MOTIFS{$b}->{"Pval"} <=> $MOTIFS{$a}->{"Pval"}} (keys (%MOTIFS))){
		my $i_motif_orf=1; my %orf_prot_cl=(); $i_orf=1; my %all_dom_cl=(); 
		foreach $ORF (@ORFS){
			$i_orf++; my %dom_cl=();
			if(exists $PROTEOME->{$ORF}){
				$seq = $PROTEOME->{$ORF}; my %orf_dom_cl=(); 
				if( $seq =~ /$motif/){
					my ($start,$end)=match_positions($motif,$seq); $orf_prot_cl{$ORF}->{"start"}=$start; $orf_prot_cl{$ORF}->{"end"}=$end;
					$Gene_name=$ORF; if(exists $Data_orf{$ORF}){ $Gene_name=$Data_orf{$ORF}->{"Gene"};	}
					$orf_prot_cl{$ORF}->{"Gene_name"} = $Gene_name;   $orf_prot_cl{$ORF}->{"Length"} = length($seq);
					my @max_end;
					if(exists $Data_domain{$ORF}){
						foreach my $id ( keys (%{$Data_domain{$ORF}})){
							foreach my $cpt ( keys (%{$Data_domain{$ORF}->{$id}})){
								my $dom_start = $Data_domain{$ORF}->{$id}->{$cpt}->{"ali_start"}; my $dom_end = $Data_domain{$ORF}->{$id}->{$cpt}->{"ali_end"};
								if(exists $orf_dom_cl{$id}){
									$orf_dom_cl{$id}->{"Count"}++;
									$orf_dom_cl{$id}->{"ORF"}.=" $ORF";
								}else{
									$orf_dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
									$orf_dom_cl{$id}->{"Count"}=1;
									$orf_dom_cl{$id}->{"ORF"}=" $ORF";
									$orf_dom_cl{$id}->{"Type"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Type"};
									$orf_dom_cl{$id}->{"Start"}=$dom_start;
									$orf_dom_cl{$id}->{"End"}=$dom_end;
									$orf_dom_cl{$id}->{"E-value"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"E-value"};
								}
								if(exists $all_dom_cl{$id}){
									$all_dom_cl{$id}->{"Count"}++;
									$all_dom_cl{$id}->{"ORF"}.=" $ORF";
								}else{
									$all_dom_cl{$id}->{"Count"}=1;
									$all_dom_cl{$id}->{"ORF"}=" $ORF";
									$all_dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
								}
								if($start >= $dom_start and $end <= $dom_end){
									if(exists $dom_cl{$id}){
										$dom_cl{$id}->{"Count"}++;
									}else{
										$dom_cl{$id}->{"Type"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Type"};
										$dom_cl{$id}->{"Start"}=$dom_start;
										$dom_cl{$id}->{"End"}=$dom_end;
										$dom_cl{$id}->{"E-value"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"E-value"};
										$dom_cl{$id}->{"Count"}=1;
										$dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
										$dom_cl{$id}->{"Length"}=$dom_end-$dom_start;
										push (@max_end,$dom_end);
									}
								}
							}
						}
					}
					$orf_prot_cl{$ORF}->{"Domain_Prot"}= {%dom_cl};
					$orf_prot_cl{$ORF}->{"Max_end_Domain_Prot"}= max(@max_end);
					$orf_prot_cl{$ORF}->{"Domains_list_Orf_prot"}= {%orf_dom_cl};
				}
			}
		}
		$MOTIFS{$motif}->{"ORF_Prot"} = {%orf_prot_cl};
		$MOTIFS{$motif}->{"ALL_Domains_Prot"} = {%all_dom_cl};
		$MOTIFS{$motif}->{"Max_Length"} = 0;
		$i_motif++;
	}

	$i_motif=1; 
	# FIND THE MOTIF IN EVERY SEQUENCES IN THE SET AND SEARCH FOR THE DOMAINS INSIDE EACH ORF AND FOR THOSE WHICH CONTAIN THE MOTIF
	@ORFS = keys %{$SEQUENCES};
	foreach my $motif (sort {$MOTIFS{$b}->{"Pval"} <=> $MOTIFS{$a}->{"Pval"}} (keys (%MOTIFS))){
		$i_orf=0; my $i_motif_orf=1; my %orf_set_cl=(); my %all_dom_cl=(); 
		foreach my $ORF (@ORFS){
			$i_orf++; my %dom_cl=(); my %orf_dom_cl=();
			my $seq = $SEQUENCES->{$ORF};
			if( $seq =~ /$motif/){
				my ($start,$end)=match_positions($motif,$seq); $orf_set_cl{$ORF}->{"start"}=$start; $orf_set_cl{$ORF}->{"end"}=$end;
				$Gene_name=$ORF; if(exists $Data_orf{$ORF}){ $Gene_name=$Data_orf{$ORF}->{"Gene"};	}
				$orf_set_cl{$ORF}->{"Gene_name"} = $Gene_name;
				$orf_set_cl{$ORF}->{"Length"} = length($seq);
				my @max_end;
				if(exists $Data_domain{$ORF}){
					foreach my $id ( keys (%{$Data_domain{$ORF}})){
						foreach my $cpt ( keys (%{$Data_domain{$ORF}->{$id}})){
							my $dom_start = $Data_domain{$ORF}->{$id}->{$cpt}->{"ali_start"}; my $dom_end = $Data_domain{$ORF}->{$id}->{$cpt}->{"ali_end"};
							if(exists $orf_dom_cl{$id}){
								$orf_dom_cl{$id}->{"Count"}++;
								$orf_dom_cl{$id}->{"ORF"}.=" $ORF";
								$all_dom_cl{$id}->{"Count"}++;
								$all_dom_cl{$id}->{"ORF"}.=" $ORF";
							}else{
								$orf_dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
								$orf_dom_cl{$id}->{"Count"}=1;
								$orf_dom_cl{$id}->{"ORF"}=" $ORF";
								$orf_dom_cl{$id}->{"Type"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Type"};
								$orf_dom_cl{$id}->{"Start"}=$dom_start;
								$orf_dom_cl{$id}->{"End"}=$dom_end;
								$orf_dom_cl{$id}->{"E-value"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"E-value"};
								$all_dom_cl{$id}->{"Count"}=1;
								$all_dom_cl{$id}->{"ORF"}=" $ORF";
								$all_dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
							}
							if($start >= $dom_start and $end <= $dom_end){
								if(exists $dom_cl{$id}){
									$dom_cl{$id}->{"Count"}++;
								}else{
									$dom_cl{$id}->{"Type"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Type"};
									$dom_cl{$id}->{"Start"}=$dom_start;
									$dom_cl{$id}->{"End"}=$dom_end;
									$dom_cl{$id}->{"E-value"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"E-value"};
									$dom_cl{$id}->{"Count"}=1;
									$dom_cl{$id}->{"Desc"}=$Data_domain{$ORF}->{$id}->{$cpt}->{"Annot"};
									$dom_cl{$id}->{"Length"}=$dom_end-$dom_start;
									push (@max_end,$dom_end);
								}
							}
						}
					}
				}
				$orf_set_cl{$ORF}->{"Domain_Set"}= {%dom_cl};
				$orf_set_cl{$ORF}->{"Max_end_Domain_Set"}= max(@max_end);
				$orf_set_cl{$ORF}->{"Domains_list_Orf_set"}= {%orf_dom_cl};
			}
		}
		$MOTIFS{$motif}->{"ORF_Set"} = {%orf_set_cl};
		$MOTIFS{$motif}->{"ALL_Domains_Set"} = {%all_dom_cl};
		$i_motif++;
	}


	# FILTER FOR SHARED DOMAINS
	$i_motif=1; my $left_motif2=0; my $i_dom=1;
	foreach my $mot (sort {$MOTIFS{$b}->{"Pval"} <=> $MOTIFS{$a}->{"Pval"}} (keys (%MOTIFS))){
		my @DOMS = keys (%{$MOTIFS{$mot}->{"ALL_Domains_Prot"}});	$i_dom=1;
		my %orf_ban;
		foreach my $dom (@DOMS){
			my @orfs = sort(split(" ",$MOTIFS{$mot}->{"ALL_Domains_Prot"}->{$dom}->{"ORF"}));
				my @count_orfs_sharing_domains=();
				foreach my $orf (@orfs){
					if(exists $Data_sharing_domain{$orf}){
						push(@count_orfs_sharing_domains,$Data_sharing_domain{$orf}->{"#ORF"});
					}else{
						push(@count_orfs_sharing_domains,0);
					}
				}
				my $min_sharing_orfs = min(@count_orfs_sharing_domains);
				my @list=which(\@count_orfs_sharing_domains,$min_sharing_orfs);
				my $rand=rand($#list+1); my $sel=$list[$rand];
				foreach (my $j=0;$j<($#orfs+1);$j++){
					if($j != $sel){
						$orf_ban{$orfs[$j]}=1;
					}
				}
			$i_dom++;
		}
		$MOTIFS{$mot}->{"ORF_Banned"}={%orf_ban};
	}

	my $main = MainWindow->new;
	$main->title("List of Motifs");
	my $xtvar; my $xtvar2; 
	$xtvar->{"0,0"} = "No";$xtvar->{"0,1"} = "MOTIF"; $xtvar->{"0,2"} = "Pval"; $xtvar->{"0,3"} = "#Set"; $xtvar->{"0,4"} = "#Prot";
	my $row = 1;
	foreach my $mot (sort { $MOTIFS{$b}->{"Pval"} <=> $MOTIFS{$a}->{"Pval"} } (keys (%MOTIFS))){
		foreach my $col (0 .. 5){
			if ($col == 0){
				$xtvar->{"$row,$col"} = "$row";
				$xtvar2->{"$row,$col"}->{"Mot"} = $mot;
			}elsif($col == 1){
				$xtvar->{"$row,$col"} = "$mot";
			}elsif($col == 2){
				$xtvar->{"$row,$col"} = $MOTIFS{$mot}->{"Pval"};
			}elsif($col == 3){
				$xtvar->{"$row,$col"} = $MOTIFS{$mot}->{"S_set"};
			}elsif($col == 4){
				$xtvar->{"$row,$col"} = $MOTIFS{$mot}->{"S_pop"};
			}
		}
		$row++;
	}

	my @motifs = keys ( %MOTIFS);
	my $xtable = $main->Scrolled(
		'TableMatrix',
		-rows          => ($#motifs+2),
		-cols          => 6,
		-ipadx         => 3,
		-variable      => $xtvar,
		-selectmode    => 'single',
		-selecttype    => 'row',
		-resizeborders => 'none',
		-state         => 'disabled',
		-cursor        => 'top_left_arrow',
		-bg            => 'white',
		-scrollbars    => 'ose',
	)->pack;

	# Clean up if mouse leaves the widget
	$xtable->bind(
		'<FocusOut>',
		sub {
			my $w = shift;
			$w->selectionClear('all');
		}
	);

	# Highlight the cell under the mouse
	$xtable->bind(
		'<Motion>',
		sub {
			my $w  = shift;
			my $Ev = $w->XEvent;
			if ( $w->selectionIncludes( '@' . $Ev->x . "," . $Ev->y ) ) {
				Tk->break;
			}
			$w->selectionClear('all');
			$w->selectionSet( '@' . $Ev->x . "," . $Ev->y );
			Tk->break;
		}
	);

	# MouseButton 1 event
	$xtable->bind(
		'<1>',
		sub {
			my $w = shift;
			$w->focus;
			my ($rc) = @{ $w->curselection };
			if(exists $xtvar2->{$rc}->{"Mot"}){
				my $window = MainWindow->new;
				my $mot = $xtvar2->{$rc}->{"Mot"};
				my @orfs = sort keys (%{$MOTIFS{$mot}->{"ORF_Set"}});
				my @orfs_prot = sort keys (%{$MOTIFS{$mot}->{"ORF_Prot"}});
				my @orfs_banned = keys (%{$MOTIFS{$mot}->{"ORF_Banned"}});
				my $adj = ($#orfs_prot+1) - $MOTIFS{$mot}->{"S_pop"};
				my $dom_space=0;
#				print join(" ",@orfs_banned)."\n";
#				print join(" ",@orfs_prot)."\n";
#				print join(" ",@orfs)."\n";
				$window->title("Motif $mot");
#				$window->Label(-text => "$id $function\n PP-Value: ".$MOTIFS{$mot}->{"Pval"}."\t| P-value: ".$MOTIFS{$mot}->{"Pval"}."\n # Set:".$MOTIFS{$mot}->{"S_set"}."\t# Proteome:".$MOTIFS{$mot}->{"S_pop"}."\n# Banned ORF (Shared Domains):".($#orfs_banned+1-$adj),-font=>"fixed 10 bold")->pack;
				my $drawing = $window->Scrolled('Canvas',-width => 1400, -height => 800,-background=>"white", -scrollregion=>[0,0,5000,5000]) -> pack;
				$i_orf=1; my $xside = 50; my $yside = 160; my $space = 15; my $header_size = 150;


				# INFORMATIONS ABOUT THE MOTIFS : ID | Feature targeted | # ORFS | Pval | Pvalue | Banned ORFS | Occurences in set | Occurences in Proteome
#				$drawing->createText(700,20,-text => "$id $function (".$MOTIFS{$mot}->{"N_set"}." ORFs)",-font=>['fixed','10','bold'],-anchor=>"center");
				$drawing->createText(700,35,-text => " P-value: ".$MOTIFS{$mot}->{"Pval"}."\t",-font=>['fixed','10','bold'],-anchor=>"center");
				$drawing->createText(700,50,-text => " # Banned ORF (Shared Domains):".($#orfs_banned+1-$adj),-font=>['fixed','10','bold'],-anchor=>"center");
				$drawing->createText(700,65,-text => " # Set:".$MOTIFS{$mot}->{"S_set"}."\t# Proteome:".$MOTIFS{$mot}->{"S_pop"},-font=>['fixed','10','bold'],-anchor=>"center");

				# LEGEND FOR THE GRAPHIC INTERFACE
				$drawing->createText(30,110,-text=>"Legend :",-anchor=>"center");
				$drawing->createText(240,130,-text=>"# of the current ORF / Total # of ORFs",-anchor=>"e",-font=> ['fax','10','italic']);
				$drawing->createRectangle(250,125,260,135,-fill=>"black",-outline=>"grey");
#				$drawing->createText(275,130,-text=>"Sequences corresponding to the ID ($id)",-anchor=>"w",-font=> ['fax','10','italic']);
				$drawing->createRectangle(250,150,260,160,-fill=>"white",-outline=>"grey");
#				$drawing->createText(275,155,-text=>"Sequences from the proteome ( without the sequences corresponding to $id)",-anchor=>"w",-font=> ['fax','10','italic']);

				# DRAW MOTIFS IN ORFS WITHIN THE SET
				foreach my $orf (@orfs){
					if(!(exists $MOTIFS{$mot}->{"ORF_Banned"}->{$orf})){
						my $Length = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Length"};
						if ($MOTIFS{$mot}->{"Max_Length"}<$Length){ $MOTIFS{$mot}->{"Max_Length"}=$Length; }
						my $motif_start = ($MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"start"});
						my $motif_end = ($MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"end"});
						my $name = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Gene_name"};
						# WRITING ORF ID, GENE NAME
						$drawing->createText(35,$header_size+$i_orf*$yside-(3.5*$space),-text=>"$i_orf/".(($#orfs_prot+1)-($#orfs_banned+1)),-anchor=>"e");
						$drawing->createRectangle(40,$header_size+$i_orf*$yside-(3.5*$space)-6,40+10,$header_size+$i_orf*$yside-(3.5*$space)+4,-fill=>"black",-outline=>"grey");
						$drawing->createText(60,$header_size+$i_orf*$yside-(3.5*$space),-text=>"$orf ",-anchor=>"w");
						$drawing->createText(120+length("$orf ")+3,$header_size+$i_orf*$yside-(3.5*$space),-text=>" $name ",-anchor=>"w",-fill=>"#3399FF");
						# DRAWING THE SEQUENCE WITH THE LENGTH IN AAs FOR THIS CURRENT SEQUENCE
						$drawing->createLine($xside,$header_size+($i_orf*$yside),($Length+$xside),$header_size+($i_orf*$yside),-capstyle=>"round",-width=>"3",-fill=>"red"); 
						$drawing->createText($Length+$xside+$space,$header_size+$i_orf*$yside,-text=>" ".length($PROTEOME->{$orf})." AAs",-anchor=>"w",-fill=>"black",-font => ['Courier', '10', 'bold']);
						# DRAWING THE MOTIF AND WRITING THE POSITIONS AND THE PATTERN OF THE MOTIF 
						$drawing->createOval($xside+$motif_start-5,$header_size+($i_orf*$yside)-7,$xside+$motif_end+5,$header_size+($i_orf*$yside)+7,-fill =>"purple",-outline=>"black",-width=>"2"); 
						$drawing->createText($xside+$motif_start,$header_size+$i_orf*$yside-$space,-text=>"$mot",-anchor=>"center",-fill=>"purple");
						$drawing->createText($xside+$motif_start,$header_size+$i_orf*$yside-(2*$space),-text=>"$motif_start-$motif_end ",-anchor=>"center",-fill=>"purple");
						$i_dom=1;
						# DRAW DOMAINS IN ORFS WITHIN THE SET
						foreach my $id (sort {$MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$a}->{"Start"} <=> $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$b}->{"Start"} || $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$b}->{"Length"} <=> $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$a}->{"Length"}} (keys (%{$MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}}))){
							$dom_space=$space*$i_dom;
							my $Type = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$id}->{"Type"};
							my $Annot = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$id}->{"Desc"}; 
							my $Start_dom = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$id}->{"Start"};
							my $End_dom = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$id}->{"End"};
							my $Max_end = $MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domain_Set"}->{$id}->{"Max_end_Domain_Set"};
							my $Eval=$MOTIFS{$mot}->{"ORF_Set"}->{$orf}->{"Domains_list_Orf_set"}->{$id}->{"E-value"};
							$drawing->createLine($Start_dom+$xside,$header_size+($i_orf*$yside+$dom_space),$End_dom+$xside,$header_size+($i_orf*$yside+$dom_space),-fill=>"blue",-width=>"2",-arrow=>"both");
							$drawing->createText($Start_dom-3+$xside,$header_size+($i_orf*$yside+$dom_space),-text=>"$Start_dom",-anchor=>"e");
							$drawing->createText($End_dom+$xside+1,$header_size+($i_orf*$yside+$dom_space),-text=>"$End_dom",-anchor=>"w");
							$drawing->createText($End_dom+$xside+(length("$End_dom")*8),$header_size+($i_orf*$yside+$dom_space),-text=>"$Type $Annot",-anchor=>"w",-fill=>"blue");
							$drawing->createText($End_dom+$xside+(length("$End_dom")+length("$Type $Annot"))*8,$header_size+($i_orf*$yside+$dom_space),-text=>"".sprintf("%.1e",$Eval),-anchor=>"w",-fill=>"red");
							$i_dom++;
						}
						# DESCRIPTION OF CURRENT ORF
						if( exists $Data_orf{$orf} ){
							$drawing->createRectangle(200+5,$header_size+$i_orf*$yside-(3.5*$space)-10,200+length($Data_orf{$orf}->{"Desc"})*5+20,$header_size+$i_orf*$yside-(3.5*$space)+10,-outline=>"black",-width=>"2");
							$drawing->createText(200+10,$header_size+$i_orf*$yside-(3.5*$space),-text=>$Data_orf{$orf}->{"Desc"},-anchor=>"w",-fill=>"grey35",-font => ['Arial','10','bold']);
						}
						$i_orf++;
					}
				}
				# DRAW MOTIFS IN ORFS WITHIN THE PROTEOME BUT NOT WITHIN THE SET
				foreach my $orf (@orfs_prot){
					if(!(exists $MOTIFS{$mot}->{"ORF_Banned"}->{$orf}) and !(exists $MOTIFS{$mot}->{"ORF_Set"}->{$orf})){
						my $Length = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Length"};
						if ($MOTIFS{$mot}->{"Max_Length"}<$Length){ $MOTIFS{$mot}->{"Max_Length"}=$Length; }
						my $motif_start = ($MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"start"});
						my $motif_end = ($MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"end"});
						my $name = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Gene_name"};
						# WRITING ORF ID, GENE NAME
						$drawing->createText(35,$header_size+$i_orf*$yside-(3.5*$space),-text=>"$i_orf/".(($#orfs_prot+1)-($#orfs_banned+1)),-anchor=>"e");
						$drawing->createRectangle(40,$header_size+$i_orf*$yside-(3.5*$space)-6,40+10,$header_size+$i_orf*$yside-(3.5*$space)+4,-fill=>"white",-outline=>"grey");
						$drawing->createText(60,$header_size+$i_orf*$yside-(3.5*$space),-text=>"$orf ",-anchor=>"w");
						$drawing->createText(120+length("$orf ")+3,$header_size+$i_orf*$yside-(3.5*$space),-text=>" $name ",-anchor=>"w",-fill=>"#3399FF");
						# DRAWING THE SEQUENCE WITH THE LENGTH IN AAs FOR THIS CURRENT SEQUENCE
						$drawing->createLine($xside,$header_size+($i_orf*$yside),($Length+$xside),$header_size+($i_orf*$yside),-capstyle=>"round",-width=>"3",-fill=>"red"); 
						$drawing->createText($Length+$xside+$space,$header_size+$i_orf*$yside,-text=>" ".length($PROTEOME->{$orf})." AAs",-anchor=>"w",-fill=>"black",-font => ['Courier', '10', 'bold']);
						# DRAWING THE MOTIF AND WRITING THE POSITIONS AND THE PATTERN OF THE MOTIF 
						$drawing->createOval($xside+$motif_start-5,$header_size+($i_orf*$yside)-7,$xside+$motif_end+5,$header_size+($i_orf*$yside)+7,-fill =>"purple",-outline=>"black",-width=>"2"); 
						$drawing->createText($xside+$motif_start,$header_size+$i_orf*$yside-$space,-text=>"$mot",-anchor=>"center",-fill=>"purple");
						$drawing->createText($xside+$motif_start,$header_size+$i_orf*$yside-(2*$space),-text=>"$motif_start-$motif_end ",-anchor=>"center",-fill=>"purple");
						$i_dom=1;
						# DRAW DOMAINS IN ORFS WITHIN THE PROTEOME BUT NOT WITHIN THE SET
						foreach my $id (sort {$MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$a}->{"Start"} <=> $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$b}->{"Start"} || $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$b}->{"Length"} <=> $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$a}->{"Length"}} (keys (%{$MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}}))){
							$dom_space=$space*$i_dom;
							my $Type = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$id}->{"Type"};
							my $Annot = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$id}->{"Desc"}; 
							my $Start_dom = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$id}->{"Start"};
							my $End_dom = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$id}->{"End"};
							my $Max_end = $MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domain_Prot"}->{$id}->{"Max_end_Domain_Prot"};
							my $Eval=$MOTIFS{$mot}->{"ORF_Prot"}->{$orf}->{"Domains_list_Orf_prot"}->{$id}->{"E-value"};
							$drawing->createLine($Start_dom+$xside,$header_size+($i_orf*$yside+$dom_space),$End_dom+$xside,$header_size+($i_orf*$yside+$dom_space),-fill=>"blue",-width=>"2",-arrow=>"both");
							$drawing->createText($Start_dom-3+$xside,$header_size+($i_orf*$yside+$dom_space),-text=>"$Start_dom",-anchor=>"e");
							$drawing->createText($End_dom+$xside+1,$header_size+($i_orf*$yside+$dom_space),-text=>"$End_dom",-anchor=>"w");
							$drawing->createText($End_dom+$xside+(length("$End_dom")*8),$header_size+($i_orf*$yside+$dom_space),-text=>"$Type $Annot",-anchor=>"w",-fill=>"blue");
							$drawing->createText($End_dom+$xside+(length("$End_dom")+length("$Type $Annot"))*8,$header_size+($i_orf*$yside+$dom_space),-text=>sprintf("%.1e",$Eval),-anchor=>"w",-fill=>"red");
							$i_dom++;
						}
						# DESCRIPTION OF CURRENT ORF
						if( exists $Data_orf{$orf} ){
							$drawing->createRectangle(200+5,$header_size+$i_orf*$yside-(3.5*$space)-10,200+length("".$Data_orf{$orf}->{"Desc"})*5+20,$header_size+$i_orf*$yside-(3.5*$space)+10,-outline=>"black",-width=>"2");
							$drawing->createText(200+10,$header_size+$i_orf*$yside-(3.5*$space),-text=>$Data_orf{$orf}->{"Desc"},-anchor=>"w",-fill=>"grey35",-font => ['Arial','10','bold']);
						}
					$i_orf++;
					}
				}
				$i++;
				my $pattern = $mot; $pattern =~ tr/./X/;
				# OUTPUT IN A POSTSCRIPT FILE
				$drawing->postscript(-file => "$resfile"."motif_".substr($rc,0,1)."_".$pattern.".ps",-pageheight=>5000,-pagewidth=>5000,-height=>$i_orf*$yside+50+$header_size,-width=>$xside+$space+$MOTIFS{$mot}->{"Max_Length"}+100,-colormode => "color");
			}
		}
	);
	MainLoop;
}else{
	my @path_to_Graphical = fileparse($0);
	print "\n".format_string(70,"#","#","#",">");
	print "\n".format_string(70,"#    For more details, consult the help to launch the program :"," ","#",">");
	print "\n".format_string(70,"#        $0       --help"," ","#",">");
	print "\n".format_string(70,"#","#","#",">")."\n\n";
	exit(0);
}
