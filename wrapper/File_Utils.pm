# --- fichier  File_Utils.pm ---
package File_Utils;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(&header_over &footer_over &header_enum &footer_enum  &possibly_listIDs &possibly_fasta &Directory_verification &File_verification &File_or_Directory &Load_data_BLAST &Load_FASTA &Load_data_DISO &guess &Extract_Last_Line &load_Gene_ID_desc &load_motifs &load_Sharing_PFAM_Domains &load_PFAM_annot &load_data_domain_SUPERFAMILY &Load_Sequence &load_IDS &load_TERMS &Count_cdef &write_IDS &List_IDS_to_Fasta  &Write_Seq_to_FASTA &compath);
use strict;
use File::Basename;
use autodie qw(open close); # open/close succeed or die
use MIME::Base64;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use List::Util qw(max reduce);


sub compath {
# Extract the longest common path
# returns the longest common path between a list of paths

	my ($sep, @paths, %hash) = @_;
	# Tokenize and tally subpaths
	foreach (@paths) {
		my @tok = split $sep, substr($_,1);
		++$hash{join $sep, @tok[0..$_]} for (0..$#tok); }
	# Return max length subpath or null
	my $max = max values %hash;
	return '' unless $max == @paths;
	my @res =  grep {$hash{$_} == $max} keys %hash;
	return $sep . reduce { length $a > length $b ? $a : $b } @res;
}


sub possibly_listIDs {
# Check if the content of a file is a list of IDS
# returns a text message
# returns 2 if the file corresponds to a list of IDs and -1 if the file is not corresponding to a list of IDs
# returns the number of sequences in the file

	my ($path) = @_;
	my $lineno=0; my $RC=0; my $ERROR_MSG="";
	open(FILEIN,"$path");
	while (my $line =<FILEIN>) {
		chop ($line); $RC+=($line =~ /^[^>]\S+/); $lineno++;
		if($RC!=$lineno){
			close(FILEIN);
			$ERROR_MSG.="------------------------->	   Not a list of IDs\n";
			$ERROR_MSG.="------------------------->	      line $lineno : $line\n";
			$ERROR_MSG.="------------------------->	      This line does not contain any ID Id.\n";
			return($ERROR_MSG,-2,0);
		}
	}
	close(FILEIN);
#	if($lineno < 2){
#		$ERROR_MSG.="------------------------->	Empty File or One-Line-File ($lineno lines) => Aborted\n";
#		return ($ERROR_MSG,0,$lineno);
#	}else{
		$ERROR_MSG.="------------------------->	List of IDs...OK\n";
		return($ERROR_MSG,2,$lineno);
#	}
}

sub possibly_fasta {
# Check the FASTA format of a file
# returns a text message
# returns 1 if the file is fasta formatted and -1 if the file is not fasta formatted
# returns the number of sequences in the file
	my ($path) = @_;
	my $lineno=0; my $RC=0; my $Nbseq=0;my $ERROR_MSG="";

	open(FILEIN,"$path");
	while (my $line =<FILEIN>) {
		chop ($line);
		$RC+=(($lineno != 0 && $line =~ /^[A-IK-NP-Z]+$/i) || $line =~ /^>\S+/);
		$Nbseq+=grep(/^>/, $line);
		$lineno++;

		if($RC!=$lineno){
			close(FILEIN);
			$ERROR_MSG.="------------------------->	   Not a FASTA File\n";
			$ERROR_MSG.="------------------------->	      $lineno : $line\n";
			$ERROR_MSG.="------------------------->	      This line does not correspond to FASTA format.\n";
			return ($ERROR_MSG,-1,0);
		}
	}
	close(FILEIN);
#	if($lineno < 2){
#		$ERROR_MSG.="------------------------->	Empty File or One-Line-File ($lineno lines) => Aborted\n";
#		return ($ERROR_MSG,0,$Nbseq);
#	}else{
		$ERROR_MSG="------------------------->	FASTA formatted...OK\n";
		return($ERROR_MSG,1,$Nbseq);
#	}
}

sub Directory_verification {
# Test wheter the path given corresponds to an existing directory
# => returns 0 if the path is an existing directory 
# => returns -1 if the path is not a directory or if it doesn't exist

# File Test Operators
#-e  file exists
#-d  file is a directory
	my ($path) = @_;
	
	if(-e $path and -d $path){ 
		return (1,"OK [Directory]\n"); 
	}elsif(-e $path and !(-d $path)){ 
		return (-1,"Not a directory.\n");
	}else{
		return (-1,"Path does not exist and/or does not correspond to a an existing path.\n");
	}
}

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

sub File_or_Directory {
# Check if the path given is a file or a directory 
# => returns 1 for a directory
# => returns 0 for a file
	my ($path) = @_;

	my $Path_rc=""; my $TEXT=""; my $TEXT2="";
	($Path_rc,$TEXT)=File_verification($path);
	if($Path_rc != 0){
		($Path_rc,$TEXT2)=Directory_verification($path);
		print $TEXT2;
		if($Path_rc != 1){ print "# The path \"$path\" does not exist ! #\n"; exit ($Path_rc); }
		return($Path_rc);
	}else{ 
		print $TEXT;
		return($Path_rc);
	}
}

sub List_IDS_to_Fasta {

	my ($file,$Reference_Proteome,$IDMapping,$Logdir,$Scriptsdir)=@_;
	my ($filename,$filedir,$filext) = fileparse($file,qr{\..*});
#	my $str = encode_base64($filename);
	my $file_newname=$filedir.$filename.".tmp";
#	print STDERR "$Scriptsdir/Select_IDS.pl --IDS $file --Fastafile $Reference_Proteome --out $file_newname --logdir $Logdir --m $IDMapping --nov";
	`$Scriptsdir/Select_IDS.pl --IDS $file --Fastafile $Reference_Proteome --out $file_newname --logdir $Logdir --m $IDMapping --nov`;
	print "- New Set of Sequences  :	\"".$file_newname."\" => "; my ($Path_rc,$newname,$nseq)=guess($file_newname,$Reference_Proteome,$Scriptsdir,$Logdir);
	return($Path_rc,$newname,$nseq);
}

sub guess {

	my ($path)=@_;
	my $Path_rc="";
	# This function returns 1 if the path corresponds to a directory, 0 if the path corresponds to a file
	$Path_rc = File_or_Directory($path); 
	if($Path_rc == 0){ # If the path is a file
		my ($TEXT,$rc,$nseq)=possibly_fasta($path);
		my ($TEXT2,$rc2,$nseq2)=possibly_listIDs($path);
#		if($rc==0 or $rc2==0){ # File is empty
#			print $TEXT; exit(-1);
		if($rc==1 and $rc2 < 0){ # File is a sequence file Fasta Formatted
			print $TEXT; return($rc,$path,$nseq);
		}elsif($rc2==2 and $rc < 0){ # File is a list of IDs
			print $TEXT2; return($rc2,$path,$nseq2);
		}elsif(($rc+$rc2)==3){ # File format is ambiguous
			print STDERR "# Impossible to define the type of the file ! (Ambiguous Format) # (file = \"$path\")\n";
			print STDERR "# Check the format of the input sequence file, and try again. #\n\n";
			print STDERR "NB  : A list of IDs should contain uniques IDs (in uppercase) followed by a carriage return (\"\\n\") in each line.\n";
			print STDERR "ex.\nYFR044C\nYFR045C\nYWC125C\n...\n";
			print STDERR "\n\n";
			print STDERR "NB2 : A FASTA formatted file should contained in each line either : \n";
			print STDERR " 	-an unique ID preceeded by a \">\" and followed by a carriage return (\"\\n\")\n";
			print STDERR "                                    OR                                           \n";
			print STDERR " 	-a sequence composed of uppercase letters which should correspond to the Amino Acid Alphabet [A-IK-NP-Z]\n";
			print STDERR "ex.\n>YFR044C\nLSKNCQPTKFHEMLNEFARDGR...\n>P17480\nMFEAINOIADCPPCAE...\n...\n";
			print STDERR "\n\n";
			exit(-2);
		}elsif(($rc+$rc2)==-3){ # File format is unrecognized
			print "------------------------------>	The input file is neither a list of IDs nor a FASTA formatted file.\n";
			print $TEXT.$TEXT2;
			print STDERR "# Check the format of the input sequence file and try again.                                             #\n";
			print STDERR "# (Don't forget that the Set sequence file must be FASTA formatted or must correspond to a list of IDs) #\n";
			exit(-3);
		}
	}
	return($Path_rc,$path,0);
}

sub header_enum {
# Write a header for enumeration result file
	my ($fname,$File_old,$File,$Kocc,$nseq,$Ref_Filters,$dmin,$dmax)=@_;
	my @Filters = @{$Ref_Filters};
	open(HEADER,">","$fname") or die "$fname";
	if($File_old ne $File){
		print HEADER "# Original File          : $File_old\n";
	}
	print HEADER "# Fasta                  : $File\n";
	print HEADER "# nseq                   : $nseq\n";
	print HEADER "# K-occ                  : $Kocc\n";
	print HEADER "# Filters                : -Homology=$Filters[0]    -Disorder=$Filters[1]   -Overlap=$Filters[2]\n";
	print HEADER "# Non-wildcard positions : -Min=$dmin -Max=$dmax\n";
	close(HEADER);
}

sub footer_enum {
# Write a footer for enumeration result file
	my ($fname)=@_;
	no warnings;
	open(FOOTER,">>","$fname") or die "$fname";
	print FOOTER "# END OF ENUMERATION #\n";
	close(FOOTER);
	use warnings;
}

sub header_over {
# Write a header for Over representation result file
	my ($fname,$Old_Set,$Set,$Old_Prot,$Prot,$Kocc,$Nseq,$nseq,$Ref_Filters,$dmin,$dmax)=@_;
	my @Filters = @{$Ref_Filters};
	open(HEADER,">","$fname") or die "$fname";
	if($Old_Set ne $Set){ print HEADER "# Original Set                : $Old_Set\n"; }
	if($Old_Prot ne $Prot){ print HEADER "# Original Proteome           : $Old_Prot\n"; }
	print HEADER "# Proteome ($Nseq sequences)  : $Prot\n";
	print HEADER "# Set      ($nseq sequences)  : $Set\n";
	print HEADER "# K-occ                       : $Kocc\n";
	print HEADER "# Filters                     : -Homology=$Filters[0]    -Disorder=$Filters[1]   -Overlap=$Filters[2]\n";
	print HEADER "# Non-wildcard positions      : -Min=$dmin -Max=$dmax\n";
	close(HEADER);
}

sub footer_over {
# Write a footer for Over representation result file
	my ($fname)=@_;
	no warnings;
	open(FOOTER,">>","$fname") or die "$fname";
	print FOOTER "# END OF OVER REPRESENTATION #\n";
	close(FOOTER);
	use warnings;
}

sub Extract_Last_Line {
# Extract and return Last Line of a file 
	my ($File) = @_;
	my $line="";
	no warnings;
	open(FIN,"<",$File) or die "Unable to open $File : $!\n";
	foreach my $tmp (reverse(<FIN> )) {
		$line=$tmp;
		last;
	}
	close(FIN);
	use warnings;
	chop($line);
	return ($line);
}

sub load_IDS{
# Returns a list of IDS from a FASTA formatted file 
	my ($file) = (@_);
	open(IN,"<$file") or print("Could not open $file\n") ;
	my $line=""; my @IDS;
#	print STDERR "$file\n";
	while( $line = <IN>){
#		print STDERR "$line\n";
		if( $line =~ /(^>?)(\S+)/){
#			print STDERR "$2\n";
			push(@IDS,$2);
		}
	}
	close(IN);
	return(\@IDS);
}

sub Load_Sequence {
# Load Fasta Sequence from a file input
# Discard sequences containing "*" (symbol corresponding to STOP codons)
	my ($fin) = @_;
	my @temp_ID;my $i=0;my $cpt=0;my $temp_seq=""; my @sel_seq;my @sel_ID;my %S;my $flag_star = 0;
	open(YEAST,$fin);
	while(my $line = <YEAST>){
		if($line =~ /^>/){
			$temp_ID[$i] = $line ;
			if($flag_star == 0){
				if($i!=0){
					$sel_ID[$cpt] = $temp_ID[$i-1]; 
					chop($temp_ID[$i-1]);
					my $ID = substr($temp_ID[$i-1],1,9);
					$ID =~ s/\s+//;
#					print $ID."\n";
					my %cl_temp ; 
					$cl_temp{"Seq"}=$temp_seq;
					$S{"$ID"}={%cl_temp}; 
					$cpt++;
				}
			}
			$flag_star = 0;
			$i++; $temp_seq="";
		}else{
			if($flag_star == 0 and $line !~ /\*/){  $temp_seq.=$line;  }else{  $flag_star = 1;  }
		}
	}

	$sel_ID[$cpt] = $temp_ID[$i-1];
	my $ID = substr($temp_ID[$i-1],1);    chop ($ID);
	my %cl_temp ; $cl_temp{"Seq"}=$temp_seq;
	$S{"$ID"}={%cl_temp};
	$cpt++;

	close(YEAST);
	return (\%S);
}

sub Load_FASTA($){
# Load sequences from a file to a Hash (path of the file given in arguments)
	my ($file) = (@_);
	my($ID, $pdb, $seq) = ("None","None", "None");
	my (%hash);

	open(IN, "<$file") or die("Could not open $file\n") ;
	while(<IN>){
		if( $_ =~ /^[^>].* .*/){
			chomp($_);
			print $_." not considered\n";
		} else {
			if($_ =~ /^>([^ ]+).*$/){
				$hash{$ID}->{"Seq"} = $seq ;
				$ID = $1;
				chomp($ID);
				$ID =~ s/ +$//;
				$seq = "" ;
			} else {
				chomp($_);
				$seq .= $_ ;
				$seq =~ s/\*$//;
			}
		}
	}
	$hash{$ID}->{"Seq"} = $seq ;
	delete($hash{"None"});
	close(IN);
	return (\%hash) ;
}

sub write_IDS{
# Write a file from a Hash, corresponding to a list of IDS sorted alphabetically 
	my ($Ref_IDS,$fname) = (@_);
	open(OUT, ">$fname") or print("Could not open $fname for writing.\n") ;
	foreach my $ID ( sort { $a cmp $b } (keys %{$Ref_IDS})){
#		print STDERR "$ID\n";
		print OUT "$ID\n";
	}
	close(OUT);
	return (0);
}


sub Load_data_DISO {
# Load the results of Disorder Predictions file and store them into a hash
	my ($file) = (@_);

	my %H=(); my $i=0;
	open(IN, "<$file") or die "Could not open $file\n";

	my @lines=<IN>;
	for($i=1;$i<($#lines+1);$i++){
		if($i%100000 == 0){
				print "Reading Disorder Predictions File...$i/".($#lines+1)." lines.\r";
		}
		my @c = split("\t",$lines[$i]);
		if($c[4] eq '.'){ $c[4] = 0; }else{ $c[4] = 1; }
		
		$H{$c[1]}->{"Seq"} 		.= $c[2];
		$H{$c[1]}->{"Dis"} 		.= $c[4];
		$H{$c[1]}->{"Dis_val"}	.= $c[5];
		$H{$c[1]}->{"CF_Dis"}	.= $c[6];
	}
	print "--- Reading Disorder Predictions File : Finished ---\n";
	close(IN);
	return (\%H);
}

sub Load_data_BLAST {
# Load the results of Basic Local Alignment when the percentage of identity is superior to the treshold 
# Returns a Hash of homologuous regions in a proteome
	my ($file) = (@_);

	my %H=();
	my $line="";
	my $cpt = 1;

	open(IN, "<$file") or print("Could not open $file\n") ;
	while($line=<IN>){
		my @c = split("\t",$line);
		if(exists $H{$c[0]}->{$c[1]}->{$cpt}){   $cpt++;   }else{   $cpt=1;   }
		$H{$c[0]}->{$c[1]}->{$cpt}->{"percId"}    = $c[2];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"alnLen"}    = $c[3];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"ID1_start"} = $c[6];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"ID1_end"}   = $c[7];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"ID2_start"} = $c[8];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"ID2_end"}   = $c[9];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"Evalue"}    = $c[10];
		$H{$c[0]}->{$c[1]}->{$cpt}->{"bit_score"} = $c[11];
#		if($i % 10000 == 0){ print $i."\n"; }
#		print $i."\n";
	}
#	print "\n";
	close(IN);

	return (\%H);
}

sub load_data_domain_SUPERFAMILY($){
# Store in a hash, informations about domains ( obtained from Superfamily database)
	my ($file) = (@_);

	my $cpt=0;
	my $motif=""; my $line="";
	my @c;
	my %H=();
	open(IN, "<$file") or print("Could not open $file\n") ;
	while( $line = <IN>){
		my %cl_tmp=();
		@c = split("\t",$line);
		$cl_tmp{$c[2]}->{"Model ID"} = $c[2];
		my @multipos=split(",",$c[3]);
		my $n=1; my @Positions;
		foreach my $p (@multipos){
			my @pos=split("-",$p);
			push(@Positions,$pos[0]);
			push(@Positions,$pos[1]);
		}
		$cl_tmp{$c[2]}->{"Positions"} = \@Positions;
		$cl_tmp{$c[2]}->{"E-val"} = $c[4];
		$cl_tmp{$c[2]}->{"Superfamily_id"} = $c[5];
		$cl_tmp{$c[2]}->{"Superfamily_desc"} = $c[6];
		$cl_tmp{$c[2]}->{"Family E-val"} = $c[7];
		$cl_tmp{$c[2]}->{"Family_id"} = $c[8];
		$cl_tmp{$c[2]}->{"Family_desc"} = $c[9];
		$cl_tmp{$c[2]}->{"Most_similar_struct"} = $c[10];
		$H{$c[1]}={%cl_tmp};
	}
	close(IN);
	return (\%H);
}

sub load_PFAM_annot($){
# Store in a hash, informations about domains ( obtained from PFAM database)
	my ($file) = (@_);

	my $cpt=1; my $i=0;
	my %H=();
	my $line="";
	open(IN, "<$file") or print("Could not open $file\n") ;

	while($line=<IN>){
		my @c = split(" ",$line);
#		print STDERR "@c\n";
		if(exists $H{$c[0]}->{$c[5]}){  $cpt++;  }else{  $cpt=1;  }
		$H{$c[0]}->{$c[5]}->{$cpt}->{"ali_start"}    = $c[1]; 
		$H{$c[0]}->{$c[5]}->{$cpt}->{"ali_end"}      = $c[2];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"env_start"}    = $c[3]; 
		$H{$c[0]}->{$c[5]}->{$cpt}->{"env_end"}      = $c[4];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"Annot"}        = $c[6]; 
		$H{$c[0]}->{$c[5]}->{$cpt}->{"Type"}         = $c[7];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"hmm_start"}    = $c[8]; 
		$H{$c[0]}->{$c[5]}->{$cpt}->{"hmm_end"}      = $c[9]; 
		$H{$c[0]}->{$c[5]}->{$cpt}->{"hmm_length"}   = $c[10];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"bit_score"}    = $c[11];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"E-value"}      = $c[12];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"significance"} = $c[13];
		$H{$c[0]}->{$c[5]}->{$cpt}->{"clan"}         = $c[14];
#		if($i % 10000 == 0){ print $i."\n"; }
#		print $i."\n";
	}
	close(IN);
	return (\%H);
}

sub load_Sharing_PFAM_Domains($$$){
# Store in a hash, informations about domains ( obtained from PFAM database)
	my ($file,$Pfam_annot,$Scriptsdir) = (@_);
	
	if(!(-f $file)){
		$file = "IDs_sharing_same_Pfam_domains.txt";
		`perl $Scriptsdir/IDs_Sharing_domains.pl --pfam --Pfam_annot $Pfam_annot --out $file --nov`;
	}
	
	if(-f $file){
		my $line=""; my @c; my %H=();
		open(IN, "<$file") or print("Could not open $file\n") ;
		while( $line = <IN>){
			@c = split(" : ",$line);
			if($#c != 0){ 
				$H{$c[0]}->{"ID_sharing_same_dom"} = $c[1];
			}else{
				$H{$c[0]}->{"ID_sharing_same_dom"} = " ";
			}
			$H{$c[0]}->{"#ID"} = 0;
			if($#c > 0){
				my @c2 = split(" ",$c[1]);
				$H{$c[0]}->{"#ID"} = ($#c2+1);
			}
		}
		close(IN);
		return (\%H);
	}
}

sub Count_cdef {
# Count fixed positions (non wildcards) in a string 
	my ($motif) = @_;
	my $flag_cdef=0;
	my @tab = split("",$motif);

	for(my $i=0;$i<($#tab+1);$i++){
		if ($tab[$i] ne "." and $tab[$i] ne "X"){
			$flag_cdef++;
		}
	}
	return ($flag_cdef);
}

sub load_motifs($){
# Load the content of over representation results file
	my ($file) = (@_);

	my $cpt=0;
	my $motif=""; my $line="";
	my @c;
	my %H=();
	open(IN, "<$file") or print("Could not open $file\n") ;
	while( $line = <IN>){
		if( $line =~ /^[^#]/){
			chop($line);
			@c = split("\t",$line);
			$cpt++;
				$H{$c[0]}->{"Len"} = length($c[0]);
				$H{$c[0]}->{"Cdef"} = Count_cdef($c[0]);
				$H{$c[0]}->{"S_set"} = $c[1];
				$H{$c[0]}->{"S_pop"} = $c[2];
				$H{$c[0]}->{"N_set"} = $c[3];
				$H{$c[0]}->{"N_pop"} = $c[4];
				$H{$c[0]}->{"Pval"} = $c[5];
#				$H{$c[0]}->{"nbdegen"} = $c[6];
#				$H{$c[0]}->{"Pval_adj"} = $c[7];
#				$H{$c[0]}->{"nwild"} = $c[8];
		}
	}
	close(IN);
	return (\%H,$cpt);
}

sub load_Gene_ID_desc($){
# Load Gene name and description from IDS
	my ($file,$link) = (@_);
	my $cpt=0;
	my $motif=""; my $line="";
	my @c;
	my %H=();

#	print STDERR "$file\n";
#	print STDERR "$link\n";

	open(IN, "<$file") or print("Could not open $file\n") ;
	while( $line = <IN>){
		chop($line);
		@c = split("\t",$line);
#		print STDERR "@c\n";
		$H{$c[0]}->{"Gene"} = $c[1];
		if(length($c[2]) > 0 ){
			$H{$c[0]}->{"Desc"} = $c[2];
		}else{
			if($link ne "" ){
				$H{$c[0]}->{"Desc"} = $link.$c[0] ;
			}else{
				$H{$c[0]}->{"Desc"} = '-' ;
			}
		}
	}
	close(IN);
	return (\%H);
}

sub load_TERMS($){

	my ($file) = (@_);
	my $line=""; my @c; my %H=();

	open(IN, "<$file") or print("Could not open $file\n") ;
	while( $line = <IN>){
		@c = split("\t",$line);
		$H{$c[0]}->{"Desc"} = $c[1];
	}
	close(IN);
	return (\%H);
}


sub fasta_format($){
# Takes a string, inserts a "\n" every 70 character. 
	my ($seq) = (@_);
	my $i=0;
	my $res="";

	while($i < length($seq) ){
		$res .= substr($seq,$i,70)."\n";
		$i+=70;
	}
	return $res;
}

sub Write_Seq_to_FASTA {
# Write Sequences from a Hash, to FASTA files. (Possibility to sort IDs in alphabetical order before writing)
	my ($Ref_Seq,$foutpath,$ordered) = @_;
	my %Seq = %{$Ref_Seq};
	my @IDS = (keys %Seq);

	if($ordered){ @IDS = sort { $a cmp $b } (keys %Seq); }
	open(OUT,">$foutpath");

	foreach my $ID (@IDS){
		if ( exists $Seq{$ID}->{"Seq"} ){
			print OUT ">$ID\n";
			print OUT fasta_format($Seq{$ID}->{"Seq"});
		}else{
			print STDERR "# No sequence found for $ID. #";
		}
	}
	close(OUT);
}

1;
