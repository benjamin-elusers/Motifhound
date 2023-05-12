# --- fichier  String_Print_Utils.pm ---
package String_Print_Utils;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(&my_print &format_string &fasta_format &print_Time_usage &write_Time_usage &match_first_positions &match_all_start_positions &match_all_end_positions &Write_sel_seq);
use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;

sub my_print {
# Build a new character chain from a string by concatenating to itself until it reaches the length "$len"
	my ($str,$len)=@_;
	my $new_str="";
	for(my $i=0;$i<$len;$i++){
		$new_str.=$str;
	}
	return($new_str);
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

sub write_Time_usage{
# Print Time Usage through the Filehandler <$fh>
	my ($fh,$Ref_List)=@_;
	my @list2=@{$Ref_List};
	my @list=times;
	my $prog = $list[0]-$list2[0]; 	my $progsys = $list[1]-$list2[1];
	my $child = $list[2]-$list2[2]; 	my $childsys = $list[3]-$list2[3];
	print $fh "Time used by the program (system) : $prog  s ($progsys  s )\n";
	print $fh "Time used by the child   (system) : $child s ($childsys s )\n";
	return(\@list);
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

sub match_first_positions {
	my ($regex, $string) = @_;
	return if not $string =~ /$regex/;
	return ($-[0], $+[0]);
}


sub match_all_start_positions {
	my ($regex, $string) = @_;
	return if not $string =~ /$regex/;
	return (\@-);
}

sub match_all_end_positions {
	my ($regex, $string) = @_;
	return if not $string =~ /$regex/;
	return (\@+);
}


sub Count_cdef {
# Count fixed positions (non wildcards) in a string 
	my ($motif) = @_;
	my $flag_cdef=0;
	my @tab = split("",$motif);

	for(my $i=0;$i<($#tab+1);$i++){
		if ($tab[$i] ne "." or $tab[$i] ne "X"){
			$flag_cdef++;
		}
	}
	return ($flag_cdef);
}

1;
