#T LEFEBURE 16 juin 2011

package Pseudogene::StructAln;
use Bio::SimpleAlign;
use Bio::Seq;
use warnings;
use strict;
use base 'Exporter';

our @EXPORT = qw(get_struct hash_codon_struct get_intron_coord get_seq_ends count_stops_and_first_position merge_consensus typeCodon compareSeq mask_introns makeStructureLine convert_to_codon get_all_stop_position); 


#get a Seq object and returns a Seq object that contains the structure line
sub makeStructureLine {
    my ($seq, $intronCoord, $verbose, $nox, $convertion) = @_;
    my $string  = $seq->seq;
#     while($string =~ /\w-*\w-*\w-*/g) #sounds more difficult to use a regexp
    my @bases = split '', $string;
    my $n = 0; #codon position
    my $j = 0; #sequence position
    my $outstring = '';
    foreach my $b (@bases) {
	#if a gap
	if($b eq '-') {
	    #test if in an intron
	    if(isInIntron($j, $intronCoord)) {
		print "position $j is in an intron\n" if $verbose;
		$outstring .= 'i';
	    }
	    else {
		if($n < 1 or $n == 3) {
		    #extra codon insertion
		    $outstring .= '-';
		}
		elsif($n >= 1 and $n < 3) {
		    #intra codon gap
		    if($nox) {
			$outstring .= '-';
		    }
		    else { $outstring .= 'x' }
		}
		else { die "WTF: not possible!\n"; }
	    }
	}
	#if a base
	elsif ($b =~ /[A-Z]/i) {
	    ++$j;
	    #test if in an intron
	    if(isInIntron($j, $intronCoord)) {
		print "position $j is inside an intron\n" if $verbose;
		$outstring .= 'i';
	    }
	    else {
		++$n;
		$n = 1 if $n == 4; #reset the counter
		#recode the numbers into letters to fit into a protein code
		$outstring .= $$convertion{$n};
	    }
	}
	print "Position $j, codon $n, base $b\n" if $verbose;
    }
    (my $temp = $outstring) =~ s/-//g;
    my $length = length($temp );
#     print "$length\n";
    my $strucSeq = Bio::LocatableSeq->new(
	    -id => 'Struct',
	    -seq => $outstring,
	    -alphabet => 'protein',
	    -start => 1,
	    -end => $length,
    );
    return $strucSeq;
}

#test if a base is within an intron
sub isInIntron {
    my ($n, $intronCoord) = @_;
#     print "$n\n";
    foreach my $intr (keys %$intronCoord) {
	my $start = $intronCoord->{$intr}->{'start'};
	my $end = $intronCoord->{$intr}->{'end'};
# 	print "$n is $start -- $end ?\n";
	if($n < $start or $n > $end) {
	    return 0;
	}
	elsif ($n >= $start and $n <= $end) {
	    return 1;
	}
	else { die "WTF: this is not possible\n" }
    }
}

#mask an intron with Ns
sub mask_introns {
    my ($string, $introns, $seqstart, $verbose) = @_;
    my $masked_string = $string;
    foreach my $nintron (keys %$introns) {
	my $s = $introns->{$nintron}->{start};
	my $l = $introns->{$nintron}->{length};
	#only mask if the sequence has started before the end of the intron!
	next if( ($s + $l) < $seqstart);
	print "\tMasking intronic region: position $s to $introns->{$nintron}->{end}\n" if $verbose;
	substr($masked_string, $s, $l) = 'N' x $l;
    }
    return $masked_string;
}


sub merge_consensus {
    my ($left, $right) = @_;
    my $new = $left . $right;
    return $new;
}

#in this first round of codon typing, we get:
# noEvent
# del3
# delfs
# ins3_inter
# ins3_intra
# insfs_intra
# insfs_inter
# indel
#
#the deletion with be futher splitted after the merging 

sub typeCodon {
    my ($consensus, $verbose) = @_;
    #removes the g's (the paired gaps)
    $consensus =~ s/g//g;
    my $type = '';
    if($consensus =~ /^m+$/) { $type = 'noEvent' }    #mmm
    elsif($consensus =~ /^d+$/) { $type = 'del3' }    #ddd
    elsif($consensus =~ /^(i+)m+$/ or $consensus =~ /^m+(i+)$/) {
		if((length($1) % 3) == 0) { $type = 'ins3_inter' } #mmmiii
		else { $type = 'insfs_inter' } #immm
    }
    elsif($consensus =~ /^m+(i+)m+$/)  {
		if((length($1) % 3) == 0) { $type = 'ins3_intra' } #miiimm
		else { $type = 'insfs_intra' } #miimm
    }
    elsif($consensus =~ /^m*d+m*$/)  { $type = 'delfs' } #mmd   mmdmm
    elsif($consensus =~ /i/ and $consensus =~ /d/)  { $type = 'indel' } #mdiim
    elsif($consensus =~ /m*d+m+d+m*/) { 
	my $nd = () = $consensus =~ /d/g;
	if($nd % 3 == 0) { $type = 'del3' } # 12------3 / -------- -> ddmmmmmmmd
	else { $type = 'delfs' } # 12---3 / T----- -> mdmmmd
    }
    else{ print "I don't know what to do with this consensus: $consensus\nThis script does not handle multiple events per codon\n" }
	#exemple of missing match: 
	    #dmd -> 2 deletions
	    #mmiimii -> 2 insertions
	
    print "\t\ttype: $consensus --> $type\n" if $verbose;
    return $type;
}


#compare two strings and return some kind of consensus
sub compareSeq {
    my ($ref, $seq, $verbose) = @_;
    my $consensus = '';
    my @r = split '', $ref;
    my @s = split '', $seq;
    print "\t\tlength: ", length($ref), "\n" if $verbose;
    for (my $i = 0; $i <= length($ref)-1; ++$i) {
# 	print "$i: Compare $r[$i] versus $s[$i]\n" if $verbose;
		my $desc = compareBase($r[$i], $s[$i]);
		$consensus .= $desc;	
    }
    return $consensus;
}

#transform into m|i|d
sub compareBase {
    my ($b1, $b2) = @_;
    #transform x in - in the seq
    if($b2 =~ /x/i) { $b2 = '-' }
    if($b1 =~ /\d/ and $b2 =~ /\w/) { return 'm' }
    elsif($b1 =~ /\d/ and $b2 =~ /-/) { return 'd' }
    elsif($b1 =~ /[-x]/i and $b2 =~ /\w/) { return 'i' }
    elsif($b1 =~ /[-x]/i and $b2 =~ /-/) { return 'g' }
    else { die "WTF: comparision of $b1 and $b2 doesn't make sense!\n" }
}

#find the struct sequence and return it as a Seq object
sub get_struct {
	my ($aln, $structName) = @_;
	if( !$aln->isa("Bio::SimpleAlign") ) {
	     die("Expecting an alignment object not a $aln");
	}
	my @array = $aln->each_seq_with_id($structName);
	if(scalar(@array) > 1 ) { die "WTF: there is more than one structure track!\n"; }
	my $structSeq = $array[0];
	return $structSeq;
}

#Hash the the codon structure
sub hash_codon_struct {
    my %args = @_;	
    my $structSeq = $args{-STRUCTSEQ};
    if(!$structSeq->isa("Bio::Seq") and !$structSeq->isa("Bio::LocatableSeq") ) {
	     die("Expecting a sequence object not $structSeq"); 
	}
    my %codonPos;
    my %site2codon;
    my $strucString = $structSeq->seq;
    my $n = 0;
    while ($strucString =~ /1.+?3-*/go) {
	#this regex is used so that insertions
	#between codon are assigned at the end of a codon struct (ATT----)
	#introns are skipped
	++$n;
	my $end = pos($strucString)-1;
	my $start = $end - length($&)+1;
	my $length = length($&);
	print "Codon $n: $&, $start -> $end\n" if $args{-VERBOSE};
	$codonPos{$n} = {
	    start => $start,
	    end => $end,
	    seq => $&,
	    length => length($&),
	};
	foreach ($start..$end) {
	    $site2codon{$_} = $n
	}
    }

    if($strucString =~ /(12*)$/) {
	#add a last incomplete codon
	++$n;
	my $start = length($strucString) - length($1);
	my $end = length($strucString) -1;
	$codonPos{$n} = {
	    start => $start,
	    end =>  $end,
	    seq => $1,
	    length => length($1),
	};
	foreach ($start..$end) {
	    $site2codon{$_} = $n
	}
    }
    my $number_of_codons = $n;
    return (\%codonPos, \%site2codon, $number_of_codons);
}


#Get the intron coordinates from the struct seq
sub get_intron_coord {
    my ($structSeq) = @_;
    my %introns;
    my $n = 0;
    my $strucString = $structSeq->seq;
    while ($strucString =~ /i+/ig) {
	++$n;
	my $end = pos($strucString)-1;
	my $start = $end - length($&)+1;
	my $length = length($&);
	$introns{$n} = {
	    start => $start,
	    end => $end,
	    length => length($&),
	};
    }
    return(\%introns);
}

#Get sequence start and stop, in base and codon position
#it positions starts at 0
sub get_seq_ends {
    my ($string, $introns, $site2codonRef, $id, $verbose) = @_;
    #start
    my $seqstart = 0;
    if($string =~ /^([-X\?]+)\w/) { $seqstart = length($1) }

    #stop
    my $seqstop = length($string) -1;
    if($string =~ /\w[-X\?]+$/) {
	   $seqstop = length($string) - length($&)
    }
    #if in the intron, make the seq start|stop after the intron
    foreach my $nintron (keys %$introns) {
		my $s = $introns->{$nintron}->{start};
		my $l = $introns->{$nintron}->{length};
		my $e = $introns->{$nintron}->{end};
		if ($seqstart >= $s and $seqstart <= $e) {
			#the seq starts in this intron, move 
			#the start to the first site downstream of the intron
			$seqstart = $e+1;
		}
		if ($seqstop >= $s and $seqstop <= $e) {
			#the seq starts in this intron, move 
			#the start to the first site downstream of the intron
			$seqstop = $e+1;
		}
    }
    #convert into codons
    my $seq_start_codon = $site2codonRef->{$seqstart} or print "Problem: $id, seqtstart: $seqstart\n";
    my $seq_stop_codon = $site2codonRef->{$seqstop};
    print "Seq starts at codon $seq_start_codon (pos $seqstart) and stops at codon $seq_stop_codon (pos $seqstop)\n" if $verbose;

    return($seqstart, $seqstop, $seq_start_codon, $seq_stop_codon);
}

#Count STOPS, and find the position of the first one
sub count_stops_and_first_position {
    my ($string, $readingFrame, $introns, $seq_start_codon) = @_;
    my $intron_seq = $string;
    #remove the introns
    foreach my $nintron (keys %$introns) {
		my $s = $introns->{$nintron}->{start};
		my $l = $introns->{$nintron}->{length};
		my $e = $introns->{$nintron}->{end};
		substr($intron_seq, $s, $l) = '';
    }
    #remove X and gaps
    $intron_seq =~ s/[-X]//g;
    my $newseq = Bio::Seq->new(
		    -id => 'temp',
		    -seq => $intron_seq,
		    -alphabet => 'dna');
    my $protseq = $newseq->translate(-frame => $readingFrame);
    my $protstring = $protseq->seq;
    #count the stops
    my $nstops = 0;
    ++$nstops while( $protstring =~ /\*/g );
  
    my $first_stop = 0;
    #position of the first STOP
    if($protstring =~ /\*/g) {
	$first_stop = pos($protstring) + $seq_start_codon -1;
	$first_stop = $first_stop . '(' . pos($protstring). ')';
    }

    return($nstops, $first_stop);
}

#Find STOPS, and return their position
#13 nov 2012
sub get_all_stop_position {
    my ($string, $readingFrame, $introns, $seq_start_codon) = @_;
    my $intron_seq = $string;
    #remove the introns
    foreach my $nintron (keys %$introns) {
		my $s = $introns->{$nintron}->{start};
		my $l = $introns->{$nintron}->{length};
		my $e = $introns->{$nintron}->{end};
		substr($intron_seq, $s, $l) = '';
    }
    #remove X and gaps
    $intron_seq =~ s/[-X]//g;
    my $newseq = Bio::Seq->new(
		    -id => 'temp',
		    -seq => $intron_seq,
		    -alphabet => 'dna');
    my $protseq = $newseq->translate(-frame => $readingFrame);
    my $protstring = $protseq->seq;
  
    my @stoppos;
    #position of the first STOP
    while($protstring =~ /\*/g) {
		my $stopp = pos($protstring) + $seq_start_codon -1;
		push @stoppos, $stopp;
    }
    return(@stoppos);
}




#convert a struct codon to a 3 base codon
#by removing the deletions and insertions
sub convert_to_codon {
    my (%args) = @_;
#     my ($seq, $ref) = @_;
    my $newseq = '';
    my @r = split '', $args{-REF};
    my @s = split '', $args{-SEQ};
    my $indel;
    for (my $i = 0; $i <= length($args{-REF})-1; ++$i) {
# 	print "$i: Compare $r[$i] versus $s[$i]\n" if $verbose;
	my $desc = compareBase($r[$i], $s[$i]);
	if($desc eq 'm') { $newseq .= $s[$i] }
	elsif($desc eq 'd') { 
	    $newseq .= 'N';
	    $indel = 1;
	}
	elsif($desc eq 'i') {
	    $indel = 1;
	    next;
	}
	elsif ($desc eq 'g') { next }
    }
    if($indel and $args{-MASK}) { $newseq =~ s/[-\w]/N/g }
    return $newseq;
}





1;