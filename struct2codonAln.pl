#!/usr/bin/perl

#21Jul2011, T Lefebure
#- mask stop codon by default


use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Seq;
use lib '/home/tristan/src/perl'; 
use Pseudogene::StructAln;

my $mask;
my $outformat = 'fasta';
my $informat = 'fasta';
my $verbose;
my $refname = 'Ref';
my $structName = 'Struct';
my $addstruct;
my $stop;

GetOptions( 
    'stop' => \$stop,
    'mask' => \$mask, 
    'outformat:s' => \$outformat, 
    'informat:s' => \$informat,
    'verbose' => \$verbose,
    'addstruct' => \$addstruct );

my $usage = <<USE;
Usage: $0 <struct Alignment> <output>

This script transforms a structure alignment, as produced
by addCodonStructureToAln.pl, into a codon alignment.
It does so by removing any insertions and filling deletions
with Ns.

Options:
    -informat <> [$informat]
    -outformat <> [$outformat]
    -mask, will completely mask codons where insertions or deletions are found
    -stop, will keep stop codon, instead of masking them
    -addstruct, add the cleaned struct track
    -verbose
USE

if($#ARGV < 1) { 
    print $usage;
    exit;
}

my $in = Bio::AlignIO->new(
	-file => $ARGV[0],
	-format => $informat);
my $aln = $in->next_aln;


my $out = Bio::AlignIO->new(
	-file => ">$ARGV[1]",
	-format => $outformat);
my $newAlign = Bio::SimpleAlign->new();


######get the Struct track
print "Getting the track ($structName)\n" if $verbose;
my $structSeq = get_struct($aln, $structName);

#####hash the codon structure
my ($codonPos, $site2codon, $number_of_codons) 
	= hash_codon_struct(
		-STRUCTSEQ => $structSeq,
		-VERBOSE => $verbose );
print "Ref contains $number_of_codons codons\n";

####get the intron coordinates
print "Getting the intron coordinates\n" if $verbose;
my $introns = get_intron_coord($structSeq);

####foreach sequence
######for each codons
########get the struct seq and the codon sequence
########remove the insertions
########mask deletions

foreach my $seq ($aln->each_seq) {
    my $id = $seq->id;
    print "###Sequence $id\n" if $verbose;

    #skip the ref and struct seq
#     if(($id eq $refname) or ($id eq $structName)) { next }
    #only skyp the struct seq
    if ($id eq $structName) { next }

    my $string = $seq->seq;

    #find where the sequence really starts and stops,
    print "Getting the sequences ends ($id)\n" if $verbose;
    my ($seqstart, $seqstop, $seq_start_codon, $seq_stop_codon) = get_seq_ends($string, $introns, $site2codon, $id, $verbose);
    
    #foreach codon transform into a 3 base codon
    my $newSeqString = '';
    foreach my $n (sort { $a <=> $b } keys %$codonPos) {
		my $newcodon = '';
		my $codonSeq = substr($string, $codonPos->{$n}->{start}, $codonPos->{$n}->{length});
		my $ref = $codonPos->{$n}->{seq};
		#the following is used for incomplete codon
		my $codonSize = 0;
		++$codonSize while $ref =~ /\d/g;

		print "\t\t$ref versus $codonSeq\n" if $verbose;

		#if sequence did not start, or stopped, skip
		if ($n < $seq_start_codon) { $newcodon = '-' x $codonSize}
		elsif ($n > $seq_stop_codon) { $newcodon = '-' x $codonSize }
		else {
	# 	   $newcodon = convert_to_codon($codonSeq, $ref);
		   $newcodon = convert_to_codon(
			-SEQ => $codonSeq,
			-REF => $ref,
			-MASK => $mask );
		   print "\t\t--> $newcodon\t\t\n" if $verbose;

		   #mask the stop codons, unless requested
		  unless ($stop) {
		    $newcodon = 'NNN' if($newcodon =~ /TAA|TAG|TGA/)
		  }  
		}
		$newSeqString .= $newcodon;
    }
    print "\tNew seq is now $newSeqString\n" if $verbose;
    my $newSeq = Bio::LocatableSeq->new(
		    -seq => $newSeqString,
		    -id => $id);
    $newAlign->add_seq($newSeq);
}

if($addstruct) {
    my $structstring = $structSeq->seq;
    (my $newStructString = $structstring) =~ s/[-XI]//gi;
    my $newSeq = Bio::LocatableSeq->new(
		    -seq => $newStructString,
		    -id => $structName);
    $newAlign->add_seq($newSeq);
}

$newAlign->set_displayname_flat();
$out->write_aln($newAlign);

    

