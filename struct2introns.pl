#!/usr/bin/perl

#5Sep2011: fix the intron coordinates problem

use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Seq;
use lib '/home/tristan/src/perl'; 
use Pseudogene::StructAln;

my $outformat = 'fasta';
my $informat = 'fasta';
my $verbose;
my $refname = 'Ref';
my $structName = 'Struct';


GetOptions( 
    'outformat:s' => \$outformat, 
    'informat:s' => \$informat,
    'verbose' => \$verbose,
);

my $usage = <<USE;
Usage: $0 <struct Alignment> <output>

This script transforms a structure alignment, as produced
by addCodonStructureToAln.pl, into an intron alignment.
Only supports one intron.


Options:
    -informat <> [$informat]
    -outformat <> [$outformat]
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
my $structSeq = get_struct($aln, $structName);

#####hash the codon structure
# my ($codonPos, $site2codon, $number_of_codons) 
# 	= hash_codon_struct(
# 		-STRUCTSEQ => $structSeq,
# 		-VERBOSE => $verbose );
# print "Ref contains $number_of_codons codons\n";

####get the intron coordinates
my $introns = get_intron_coord($structSeq);

####slice the alignment given the intron coordinates
my $n = 0;
foreach my $nintron (keys %$introns) {
	++$n;
	if($n == 2) {
		print "Watch out: more than one intron! Only the first one has been removed\n";
		last;
	}
	my $s = $introns->{$nintron}->{start};
# 	my $l = $introns->{$nintron}->{length};
	my $e = $introns->{$nintron}->{end};
	$newAlign = $aln->slice($s + 1, $e + 1);
}

#remove the struc track
my $seqStruct = $newAlign->get_seq_by_id($structName);
$newAlign->remove_seq($seqStruct);

$newAlign->set_displayname_flat();
$out->write_aln($newAlign);

 
