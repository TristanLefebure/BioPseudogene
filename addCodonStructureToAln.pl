#!/usr/bin/perl

#23Mai2011 Tristan Lefebure

use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Seq;
use lib '/home/tristan/src/perl'; 
use Pseudogene::StructAln;

my $informat = 'mase';
my $outformat = 'fasta';
# my $update;
my $refname = 'Ref';
my $first;
# my $borders;
my @introns;
my $fst;
my $nox;
my $verbose;

GetOptions(
    'informat:s' => \$informat,
    'outformat:s' => \$outformat,
#     'update' => \$update,
    'refname:s' => \$refname,
    'first' => \$first,
#     'borders' => \$borders,
    'intron=s' => \@introns,
    'fst' => \$fst,
    'nox' => \$nox,
    'verbose' => \$verbose,
);

my $usage =
"Usage: $0 <input alignment> <output alignment>
Options:
    -informat <> [$informat]
    -outformat <> [$outformat]
    -intron <xxx,xxx>, gives intron position on the ref
     sequence
    -refname <>, set the reference sequence name [$refname]
    -first, instead of looking for the refname, uses 
     the first sequence, whatever its name
    -fst, uses First-Second-Third instead of 123
    -nox, don't mark intra-codon frameshift with x
    -verbose\n";

if($#ARGV < 1) {
    print $usage;
    exit;
}

#if conversion is needed:
my $convertion = {
	1 => '1',
	2 => '2',
	3 => '3',
};
if($fst) {
    $convertion = {
	    1 => 'F',
	    2 => 'S',
	    3 => 'T',
    };
}

#format the introns
my $intronCoord = {}; #hashref!
my $i = 0;
foreach (@introns) {
    ++$i;
    my @f = split ',';
    $intronCoord->{$i} = {
		start => $f[0],
		end => $f[1],
    };
#     print "Intron start position: $intronCoord->{$i}->{start}\n";
}

#import/export
my $in = Bio::AlignIO->new(
	-file => $ARGV[0],
	-format => $informat,
);

my $out = Bio::AlignIO->new(
	-file => ">$ARGV[1]",
	-format => $outformat,
);


my $aln = $in->next_aln();
$aln->set_displayname_flat();

my $n = 0;
foreach my $seq ($aln->each_seq) {
    ++$n;
    my $id = $seq->id;
    if(($id eq $refname) or ($n == 1 and $first)) {
	my $structSeq = makeStructureLine($seq, $intronCoord, $verbose, $nox, $convertion);
	$aln->add_seq(-SEQ=>$structSeq, -ORDER=>1);
	$aln->set_displayname_flat();
	$out->write_aln($aln);
	last;
    }
}

