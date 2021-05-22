#!/usr/bin/env perl 
use strict;
use warnings;

use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempfile);

my $cdsfolder = 'cds';
my $pepfolder = 'pep';
my $cdsext    = 'cds.fasta';
my $pepext    = 'aa.fasta';
my $outfolder = 'pairs';
my $temp      = '/tmp';
my $gzip      = 'pigz';
my $file;
GetOptions(
    'c|cds:s'     => \$cdsfolder,
    'p|pep:s'     => \$pepfolder,
    'o|out:s'     => \$outfolder,
    'i|in|file:s' => \$file,
    'pe|pepext:s' => \$pepext,
    'ce|cdsext:s' => \$cdsext,
    't|tmp:s'     => \$temp,
    'z|gzip:s'    => \$gzip,
    );
$file ||= shift || die "need a file";

open(my $fh => $file) || die "cannot open $file: $!";

mkdir($outfolder) unless -d $outfolder;
my (undef,undef,$filename) = File::Spec->splitpath($file);
my $basename;
if ( $filename =~ /(\S+)\.FASTA.tab$/ ) {
    $basename = $1;
} else {
    die "input report unexpected name: $filename"
}

my $pairct = 0;
my $pepfile = File::Spec->catfile($pepfolder,sprintf("%s.%s",$basename,$pepext));
my $pepdb = Bio::DB::Fasta->new($pepfile);
my $cdsfile = File::Spec->catfile($cdsfolder,sprintf("%s.%s",$basename,$cdsext));
my $tdir;
if ( ! -f $cdsfile && -f "$cdsfile.gz" ) {
    $tdir = File::Temp->newdir(TEMPLATE => 'TMPXXXX', DIR=> $temp);
    my $tfilename = File::Spec->catfile($tdir,sprintf("%s.%s",$basename,$cdsext));
    `$gzip -dc $cdsfile.gz > $tfilename`;
    $cdsfile = $tfilename;
}

my $cdsdb = Bio::DB::Fasta->new($cdsfile);

my $min_aligned = 0.60;
my %seen;
while(<$fh>) {
    next if /^\#/;
    my ($q,$t,$pid,$alnlen) = split;
    next if $q eq $t;		# skip self-vs-self
    my $sp;
    if ( $q !~ /\|/ ) {
	warn("expected a species name to be encoded in $q as >species|ACCESSION, will use $basename");
	$sp = $basename;
    } else {
	($sp) = split(/\|/,$q);
    }
    my $outfile = File::Spec->catdir($outfolder,$sp);
    mkdir($outfile);
    my $qseq = $pepdb->get_Seq_by_acc($q);
    my $tseq = $pepdb->get_Seq_by_acc($t); 
    if ( ! $qseq ) {
	warn("no qseq $q");
	next;
    }
    if ( ! $tseq ) {
	warn("no tseq $t");
	next;
    }
    # this is 1 selection criteria - make sure is at least (min_aligned)
    # aligned
    if(( $alnlen / $qseq->length)  < $min_aligned ) {
	next;
    }
    # only dump a gene/protein in a single pair, if it participated
    # as a 'query' or 'target' in any other results, skip it
    # This could have the effect of picking 2nd best hit for something
    # not sure should further restrict how this allows a pair to be written
    next if $seen{$q}++ || $seen{$t}++;
    $outfile = File::Spec->catfile($outfile,sprintf("p%06d",$pairct++));
    my $cdsout = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile.cds");
    my $pepout = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile.pep");
    my $qcds = $cdsdb->get_Seq_by_acc($q);
    my $tcds = $cdsdb->get_Seq_by_acc($t);

# deal with potential stop codons in these sequences for various reasons
# JGI doesn't produce clean CDS files, not sure why, recode these
# positions as 'X' then encode these as gaps after the alignment is performed
# see the mr_trans.pl script in this folder

    if( $qseq->seq =~ /\*/ ) { 
	my $s = $qseq->seq;
	my $scds = $qcds->seq;
	if( $s =~ s/\*$// ) {	# drop trailing stop codon completely
	    $s = substr($s,0,$qseq->length - 1);
	    $scds = substr($scds,0,$qcds->length - 3)
	}
	# any other stop codon should be recoded as an X to be handled by
	# the aligner and then swapped out by mrtrans
	$s =~ s/\*/X/g;
	# have to rewrite the seq object because Bio::DB::Fasta essentially produces
	# read-only objects , otherwise would just update ->seq() in place   
	$qseq = Bio::PrimarySeq->new(-id => $qseq->id, -seq => $s, 
				     -desc => $qseq->description);
	$qcds = Bio::PrimarySeq->new(-id => $qcds->id, -seq => $scds, 
				     -desc => $qcds->description);
    } 

    if( $tseq->seq =~ /\*/ ) {
	my $s = $tseq->seq;	# the protein
	my $scds = $tcds->seq;	# the CDS

	if( $s =~ s/\*$// ) {	# drop trailing stop codon completely
	    $s = substr($s,0,$tseq->length - 1);
	    $scds = substr($scds,0,$tcds->length - 3);
	}
	# any other stop codon should be recoded as an X to be handled by  
	# the aligner and then swapped out by mrtrans
	$s =~ s/\*/X/g;
	# have to rewrite the seq object because Bio::DB::Fasta essentially produces
	# read-only objects , otherwise would just update ->seq() in place   

	$tseq = Bio::PrimarySeq->new(-id   => $tseq->id, 
				     -seq  => $s, 
				     -desc => $tseq->description);
	$tcds = Bio::PrimarySeq->new(-id   => $tcds->id,
				     -seq  => $scds, 
				     -desc => $tcds->description);
    }
    $cdsout->write_seq($qcds,$tcds);
    $pepout->write_seq($qseq,$tseq);
}
