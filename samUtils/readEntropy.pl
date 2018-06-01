#!/usr/bin/env perl

=head1 NAME

readEntropy.pl - Calculate DNA sequence entropy and/or redundancy per sequence

=head1 SYNOPSIS

perl readEntropy.pl -i sequencefile --fmt [sam|fasta|fastq] -k 8 --out output.tsv

=head1 DESCRIPTION

Calculate DNA sequence entropy and redundancy for kmers each sequence in a file
(SAM, Fasta, or Fastq formats). Entropy and redundancy are reported in bits.

=head2 Entropy

=over 4

p_i = (occurrences kmer i) / (total kmer count)

H = - \sum_i^{N} p_i log_2 p_i

for all N kmers

=back

=head2 Redundancy

=over 4

A = (number of possible kmers for a sequence of length L)

R = 1 - \log_2 (H / A)

=back

=cut

use strict;
use warnings;
use diagnostics;
use 5.010;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my %hash;

my $infile;
my $infmt = 'sam';
my $ignoreunmapped = 1;
my $out = 'test.entropy';
my $calc_redundancy = 1;
my $kmer_length = 5;
my $reportsequence = 0;

pod2usage(-verbose=>1) if !@ARGV;

GetOptions ("in|i=s" => \$infile,
            "fmt|f=s" => \$infmt,
            "out|o=s" => \$out,
            "k=i" => \$kmer_length,
            "redundancy!" => \$calc_redundancy,
            "ignoreunmapped!" => \$ignoreunmapped,
            "reportsequence!" => \$reportsequence,
            "help|h" => sub { pod2usage(-verbose=>1); },
            "man|m" => sub { pod2usage(-verbose=>2); },
            ) or pod2usage(-verbose=>1);

=head1 ARGUMENTS

=over 8

=item --in|-i F<FILE>

Input file with sequences

=item --fmt|-f I<STRING>

Format of input file. Either 'sam', 'fasta', or 'fastq'

Default: 'sam'

=item --out|-o I<STRING>

Name for output file

Default: 'test.entropy'

=item --k I<INT>

Length of kmer to use to calculate entropy

Default: 5

=item --ignoreunmapped

For SAM files - ignore alignment records with bit FLAG 0x4 - 'segment unmapped'. 
This is in case there are read pairs where one segment maps but the other does
not.

Turn off with I<--noignoreunmapped>

Default: Yes

=item --help|-h

Help message

=item --man|-h

Full manual page

=back

=head1 OUTPUT

Output is a tab-separated table with three columns: sequence ID, entropy (in
bits), and redundancy (in bits). Header line is prefixed wtih "#" character.

=cut 

# Read in the SAM files
if ($infmt eq 'sam') {
    hash_sam($infile,\%hash);
} elsif ($infmt eq 'fastq') {
    hash_fastq ($infile, \%hash);
} elsif ($infmt eq 'fasta') {
    hash_fasta ($infile, \%hash);
} else {
    say STDERR "Please specify valid input file format: sam, fasta, or fastq";
    pod2usage(-verbose=>1);
    exit;
}

my @allid = keys %hash;
my $entropy_aref = report_entropy (\@allid, \%hash, $kmer_length);
arr2file($entropy_aref, $out);

## SUBS ########################################################################

sub arr2file {
    my ($aref,
        $file,
       ) = @_;
    open(my $fh, ">", $file) or die ("$!");
    foreach my $line (@$aref) {
        print $fh $line;
        print $fh "\n";
    }
    close ($fh);
}

sub report_entropy {
    my ($ids_aref,
        $sam_href,
        $k
        ) = @_;
    my @outarr;
    my @outheader = ('#ID', 'entropy');
    push @outheader, 'redundancy' if $calc_redundancy == 1;
    push @outheader, 'sequence' if $reportsequence == 1;
    push @outarr, join("\t", @outheader); # Header line
    foreach my $id (@$ids_aref) {
        my ($entropy,$redundancy) = calc_entropy($sam_href->{$id}{'seq'},
                                                 $k,
                                                 $calc_redundancy);
        my @outfields = ($id, $entropy);
        push @outfields, $redundancy if defined $redundancy;
        push @outfields, $sam_href->{$id}{'seq'} if $reportsequence == 1;
        push @outarr, join("\t", @outfields);
        #push @outarr, $redundancy;
        #push @outarr, $entropy;
    }
    return \@outarr;
}

sub mean_arr {
    my $aref = shift;
    my $total = sum_arr ($aref);
    my $length = scalar @$aref;
    return ($total / $length);
}

sub report_sam {
    my ($aref,
        $href,
        $outfile) = @_;
    open(my $fh_out, ">", $outfile) or die ("$!");
    foreach my $id (sort @$aref) {
        foreach my $line (@{$href->{$id}{'sam'}}) {
            say $fh_out $line;
        }
    }
    close($fh_out);
}

sub calc_entropy {
    my ($seq, # DNA sequence
        $k,   # Kmer length
        $report_redundancy, # Flag to report redundancy instead of entropy
       ) = @_;
    
    my %hash;
    for (my $i=0;$i<=length($seq)-$k;$i++) {
        # Sliding window of width k across the sequence, count kmer occurrences
        my $word = substr($seq,$i,$k);
        $hash{$word} ++ unless $word =~ m/N/; # Ignore sequences with Ns
    }
    my @vals = values (%hash);
    my $total = sum_arr(\@vals);
    my @p = map {$_ / $total} @vals;
    my @h = map { 0 - $_ * log2($_) } @p;
    
    my $H =sum_arr (\@h);   # Entropy
    my $alphabet = 4**$k;   # Number possible kmers of length k for DNA seq (4 bases)
    my $max_states;
    if ($alphabet < length($seq)-$k + 1) {
        $max_states = $alphabet;       # Max states is limited by length of kmer
    } else {
        $max_states = length($seq)-$k + 1; # Max states is limited by length of sequence
    }
    my $redundancy = 1 - ($H / log2($max_states)) if $report_redundancy == 1;    # Redundancy
    #say "Entropy $H - kmer length $k - alphabet size $alphabet - redundancy $redundancy";
    
    if ($report_redundancy == 1) {
        return ($H, $redundancy);
    } else {
        return ($H);
    }
}

sub sum_arr {
    my $aref = shift;
    my $total = 0;
    foreach my $val (@$aref) {
        $total += $val;
    }
    return ($total);
}

sub log2 {
    my $n = shift;
    return log ($n) / log(2);
}

sub hash_sam {
    my ($file,
        $href,
        ) = @_;
    open(my $fh_in, "<", $file) or die ("$!");
    while (my $line = <$fh_in>) {
        chomp $line;
        next if $line =~ m/^@/; # skip header
        my ($id, @discard) = split /\s/, $line;
        my @splitsam = split /\t/, $line;
        if ($splitsam[1] & 0x40) {
            $id .= '_1';
        } elsif ($splitsam[1] & 0x80) {
            $id .= '_2';
        }
        unless ($ignoreunmapped == 1 && $splitsam[1] & 0x4) {
            $href->{$id}{'seq'} = $splitsam[9] unless defined $href->{$id}{'seq'};
            push @{$href->{$id}{'sam'}}, $line;
        }
    }
    close($fh_in);
}

sub hash_fasta {
    my ($file,
        $href,
        ) = @_;
    my $curr_seq;
    my $curr_id;
    open(my $fh_in, "<", $file) or die ("$!");
    while (my $line = <$fh_in>) {
        chomp $line;
        if ($line =~ m/^>(.+)/) {
            my $id = $1;
            # Header line
            if (defined $curr_id) {
                # Not the first entry in file
                $href->{$curr_id}{'seq'} = $curr_seq;
                $curr_seq = undef;
                $curr_id = $id;
            } else {
                # First entry in file
                $curr_id = $id;
            }
        } else {
            # Sequence line
            $curr_seq .= $line;
        }
    }
    # Scoop up the last entry
    $href->{$curr_id}{'seq'} = $curr_seq;
    close ($fh_in);
}

sub hash_fastq {
    my ($file,
        $href,
       ) = @_;
    my $curr_seq;
    my $curr_id;
    my $curr_qual;
    my $counter = 0;
    open(my $fh, "<", $file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        if ($counter % 4 == 0 && $line =~ m/^@(.+)/) {
            my $id = $1;
            # Header line
            if (defined $curr_id) {
                # Not the first header in file
                $href->{$curr_id}{'seq'} = $curr_seq;
                $href->{$curr_id}{'qual'} = $curr_qual;
                $curr_seq = undef;
                $curr_qual = undef;
                $curr_id = $id;
            } else {
                # First line in file
                $curr_id = $id;
            }
        } elsif ($counter % 4 == 1 ) {
            # Sequence line
            $curr_seq = $line;
        } elsif ($counter % 4 == 3) {
            # Quality line
            $curr_qual = $line;
        }
        # Update counter
        $counter ++;
    }
    # Scoop up last entry
    $href->{$curr_id}{'seq'} = $curr_seq;
    $href->{$curr_id}{'qual'} = $curr_qual;
    close ($fh);
    
}