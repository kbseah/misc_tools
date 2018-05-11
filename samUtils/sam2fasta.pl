#!/usr/bin/env perl

=head1 NAME

sam2fastx.pl - Convert SAM to Fasta or Fastq file

=head1 SYNOPSIS

perl -i file.sam > out.fasta

perl -i file.sam --outfmt fastq > out.fastq

=head1 DESCRIPTION 

Convert SAM file to Fasta or Fastq file using sequence names found in SAM file.
Does so naively, i.e. does not distinguish paired reads if recorded under the
same name.

Circumvents the requirement of reformat.sh from BBmap suite, which requires
header information in SAM file for conversion.

=cut

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;

my %hash;

my $insam;
my $out;
my $outfmt = 'fasta';
my $limit;

pod2usage(-verbose=>1) if (!@ARGV);

GetOptions ("--in|i=s" => \$insam,
            "--out|o=s" => \$out,
            "--outfmt=s" => \$outfmt,
            "--limit=i" => \$limit,
            "--help|h" => sub { pod2usage (-verbose=>1) },
            "--man|m" => sub { pod2usage (-verbose=>2) },
            ) or pod2uage(-verbose=>1);

=head1 ARGUMENTS

=item --in|-i F<FILE>

SAM file input

=item --out|-o F<FILE>

Output filename

Default: STDOUT

=item --outfmt I<STRING>

Format for output file. Either 'fasta' or 'fastq'

Default: 'fasta'

=item --limit I<INT>

Stop processing after N reads

Default: None

=cut

open(my $fh, "<", $insam) or die ("$!");
my $counter = 0;
while (my $line = <$fh>) {
    next if $line =~ m/^@/; # Skip headers
    if (defined $limit && $counter > $limit) {
        last; # Stop processing after $limit reads, if a limit is defiend
    }
    my @split = split /\t/, $line;
    $hash{$split[0]}{'seq'} = $split[9] unless defined $hash{$split[0]}{'seq'};
    $hash{$split[0]}{'qual'} = $split[10] if $outfmt eq 'fastq';
}
close($fh);


my $fhout;
if (defined $out) {
    open($fhout, ">", $out) or die ("$!");
} else {
    $fhout = *STDOUT; # Write to STDOUT by default
}
if ($outfmt eq 'fastq') {
    write_fastq ($fhout, \%hash);
} else {
    write_fasta ($fhout, \%hash);
}
close ($fhout) if defined $out;

## SUBS ########################################################################

sub usage {
    say STDERR "Usage: perl sam2fasta.pl file.sam > out.fasta";
    exit;
}

sub write_fasta {
    my ($fh,
        $href) = @_;
    foreach my $id (sort keys %$href) {
        say $fh ">$id";
        say $fh $href->{$id}{'seq'};
    }
}

sub write_fastq {
    my ($fh,
        $href) = @_;
    foreach my $id (sort keys %$href) {
        say $fh "@$id";
        say $fh $href->{$id}{'seq'};
        say $fh "+";
        say $fh $href->{$id}{'qual'};
    }
}