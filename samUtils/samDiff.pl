#!/usr/bin/env perl

=head1 NAME

samDiff.pl - Compare two SAM files and report subset of reads that are in one but not the other

=head1 SYNOPSIS

perl samDiff.pl --1 file1.sam --2 file2.sam --out prefix

=cut

use strict;
use warnings;
use 5.010;

#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my ($sam1, $sam2);
my $out = "test";
my $trd;
my $writesam = 1;
my $verbose = 0;
my $header = 1;

pod2usage(-verbose=>1) if !@ARGV;

GetOptions ("1=s" => \$sam1,
            "2=s" => \$sam2,
            "trimreaddescription|trd" => \$trd,
            "writesam!" => \$writesam,
            "verbose!" => \$verbose,
            "header!" => \$header,
            "out|o=s" => \$out,
            "help|h" => sub {pod2usage(-verbose=>1); },
            "man|m" => sub {pod2usage (-verbose=>2); },
            );

=head1 ARGUMENTS

=over 8

=item --1 F<FILE>

First SAM file to compare

=item --2 F<FILE>

Second SAM file to compare

=item --out I<STRING>

Prefix for output file names

Default: 'test'

=item --writesam

Write SAM output files?

Default: Yes (turn off with --nowritesam)

=item --trimreaddescription

Truncate QNAME and RNAME fields on first whitespace.

=item --verbose

Verbosely report statistics? Otherwise report stats as a tab-separated list

Default: No (--noverbose)

=back

=head1 OUTPUT

=over 8

=item [OUT].in1not2.sam

SAM fields for reads present in SAM file 1 but not file 2

=item [OUT].in2not1.sam

SAM fields for reads present in SAM file 2 but not file 1

=back

=cut

my %hash1;
my %hash2;

# Read in the SAM files
my ($sam1_totalaln, $sam1_totalreadpairs) = hash_sam($sam1, \%hash1);
my ($sam2_totalaln, $sam2_totalreadpairs) = hash_sam($sam2, \%hash2);
# Find sequences extracted by one pipeline but not the other
my @in1not2;
my @in2not1;
my @inboth;
foreach my $id (keys %hash1) {
    if (!defined $hash2{$id}) {
        push @in1not2, $id;
    } else {
        push @inboth, $id;
    }
}
foreach my $id (keys %hash2) {
    if (!defined $hash1{$id}) {
        push @in2not1, $id;
    }
}

my $in1not2_count = scalar (@in1not2);
my $in2not1_count = scalar (@in2not1);
my $inboth_count = scalar (@inboth);

if ($verbose == 1) {
    say "File1 - Detected $sam1_totalaln alignments";
    say "        representing $sam1_totalreadpairs read pairs";
    say "File2 - Detected $sam2_totalaln alignments";
    say "        representing $sam2_totalreadpairs read pairs";
    say "In file 1 but not file 2: $in1not2_count read pairs";
    say "In file 2 but not file 1: $in2not1_count read pairs";
    say "In both SAM files: $inboth_count read pairs";
} else {
    my @outarr = ($sam1,
                  $sam1_totalaln,
                  $sam1_totalreadpairs,
                  $sam2,
                  $sam2_totalaln,
                  $sam2_totalreadpairs,
                  $in1not2_count,
                  $in2not1_count,
                  $inboth_count);
    say join "\t", qw(file1 file1_totalaln file1_totalPE file2 file2_totalaln file2_totalPE infile1notfile2 infile2notfile1 inboth) if $header;
    say join "\t", @outarr;
}

if ($writesam == 1) {
    report_sam (\@in1not2,\%hash1,"$out.in1not2.sam");
    report_sam (\@in2not1,\%hash2,"$out.in2not1.sam");
    report_sam (\@inboth, \%hash1,"$out.inboth.from1.sam");
    report_sam (\@inboth, \%hash2,"$out.inboth.from2.sam");
}

## SUBS ########################################################################


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

sub hash_sam {
    my ($file,
        $href,
        ) = @_;
    my $counter = 0;
    open(my $fh_in, "<", $file) or die ("$!");
    while (my $line = <$fh_in>) {
        chomp $line;
        next if $line =~ m/^@/; # skip header
        $counter ++;
        my $id;
        my @splitsam = split /\t/, $line;
        if ($trd) {
            my @discard;
            ($id,@discard) = split /\s+/, $splitsam[0]; # QNAME
            ($splitsam[2],@discard) =  split /\s+/, $splitsam[2]; # RNAME
        } else {
            $id = $splitsam[0];
        }
        $href->{$id}{'seq'} = $splitsam[9] unless defined $href->{$id}{'seq'};
        push @{$href->{$id}{'sam'}}, $line;
    }
    close($fh_in);
    my $num_readpairs = scalar(keys %$href);
    return ($counter,       # Total number of alignments in SAM file
            $num_readpairs  # Total number of reads (or read pairs) in SAM file
            );
    #say "Total alignments in SAM file $file: $counter";
    #say "Total read pairs aligned in SAM file $file: ".scalar(keys %$href);
}

