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

pod2usage(-verbose=>1) if !@ARGV;

GetOptions ("1=s" => \$sam1,
            "2=s" => \$sam2,
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
hash_sam($sam1, \%hash1);
hash_sam($sam2, \%hash2);
# Find sequences extracted by one pipeline but not the other
my @in1not2;
my @in2not1;
foreach my $id (keys %hash1) {
    if (!defined $hash2{$id}) {
        push @in1not2, $id;
    }
}
foreach my $id (keys %hash2) {
    if (!defined $hash1{$id}) {
        push @in2not1, $id;
    }
}
say "In file 1 but not file 2: ".scalar (@in1not2)." read pairs";
say "In file 2 but not file 1: ".scalar (@in2not1)." read pairs";

report_sam (\@in1not2,\%hash1,"$out.in1not2.sam");
report_sam (\@in2not1,\%hash2,"$out.in2not1.sam");


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
    open(my $fh_in, "<", $file) or die ("$!");
    while (my $line = <$fh_in>) {
        chomp $line;
        next if $line =~ m/^@/; # skip header
        my ($id, @discard) = split /\s/, $line;
        my @splitsam = split /\t/, $line;
        $href->{$id}{'seq'} = $splitsam[9] unless defined $href->{$id}{'seq'};
        push @{$href->{$id}{'sam'}}, $line;
    }
    close($fh_in);
}

