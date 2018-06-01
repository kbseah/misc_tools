#!/usr/bin/env perl

=head1 NAME

samBitflagExplain.pl - Explain meaning of bitflags in SAM file

=head1 SYNOPSIS

perl samBitflagExplain.pl --sam file.sam > output

perl samBitflagExplain.pl --help

=cut

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Pod::Usage;

my $samfile;
my $out;
my $limit = 10;

pod2usage(-verbose=>0) if !@ARGV;
GetOptions ("sam=s"=> \$samfile,
            "out|o=s" => \$out,
            "limit=i" => \$limit,
            "help|h" => sub {pod2usage(-verbose=>1); },
            "man|m" => sub{ pod2usage(-verbose=>2); },
            ) or pod2usage(-verbose=>1);

=head1 ARGUMENTS

=over 8

=item --sam F<FILE>

SAM file to parse

=item --out F<FILE>

Path to write output.

Default: STDOUT

=item --limit I<INT>

Stop parsing after N entries. Specify -1 for no limit.

Default: 10

=item --help

Help message

=item --man

Full manual page

=back

=cut

if (!defined $samfile) {
    say STDERR "SAM file not specified";
    pod2usage(-verbose=>1);
    exit;
}


my $fhout;
if (defined $out) {
    open($fhout, ">", $out) or die ("$!");
} else {
    $fhout = *STDOUT;
}

open(my $fhin, "<", $samfile) or die ("$!");
my $counter = 0;
while (my $samline = <$fhin>) {
    next if $samline =~ m/^@/; # Skip headers
    last if $limit > 0 && $counter >= $limit;
    my ($read, $bitflag, @discard) = split /\t/, $samline;
    my $binary = sprintf "%014b", $bitflag;

    # Using definitions of bitflags from SAM v1 specification 2017-05-10
    # https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    say "Read $read\tBitflag $bitflag\tBinary $binary";
    if ($bitflag & 0x1) {
        say $fhout "Bit 0x1\ttemplate having multiple segments in sequencing";
    }
    if ($bitflag & 0x2) {
        say $fhout "Bit 0x2\teach segment properly aligned according to the aligner";
    }
    if ($bitflag & 0x4) {
        say $fhout "Bit 0x4\tsegment unmapped";
    }
    if ($bitflag & 0x8) {
        say $fhout "Bit 0x8\tnext segment in the template unmapped";
    }
    if ($bitflag & 0x10) {
        say $fhout "Bit 0x10\tSEQ being reverse complemented";
    }
    if ($bitflag & 0x20) {
        say $fhout "Bit 0x20\tSEQ of the next segment in the template being reverse complemented";
    }
    if ($bitflag & 0x40) {
        say $fhout "Bit 0x40\tthe first segment in the template";
    }
    if ($bitflag & 0x80) {
        say $fhout "Bit 0x80\tthe last segment in the template";
    }
    if ($bitflag & 0x100) {
        say $fhout "Bit 0x100\tsecondary alignment";
    }
    if ($bitflag & 0x200) {
        say $fhout "Bit 0x200\tnot passing filters, such as platform/vendor quality controls";
    }
    if ($bitflag & 0x400) {
        say $fhout "Bit 0x400\tPCR or optical duplicate";
    }
    if ($bitflag & 0x800) {
        say $fhout "Bit 0x800\tsupplementary alignment";
    }
    say $fhout "";
    $counter ++;
}
close($fhin);
close ($fhout) if defined $out;