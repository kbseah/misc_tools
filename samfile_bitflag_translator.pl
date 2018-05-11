#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

if (!@ARGV) {
    say STDERR "Usage: cat file.sam | perl samfile_bitflag_translator.pl > output";
    say STDERR "(Ctrl-c) to stop";
}

while (my $samline = <>) {
    my ($read, $bitflag, @discard) = split /\t/, $samline;

    # Using definitions of bitflags from SAM v1 specification 2017-05-10
    # https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    say "Read $read \t Bitflag $bitflag";
    if ($bitflag & 0x1) {
        say "Bit 0x1\ttemplate having multiple segments in sequencing";
    }
    if ($bitflag & 0x2) {
        say "Bit 0x2\teach segment properly aligned according to the aligner";
    }
    if ($bitflag & 0x4) {
        say "Bit 0x4\tsegment unmapped";
    }
    if ($bitflag & 0x8) {
        say "Bit 0x8\tnext segment in the template unmapped";
    }
    if ($bitflag & 0x10) {
        say "Bit 0x10\tSEQ being reverse complemented";
    }
    if ($bitflag & 0x20) {
        say "Bit 0x20\tSEQ of the next segment in the template being reverse complemented";
    }
    if ($bitflag & 0x40) {
        say "Bit 0x40\tthe first segment in the template";
    }
    if ($bitflag & 0x80) {
        say "Bit 0x80\tthe last segment in the template";
    }
    if ($bitflag & 0x100) {
        say "Bit 0x100\tsecondary alignment";
    }
    if ($bitflag & 0x200) {
        say "Bit 0x200\tnot passing filters, such as platform/vendor quality controls";
    }
    if ($bitflag & 0x400) {
        say "Bit 0x400\tPCR or optical duplicate";
    }
    if ($bitflag & 0x800) {
        say "Bit 0x800\tsupplementary alignment";
    }
    say "";
}
