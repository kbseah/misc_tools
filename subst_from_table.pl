#!/usr/bin/perl

# Given a tab-separated table of string pairs, substitute one for the other

use strict;
use warnings;
use Getopt::Long;

if (! @ARGV) {
    print "\n";
    print STDERR "\tperl subst_from_table.pl --input file --table substitution_table [--flip] \n\n";
    exit;
}

my $file_in;
my $table_in;
my $flip=0; # If activated, column 2 substituted by column 1 instead of default 
my %subst_hash;

GetOptions('input|i=s'=>\$file_in,
           'table|t=s'=>\$table_in,
           'flip|f'=>\$flip) or die ("$!\n");

open(TAB, "<", $table_in) or die ("$!\n");
while (<TAB>) {
    chomp;
    my @theline=split"\t";
    if ($flip == 0) { $subst_hash{$theline[0]} = $theline[1]; }
    elsif ($flip == 1) { $subst_hash{$theline[1]} = $theline[0]; }
}
close(TAB);

open(FILE, "<", $file_in) or die ("$!\n");
while (<FILE>) {
    chomp;
    my $theline = $_;
    foreach my $thekey (keys %subst_hash) {
        $theline =~ s/$thekey/$subst_hash{$thekey}/;
    }
    print STDOUT $theline."\n";
}
close(FILE);
