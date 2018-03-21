#!/usr/bin/env perl
# Tabulate KO numbers from IMG xls file

=head1 NAME

img_xls_extract_KO.pl - Extract KO numbers vs gene ids from IMG annotation table

=head1 SYNOPSIS

perl img_xls_extract_KO.pl --input <table>

perl img_xls_extract_KO.pl --help

perl img_xls_extract_KO.pl --man

=head1 DESCRIPTION

Extracts a list of KO numbers and corresponding gene IDs from an IMG annotation
table in tab-separated format (named IMG_....xls in download bundle).

The output are tab-separated tables, one per genome, which can be fed to the
KEGG Mapper Reconstruct Modules tool.

=cut

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

my $img_tab_file; # TSV table of genome shortnames and paths to annotation files
my %names_file_hash; # Hash for genome shortnames and paths to annotation files

pod2usage(verbose=>0,exit_status=>1) if !@ARGV;

GetOptions("input=s"=>\$img_tab_file,
           "help|h"=> sub{pod2usage(verbose=>1,exit_status=>1);},
           "man|m"=> sub{pod2usage(verbose=>2,exit_status=>1);},
           ) or pod2usage(verbose=>1,exit_status=>1);

=head1 ARGUMENTS

=over 8

=item --input <file>

Tab-separated table of genome names (column 1) and paths to IMG annotation table
files (column 2) for each genome. The genome names will be used to name the
output files.

=item --help|-h

Concise help message

=item --man|-m

Full manual page

=back

=head1 OUTPUT

Tab-separated tables, named KO_list_[genome name].tsv, where [genome name] is
read from column 1 of the input table. One TSV file per genome annotation.

The tables can be uploaded to the KEGG Mapper tool for module prediction.

=cut


%names_file_hash = %{hashTSV_KV($img_tab_file)};

foreach my $genome (sort {$a cmp $b} keys %names_file_hash) {
    my $outfile = "KO_list_$genome.tsv";
    my $infile = $names_file_hash{$genome};
    make_outfile ($infile, $outfile);
}

sub make_outfile {
    my ($infile, $outfile) = @_;
    open(IN, "<", $infile) or die ("$!");
    open(OUT, ">", $outfile) or die ("$!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t", $_;
        if (defined $splitline[2] && $splitline[2] =~ m/KO:(K\d{5})/) {
            print OUT "$splitline[1]\t$1\n";
        }
    }
    close (IN);
    close (OUT);
}

sub hashTSV_KV {
# Open TSV file, split cols by tab, and hash with col1 as key and col2 as val
    my ($file) = @_;
    my %hash;
    open(IN, "<", $file) or die ("File $file not found: $!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        $hash{$splitline[0]} = $splitline[1];
    }
    close(IN);
    return (\%hash);
}