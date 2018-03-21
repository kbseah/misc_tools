#!/usr/bin/env perl

=head1 NAME

kegg_module_reconstructions_tabulate.pl - Tabulate KEGG module reconstructions

=head1 SYNOPSIS

perl kegg_module_reconstructions_tabulate.pl --table [table of KEGG module output files]

perl kegg_module_reconstructions_tabulate.pl --help

perl kegg_module_reconstructions_tabulate.pl --man

=head1 DESCRIPTION

The KEGG Mapper - Reconstruct Module tool (http://www.kegg.jp/kegg/tool/map_module.html)
predicts the completeness of KEGG metabolic modules, given a list of KEGG
orthology (KO) numbers for a given organism. The tool is available online and
reports the presence/completness of each module in a HTML-formatted table.
However this output is not easy to parse directly, e.g. for synoptic comparison
of several genomes in a single table.

This script takes the output from the Reconstruct Module tool, for several
genomes, and reports a table of modules vs. genomes and their respective
completness. The Reconstruct Module tool will have to be run manually for each
genome, using a list of KO numbers as input. When the results are displayed in
the web browser, uncollapse all the hidden text by clicking on "Show all objects"
and then copy-paste the page contents to a text file, one file per genome.

These results are then tabulated with this script into a single table for
downstream analysis.

=cut

use strict;
use warnings;
#use diagnostics;

use Getopt::Long;
use Pod::Usage;

my $intab;
my $outfile = "kegg_module_reconstruction_tabulated.tsv";
my $keggmodules = "ko00002.keg"; # KEGG modules brite hierarchy
my %genome_file_hash;
my @genome_arr;
my %module_name_hash;
my @module_sort_arr;
my %module_genome_completeness_hash;

pod2usage(verbose=>0,exit_status=>1) if !@ARGV;

GetOptions ("table=s"=>\$intab,
            "module=s"=>\$keggmodules,
            "out=s"=>\$outfile,
            "help|h"=>sub{pod2usage(verbose=>1,exit_status=>1);},
            "man|m"=>sub{pod2usage(verbose=>2,exit_status=>1) },
            ) or pod2usage(verbose=>1,exit_status=>1) ;

=head1 ARGUMENTS

=over 8

=item --table <file>

Tab-separated table, listing all the genome names (column 1) and the text files
containing the copied results from the KEGG Reconstruct Modules web tool (column
2), one line per genome.

=item --module <file>

Path to the KEGG modules BRITE hierarchy file. Should be named ko000002.keg,
available from: http://www.kegg.jp/kegg-bin/get_htext?ko00002.keg

(Download the version in htext format)

=item --out <filename>

Name for the output file, in tab-separated format.

=item --help|-h

Concise help message.

=item --man|-m

Full manual page.

=back

=head1 OUTPUT

Tab separated table of genomes (columns) vs. KEGG modules (rows), with each
indicated as complete ("C"), missing one block ("M"), missing two blocks ("MM"),
or incomplete ("I"). Blanks represent modules that are completely absent.

=cut

## MAIN #######################################################################

# Hash the genomes
%genome_file_hash = %{hashTSV_KV($intab)};
# Get and sort list of genome names
@genome_arr = sort {$a cmp $b} keys %genome_file_hash;
# Hash module names
my ($ref1, $ref2) = hash_module_names($keggmodules);
%module_name_hash = %$ref1;
@module_sort_arr = @$ref2;
# For each genome, note completeness status of each module
foreach my $genome (@genome_arr) {
    hash_module_genome_counts($genome,$genome_file_hash{$genome},\%module_genome_completeness_hash);
}

# Print results
open(OUT, ">", $outfile) or die ("Cannot open $outfile: $!");
my @outheader = ("Module","Name",@genome_arr);
print OUT join "\t", @outheader;
print OUT "\n";
foreach my $module (@module_sort_arr) {
    my $modulename = $module_name_hash{$module};
    my @outline = ($module,$modulename);
    my $count_hits = 0;
    foreach my $genome (@genome_arr) {
        if (defined $module_genome_completeness_hash{$module}{$genome}) {
            push @outline, $module_genome_completeness_hash{$module}{$genome};
            $count_hits++;
        } else {
            push @outline, "";
        }
    }
    if ($count_hits > 0) { # Print output for this module only if module is not completely missing
        print OUT join "\t", @outline;
        print OUT "\n";
    }
}
close(OUT);


## SUBROUTINES ################################################################

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

sub hash_module_names {
    my ($file) = @_;
    my %hash;
    my @arr;
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        if (m/^\w\s+(M\d{5})  (.*)$/) {
            $hash{$1} = $2;
            push @arr, $1;
        }
    }
    close(IN);
    return \%hash, \@arr;
}

sub hash_module_genome_counts {
    my ($genome, $file, $hashref) = @_;
    my %statushash = ("incomplete" => "I",
                      "complete" => "C",
                      "1 block missing" => "M",
                      "2 blocks missing" => "MM");
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        if (m/^\s+(M\d{5}).*\((.+)\)$/) {
            my ($module, $status) = ($1, $2);
            ${hashref}->{$module}{$genome} = $statushash{$status};
            #print "$genome\t$module\t".$statushash{$status}."\n";
        }
    }
    close(IN);
}