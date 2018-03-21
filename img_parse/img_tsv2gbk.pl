#!/usr/bin/env perl

=head1 NAME

img_tsv2gbk.pl - Convert IMG annotation table to Genbank format

=head1 SYNOPSIS

perl img_tsv2gbk.pl --input [input_file_table]

perl img_tsv2gbk.pl --help # Short help

perl img_tsv2gbk.pl --man # Full manual

=head1 DESCRIPTION

The JGI IMG annotation pipeline produces a "download bundle" with the raw
annotation results spread across several files. The final product name
assignment is in an annotation table (tab-separated format), but other database
cross-references are also useful for downstream analyses.

However, they do not provide a Genbank-formatted feature table "flatfile", which
is often required for other analysis software, e.g. Pathway Tools. The Genbank-
formatted file in the IMG download bundle contains only feature locations but
not product names or functional annotations.

This script combines the annotations in the IMG annotation table, the Genbank
file from the download bundle, and the table of Interpro mappings, in order to
produce a Genbank file that can be used for downstream analysis that require
Genbank format. This file is not necessarily compliant with the full Genbank
standard, and will not be suitable for submission to INSDC (for sequence
submission, see script "tidy_img_gff3.pl").

=cut

use strict;
use warnings;

# Parse IMG spreadsheet export (TSV) and combine with GBK file
# to have single GBK file with structural and functional annotations
# (e.g. for import to Pathway Tools)

use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);

# Input files
my $input_list; # Table of paths to input files for each genome
my $add_sym_file = "$Bin/product_synonyms.tsv"; # TSV file with additional product name synonyms

# Hashes
my %input_files_hash;
my %func_annot_hash; # Keyed by locus tag
my %oid_locus_hash; # Key - OID, value - locus tag
my %synonym_hash;

pod2usage(verbose=>0,exit_status=>1) if !@ARGV;

GetOptions ("input=s"=>\$input_list,
            "synonym=s"=>\$add_sym_file,
            "help|h"=>sub{pod2usage(verbose=>1,exit_status=>1);},
            "man|m"=>sub{pod2usage(verbose=>2,exit_status=>1);},
            ) or pod2usage(verbose=>1,exit_status=>1);

=head1 ARGUMENTS

=over 8

=item --input <file>

Tab-separated table of IMG annotation files from the download bundle to be
combined into new Genbank files, each line represents annotation files for one
genome.

Columns:
(1) Genome name abbreviation (used to name the output files), (2) IMG annotation
table in tab-separated format (named IMG....info.xls), (3) IMG feature table in
Genbank format, without product names, etc. (named IMG....assembled.gbk), (4)
IMG table of mappings to Interpro, to get cross-references to GO numbers (named
IMG...ipr.tab.txt).

Lines beginning with hash (#) are interpreted as comments and ignored.

=item --synonym <file>

Tab-separated table of additional product name synonyms to be changed in the
output Genbank file. 

=item --help|-h

Concise help message.

=item --man|-m

Full manual page.

=back

=head1 OUTPUT

Output files are named [genome_name]_augmented.gbk, where [genome name] is from
the input table, column 1.

=cut

## MAIN #######################################################################

read_input_file_list($input_list,\%input_files_hash);

foreach my $genome (keys %input_files_hash) {
    my $tsv_file = $input_files_hash{$genome}{"tsv"} if defined $input_files_hash{$genome}{"tsv"};
    my $gbk_file = $input_files_hash{$genome}{"gbk"} if defined $input_files_hash{$genome}{"gbk"};
    my $ipr_file = $input_files_hash{$genome}{"ipr"} if defined $input_files_hash{$genome}{"ipr"};
    my $out_file = $genome."_augmented.gbk";
    hash_synonyms($add_sym_file,\%synonym_hash);
    open_hash_tsv($tsv_file, \%oid_locus_hash, \%func_annot_hash, \%synonym_hash);
    read_hash_ipr($ipr_file, \%oid_locus_hash, \%func_annot_hash) unless !defined $ipr_file;
    augment_GBK($gbk_file, \%func_annot_hash, $out_file);
}


## SUBROUTINES ################################################################

sub hash_synonyms {
    my ($file, $hashref) = @_;
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        ${hashref}->{$splitline[0]} = $splitline[1];
    }
    close(IN);
}

sub open_hash_tsv {
    my ($file, $oidlocushref, $hashref, $synhref) = @_;
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        if (defined $splitline[0] && $splitline[1]) { # Skip blank lines
            my $curr_geneoid = $splitline[0];
            my $curr_locus = $splitline[1];
            ${oidlocushref}->{$curr_geneoid} = $curr_locus;
            ${hashref}->{$curr_locus}{"gene_oid"} = $curr_geneoid;
            # Check for EC numbers
            if (defined $splitline[2] && $splitline[2] =~ m/EC:(\S+\.\S+\.\S+\.\S+)/) {
                my $curr_EC = $1;
                push @{${hashref}->{$curr_locus}{"EC_number"}}, $curr_EC;
            } elsif (defined $splitline[2] && $splitline[2] =~ m/Product_name/) {
                my $curr_product = $splitline[4];
                $curr_product =~ s/^\s+//; # Strip leading whitespace
                if (defined ${synhref}->{$curr_product}) { # Check for synonyms
                    ${hashref}->{$curr_locus}{"product"} = ${synhref}->{$curr_product};
                } else {
                    ${hashref}->{$curr_locus}{"product"} = $curr_product;
                }
            } elsif (defined $splitline[2] && $splitline[2] =~ m/(KO:\S+)/) { # KO numbers
                push @{${hashref}->{$curr_locus}{"db_xref"}}, $1;
            }
        }
    }
    close(IN);
}

sub read_hash_ipr {
    my ($file, $oidlocushref, $hashref) = @_;
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t", $_;
        my $curr_geneoid = $splitline[0];
        my $curr_locus;
        if (defined ${oidlocushref}->{$curr_geneoid}) {
            $curr_locus = ${oidlocushref}->{$curr_geneoid};
        }
        if (defined $splitline[8]) {
            # Get list of GO terms separated by | character
            my @split_go = split /\|/, $splitline[8];
            foreach my $go (@split_go) {
                push @{${hashref}->{$curr_locus}{"db_xref"}}, $go unless !defined $curr_locus;
            }
        }
    }
    close(IN);
}

sub augment_GBK {
    my ($file, $hashref, $outfile) = @_;
    open(my $outfh, ">", $outfile) or die ("Cannot open file $file: $!");
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        print $outfh $_;
        print $outfh "\n";
        if (m/^\s+\/locus_tag="(\S+)"/) {
            my $curr_locus = $1;
            print_extra_fields($curr_locus,$hashref, $outfh);
        }
    }
    close(IN);
    close ($outfh);
}

sub print_extra_fields {
    my ($locus,$hashref, $outfh) = @_;
    my $spaces = "                     ";
    # Product name
    print $outfh genbank_param_line ("product",${hashref}->{$locus}{"product"}) if defined ${hashref}->{$locus}{"product"};
    # Use locus name as gene name to make Pathologic happy **HACKY**
    print $outfh genbank_param_line("gene",$locus) if defined ${hashref}->{$locus}{"product"};
    if (defined ${hashref}->{$locus}{"EC_number"}) { # EC_number fields
        foreach my $EC (@{${hashref}->{$locus}{"EC_number"}}) {
            print $outfh genbank_param_line("EC_number",$EC);
        }
    }
    if (defined ${hashref}->{$locus}{"db_xref"}) { # db_xref fields
        foreach my $xref (@{${hashref}->{$locus}{"db_xref"}}) {
            print $outfh genbank_param_line("db_xref",$xref);
        }
    }
}

sub genbank_param_line {
    # Format attribute line for Genbank file (including newline)
    my ($param, $val) = @_;
    my $spaces = "                     ";
    return "$spaces/$param=\"$val\"\n";
}

sub check_hash { # For testing
    foreach my $locus (keys %func_annot_hash) {
        print $locus."\t";
        print $func_annot_hash{$locus}{"product"} if defined print $func_annot_hash{$locus}{"product"};
        print "\t";
        print join ",", @{$func_annot_hash{$locus}{"EC_number"}} if defined $func_annot_hash{$locus}{"EC_number"};
        print "\n";
    }
}

sub read_input_file_list {
    # Read TSV table of input files
    my ($file, $hashref) = @_;
    open(IN, "<", $file) or die ("Cannot open file $file: $!");
    while (<IN>) {
        chomp;
        unless (/^#/) { # Ignore comment lines
            my @splitline = split "\t";
            ${hashref}->{$splitline[0]}{"tsv"} = $splitline[1] if defined $splitline[1];
            ${hashref}->{$splitline[0]}{"gbk"} = $splitline[2] if defined $splitline[2];
            ${hashref}->{$splitline[0]}{"ipr"} = $splitline[3] if defined $splitline[3];
        }
    }
    close(IN);
}