#!/usr/bin/env perl

# Reformat .end file from FastOrtho into tabular TSV format

=head1 NAME

fastorthoend2tsv.pl - Convert FastOrtho .end output to tabular TSV format

=head1 SYNOPSIS

perl fastorthoend2tsv.pl -e input.end -o output.tsv -f list

perl fastorthoend2tsv.pl --help

perl fastorthoend2tsv.pl --man

=head1 DESCRIPTION

FastOrtho is a reimplementation of the OrthoMCL pipeline for finding gene
ortholog clusters from a set of genomes, given their all-vs-all pairwise Blast
results. The output from FastOrtho is an I<.end> file, which contains a list of
ortholog clusters and the genes (and parent genomes) in each cluster. This
script converts the I<.end> file into a TSV formatted table.

=cut

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

my $fastortho_result;
my $outfile = "test.out";
my $coreflag = 0;
my $num_genomes;
my $format = "list";
my $singlecopyflag = 0;
my %id_lib_ortho_hash;

pod2usage(-verbose=>0) if !@ARGV;

GetOptions("fastortho_end|e=s" => \$fastortho_result,   # .end file output from FastOrtho
           "out|o=s" => \$outfile,                      # Name for output TSV file
           "format|f=s" => \$format,                    # Input format
           "core" => \$coreflag,                        # Flag: Output only core genome members?
           "num_genomes=i" => \$num_genomes,            # Number of genomes (needed if core flag is on)
           "singlecopy" => \$singlecopyflag,            # Flag: Output only single copy core genome?
           "help|h" => sub{ pod2usage(-verbose=>1); },
           "man|m" => sub{ pod2usage(-verbose=>2); }
           ) or pod2usage(-verbose=>1);

=head1 ARGUMENTS

=over 8

=item --fastortho_end | -e I<FILE>

Path to input file. Must be .end output file from a FastOrtho run.

=item --out | -o I<FILE>

Name for output file.

Default: "test.out"

=item --format | -f I<STRING>

Output format.

"list" will report table with three columns: product (gene name), library
(genome name), and ortholog (FastOrtho-assigned ortholog name)

"matrix" will report a table of library (genome) vs. ortholog, with the counts
of orthologs per genome as the values of the matrix.

Default: "list"

=item --core

Logical: Output only core genome? I.e. genes that occur in all genomes.

Default: No

=item --num_genomes I<INTEGER>

Number of genomes in the input. This is compulsory if the I<--core> option is
supplied.

=item --singlecopy

Logical: Output only the single-copy core genome? I.e. genes that occur in all
genomes once and only once.

Default: No.

=item --help | -h

Short help message

=item --man | -m

Full manual page

=back

=cut


# Check that option flags make sense
if ($singlecopyflag > 0 && $coreflag == 0) { # Find single-copy core implies find core-genome
    $coreflag = 1;
}
if ($coreflag > 0) {
    if (!defined $num_genomes) {
        print STDERR "Please specify number of genomes analyzed with --num_genomes option \n";
        exit;
    } else {
        print STDERR "Output contains only core genes found in all $num_genomes genomes\n";
    }
} 

# Main
read_FastOrtho2();

if ($format eq "list") {
    write_tsv_list();
} elsif ($format eq "matrix") {
    write_tsv_matrix();
}




## SUBROUTINES ################################################################

sub read_FastOrtho2 {
    # Read and parse FastOrtho output for ALL clusters
    open(FASTORTHO, "<", $fastortho_result) or die ("$!\n");
        while (<FASTORTHO>) {
        chomp;
        my @theline = split "\t", $_;
        if ($theline[0] =~ /^(ORTHOMCL\d+) \((\d+) genes,(\d+) taxa\)/) {
            my $current_cluster = $1;
            my ($num_genes,$num_taxa) = ($2, $3);
            my @theproteins = split /\s+/, $theline[1];
            my $flag = 1; # Flag - hash the data or not?
            # Option flags for core genome
            if ($coreflag > 0) {
                if ($num_taxa == $num_genomes) {
                    $flag = 1;
                    if ($singlecopyflag > 0) { # Catch genes that are not single-copy
                        if ($num_genes != $num_genomes) {
                            $flag = 0;
                        } # Else $flag remains 1 
                    }
                } elsif ($num_taxa < $num_genomes) { # Not found in all genomes
                    $flag = 0;
                } elsif ($num_taxa > $num_genomes) { # Number of genomes specified by user does not match FastOrtho file
                    print STDERR "Number of genomes specified is less than maximum found in FastOrtho output file - check your input?\n";
                    exit;
                } else { # If it is some other exception
                    $flag = 0;
                }
            }
            # Parse entry and read into hash
            if ($flag == 1) {
                foreach my $theentry (@theproteins) {
                    if ($theentry =~ /^(\S+)\((\S+)\)/) {
                        my $product = $1;
                        my $lib = $2;
                        $id_lib_ortho_hash{$product}{"lib"} = $lib;
                        $id_lib_ortho_hash{$product}{"ortho"} = $current_cluster;
                    }
                }
            }
        }
    }
    close (FASTORTHO);
}

sub write_tsv_matrix {
    # Refactor the hash
    my %lib_ortho_count_hash;
    my %lib_hash;
    my %ortho_hash;
    
    foreach my $product (keys %id_lib_ortho_hash) {
        my $currlib = $id_lib_ortho_hash{$product}{"lib"};
        my $currortho = $id_lib_ortho_hash{$product}{"ortho"};
        $lib_ortho_count_hash{$currlib}{$currortho} ++;
        $lib_hash{$currlib} ++;
        $ortho_hash{$currortho} ++;
    }
    
    # Write output
    open(OUT, ">", $outfile) or die ("Cannot open output file: $!");
    # header line
    my @outhead = sort {$a cmp $b} keys %ortho_hash;
    unshift @outhead, "library";
    print OUT join "\t", @outhead;
    print OUT "\n";
    foreach my $lib (sort {$a cmp $b} keys %lib_ortho_count_hash) {
        my @out = $lib;
        foreach my $ortho (sort {$a cmp $b} keys %ortho_hash) {
            # If nothing counted then print 0
            my $count = defined $lib_ortho_count_hash{$lib}{$ortho} ? $lib_ortho_count_hash{$lib}{$ortho} : 0;
            push @out, $count;
        }
        print OUT join "\t", @out;
        print OUT "\n";
    }
    close(OUT);
}

sub write_tsv_list {
    open(OUT, ">", $outfile) or die ("Cannot open output file: $!"); 
    # Header line
    print OUT join "\t", qw(product library ortholog);
    print OUT "\n";
    foreach my $prod (sort {$a cmp $b} keys %id_lib_ortho_hash) {
        my @out = ($prod,
                   $id_lib_ortho_hash{$prod}{"lib"},
                   $id_lib_ortho_hash{$prod}{"ortho"});
        print OUT join "\t", @out;
        print OUT "\n";
    }
    close (OUT);
}