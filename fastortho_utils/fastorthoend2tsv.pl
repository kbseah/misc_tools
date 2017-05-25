#!/usr/bin/env perl

# Reformat .end file from FastOrtho into tabular TSV format

use strict;
use warnings;

use Getopt::Long;

my $fastortho_result;
my $outfile = "test.out";
my $coreflag = 0;
my $num_genomes;
my $singlecopyflag = 0;
my %id_lib_ortho_hash;

GetOptions("fastortho_end|e=s"=>\$fastortho_result, # .end file output from FastOrtho
           "out|o=s"=>\$outfile, # Name for output TSV file
           "core"=>\$coreflag, # Flag: Output only core genome members?
           "num_genomes=i"=>\$num_genomes, # Number of genomes (needed if core flag is on)
           "singlecopy"=>\$singlecopyflag, # Flag: Output only single copy core genome?
           ) or die ("Invalid options: $!");


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
write_tsv();


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

sub write_tsv {
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