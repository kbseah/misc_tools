#!/usr/bin/env perl

=head1 NAME

aa_composition_from_fasta -- Tabulate AA composition from a Fasta file

=head1 SYNOPSIS

perl aa_composition_from_fasta.pl -i input_table

perl aa_composition_from_fasta.pl --help

perl aa_composition_from_fasta.pl --man

=head1 DESCRIPTION

Tabulate amino acid composition per sequence from a set of Fasta files.
Assume that each fasta file represents one gene, and each sequence one taxon.

=cut

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;

# Tabulate amino acid composition per sequence per fasta file
# Assume each fasta file represents one gene and each sequence one taxon

my @aa = qw(R H K D E S T N Q C U G P A V I L M F Y W X); # List of amino acid characters
my @squareplot_chars = qw(AA AG GA GG); # List of squarelpot states
my @dayhoff6_states = qw(1 2 3 4 5 6); # List of dayhoff6 states
my %squareplot = ( # Hash to convert amino acids to Squareplot states
    "K" => "AA",
    "M" => "AA",
    "N" => "AA",
    "I" => "AA",
    "L" => "AA",
    "Y" => "AA",
    "F" => "AA",
    "S" => "AG",
    "T" => "AG",
    "W" => "AG",
    "C" => "AG",
    "E" => "GA",
    "D" => "GA",
    "V" => "GA",
    "Q" => "GA",
    "H" => "GA",
    "L" => "GA",
    "G" => "GG",
    "R" => "GG",
    "A" => "GG",
    "P" => "GG",
);
my %dayhoff6 = ( # Hash to convert amino acids to dayhoff6 states
    "A" => "1",
    "G" => "1",
    "P" => "1",
    "S" => "1",
    "T" => "1",
    "D" => "2",
    "E" => "2",
    "N" => "2",
    "Q" => "2",
    "H" => "3",
    "K" => "3",
    "R" => "3",
    "I" => "4",
    "L" => "4",
    "M" => "4",
    "V" => "4",
    "F" => "5",
    "W" => "5",
    "Y" => "5",
    "C" => "6",
);

my $input_table;
my $recoding = 'aa';

pod2usage (-exitstatus=>2, -verbose=>1) if (!@ARGV); 

GetOptions("input|i=s" => \$input_table,
           "recoding|r=s" => \$recoding,
           "help|h" => sub { pod2usage(-exitstatus=>2, -verbose=>2); },
           "man|m" => sub { pod2usage(-exitstatus=>0, -verbose=>2); },
           ) or pod2usage(-exitstatus=>2, -verbose=>2);

=head1 ARGUMENTS

=over 8

=item --input|-i <file>

Tab-separated table of Fasta files; column 1 is filename, col 2 is the gene
name or abbreviation.

=item --recoding|-r <string>

Recoding to use for amino acid sequences. Possible values 'dayhoff' (Dayhoff
6-state) or 'squareplot' (Foster 1997) (Default: none)

=item --help|-h

Help message

=item --man|-m

This manual page.

=back

=cut


## MAIN ########################################################################

my $href = read_input_table($input_table);
my %fasta_by_gene = %$href;
my %counts;
my $char_aref;
if ($recoding eq 'dayhoff') {
    $char_aref = \@dayhoff6_states;
} elsif ($recoding eq 'squareplot') {
    $char_aref = \@squareplot_chars;
} else {
    $char_aref = \@aa;
}

foreach my $gene (keys %fasta_by_gene) {
    my $fasta_file = $fasta_by_gene{$gene};
    fasta_to_counts($fasta_file, $gene, \%counts);
}

# Print output in "long" format
say join "\t", qw (gene species state count);
foreach my $gene (sort {$a cmp $b} keys %counts) {
    foreach my $species (sort {$a cmp $b} keys %{$counts{$gene}}) {
        foreach my $char (@$char_aref) {
            my $thecount;
            if (defined $counts{$gene}{$species}{$char}) {
                $thecount = $counts{$gene}{$species}{$char};
            } else {
                $thecount = 0;
            }
            say join "\t", ($gene, $species, $char, $thecount);
        }
    }
}

## SUBS ########################################################################

sub read_input_table {
    my ($file) = @_;
    my %hash;
    my $fh;
    open($fh, "<", $file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        my ($fasta_file, $gene) = split /\t/, $line;
        $hash{$gene} = $fasta_file;
    }
    close($fh);
    return (\%hash);
}

sub fasta_to_counts {
    my ($fasta_file, $gene, $href) = @_;
    my $current_species;
    my $current_seq;
    open(my $fh, "<", $fasta_file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^>(.*)/) {
            my $new_species = $1;
            #$href->{$gene}{$current_species}{'string'} = $current_seq;
            count_and_hash_aa_sequence ($current_seq, $current_species, $gene, $href) unless !defined $current_seq;
            $current_species = $new_species;
            $current_seq = "";
        } else {
            $current_seq .= $line;
        }
    }
    # Sweep up last sequence
    #$href->{$gene}{$current_species}{'string'} = $current_seq;
    count_and_hash_aa_sequence ($current_seq, $current_species, $gene, $href);
    close($fh);
}

sub count_and_hash_aa_sequence {
    my ($seq, $species, $gene, $href) = @_;
    chomp $seq;
    my @seqsplit = split //, $seq;
    foreach my $rawchar (@seqsplit) {
        my $char;
        # Recoding to Dayhoff or Squareplot states
        if ($recoding eq 'dayhoff') {
            if (defined $dayhoff6{$rawchar}) {
                $char = $dayhoff6{$rawchar};
            }
        } elsif ($recoding eq 'squareplot') {
            if (defined $squareplot{$rawchar}) {
                $char = $squareplot{$rawchar};
            }
        } else {
            $char = $rawchar;
        }
        # Hash the character
        $href->{$gene}{$species}{$char} ++ if defined $char;
    }
}



=head1 COPYRIGHT AND LICENSE

Copyright 2018, Brandon Seah (kbseah@mpi-bremen.de)

LICENSE
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut