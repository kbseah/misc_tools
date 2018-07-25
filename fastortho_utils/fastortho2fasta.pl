#!/usr/bin/perl

=head1 NAME

fastortho2fasta.pl - Convert FastOrtho orthologs to Fasta-formatted alignment

=head1 SYNOPSIS

perl fastortho2fasta.pl -e fastortho.end -f fastortho.faa -n [num taxa] -c concat.fasta

perl fastortho2fasta.pl --help

perl fastortho2fasta.pl --man

=head1 DESCRIPTION

Parse FastOrtho results file in I<.end> format to get the list of ortholog
clusters found in all genomes only once (i.e. single-copy markers). Align each
of these clusters with Muscle (should be in path), and concatenate the alignment
for phylogenetic analysis.

=cut

use strict;
use warnings;

use Bio::Align::Utilities qw(:all);
use Bio::LocatableSeq;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;

use Getopt::Long;
use Pod::Usage;

my $fastortho_result;
my $num_taxa;
my $fasta_db;
my $concat_file;

pod2usage(-verbose=>0) if !@ARGV;

GetOptions ('fastortho_end|e=s'=>\$fastortho_result,
            'num_taxa|n=i'=>\$num_taxa,
            'fasta|f=s'=>\$fasta_db,
            'concat|c=s'=>\$concat_file,
            'help|h' => sub {pod2usage(-verbose=>1);},
            'man|m' => sub {pod2usage(-verbose=>2);}
            ) or pod2usage(-verbose=>1);

=head1 ARGUMENTS

=over 8

=item --fastortho_end | -e I<FILE>

Output file in I<.end> format from FastOrtho.

=item --num_taxa | -n I<INTEGER>

Number of genomes/taxa in the ortholog clustering analysis. If the number 0 is
specified, then a separate Fasta file and alignment will be produced for each
ortholog cluster. Otherwise, only those taxa appearing once in all genomes will
be extracted and aligned.

=item --fasta | -f I<FILE>

Fasta file of protein sequences used as input for FastOrtho clustering.

=item --concat | -c I<FILE>

If specified, concatenated alignment of the extracted single-copy markers will
be written to this file.

Default: None

=item --help | -h

Short help message

=item --man | -m

Full manual page

=back

=cut

## MAIN ########################################################################

my %cluster_hash;   # Store the protein IDs that correspond to each cluster
my %clusters_with_paralogs; # Clusters that contain members from same genome (i.e. putative paralogs)
my %alignment_hash;
my %genomes_all;   # Store the names of all the genomes, as read from the FastOrtho.end file

my $concat_aln;

if ($num_taxa == 0) {
    read_FastOrtho2();
}
elsif ($num_taxa > 0 ) {
    read_FastOrtho();
}

print STDERR "Clusters with putative paralogs: ".scalar (keys %clusters_with_paralogs)."\n";
#test_hash();
read_write_Fasta();
align_Fasta();
if (defined $concat_file) { concat_alignment(); }

## SUBS ########################################################################

sub read_FastOrtho {        # Read and parse FastOrtho output, given a minimum number of taxa that must occur in the cluster
open(FASTORTHO, "<", $fastortho_result) or die ("$!\n");
while (<FASTORTHO>) {
    chomp;
    my @theline = split "\t", $_;
    if ($theline[0] =~ /^(ORTHOMCL\d+) \($num_taxa genes,$num_taxa taxa\)/) {
        my $current_cluster = $1;
        #print $current_cluster."\n";
        my @theproteins = split /\s+/, $theline[1];
        foreach my $theentry (@theproteins) {
            if ($theentry =~ /^(\S+)\((\S+)\)/) {
                #print $1 ."\t". $2. "\t";
                if (defined ${$cluster_hash{$current_cluster}}{$2}) {   # Catch clusters which contain paralogs - not nice for phylogenetics!
                    print STDERR "$2 appears more than once in $current_cluster\n";
                    $clusters_with_paralogs{$current_cluster}++;
                }
                ${$cluster_hash{$current_cluster}}{$2} = $1;    # Hash-of-hash - store product ID as value of hash with genome name as key
                $genomes_all{$2}=1;     # Record all the genomes that occur in these clusters
            }
        }
        #print "\n";
    }
}
close(FASTORTHO);
}

sub read_FastOrtho2 {       # Read and parse FastOrtho output for ALL clusters
open(FASTORTHO, "<", $fastortho_result) or die ("$!\n");
    while (<FASTORTHO>) {
    chomp;
    my @theline = split "\t", $_;
    if ($theline[0] =~ /^(ORTHOMCL\d+) /) {
        my $current_cluster = $1;
        #print $current_cluster."\n";
        my @theproteins = split /\s+/, $theline[1];
        foreach my $theentry (@theproteins) {
            if ($theentry =~ /^(\S+)\((\S+)\)/) {
                #print $1 ."\t". $2. "\t";
                if (defined ${$cluster_hash{$current_cluster}}{$2}) {   # Catch clusters which contain paralogs - not nice for phylogenetics!
                    print STDERR "$2 appears more than once in $current_cluster\n";
                    $clusters_with_paralogs{$current_cluster}++;
                }
                ${$cluster_hash{$current_cluster}}{$2} = $1;    # Hash-of-hash - store product ID as value of hash with genome name as key
                $genomes_all{$2}=1;         # Record all the genomes that occur in these clusters
            }
        }
        #print "\n";
    }
}
close (FASTORTHO);
}


sub test_hash {
foreach my $thekey (sort keys %cluster_hash) {
    print $thekey."\n";
    for my $thekey2 (sort keys %{$cluster_hash{$thekey}}) {
        print $thekey2."\t".$cluster_hash{$thekey}{$thekey2}."\t";
    }
    print "\n";
}
}

sub read_write_Fasta {      # For each cluster, write a separate Fasta file, with the genome names as headers 
    my $db = Bio::DB::Fasta->new($fasta_db);
    my @ids = $db->get_all_primary_ids;
    print STDERR scalar @ids."\t sequences read from Fasta file\n";    # Report number of sequences
    foreach my $thecluster (sort keys %cluster_hash) {
        #my @seqarray;
        my $seqio_obj = Bio::SeqIO->new(-file=>">$thecluster.fasta",-format=>'fasta');
        for my $thegenome (sort keys %{$cluster_hash{$thecluster}}) {
            my $theprotein = $cluster_hash{$thecluster}{$thegenome};
            my $seq = $db->get_Seq_by_id($theprotein);
            my $newseq=Bio::Seq->new(-seq=>$seq->seq(), -id=>$thegenome);
            $seqio_obj->write_seq($newseq);
        }
    }
}

sub align_Fasta {       # For each cluster's Fasta file, align with Muscle
    foreach my $thecluster (sort keys %cluster_hash) {
        system ("muscle -in $thecluster.fasta -out $thecluster.muscle.aln");
        my $in=Bio::AlignIO->new(-file=>"$thecluster.muscle.aln",-format=>'fasta');
        my $thein=$in->next_aln();
        $alignment_hash{$thecluster}=$thein;
        #system ("rm $thecluster.muscle.aln");
    }
}

sub concat_alignment {  # Concatenate all the aligned clusters into a single alignment, for phylogenetic analysis
    my @cluster_array;
    foreach my $thekey (sort keys %cluster_hash) {      # Make a list of all the clusters to concatenate
        push @cluster_array, $thekey unless (defined $clusters_with_paralogs{$thekey}); # Skip those clusters which contain paralogs
    }
    #my $current_cluster = $cluster_array[0];
    my $in=Bio::AlignIO->new(-file=>"$cluster_array[0].muscle.aln",-format=>'fasta');   # Read in the alignment file produced in align_Fasta()
    my $thein=$in->next_aln();
    
    my %genomes_in_alignment;
    my $aln_length= $thein->length();
    my $dummy_seq = '-' x $aln_length;
    foreach my $seq ($thein->each_seq()) { 
        $genomes_in_alignment{$seq->id()} = 1; # Record which genomes are represented in the alignment
    }
    foreach my $thegenome (keys %genomes_all) { # For each genome that is not represented in this cluster, add a dummy sequence comprising only gaps
        if (!defined $genomes_in_alignment{$thegenome}) {
            my $newseq = Bio::LocatableSeq->new(-seq=>$dummy_seq,-id=>$thegenome);
            $thein->add_seq($newseq);
        }
    }

    $concat_aln = cat($thein);  # Cat the first alignment to the concatenated alignment (cannot cat to an empty alignment, for some reason)
    
    for (my $x=1; $x <=(scalar @cluster_array - 1); $x++ ) {    # For the rest of the alignments, do the same as above, but in a loop
        my $in = Bio::AlignIO->new(-file=>"$cluster_array[$x].muscle.aln",-format=>'fasta');
        my $thein=$in->next_aln();
        my %genomes_in_alignment;
        my $aln_length= $thein->length();
        my $dummy_seq = '-' x $aln_length;
        foreach my $seq ($thein->each_seq()) { 
            $genomes_in_alignment{$seq->id()} = 1; # Record which genomes are represented in the alignment
        }
        foreach my $thegenome (keys %genomes_all) {
            if (!defined $genomes_in_alignment{$thegenome}) {
                my $newseq = Bio::LocatableSeq->new(-seq=>$dummy_seq,-id=>$thegenome);
                $thein->add_seq($newseq);
            }
        }
        $concat_aln = cat ($concat_aln,$thein);
    }
    print STDERR "Number of concatenated alignments:\t". scalar @cluster_array."\n";
    print STDERR "Length of the concatenated alignment is\t". $concat_aln->length()."\n";

    my $out_fasta = Bio::AlignIO->new(-file=>">$concat_file",-format=>'fasta');	# define the Multifasta output file
    $out_fasta -> write_aln($concat_aln);	# Write the alignments to file defined above
}


=head1 COPYRIGHT AND LICENSE

Copyright 2016-2018, Brandon Seah (kbseah@mpi-bremen.de)

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