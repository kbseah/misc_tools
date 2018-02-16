#!/usr/bin/env perl

# Concatenate multiple species alignments
# If some species are missing from some alignments (incomplete data),
# fill in those spaces with gap character "-"

=head1 NAME

concat_alignments_with_gaps.pl - Concatenate Fasta alignments

=head1 SYNOPSIS

    perl concat_alignments_with_gaps.pl -a <table of alignment file names> \
                                        -s <list of species names> \
                                        -w <path to working folder> \
                                        -o <prefix for output files> \
                                        -m <find best AA substitution model with RAxML?>
                                        -t <num threads for RAxML>

    perl concat_alignments_with_gaps.pl --help [Full help message]

    perl concat_alignments_with_gaps.pl --man [Manual page]


=head1 DESCRIPTION

Concatenate several multisequence alignments in Fasta format. Given individual
alignment files, each alignment representing a separate gene, and the sequence
names corresponding to OTUs or species, concatenate them into a single super-
gene alignment. The sequence names in the alignment files must correspond to
the names in the list of species names. Missing data (species not represented
in a given alignment) will be filled in with gap character "-".

=head1 ARGUMENTS

=over 8

=item --alignment|-a <string>

File of tab-separated table of alignment file names and their corresponding
marker names. Marker names will be used to report partition positions in the
output.partition file.

=item --species|-s <string>

File containing list of species names. Species names must match the
sequence names used in the alignment files. If a sequence name in the
alignment file is not represented in this list, then those sequences will be
excluded from the final concatenated alignment.

=item --wd|-w <path>

Path to the folder containing input files; output and intermediate files will
also be written to this folder.
Default: current folder

=item --output|-o <string>

Prefix for output file names.
Default: "test"

=item --model|-m <logical>

Flag - Use RAxML to find best-fitting AA substitution model for each partition?
Default: no

=item --threads|-t <integer>

Number of threads for RAxML
Default: 8

=item --help|-h

Full help message

=item --man|-m

Manual page

=back

=head1 COPYRIGHT AND LICENSE

Copyright 2016, Brandon Seah (kbseah@mpi-bremen.de)

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

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::Align::Utilities qw(:all);
use Pod::Usage;

my $alnTableFile;       # Table of alignment files and marker names (tab-separated)
my %alnTable;
my @markers;
my $speciesTableFile;   # List of species names
my %speciesTable;
my $wdPath=".";         # Path to working folder
my $outfix="test";      # Output prefix
my %concatHash;         # Hash to store alignments for concatenation
my $concat_aln;         # Concatenated alignment
my %genePos;
my $raxmlExecutable="/usr/local/bin/raxmlHPC-PTHREADS-AVX";
my $numThreads = 8;
my %best_model_choice;
my $model_choose = 0;
my $use_mask = 0;
my $concat_len;

GetOptions (
    "alignments|a=s" => \$alnTableFile,
    "species|s=s" => \$speciesTableFile,
    "wd|w=s" => \$wdPath,
    "output|o=s" => \$outfix,
    "model|m" => \$model_choose,
    "zorro|z" => \$use_mask,
    "threads|t=i" => \$numThreads,
    'help|h'=> sub { pod2usage( -exitstatus => 2, -verbose => 2); },
    'man|m'=> sub { pod2usage ( -exitstatus => 0, -verbose => 2) },
);

if (!defined $alnTableFile|| !defined $speciesTableFile) {
    pod2usage(-message => "Insufficient options were supplied", -existatus => 2);
}


## MAIN #######################################################################

read_alntable();
read_specieslist();

@markers = sort {$a cmp $b} keys %alnTable;
$concat_len = scalar @markers;
for my $marker (@markers) {
    # Read alignment file
    my $in = Bio::AlignIO->new(-file=>"$wdPath/$alnTable{$marker}",
                               -format=>"fasta"
                               );
    my $aln = $in->next_aln();
    
    # Make dummy sequence matching alignment length
    my $alnLen = $aln->length();
    my $dummyseq = "-" x $alnLen;
    my $alphabet;
    
    # Check which species are in this alignment
    foreach my $seq ($aln->each_seq()) {
        my $theid = $seq->id();
        $alphabet = $seq->alphabet();
        if (defined $speciesTable{$theid}) {
            $speciesTable{$theid}++;
        } else {
            # Remove sequences not in the speciesTable list
            $aln->remove_seq($seq);
        }
    }
    
    # Write alignment to phylip format for RAxML model testing
    my $outfile = Bio::AlignIO->new(-file=>">$wdPath/$outfix.$marker.phy",
                                    -format=>"phylip");
    $outfile -> write_aln($aln);
    
    # For species that are absent from this alignment, make dummy sequence
    # comprising gap characters
    foreach my $species (keys %speciesTable) {
        if ($speciesTable{$species} == 0) {
            print STDERR "Species $species not found in $marker alignment\n";
            my $newseq = Bio::LocatableSeq->new(-seq=>$dummyseq,
                                             -id => $species,
                                             -alphabet => $alphabet
                                             );
            # Add dummy sequence to alignment
            $aln->add_seq($newseq);
        }
    }

    # Add alignment with dummies to hash
    $concatHash{$marker} = $aln;
    
    # Reset counter for species presence/absence
    foreach my $species (keys %speciesTable) {
        $speciesTable{$species} = 0;
    }
}

if ($model_choose == 1) {
    model_choice_wrapper();

    foreach my $marker (keys %best_model_choice) {
        print $marker."\t".$best_model_choice{$marker}."\n";
    }
}

my $concat_fasta= Bio::AlignIO->new(-file=>">$wdPath/$outfix.concat.fasta",-format=>'fasta');	# define the Multifasta output file
my $concat_phylip = Bio::AlignIO->new(-file=>">$wdPath/$outfix.concat.phy",-format=>'phylip');	# define the Phylip output file
write_concat_alignments();

$concat_fasta->write_aln($concat_aln);
$concat_phylip->write_aln($concat_aln);

###############################################################################

sub read_alntable {
    open(IN, "<", $alnTableFile) or die ("$!\n");
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        # Column 1 - filename; Column 2 - marker name
        $alnTable{$splitline[1]} = $splitline[0];
    }
    close(IN);
}

sub read_specieslist {
    open(IN, "<", $speciesTableFile) or die ("$!\n");
    while (<IN>) {
        chomp;
        $speciesTable{$_} = 0;
    }
    close(IN);
}

sub getLH {
    # Subroutine required for model test script by Alexis Stamatakis
    my $fileID = $_[0];  
    open(CPF, $fileID);
    my @lines = <CPF>;	
    close(CPF);	
    my $numIT = @lines;   	
    my $lastLH = pop(@lines);  
    my $k = index($lastLH, '-');   
    my $LH = substr($lastLH, $k);     
    return $LH;
}

sub choose_prot_model {
    # Adapted from script by Alexis Stamatakis
    my $alignmentName = $_[0]; # Read the MARKER name
    my $UNLIKELY = -1.0E300;
    my @lh;

    my @AA_Models = ("DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", 
        "LG", "MTART", "MTZOA", "PMB", "HIVB", "HIVW", "JTTDCMUT", "FLU",
        "DAYHOFFF", "DCMUTF", "JTTF", "MTREVF", "WAGF", "RTREVF", "CPREVF", "VTF", "BLOSUM62F", 
        "MTMAMF", "LGF", "MTARTF", "MTZOAF", "PMBF", "HIVBF", "HIVWF", "JTTDCMUTF", "FLUF");
    #print "Determining AA model data\n";
    #print "Computing randomized stepwise addition starting tree number :".$i."\n";
    my $cmd = $raxmlExecutable." -T $numThreads -y -p 12345 -m PROTCATJTT -s $wdPath/$outfix.$alignmentName.phy -n ST_".$alignmentName." \> ST_".$alignmentName."_out";
    system($cmd);
    my $numberOfModels = @AA_Models;
    for(my $i = 0; $i < $numberOfModels; $i++) {
        my $aa = "PROTGAMMA".$AA_Models[$i];
        my $cmd = $raxmlExecutable." -T $numThreads -f e -m ".$aa." -s $wdPath/$outfix.$alignmentName\.phy -t RAxML_parsimonyTree.ST_".$alignmentName." -n ".$AA_Models[$i]."_".$alignmentName."_EVAL \> ".$AA_Models[$i]."_".$alignmentName."_EVAL.out\n";  
        system($cmd);
    }

    for(my $i = 0; $i < $numberOfModels; $i++) {
        my $logFileName = "RAxML_log.".$AA_Models[$i]."_".$alignmentName."_EVAL";
        $lh[$i] = getLH($logFileName);
    }

    my $bestLH = $UNLIKELY;
    my $bestI = -1;

    for(my $i = 0; $i < $numberOfModels; $i++) {
        if($lh[$i] > $bestLH) {
            $bestLH = $lh[$i];
            $bestI = $i;
        }
    }

    $best_model_choice{$alignmentName} = $AA_Models[$bestI];	# Write best model chosen to the hash of best model choices
    
    # cleanup RAxML files
    my $cleanupcmd = "rm $wdPath/RAxML\*$alignmentName\_EVAL $wdPath/RAxML\*\_$alignmentName $wdPath/\*EVAL.out $wdPath/ST_\*out";
    system ($cleanupcmd);
}

sub model_choice_wrapper {		# Wrapper routine to perform model test for each marker gene in the list
    foreach my $marker (@markers) {
	choose_prot_model($marker);
    }
}

sub write_concat_alignments {		# Concatenate the alignments into new alignment file
    $concat_aln = cat($concatHash{$markers[0]});			# Can't recursively cat on empty string, so must start from first element
    $genePos{$markers[0]}= $concat_aln->length();
    
    # Initialize a file to contain the list of marker genes and their positions in the alignment
    open (my $PARTITION, "> $wdPath/$outfix\.partitions") || die ("Cannot write partitions file: $!");
    if ( $model_choose == 0 ) {
        print $PARTITION join("" , "WAG" , ", " , $markers[0] , " = ", join("-","1",$genePos{$markers[0]})), "\n"; }
    elsif ( $model_choose == 1 ) {
        print $PARTITION join("" , $best_model_choice{$markers[0]} , ", ", $markers[0], " = ", join("-","1",$genePos{$markers[0]})), "\n"; }
    close ($PARTITION);

    if ( $use_mask == 1 ) {
        system ("cat $wdPath/$outfix.$markers[0]\.mask >> $wdPath/$outfix\.mask");
    } 
    
    for (my $x = 1; $x <= ($concat_len - 1); $x++){			# Iterate through the hash and recursively concatenate the alignments
        $concat_aln = cat($concat_aln, $concatHash{$markers[$x]});
        $genePos{$markers[$x]}= $concat_aln->length();	# store the end position of this marker gene in the alignment
        # Add the marker and its positions to the partitions file for later use with phylogenetic software
        open (my $PARTITION, ">> $wdPath/$outfix\.partitions") || die ("Cannot write partitions file: $!");
        if ( $model_choose == 0 ) {
            print $PARTITION join("", "WAG", ", " , $markers[$x], " = ", join("-",($genePos{$markers[$x-1]} + 1),$genePos{$markers[$x]}) ), "\n"; }
        elsif ( $model_choose == 1 ) {
            print $PARTITION join("", $best_model_choice{$markers[$x]}, ", " , $markers[$x], " = ", join("-",($genePos{$markers[$x-1]} + 1),$genePos{$markers[$x]}) ), "\n"; }
        close ($PARTITION);
        if ( $use_mask == 1 ) {
            system ("cat $wdPath/$outfix.$markers[$x]\.mask >> $wdPath/$outfix\.mask");
        } 
    }
}
