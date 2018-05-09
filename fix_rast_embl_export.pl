#!/usr/bin/env perl

=head1 NAME

fix_rast_embl_export.pl - Fix RAST EMBL-flatfile export for submission to ENA

=head1 SYNOPSIS

perl fix_rast_embl_export.pl --embl [flatfile] --locus [LOCUS] --os [organism] > stdout

perl fix_rast_embl_export.pl --help

=head1 DESCRIPTION

Tidy RAST EMBL-format feature table export for validation with the ENA flatfile
validator tool.

Assumptions:
 * Prokaryotic genome 
 * Only source, CDS, tRNA and rRNA features
 * CDS elements all have codon_start=1

What it does:
 * Add locus tags with supplied prefix, and numbers parsed from RAST SEED IDs
 * Add codon_start and transl_table attributes to CDSs
 * Remove translation field (ENA validator adds translation based on transl_table)
 * Unwrap product names that span multiple lines and add spaces

=cut 

# To do:
# Running contig number counter

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Pod::Usage;
use Time::Piece;

my $t = localtime();
my $date = uc($t->strftime ("%d-%b-%Y"));

my ($embl_in, $classification);
my $locus_prefix = "TESTLOCUS";
my $transl_table = "11";
my $organism = "TEST ORGANISM";
my $strain = "TEST STRAIN";
my $spacing = 20; # Spacing between "FT" label and beginning of annotation
my $project_id = 'XXX';
my $ra = 'XXX';
my $rg = 'XXX';


pod2usage (-verbose=>0) if (!@ARGV);

GetOptions("embl=s"=>\$embl_in,
           "locus=s"=>\$locus_prefix,
           "project=s"=>\$project_id,
           "os=s" =>\$organism,
           "oc=s" =>\$classification,
           "ra=s" => \$ra,
           "rg=s" => \$rg,
           "strain=s"=>\$strain,
           "transl=i"=>\$transl_table,
           "help|h" => sub{ pod2usage(-verbose=>1,-exitstatus=>1); },
           "man" => sub{ pod2usage(-verbose=>2,-exitstatus=>0); },
           ) or die ("$!");

=head1 ARGUMENTS

=over 8

=item --embl I<FILE>

EMBL feature table flatfile exported from RAST

=item --locus I<STRING>

Locus tag prefix

Default: TESTLOCUS

=item --os I<STRING>

Organism name, should match entry in NCBI taxonomy database

Default: "TEST ORGANISM"

=item --oc I<STRING>

Organism classification, should match entry in NCBI taxonomy database

=item --strain I<STRING>

Strain name. Required for prokaryotic sequences.

Default: "TEST STRAIN"

=item --project I<STRING>

Project ID.

Default: "XXX"

=item --ra I<STRING>

Reference author(s) in format "Lastname Initial, Lastname Initial". No comma
between lastname and initial, period after intials.

Default: "XXX"

=item --rg I<STRING>

Reference group (work group/institution of author)

Default: "XXX"

=item --transl I<INTEGER>

Translation table number. Default: 11

=item --help

Help message

=item --man

Full manual page

=back

=cut

open(my $fhin, "<", $embl_in) or die ("$!");
my @outlines;
my $prevline;
my $current_feature_type;
my $current_length;
while (my $line = <$fhin>) {
    chomp $line;
    
    if ($line =~ m/^ID   (.+)/) {
        # Reformat ID line
        say STDERR "Replacing ID line with defaults assuming prokaryotic genome";
        my @fields = split /;/, $1;
        $fields[0] = 'XXX';
        $fields[1] = ' SV XXX';
        $fields[3] = ' genomic DNA';
        $fields[5] = ' PRO';
        if ($fields[6] =~ m/(\d+) BP/) {
            $current_length = $1;
        }
        
        push @outlines, "ID   ".join(';', @fields);
        push @outlines, 'XX';
    } elsif ($line =~ m/^AC   /) {
        say STDERR "Replacing AC (accession number) line with XXX";
        say STDERR "Adding PR (project identifier) line with accession $project_id";
        push @outlines, 'AC   ;';
        push @outlines, 'XX';
        push @outlines, "PR   Project:$project_id" if defined $project_id;
        push @outlines, 'XX';
        #$line = join "\n", @newlines;
    } elsif ($line =~ m/^DE/) {
        say STDERR "Replacing DE (description) field with XXX";
        say STDERR "Adding KW (keyword) field with blank";
        push @outlines, 'DE   XXX';
        push @outlines, 'XX';
        push @outlines, 'KW   .';
        push @outlines, 'XX';
        say STDERR "Adding OS (organism species) line with $organism" if defined $organism;
        say STDERR "Adding OC (organism classification) line with $classification" if defined $classification;
        push @outlines, "OS   $organism" if defined $organism;
        push @outlines, "OC   $classification" if defined $classification;
        push @outlines, 'XX';
        say STDERR "Adding RN, RP, RG, RA, RL lines with dummy citation";
        push @outlines, 'RN   [1]';
        push @outlines, "RP   1-$current_length";
        push @outlines, "RG   $rg";
        push @outlines, "RA   $ra;";
        push @outlines, 'RT   ;';
        push @outlines, "RL   Submitted ($date) to the INSDC.";
    } else {
        # Count length of spacing between data class line def and the entry
        # Will be used later...
        if ($line =~ m/^FH(\s+Key\s+)Location/) {
            $spacing = length $1;
        }
        
        if ($line =~ m/^FT   (\S+)\s+\S+/) {
            # Get current feature type
            $current_feature_type = $1;
        }
        
        if (defined $current_feature_type && $current_feature_type eq 'source') {
            # Skip source attributes that are RAST-specific
            next if ($line =~ m/\/db_xref=/);
            next if ($line =~ m/\/genome_md5=/);
            next if ($line =~ m/\/project=/);
            next if ($line =~ m/\/genome_id=/);
            
            # Substitute RAST organism name with the registered name
            if ($line =~ m/\/organism=/) {
                $line = "FT                   /organism=\"$organism\"";
                # Add strain name if defined (must have strain name for prokaryote)
                $line .= "\nFT                   /strain=\"$strain\"";
            }
            
        }
        
        # Add locus tags based on RAST PEG and RNA IDs
        if ($line =~ m/\/db_xref="SEED:fig\|(.+)\.peg\.(\d+)/) {
            # Check for SEED peg ID and add locus tag with given prefix
            my $locus_tag = "$locus_prefix\_PEG_$2";
            $line .= "\nFT                   /locus_tag=\"$locus_tag\"";
            # Add codon start and transl_table fields
            $line .= "\nFT                   /codon_start=1";
            $line .= "\nFT                   /transl_table=$transl_table";
        } elsif  ($line =~ m/\/db_xref="SEED:fig\|(.+)\.rna\.(\d+)/) {
            # Check for SEED RNA id and add locus tag with given prefix
            my $locus_tag = "$locus_prefix\_RNA_$2";
            $line .= "\nFT                   /locus_tag=\"$locus_tag\"";
        }
        
        # Unwrap feature qualifier values that span multiple lines
        if ($line =~ m/^FT\s{$spacing}([^\/].*)$/) {
            # If current line is a feature that doesn't start with a qualifier
            # tag, it must be a qualifier value that has been broken over multiple
            # lines. If so, add it to the previous line and move on
            my $piece = $1;
            $prevline .= " $piece" if $prevline =~ m/\/product=/; # Add space if product name
            $prevline .= "$piece" if $prevline =~ m/\/translation=/; # Do not add space if translation
        } else {
            # Print lines
            if (defined $prevline) {
                push @outlines, $prevline unless $prevline =~ m/\/translation=/;
            }
            $prevline = $line;
        }
    }
}
# Print last line
push @outlines, $prevline;
close($fhin);

foreach my $line (@outlines) {
    print STDOUT $line;
    print STDOUT "\n";
}
