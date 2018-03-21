#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

tidy_img_gff3.pl - Tidy up GFF and FNA files from IMG for conversion to EMBL

=head1 SYNOPSIS

perl tidy_img_gff3.pl --gff xxx.gff --fna xxx.fna --xls xxx.xls --convert

perl tidy_img_gff3.pl --help

perl tidy_img_gff3.pl --man

=head1 DESCRIPTION

Given IMG annotation (GFF feature flatfile, FNA sequences, XLS annotation table)
tidy up the GFF output for conversion to EMBL flatfile format with EMBLmyGFF3
(https://github.com/NBISweden/EMBLmyGFF3)

This script has been tested with EMBLmyGFF3 v1.2.1 and ENA flat file validator
v1.1.200, and output from IMG pipeline v4.

It assumes that genomes are prokaryotic. A few ad hoc fixes for the output from
EMBLmyGFF3 are also implemented. Before submitting anything to ENA you should
first check with the ENA flat file validator. Make sure to get the latest 
version because the software is constantly updated.

An internet connection is required for EMBLmyGFF3 to check the NCBI taxonomy
database.

=cut

my ($gff, $fna, $xls);
my $out = 'test';
my $do_convert;
#my ($species, $locus_tag, $project_id);
my $strain;
my $emblmygff3_args;
my %prod; # Product name hash from xls file

pod2usage (-verbose=>1) if !@ARGV;

GetOptions ("gff=s" => \$gff,   # GFF file from download bundle
            "fna=s" => \$fna,   # Genome contig FNA file from download bundle
            "xls=s" => \$xls,   # Annotation spreadsheet from IMG
            "strain=s" => \$strain, # Strain name, if desired
            "convert" => \$do_convert, # Flag - do conversion with EMBLmyGFF3? Default: no
            "emblmygff3=s"=> \$emblmygff3_args,
            #"species=s" => \$species,
            #"locus_tag=s" => \$locus_tag,
            #"project_id=s" => \$project_id,
            "out=s" => \$out,   # Output file prefix
            "help|h" => sub{ pod2usage (-verbose=>2, -exitstatus=>2) },
            "man" => sub{ pod2usage (-verbose=>2, -exitstatus=>1) },
            ) or pod2usage (-verbose=>2);

=head1 ARGUMENTS

=over 8

=item --gff <file>

GFF file containing feature coordinates. This should be from the download_bundle.tar.gz
file within the IMG zip archive. The name should be <genome_id>.gff. There are
three GFF files within the IMG zip archive - this one should contain product
names for CDSs, and the product names should be fixed (e.g. "XXX protein"
corrected to "XXX-domain containing protein" when evidence is insufficient). This
differs from the other <genome_id>.gff file in the IMG_Data folder, where the
product names are not fixed! 

=item --fna <file>

Contig nucleotide sequence in Fasta format. Contig names should match GFF file.
Should also be in download_bundle.tar.gz and be named <genome_id>.fna

=item --xls <file>

IMG annotation in spreadsheet format (tab-separated variables, but from the IMG
platform the file is named <accession>.xls).

=item --strain <string>

Strain name, if the organism is a prokaryote. Required for annotating "source"
qualifier.

=item --help|-h

Print help to screen.

=item --man

Manual page.

=item --convert

Call EMBLmyGFF3 to do conversion of GFF flatfile immediately, with the following
defaults assuming prokaryotic genome: --topology linear --molecule_type 'genomic
DNA' --transl_table 11 --taxonomy PRO

Assumes that EMBLmyGFF3 is in path.

=item --emblmygff3="[args]"

Additional arguments to be passed to EMBLmyGFF3, e.g. species name, accession,
locus tag prefix... 

=back

=cut

## MAIN ########################################################################

# Read annotation from XLS file
read_prods($xls,\%prod);

# Process FNA and GFF files
my $outgff_aref = fix_gff_file($gff);
my $outfna_aref = fix_fna_file($fna);

# Write output
print STDERR "Writing tidied GFF file to $out.gff...\n";
open(my $fhout, ">", "$out.gff") or die ("$!");
print $fhout join "\n", @$outgff_aref;
close($fhout);

print STDERR "Writing tidied FNA file to $out.fna...\n";
open(my $fhout2, ">", "$out.fna") or die ("$!");
print $fhout2 join "\n", @$outfna_aref;
close($fhout2);

# Run conversion with EMBLmyGFF3
if ($do_convert) {
    run_EMBLmyGFF3();            # Run EMBLmyGFF3
    fix_EMBLmyGFF3_flatfile();   # Fix some issues in output from EMBLmyGFF3 v1.2.1
    print STDERR "You can now run the EMBL flat file validator with the -fix option \n";
    print STDERR "to check the $out.embl flatfile for errors and to fix known issues \n";
} else {
    print STDERR "You can now run EMBLmyGFF3 to convert the GFF3 and Fasta files to \n";
    print STDERR "EMBL flatfile format \n";
}


## SUBS ########################################################################

sub fix_EMBLmyGFF3_flatfile {
    # Fix product lines broken over multiple lines by EMBLmyGFF3
    # Add strain information to source qualifiers
    print STDERR "\nFixing raw output from EMBLmyGFF3 ... \n";
    print STDERR "Writing final EMBL flatfile to $out.embl\n";
    
    my $outfile = "$out.embl";
    my $spacing = 20; # Default value
    
    my $prevline;
    open(my $fh, "<", "$out.raw.embl") or die ("$!"); # Input raw EMBLmyGFF3 output
    open(my $fhout, ">", $outfile) or die ("$!"); # Fixed EMBL flatfile
    while (my $line = <$fh>) {
        chomp $line;
        # Count length of spacing between data class line def and the entry
        # Will be used later...
        if ($line =~ m/^FH(\s+Key\s+)Location/) {
            $spacing = length $1;
        }
        # Add strain name
        if ($line =~ m/^FT\s{$spacing}\/organism/ && defined $strain) {
            $line .= "\nFT                   /strain=\"$strain\"";
        }
        # Locus tags to Uppercase
        if ($line =~ m/^(FT\s{$spacing}\/locus_tag)=\"(.+)\"$/) {
            my $locus_tag_uc = uc $2;
            $line = "$1=\"$locus_tag_uc\"";
        }
        # Update ID lines
        if ($line =~ m/^FT\s+\/note="ID:\d+"/) {
            $line =~ s/note="ID:/note="IMG_gene_oid:/;
        }
        # Unwrap feature qualifier values that span multiple lines
        if ($line =~ m/^FT\s{$spacing}([^\/].*)$/) {
            # If current line is a feature that doesn't start with a qualifier
            # tag, it must be a qualifier value that has been broken over multiple
            # lines. If so, add it to the previous line and move on
            $prevline .= $1; 
        } else {
            print $fhout $prevline."\n" if defined $prevline;
            $prevline = $line;
        }
    }
    # Print last line - important!
    print $fhout $prevline."\n";
    close ($fh);
    close ($fhout);
}

sub run_EMBLmyGFF3 {
    print STDERR "Running EMBLmyGFF3 with the following command: \n";
    #print STDERR "Using species $species, locus tag $locus_tag, and project id $project_id\n";
    my @cmd = ("EMBLmyGFF3",
               "$out.gff",
               "$out.fna",
               "--topology linear",
               "--molecule_type \'genomic DNA\'",
               "--data_class STD",
               "--taxonomy PRO",
               "--transl_table 11",
               "-o $out.raw.embl"
               );
    # Optional params
    push @cmd, $emblmygff3_args if defined $emblmygff3_args;
    #push @cmd, "--species \'$species\'" if defined $species;
    #push @cmd, "--locus_tag $locus_tag" if defined $locus_tag;
    #push @cmd, "--project_id $project_id" if defined $project_id;
    # Run command
    print STDERR join " ", @cmd;
    print STDERR "\n";
    print STDERR "...\n";
    system (join " ", @cmd);
}

sub fix_gff_file {
    # Fix various issues in GFF file and return fixed file as array (in ref)
    my ($file) = @_;
    my @outarr;
    open(my $fh, "<", $file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        my @split = split /\t/, $line;
        if (defined $split[2]) {
            next if $split[2] =~ m/CRISPR/; # Ignore CRISPR annotations
            #fix_IMG_Gene_ID(\$line);
            quote_product_names(\$line);
            add_tRNA_product_name(\$line,\%prod) if $split[2] eq 'tRNA';
            fix_RNA_entry(\$line,\%prod) if $split[2] eq 'RNA'; 
        }
        push @outarr, $line;
    }
    close($fh);
    return (\@outarr);
}

sub fix_fna_file {
    # Strip header line to contig ID only, everything after first space
    # Return fixed file in array (as ref)
    my ($file) = @_;
    my @outarr;
    open(my $fh, "<", $file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^>/) {
            my @split = split /\s/, $line;
            push @outarr, $split[0];
        } else {
            push @outarr, $line;
        }
    }
    close($fh);
    return (\@outarr);
}

sub fix_RNA_entry {
    # "RNA" is not a valid feature type. Most of these are misc_RNA features
    my ($ref, $href) = @_;
    my @split = split /\t/, $$ref;
    my ($locus_tag) = $split[8] =~ m/locus_tag=([^;]+)/;
    if ($split[8] !~ m/product=/) { # Check that product naem not already in place
        $split[2] = $href->{$locus_tag}{'Locus_type'};
        my $product = $href->{$locus_tag}{'Product_name'};
        $split[8] .= ";product=\"$product\";";
    }
    $$ref = join "\t", @split;
    
}

sub add_tRNA_product_name {
    # Add product name to tRNA
    my ($ref, $href) = @_;
    my @split = split /\t/, $$ref;
    my ($locus_tag) = $split[8] =~ m/locus_tag=([^;]+)/;
    if ($split[8] !~ m/product=/) { # Check that product name is not already in place
        my $product = $href->{$locus_tag}{'Product_name'};
        $split[8] .= ";product=\"$product\";";
    }
    $$ref = join "\t", @split;
}

sub fix_IMG_Gene_ID {
    # Note: ID field is required for correctly-formed GFF file
    # Change "ID" to "IMG_Gene_ID" which will end up in a "note" feature
    my ($ref) = @_; # Use reference so that string can be modded in-place
    my @split = split /\t/, $$ref;
    $split[8] =~ s/^ID=/IMG_Gene_ID=/;
    $$ref = join "\t", @split;
}

sub quote_product_names {
    # Quote product names so that they will not be broken over lines, unless already quoted
    my ($ref) = @_;
    my @split = split /\t/, $$ref;
    $split[8] =~ s/product=([^"][^;]+[^"][^;])$/product=\"$1\";/;
    $$ref = join "\t", @split;
}

sub read_prods {
    my ($file, $href) = @_;
    open(my $fh, "<", $file) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        my @split = split /\t/, $line;
        if (defined $split[2]) {
            if ($split[2] eq 'Locus_type') {
                $split[4] =~ s/miscRNA/misc_RNA/; # Add underscore to match feature type in EMBL format def
                $href->{$split[1]}{'Locus_type'} = $split[4];
            } elsif ($split[2] eq 'Product_name') {
                $href->{$split[1]}{'Product_name'} = $split[4];
            }
        } 
    }
    close ($fh);
}
