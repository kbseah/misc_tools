#!/usr/bin/env perl

=head1 NAME

reduce_table.pl - Reduce a table to selected rows or columns

=head1 SYNOPSIS

perl reduce_table.pl -table INPUT -list SHORTLIST -out OUTPUT

perl reduce_table.pl --help

=cut 
# Select rows or columns from a table by a list of row- or column-names

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $delim = "\t"; # Default assume TSV table
my $by = "col";
my ($infile, $listfile, $outfile);

my %list;   # Hash to store list of row or col names
my @outarr; # Array to hold output lines

GetOptions ("delim|d=s" => \$delim,     # Delimiter for table
            "by=s" => \$by,             # Subset by row or column? Default: "col"
            "table|t=s" => \$infile,    # Input table file (else STDIN)
            "list|l=s" => \$listfile,   # List of row or column names to choose (no default)
            "out|o=s" => \$outfile,     # Output table file (else STDOUT)
            "help|h" => sub {pod2usage(-verbose => 2)},         # Help message
            ) or pod2usage (-verbose => 1);

if (!@ARGV) {
    pod2usage (-verbose => 1);
}

if (!defined $listfile) {
    die ("You did not specify a list of rows or columns to be subsetted.\n");
}

=head1 OPTIONS

=over 8

=item -delim STRING

Delimiter for table. Default: \t (tab character)

=item -by STRING

Subset by "row" or "column"? Default: "col"

=item -table FILENAME

Table to be subsetted. First row is assumed to be header, first column is
assumed to be row names. Default: STDIN

=item -list FILENAME

List of row or column names to be selected from the table file. One entry per
line. Default: None (file must be specified)

=item -out FILENAME

Name of file to write output. Default: STDOUT

=item -help

Help message

=back

=cut

# Get list of row or col names 
open(my $listfh, "<", $listfile) or die ("Cannot open $listfile: $!");
while (my $line = <$listfh>) {
    chomp $line;
    $list{$line}++;     # Save the name to hash
}
close($listfh);

# Screen the input table
# Open the file
my $infh;
if (defined $infile) {
    # If input file specified, open file and assign filehandle
    open($infh, "<", $infile) or die ("Cannot open $infile for reading: $!");
} else {
    # If input file not specified, get from STDIN
    $infh = *STDIN;
}
# Get header line and split by delimiter
my $header = <$infh>;
chomp $header;
my @headsplit = split /$delim/, $header;

# If subset by column names
if ($by =~ m/^c/) {
    # Initialize array of which columns to keep; 0 is first col, row names
    my @colsarr = (0);
    # Run through each col name and check if in table
    for (my $i=1; $i <= $#headsplit; $i++) {
        push @colsarr, $i if (defined $list{$headsplit[$i]});
    }
    # Save the subsetted header line
    push @outarr, join ($delim, @headsplit[@colsarr]);
    # Process the rest of the lines
    while (my $line = <$infh>) {
        chomp $line;
        my @split = split /$delim/, $line;
        # Save only the specified columns
        push @outarr, join($delim, @split[@colsarr]);
    }

# If subset by row names
} elsif ($by =~ m/^r/) {
    # Save header automatically
    push @outarr, $header;
    # Process the remaining lines
    while (my $line = <$infh>) {
        chomp $line;
        my @split = split /$delim/, $line;
        # Save line if row name matches something in shortlist
        if (defined $list{$split[0]}) {
            push @outarr, join ($delim, @split);
        }
    }
}

# Error if neither column nor row specified
else {
    die ("Please specify either by row or col\n");
}

close($infh) if defined $infile;

# Write output
my $outfh;
if (defined $outfile) {
    open ($outfh, ">", $outfile) or die ("Cannot open $outfile: $!");
} else {
    # If no output file specified, write to STDOUT
    $outfh = *STDOUT;
}
# Write each line from saved array
foreach my $line (@outarr) {
    print $outfh $line."\n";
}

# Close filehandle if open() was called before
close($outfh) if defined $outfile;


=head1 COPYRIGHT AND LICENSE

Copyright (C) 2017 by  Brandon Seah <kbseah@mpi-bremen.de>

LICENCE

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