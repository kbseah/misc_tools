#!/usr/bin/env perl
=head1 NAME

svg_color_scale.pl - Colorize elements in SVG file by user-supplied values and color scale

=head1 SYNOPSIS

perl svg_color_scale.pl --svg <svg file>
                        --tab <color table>
                        --out <output svg file>
                        --scale <output scalebar svg file>

perl svg_color_scale.pl --help

=head1 DESCRIPTION

This tool was designed to simplify re-drawing diagrams, e.g. metabolic
pathways, which are annotated with colored symbols corresponding to some
numerical parameter, such as gene expression level.

=head1 ARGUMENTS

=over 8

=item --svg <file>

SVG file to be modified

=item --tab <file>

Comma- or tab-separated table. Column 1 - SVG object IDs, column 2 - numerical
parameter to be converted to color value. The range of the parameter is
automatically used to define color scale over that range, and objects with the
corresponding IDs have their fill colors changed according to that scale.

=item --out <file>

File name for SVG output file with changed colors

=item --scale <file>

File name for SVG output file containing a scale bar.

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
# Brandon Seah (kbseah@mpi-bremen.de
# Version 001 -- 2016-01-19

use strict;
use warnings;
use SVG;
use SVG::Parser;
use List::Util qw(min max);
use Getopt::Long;
use Pod::Usage;

#my $parser = new SVG::Parser(-debug => 1); # Debugging version
my $parser = new SVG::Parser();

my $thefile; # Specify input SVG file
my $outfile; # Output SVG filename
my $valsTab; # Table with list of SVG object IDs and values to colorify
my $scalebar; # File for scale bar output
my $userMin;
my $userMax; # user-specified min and max values for scale

if (!@ARGV) {
    pod2usage (-message=>"No arguments given", -exitstatus=>2);
}


GetOptions ("svg=s"=>\$thefile,
            "out=s"=>\$outfile,
            "tab=s"=>\$valsTab,
            "scale=s"=>\$scalebar,
            "min=i"=>\$userMin,
            "max=i"=>\$userMax,
            "help|h"=> sub { pod2usage(-exitstatus=>2, -verbose=>2); },
            "man|m"=> sub { pod2usage (-exitstatus=>0, -verbose=>2); }
            ) or pod2usage (-message=>"Invalid arguments", -exitstatus=>2) ;

my %vals_hash; # Hash to store values from table
my @minmax; # Array to store min and max values for values

## MAIN #################################################

# Catch exceptions
if ($thefile eq "") {
    pod2usage (-message=>"Invalid arguments", -exitstatus=>2);
} elsif ($valsTab eq "") {
    pod2usage (-message=>"Invalid arguments", -exitstatus=>2) ;
}

my $svg = $parser->parse_file($thefile); # Parse SVG file

%vals_hash = readValsTab($valsTab); # Read values from table

# use user-specified minmax values if given
if (defined ($userMin)) {
    $minmax[0] = $userMin;
    $minmax[1] = $userMax;
} else { # Else automatically scale from the values table
    $minmax[0] = min (values %vals_hash);
    $minmax[1] = max (values %vals_hash);
}

foreach my $theID (keys %vals_hash) {
    print STDERR "Replacing fill for element ID $theID \n";
    # Replace fill parameter for objects, given the color scale
    replaceFillGivenID($theID,
                       generateColorsFromVal($minmax[0],
                                             $minmax[1],
                                             $vals_hash{$theID},
                                             "redgreen")
                       );
}

if (! $outfile eq "") { # If output filename specified, write to file
    open(SVGOUT, ">", $outfile) or die ("$!\n");
    print SVGOUT $svg->xmlify;
    close(SVGOUT);
} else { # Else print to STDOUT
    print STDOUT $svg->xmlify;
}

if (! $scalebar eq "") { # If scale bar desired, write to file
    my $scaleSVG = makeScaleBar($minmax[0],$minmax[1],"redgreen");
    open(SCALEOUT, ">", $scalebar) or die ("$!\n");
    print SCALEOUT $scaleSVG->xmlify;
    close(SCALEOUT);
} else {
    print STDERR "No scale bar drawn\n";
}


## SUBROUTINES ############################################


sub makeScaleBar {
    my ($min, $max, $type) = @_; 
    my $incr = ($max - $min ) / 10; # Default ten steps in scale
    # Create SVG object of the scale bar
    my $svgScale = SVG->new(width=>110,height=>20);
    my %svgTagHash;
    my %svgLabelHash;
    # Draw colored boxes for scale bar
    for (my $j=0; $j <=10; $j++) {
        $svgTagHash{$j} = $svgScale->rect(x=>$j*10,
                                          y=>10,
                                          width=>10, height=>10,
                                          fill=>generateColorsFromVal($min,$max,($min + $j*$incr),$type));
        #$svgLabelHash{$j} = $svgScale->text(id=>$min+$j*$incr,
        #                                    x=>$j*10,
        #                                    y=>10,
        #                                    "font-family"=>"Sans",
        #                                    "font-size"=>4)->cdata($min+$j*$incr);
    }
    # Write text labels for min and max values only
    $svgScale->text(id=>"minLabel",
                    x=>2,
                    y=>10,
                    "font-family"=>"Sans",
                    "font-size"=>4)->cdata($min);
    $svgScale->text(id=>"maxLabel",
                    x=>102,
                    y=>10,
                    "font-family"=>"Sans",
                    "font-size"=>4)->cdata($max);
    return $svgScale;
}

sub readValsTab {
    my ($table_file) = @_;
    my %output_hash;
    open(VALTAB, "<", $table_file) or die ("$!\n");
    while (<VALTAB>) {
        chomp;
        my @line_split = split /[,;\t]/, $_; #CSV or TSV file format
        $output_hash{$line_split[0]} = $line_split[1];
    }
    close(VALTAB);
    return %output_hash;
}

sub generateColorsFromVal {
    # Generate colors given min, max, val, and palette type
    my ($min, $max, $val, $type) = @_;
    if ($type eq "redgreen") { # Red-Green color scale
        # Code adapted from http://blogs.perl.org/users/ovid/2010/12/perl101-red-to-green-gradient.html
        my $middle = ($min + $max) / 2;
        my $scale = 255 / ($middle - $min);
        return "#ff0000" if $val <= $min;
        return "#00ff00" if $val >= $max;
        if ($val < $middle) {
            return sprintf "#ff%02x00" => int ( ($val - $min) * $scale);
        } else {
            return sprintf "#%02xff00" => 255 - int ( ( $val - $middle) * $scale);
        }
    }
}

sub replaceFillGivenID {
    # Given object of given ID, replace "fill" value of "style" attribute with new value
    # Assumes that fill color is child of attribute "style" (e.g. for "path" objects)
    my ($theID, $newFill) = @_;
    my $ref = $svg->getElementByID($theID); # Generates a ref for object of given ID
    my $attr = $ref->getAttributes(); # Generates a ref for attributes of given object
    ${$attr}{"style"}{"fill"} = $newFill; # Change the "fill" value in the "style" attribute

}
