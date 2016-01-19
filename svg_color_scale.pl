#!/usr/bin/env perl

# Brandon Seah (kbseah@mpi-bremen.de
# Version 001 -- 2016-01-19

use strict;
use warnings;
use SVG;
use SVG::Parser;
use List::Util qw(min max);
use Getopt::Long;

#my $parser = new SVG::Parser(-debug => 1); # Debugging version
my $parser = new SVG::Parser();

my $thefile; # Specify input SVG file
my $outfile; # Output SVG filename
my $valsTab; # Table with list of SVG object IDs and values to colorify
my $scalebar; # File for scale bar output
my $userMin;
my $userMax; # user-specified min and max values for scale

GetOptions ("svg=s"=>\$thefile,
            "out=s"=>\$outfile,
            "tab=s"=>\$valsTab,
            "scale=s"=>\$scalebar,
            "min=i"=>\$userMin,
            "max=i"=>\$userMax) or die ("$!");

my %vals_hash; # Hash to store values from table
my @minmax; # Array to store min and max values for values

## MAIN #################################################

# Catch exceptions
if ($thefile eq "") {
    usage();
} elsif ($valsTab eq "") {
    usage();
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

sub usage {
    print STDERR "***********************************************************************\n";
    print STDERR "Colorize elements in SVG file by user-supplied values and color scale\n\n";
    print STDERR "Usage: \n";
    print STDERR "       perl $0 --svg input.svg \\ \n ";
    print STDERR "                             --tab values.csv \\ \n";
    print STDERR "                             --out output.svg \\ \n";
    print STDERR "                             --scale scalebar_output.svg\n\n";
    print STDERR "User-supplied input table (comma- or tab-separated) contains SVG object IDs\n";
    print STDERR "in first column, and numerical parameters in second column. The numerical \n";
    print STDERR "parameters are used to make a color scale over their range, and the objects\n";
    print STDERR "with the corresponding IDs have their fill colors changed according to that \n";
    print STDERR "color scale. \n\n";
    print STDERR "This tool was designed to simplify re-drawing diagrams, e.g. metabolic \n";
    print STDERR "pathways, which are annotated with colored symbols corresponding to some \n";
    print STDERR "numerical parameter, such as gene expression level.\n\n";
    print STDERR "***********************************************************************\n";
    exit;
}

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
