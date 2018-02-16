#!/usr/bin/env perl

=head1 NAME

renumber_partitions_from_trimmed_aln.pl - Renumber alignment partitions after trimming

=head1 SYNOPSIS

perl renumber_partitions_from_trimmed_aln.pl -p partitions_block -t position_list

perl renumber_partitions_from_trimmed_aln.pl --help

perl renumber_partitions_from_trimmed_aln.pl --man

=head1 DESCRIPTION

Concatenated sequence alignments can be partitioned into different 'blocks', e.g.
representing different genes in the alignment. The partitions are usually
indicated by a partition file or by charsets in a 'sets' block in a Nexus file.

However, if an alignment has been trimmed, e.g. by BMGE, the partition blocks
will have to be renumbered.

Output is written to STDOUT

=cut 

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;

# Given a set of alignment partitions in Nexus block format and a list of
# alignment positions to be trimmed, renumber the charsets to match the new
# trimmed alignment positions
# Assume 1-based numbering in both

# Inputs
my ($parts_file, $trim_file);
my $remove = 0;
my @trimpos;
# Hash of aln positions and partition assignments
my %pos2part;
# Outputs
my @outarr;
my @charpartitionarr;

pod2usage { -exitstatus=>2, -verbose=>0 } if (!@ARGV);


GetOptions ("parts|p=s" => \$parts_file, # Partitions in Nexus sets block format, 1-based
            "trim|t=s" => \$trim_file, # List of positions for trimmed aln, whitespace-separated, 1-based 
            "remove" => \$remove, # Keep or remove positions in trim file? (default: 0, keep)
            "help|h" => sub { pod2usage(-exitstatus=>2, -verbose=>1); },
            "man|m" => sub { pod2usage(-exitstatus=>0, -verbose=>2); },
            ) or pod2usage(-exitstatus=>2, -verbose=>0);

=head1 ARGUMENTS

=over 8

=item --parts|-p <file>

Partition definition in Nexus 'sets' block format, 1-based numbering

=item --trim|-t <file>

List of positions for the trimmed alignment, separated by whitespace (e.g.
space, newline, tab), 1-based numbering

=item --remove

Logical: Do the positions in the file specified to 'trim' option represent
positions that should be kept, or removed? (Default: kept)

=back

=cut

# Read in partitions/charsets in 
open(my $partsfh, "<", $parts_file) or die ("$!");
while (my $line = <$partsfh>) {
    chomp $line;
    if ($line =~ m/charset (.*) = (\d+)\s*-\s*(\d+)/) {
        my ($set_name, $startpos, $endpos) = ($1, $2, $3);
        for (my $i = $startpos; $i <= $endpos; $i++) {
            $pos2part{$i}{'part'} = $set_name;
            $pos2part{$i}{'flag'} = 1;
        }
    } elsif ($line =~ m/charpartition/) {
        push @charpartitionarr, $line;
    } elsif ($line =~ m/begin sets;/) {
        
    } elsif ($line =~ m/end;/) {
        
    } else {
        say STDERR "Unrecognized line format: $line";
    }
}
close($partsfh);

open(my $trimfh, "<", $trim_file) or die ("$!");
while (my $line = <$trimfh>) {
    chomp $line;
    my @vals = split /\s+/, $line;
    foreach my $val (@vals) {
        if ($val =~ m/\d+/) {
            # Check that it is numeric;
            push @trimpos, $val;
        }
    }
}
close($trimfh);

my $newpos2parthref = trimpos(\%pos2part,\@trimpos);
my $charset_aref = part2pos($newpos2parthref);

push @outarr, "begin sets;";
push @outarr, @$charset_aref;
push @outarr, @charpartitionarr;
push @outarr, "end;";

print STDOUT join "\n", @outarr;
print STDOUT "\n";

## SUBS ########################################################################

sub trimpos {
    my ($href_in, $aref_trim) = @_;
    my %hash = %$href_in;
    my @trim = @$aref_trim;
    my %newhash;
    my $newpos = 0;
    foreach my $pos (@trim) {
        $hash{$pos}{'flag'} = 0;
    }
    foreach my $pos (sort {$a <=> $b} keys %hash) {
        if ($remove == 0) {
            # Keep positions with flag 0
            if ($hash{$pos}{'flag'} == 0) {
                $newpos++;
                $newhash{$newpos}{'part'} = $hash{$pos}{'part'};
            }
        } else {
            # Keep positions with flag 1
            if ($hash{$pos}{'flag'} == 1) {
                $newpos++;
                $newhash{$newpos}{'part'} = $hash{$pos}{'part'};
            }
        }
    }
    return (\%newhash);
}

sub part2pos {
    my ($href_in) = @_;
    my %hash = %$href_in;
    my %hout;
    my @outarr;
    my ($last_part, $curr_part, $pos);
    foreach my $pos (sort {$a <=> $b} keys %hash) {
        $curr_part = $hash{$pos}{'part'};
        if (!defined $last_part) {
            $last_part = $curr_part;
            $hout{$curr_part}{'firstpos'} = $pos;
            $hout{$curr_part}{'currpos'} = $pos;
            $hout{$curr_part}{'string'} .= $pos;
        } else {
            if (defined $hout{$curr_part}{'currpos'} && $pos == $hout{$curr_part}{'currpos'} + 1) {
                # Current partition is same as previous position
                # Continue the range...
                $hout{$curr_part}{'currpos'} = $pos; # Update last position
            } else {
                # Current partition differs from previous
                # Close previous range, start new range
                if ($hout{$last_part}{'currpos'} == $hout{$last_part}{'firstpos'}) {
                    $hout{$last_part}{'string'} .= " "; # Add space if not a range
                } else {
                    $hout{$last_part}{'string'} .= '-'.$hout{$last_part}{'currpos'}." ";
                }
                # Update counters
                $last_part = $curr_part;
                $hout{$curr_part}{'firstpos'} = $pos;
                $hout{$curr_part}{'currpos'} = $pos;
                $hout{$curr_part}{'string'} .= $pos;
            }
        }
    }
    # Sweep up the last element
    if ($hout{$last_part}{'currpos'} == $hout{$last_part}{'firstpos'}) {
        $hout{$last_part}{'string'} .= " "; # Add space if not a range
    } else {
        $hout{$last_part}{'string'} .= '-'.$hout{$last_part}{'currpos'}." ";
    }
    
    # Sort and return output
    foreach my $part (sort {$a cmp $b} keys %hout) {
        $hout{$part}{'string'} =~ s/\s+$//g; # Remove trailing whitespace
        push @outarr, "  charset ".$part." = ".$hout{$part}{'string'}.";";
    }
    
    return (\@outarr);
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