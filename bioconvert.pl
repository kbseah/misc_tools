#!/usr/bin/perl

=head1 NAME

bioconvert.pl - Convert biosequence formats with BioPerl

=head1 SYNOPSIS

perl bioconvert.pl -i <input> -o <output> -fi <input format> -fo <output format> [-ali]

perl bioconvert.pl --help

=head1 DESCRIPTION

Use BioPerl to convert between different file formats.

=head1 ARGUMENTS

=over 8

=item -i <file>

Input filename

=item -o <file>

Name for output file

=item -fi <string>

=item -fo <string>

Specify input and output formats. Choose from: fasta, phylip, nexus.
(Default: both "fasta")

=item -ali

Flag if input is an alignment

=item -revcomp

Flag to report reverse complement sequence (nucleotide only)

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
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;
use Pod::Usage;

my ($input_file, $output_file);
my $format_in = "fasta";
my $format_out = "fasta";
my $is_alignment = 0;
my $revcomp_flag = 0;

if (! @ARGV) {
    pod2usage(-message=>"Invalid input arguments", -exitstatus=>2);
    exit;
}

GetOptions ('input|i=s'=>\$input_file,
            'output|o=s'=>\$output_file,
            'fi=s'=>\$format_in,
            'fo=s'=>\$format_out,
            'ali'=>\$is_alignment,
            'revcomp' =>\$revcomp_flag,
            'help|h' => sub { pod2usage(-exitstatus=>2, -verbose=>2); },
            'man|m' => sub { pod2usage(-exitstatus=>0, -verbose=>2); }
            ) or pod2usage(-exitstatus=>2, -verbose=>2);



if ($is_alignment == 1) {
    convert_alignment();
}

elsif ($is_alignment == 0) {
    convert_sequence();
    }

sub convert_alignment {
    my $input = Bio::AlignIO->new(-file=>"<$input_file",-format=>$format_in);
    my $output = Bio::AlignIO->new(-file=>">$output_file",-format=>$format_out);
    while (my $aln = $input->next_aln) {
        $output->write_aln($aln);
    }
}

sub convert_sequence {
    my $input = Bio::SeqIO->new(-file=>"<$input_file",-format=>$format_in);
    my $output = Bio::SeqIO->new(-file=>">$output_file",-format=>$format_out);
    while (my $seq = $input->next_seq) {
        if ($revcomp_flag == 0) {
            $output->write_seq($seq);
        }
        elsif ($revcomp_flag == 1) {
            $output->write_seq($seq->revcom());
        }
        
    }
}
