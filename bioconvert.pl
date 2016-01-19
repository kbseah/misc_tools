#!/usr/bin/perl

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;

my ($input_file, $output_file);
my $format_in = "fasta";
my $format_out = "fasta";
my $is_alignment = 0;
my $revcomp_flag = 0;

if (!defined @ARGV) {
    usage();
    exit;
}

GetOptions ('input|i=s'=>\$input_file,
            'output|o=s'=>\$output_file,
            'fi=s'=>\$format_in,
            'fo=s'=>\$format_out,
            'ali'=>\$is_alignment,
            'revcomp' =>\$revcomp_flag) or die ("$!\n");



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

sub usage {
    print STDERR "\nConvert bio sequence formats with BioPerl\n\n";
    print STDERR "\tUsage: perl bioconvert.pl -i input_file -o output_file -fi <input_format> -fo <output_format> (-ali)\n\n";
    print STDERR "\t-i\tInput filename\n\t-o\tOutput filename\n\t-fi\tInput format, e.g. fasta, phylip, nexus\n\t\-fo\tOutput format\n\t-ali\tFlag if sequences are aligned\n\t-revcomp\tFlag if you want reverse complement\n\n";
}