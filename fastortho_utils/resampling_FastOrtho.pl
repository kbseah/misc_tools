#!/usr/bin/perl

# Parse FastOrtho results (.end) file and options file (.opt) and report singleton proteins
# Report "core" genome (clusters present in at least one copy in every genome)
# Generate accumulation curve for number of clusters as number of genomes is varied

use strict;
use warnings;

my ($fastortho_results, $options_file);

if (!defined @ARGV) {
    print "\n\t";
    print "Usage: perl resampling_FastOrtho.pl fastortho.end fastortho.opt > resampling_table\n\n";
    print "\tCaution: Performs sampling WITH replacement\n\n";
    exit;
}
else { ($fastortho_results, $options_file) = @ARGV; }

my %proteins_clusters_hash;  # Hash of all protein IDs (keys) and the clusters they belong to (value)
my %proteins_genomes_hash;  # Hash of all protein IDs (keys) and the genomes they belong to (value)
my %genome_files;   # List of single_genome_fasta files (keys) and paths (values)
my $number_of_genomes=0;    # Stores count of total number of genomes
my %clusters_genomes_hash;  # Hash-of-hash with counts of number of cluster members (value) per cluster (key1) per genome (key2)

# Read list of genome files (fasta) from options file
open(OPTIONS, "<", $options_file) or die ("$!\n");
while (<OPTIONS>) {
    chomp;
    my $filename;
    my $pathname;
    if ($_ =~ /--single_genome_fasta (\S+)/) {
        $number_of_genomes++;   # Count total number of genomes listed
        $pathname = $1;
        if ($pathname =~ /.*\/(\S+)/) {
            $filename = $1;
        }
        else {$filename = $pathname;}
    }
    $genome_files{$filename} = $pathname unless (!defined $filename);
}
close(OPTIONS);

# Read list of proteins and associate them to genome names
foreach my $thekey (keys %genome_files) {
    open(GENOME, "<", "$genome_files{$thekey}") or die ("$!\n");
    while (<GENOME>) {
        if ($_ =~ /^>(\S+)/) {
            $proteins_genomes_hash{$1} = $thekey;
            $proteins_clusters_hash{$1} = $1;  # Default cluster is the protein name itself unless later overwritten
        }
    }
    close(GENOME);
}

# Read FastOrtho results file to associate protein IDs to cluster names
open(FASTORTHO, "<", $fastortho_results) or die ("$!\n");
while (<FASTORTHO>) {
    chomp;
    my @theline = split /:\t/, $_;
    if ($theline[0] =~ /^(ORTHOMCL\d+) /) {
        my $currentcluster = $1;
        my @theproteins = split /\s/, $theline[1];
        foreach my $theprotein (@theproteins) {
            if ($theprotein =~ /^(\S+)\((\S+)\)/) {
                $proteins_clusters_hash{$1} = $currentcluster;
            }
        }
    }
}
close(FASTORTHO);

# Report each protein ID with its associated cluster and genome 
foreach my $thekey (sort keys %proteins_clusters_hash) {
    #print $proteins_clusters_hash{$thekey}."\t".$thekey."\t".$proteins_genomes_hash{$thekey}."\n";
    $clusters_genomes_hash{$proteins_clusters_hash{$thekey}}{$proteins_genomes_hash{$thekey}}++;
}

#print "cluster\t"; # Header for report counts per cluster per genome
#foreach my $filename (sort keys %genome_files) { print $filename."\t";}
#print "\n";

foreach my $thekey (sort keys %clusters_genomes_hash) {
    my $sum_counts=0;
    my $num_genomes=0;
    foreach my $thekey2 (sort keys %{$clusters_genomes_hash{$thekey}}) {
        $num_genomes++;
        $sum_counts+= $clusters_genomes_hash{$thekey}{$thekey2};
    }
    #print $thekey."\t".$sum_counts."\t".$num_genomes."\n";     # Report number of counts per cluster
    
    #print $thekey."\t";    # Report counts per cluster per genome
    #foreach my $filename (sort keys %genome_files) {
    #    if (defined $clusters_genomes_hash{$thekey}{$filename}) {
    #        print $clusters_genomes_hash{$thekey}{$filename}."\t";
    #    }
    #    else {print "0\t";}
    #}
    #print "\n";
    
    #if ($num_genomes == $number_of_genomes) { print $thekey."\n"; }    # Report clusters of "core genome"
}

my @genomes_array;
foreach my $filename (sort keys %genome_files) {
    push @genomes_array, $filename;
}

# Actual observed counts
my $num_genomes_real = scalar @genomes_array;
my $num_clusters_real = scalar (keys %clusters_genomes_hash);
my $num_core_real = 0;
my $num_singletons_real = 0;
foreach my $thecluster (sort keys %clusters_genomes_hash) {
    if (scalar (keys %{$clusters_genomes_hash{$thecluster}}) == scalar @genomes_array) {
        $num_core_real++;
    }
    elsif (scalar (keys %{$clusters_genomes_hash{$thecluster}}) == 1) { $num_singletons_real++; }
}

print "# Observed counts\n";
print "# num_genomes\tpan\tcore\tsingleton\n";
print "# ";
print join "\t", ($num_genomes_real, $num_clusters_real, $num_core_real, $num_singletons_real);
print "\n\n\n";

## Perform resampling

my $reps = 200;  # Number of pseudoreplicates
my @headerline = ("num_genomes","pan","core","singleton");  # Header line for resampling output
print join ("\t", @headerline). "\n";

# NB: SAMPLING WITH REPLACEMENT!
for (my $gen=1; $gen <= $number_of_genomes; $gen++) {    # For number of pseudoreplicates between 1 and $number_of_genomes 
    #print $gen."\n";
    for (my $x=1; $x <= $reps; $x++) {  # Resample $reps number of times
        my %resampled_genomes;
        my %resampled_clusters_genomes_hash;
        my $num_clusters=0;
        my $num_core=0;
        my $num_singletons=0;
        
        for (my $i=1; $i <= $gen; $i++) {       # Generate list of resampled genomes
            my $random_number = int (rand($number_of_genomes));
            #print $random_number." ";
            $resampled_genomes{$genomes_array[$random_number]} = 1;
        }
        
        my $true_gen = scalar (keys %resampled_genomes); # Count number of resampled genomes; because of sampling WITH replacement, this is less than or equal to $gen; affects the core count!
        
        foreach my $thecluster (sort keys %clusters_genomes_hash) {     # Generate the resampled clusters by genomes hash
            foreach my $thegenome (sort keys %resampled_genomes) {
                #print $thecluster."\t".$thegenome."\n";
                $resampled_clusters_genomes_hash{$thecluster}{$thegenome} = $clusters_genomes_hash{$thecluster}{$thegenome} unless (!defined $clusters_genomes_hash{$thecluster}{$thegenome});
            }
        }
        
        $num_clusters = scalar (keys %resampled_clusters_genomes_hash);
        foreach my $thecluster (sort keys %resampled_clusters_genomes_hash) {
            if (scalar (keys %{$resampled_clusters_genomes_hash{$thecluster}}) == $true_gen) {
                $num_core++;
            }
            elsif (scalar (keys %{$resampled_clusters_genomes_hash{$thecluster}}) == 1) {
                $num_singletons++;
            }
        }
        print $gen."\t". $num_clusters."\t".$num_core."\t".$num_singletons;
        print "\n";
    }
    #print "\n";
}


######
# Results can be imported into R and plotted:

#d <- read.table("result.out",header=T)
#plot(d$num_genomes,d$pan,pch=".",ylim=c(0,max(d$pan)))
#points(d$num_genomes,d$core,pch=".",col="blue")
#points(d$num_genomes,d$singleton,pch=".",col="red")
#points(tapply(X=d$pan,INDEX=d$num_genomes,FUN=mean),type="l")
#points(tapply(X=d$core,INDEX=d$num_genomes,FUN=mean),type="l",col="blue")
#points(tapply(X=d$singleton,INDEX=d$num_genomes,FUN=mean),type="l",col="red")
#legend(x="topleft",col=c("black","red","blue"),legend=c("pan","singleton","core"),lty=1)