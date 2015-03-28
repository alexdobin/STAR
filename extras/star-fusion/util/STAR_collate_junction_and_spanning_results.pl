#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;

my $usage = "\n\n\tusage: $0 Chimeric.out.junction.genes Chimeric.out.sam.spanning.genes\n\n";

my $junction_genes = $ARGV[0] or die $usage;
my $spanning_genes = $ARGV[1] or die $usage;

main: {

    my %fusion_to_counts;

    my %fusion_to_coords_counts;

    {
        open (my $fh, $junction_genes) or die "Error, cannot open file $junction_genes";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $fusion = $x[14];
            my $fusion_coords_text = join("\t", @x[15..$#x]);
            
            if ($fusion =~ /unknown/) { next; }
            
            $fusion_to_counts{$fusion}->{junction}++;
            
            $fusion_to_coords_counts{$fusion}->{$fusion_coords_text}++;
        }
        close $fh;
    }

    {
        open (my $fh, $spanning_genes) or die "Error, cannot open file $spanning_genes";
        while (<$fh>) {
            chomp;
            my ($read_name, $fusion_A) = split(/\t/);
            my ($geneA, $geneB) = split(/--/, $fusion_A);
            
            my $fusion_B = "$geneB--$geneA";
            
            if (exists $fusion_to_counts{$fusion_A}) {
                $fusion_to_counts{$fusion_A}->{spanning}++;
            }
            if (exists $fusion_to_counts{$fusion_B}) {
                $fusion_to_counts{$fusion_B}->{spanning}++;
            }
        }
        close $fh;


    }

    

    foreach my $fusion (keys %fusion_to_counts) {


        my ($geneA, $geneB) = split(/--/, $fusion);
        
        ## no self-fusions
        if ($geneA eq $geneB) { next; }

        my $junction_count = $fusion_to_counts{$fusion}->{junction};
        my $spanning_count = $fusion_to_counts{$fusion}->{spanning} || 0;
        
        my $fusion_coords_info_href = $fusion_to_coords_counts{$fusion};
        
        # get the coordinate junction set with the highest read support
        my @coords_info = reverse sort {$fusion_coords_info_href->{$a} <=> $fusion_coords_info_href->{$b}} keys %$fusion_coords_info_href;

        print join("\t", $fusion, $junction_count, $spanning_count, $coords_info[0]) . "\n";
    }
    
    exit(0);
}
            
