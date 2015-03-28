#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Set::IntervalTree;
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 star.spanning.sam ref_annot.gtf\n\n";

my $spans_file = $ARGV[0] or die $usage;
my $gtf_file = $ARGV[1] or die $usage;


main: {
    
    my %chr_to_interval_tree = &parse_junctions($gtf_file);
    
    my $sam_reader = new SAM_reader($spans_file);
    
    my @read_info;
    my $prev_core_read_name = "";
    
    while (my $sam_entry = $sam_reader->get_next()) {
        
        my $core_read_name = $sam_entry->get_core_read_name();
        my $full_read_name = $sam_entry->reconstruct_full_read_name();
        
        if ($core_read_name ne $prev_core_read_name) {
            &describe_read_pair_info($core_read_name, \@read_info);
            @read_info = ();
        }


        $full_read_name =~ m|/([12])$| or die "Error, cannot decipher read name $full_read_name as /1 or /2 ";
        my $read_dir = $1;
        
        my $chr = $sam_entry->get_scaffold_name();

        unless (exists $chr_to_interval_tree{$chr}) { next; } # no annotations to search overlaps of.
        
        my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();

        foreach my $exon_segment (@$genome_coords_aref) {
            my ($lend, $rend) = @$exon_segment;
            my $overlaps_aref = $chr_to_interval_tree{$chr}->fetch($lend, $rend);
            if (@$overlaps_aref) {
                &add_overlapping_genes(\@read_info, $overlaps_aref, $read_dir);                
            }
        }
        
        $prev_core_read_name = $core_read_name;
    }

    # get last one
    &describe_read_pair_info($prev_core_read_name, \@read_info);

    exit(0);
    
}

####
sub describe_read_pair_info {
    my ($read_name, $read_info_aref) = @_;

    my @read_1_genes;
    my @read_2_genes;

    my $read_1_genes_href = $read_info_aref->[1];
    if ($read_1_genes_href) {
        @read_1_genes = sort keys (%$read_1_genes_href);
    }

    my $read_2_genes_href = $read_info_aref->[2];
    if ($read_2_genes_href) {
        @read_2_genes = sort keys (%$read_2_genes_href);
    }
    
    if (@read_1_genes && @read_2_genes) {
        my $fusion_token = join(",", @read_1_genes) . "--" . join(",", @read_2_genes);
                
        print "$read_name\t$fusion_token\n";
    }
    
    return;
}


####
sub add_overlapping_genes {
    my ($read_info_aref, $overlaps_aref, $read_dir) = @_;
    
    foreach my $gene (@$overlaps_aref) {
        $read_info_aref->[$read_dir]->{$gene} = 1;
    }

    return;
}


####
sub parse_junctions {
    my ($gtf_file) = @_;

    my %chr_exons_interval_tree;

    my $counter = 0;
    
    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
                
        my @x = split(/\t/);
        my $chr = $x[0];
        my $type = $x[2];
    
        unless ($type eq 'exon') { next; }
        
        $counter++;
        print STDERR "\r[$counter]    " if $counter % 1000 == 0;
        

        my $lend = $x[3];
        my $rend = $x[4];
        
        my $orient = $x[6];


                
        my $info = $x[8];

        my $gene_id;
        if ($info =~ /gene_id \"?([^;\"]+)\"?;/) {
            $gene_id = $1;
        }

        if ($info =~ /gene_name \"?([^;\"]+)\"?;/) {
            $gene_id = $1;
        }

        unless ($gene_id) {
            confess "Error, cannot find gene_id from $info";
        }

        
        if ($rend - $lend > 1) {
            &add_chr_interval(\%chr_exons_interval_tree, $chr, $gene_id, $lend, $rend);
        }

    }

    return(%chr_exons_interval_tree);
}
        
####
sub add_chr_interval {
    my ($chr_exons_interval_tree_href, $chr, $gene_id, $lend, $rend) = @_;

    my $i_tree = $chr_exons_interval_tree_href->{$chr};
    unless ($i_tree) {
        $i_tree = $chr_exons_interval_tree_href->{$chr} = Set::IntervalTree->new;
    }

    $i_tree->insert($gene_id, $lend, $rend);

    return;
}
        
