#!/usr/bin/env perl

use strict;
use warnings;
use Carp;


my $usage = "usage: $0 star.junctions ref_annotation.gtf\n\n";

my $junctions_file = $ARGV[0] or die $usage;
my $annot_gtf_file = $ARGV[1] or die $usage;



main: {
    
    my %gene_to_junctions = &parse_junctions($annot_gtf_file);

    if ($junctions_file =~ /\.gz$/) {
        $junctions_file = "gunzip -c $junctions_file | ";
    }
    
    open (my $fh, $junctions_file) or die "Error, cannot open file $junctions_file";
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my ($chrA, $coordA) = ($x[0], $x[1]);
        my ($chrB, $coordB) = ($x[3], $x[4]);
        

        my @acceptors;
        my @donors;
        if (my $gene = $gene_to_junctions{"A:$chrA;$coordA"}) {
            push (@acceptors, "$gene\t$chrA:$coordA");
        }
        if (my $gene = $gene_to_junctions{"A:$chrB;$coordB"}) {
            push (@acceptors, "$gene\t$chrB:$coordB");
        }
        if (my $gene = $gene_to_junctions{"D:$chrA;$coordA"}) {
            push (@donors, "$gene\t$chrA:$coordA");
        }
        if (my $gene = $gene_to_junctions{"D:$chrB;$coordB"}) {
            push (@donors, "$gene\t$chrB:$coordB");
        }
        
        foreach my $donor (@donors) {
            my ($donor_gene, $donor_coords) = split(/\t/, $donor);
            foreach my $acceptor (@acceptors) {
                my ($acceptor_gene, $acceptor_coords) = split(/\t/, $acceptor);
                my $fusion = "$donor_gene--$acceptor_gene";
                print "$line\t$fusion\t$donor\t$acceptor\n";
            }
        }
        


    }

    close $fh;
    
    exit(0);
    
}

####
sub parse_junctions {
    my ($gtf_file) = @_;

    my %junctions;

    my $counter = 0;

    if ($gtf_file =~ /\.gz$/) {
        $gtf_file = "gunzip -c $gtf_file | ";
    }
    
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
        
        $lend--;  ## storing first bases of candidate splice sites.
        $rend++;


        my $orient = $x[6];

        my ($acceptor, $donor) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);


                
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

        $junctions{"A:${chr};$acceptor"} = $gene_id;
        $junctions{"D:${chr};$donor"} = $gene_id;
        
    }

    return(%junctions);
}
        
