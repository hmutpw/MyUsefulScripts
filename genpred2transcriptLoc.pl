#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $input_file;
my $output_file;

GetOptions(
    'output=s' => \$output_file,
    'help'     => sub { print_usage(); exit; },
) or die "Invalid options passed. Use --help for usage.\n";

# check the input
$input_file = shift or die "Usage: $0 [--output OUTPUT] INPUT\n";

# open genPred
open(my $genpred_fh, '<', $input_file) or die "Could not open $input_file: $!";

# open output or STD
my $output_fh;
if (defined $output_file) {
    open($output_fh, '>', $output_file) or die "Could not open $output_file: $!";
} else {
    $output_fh = \*STDOUT;
}

# output header
print $output_fh join("\t", "gene_id", "tx_id", "chr", "start", "strand", "genomic_pos", "tx_pos", "exon_number", "total_exons"), "\n";

while (my $line = <$genpred_fh>) {
    chomp $line;
    my @fields = split "\t", $line;

    # get information from genpred
    my $gene_id = $fields[0];
    my $transcript_id = $fields[1];
    my $chromosome = $fields[2];
    my $strand = $fields[3];
    my @exon_starts = split ',', $fields[9];
    my @exon_ends = split ',', $fields[10];

    # remove null element
    pop @exon_starts if $exon_starts[-1] eq '';
    pop @exon_ends if $exon_ends[-1] eq '';

    my $transcript_position = 1;
    if ($strand eq '+') {
        for (my $i = 0; $i < scalar @exon_starts; $i++) {
            my $exon_start = $exon_starts[$i];
            my $exon_end = $exon_ends[$i];
            my $exon_length = $exon_end - $exon_start;

            for (my $j = 0; $j < $exon_length; $j++) {
                my $genomic_position = $exon_start + $j + 1;  # Adjust to 1-based position
                my $exon_number = $i + 1;
                my $total_exons = scalar @exon_starts;

                # print basepair information
                print $output_fh join("\t",
                    $gene_id,
                    $transcript_id,
                    $chromosome,
                    $genomic_position,
                    $strand,
                    "$chromosome:$genomic_position:$strand",
                    $transcript_position,
                    $exon_number,
                    $total_exons
                ), "\n";

                $transcript_position++;
            }
        }
    } else {
        for (my $i = scalar @exon_starts - 1; $i >= 0; $i--) {
            my $exon_start = $exon_starts[$i];
            my $exon_end = $exon_ends[$i];
            my $exon_length = $exon_end - $exon_start;

            for (my $j = $exon_length - 1; $j >= 0; $j--) {
                my $genomic_position = $exon_start + $j + 1;  # Adjust to 1-based position
                my $exon_number = $i + 1;
                my $total_exons = scalar @exon_starts;

                # print basepair information
                print $output_fh join("\t",
                    $gene_id,
                    $transcript_id,
                    $chromosome,
                    $genomic_position,
                    $strand,
                    "$chromosome:$genomic_position:$strand",
                    $transcript_position,
                    $exon_number,
                    $total_exons
                ), "\n";

                $transcript_position++;
            }
        }
    }
}

close($genpred_fh);
close($output_fh) if defined $output_file;

sub print_usage {
    print << "USAGE";
Usage: $0 [--output OUTPUT] INPUT

Description:
    This script processes a genpred file and outputs the annotation of each basepair
    in each transcript. The output includes gene id, transcript id, chromosome,
    strand, genomic position, transcript position, exon number, and total exons.
    
    INPUT is the path to the input genpred file.
    --output OUTPUT is the path to the output file. If not provided, the output will be
    written to standard output.

Options:
    --help      Show this help message and exit.

Example:
    perl $0 --output output.txt input.gpd
    perl $0 input.gpd > output.txt

USAGE
}
