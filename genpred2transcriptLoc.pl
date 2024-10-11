#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# ????
my $input_file;
my $output_file;
my $genePredExt = 0;  # ????? genePredExt ??

# ???????
GetOptions(
    'output=s'      => \$output_file,
    'genePredExt'   => \$genePredExt,
    'help'          => sub { print_usage(); exit; },
) or die "Invalid options passed. Use --help for usage.\n";

# ??????
$input_file = shift or die "Usage: $0 [--output OUTPUT] [--genePredExt] INPUT\n";

# ?? genePred ? genePredExt ??
open(my $genpred_fh, '<', $input_file) or die "Could not open $input_file: $!";

# ????????? STDOUT
my $output_fh;
if (defined $output_file) {
    open($output_fh, '>', $output_file) or die "Could not open $output_file: $!";
} else {
    $output_fh = \*STDOUT;
}

# ????
print $output_fh join("\t", "gene_id", "tx_id", "chr", "genomic_pos", "strand", "genomic_pos_strand", "tx_pos", "exon_number", "total_exons"), "\n";

# ?????
while (my $line = <$genpred_fh>) {
    chomp $line;
    my @fields = split "\t", $line;

    # ??????
    if ($genePredExt) {
        if (@fields < 15) {
            warn "Skipping line with insufficient fields for genePredExt: $line\n";
            next;
        }
    } else {
        if (@fields < 10) {
            warn "Skipping line with insufficient fields for genePred: $line\n";
            next;
        }
    }

    # ????
    my ($gene_id, $transcript_id, $chromosome, $strand, @exon_starts, @exon_ends);

    if ($genePredExt) {
        $transcript_id = $fields[0];
        $chromosome    = $fields[1];
        $strand        = $fields[2];
        # ?? genePredExt ?????????
        @exon_starts = split ',', $fields[8];
        @exon_ends   = split ',', $fields[9];
        $gene_id      = $fields[11];
    } else {
        $transcript_id = $fields[0];
        $chromosome    = $fields[1];
        $strand        = $fields[2];
        # ?? genePred ?????????
        @exon_starts = split ',', $fields[8];
        @exon_ends   = split ',', $fields[9];
        $gene_id      = $fields[0];  # ?? genePred ? gene_id ? transcript_id ??,?????
    }

    # ????????(???)
    pop @exon_starts if defined $exon_starts[-1] && $exon_starts[-1] eq '';
    pop @exon_ends if defined $exon_ends[-1] && $exon_ends[-1] eq '';

    # ?? exon_starts ? exon_ends ????
    if (@exon_starts != @exon_ends) {
        warn "Mismatch between exon starts and ends in line: $line\n";
        next;
    }

    my $transcript_position = 1;
    my $total_exons = scalar @exon_starts;

    if ($strand eq '+') {
        for (my $i = 0; $i < $total_exons; $i++) {
            my $exon_start = $exon_starts[$i];
            my $exon_end   = $exon_ends[$i];
            my $exon_length = $exon_end - $exon_start;

            # ?? exon_length ??
            if ($exon_length <= 0) {
                warn "Invalid exon length in line: $line\n";
                next;
            }

            for (my $j = 0; $j < $exon_length; $j++) {
                my $genomic_position = $exon_start + $j + 1;  # ???1-based??
                my $exon_number = $i + 1;

                # ???????
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
    } elsif ($strand eq '-') {
        for (my $i = $total_exons - 1; $i >= 0; $i--) {
            my $exon_start = $exon_starts[$i];
            my $exon_end   = $exon_ends[$i];
            my $exon_length = $exon_end - $exon_start;

            # ?? exon_length ??
            if ($exon_length <= 0) {
                warn "Invalid exon length in line: $line\n";
                next;
            }

            for (my $j = $exon_length - 1; $j >= 0; $j--) {
                my $genomic_position = $exon_start + $j + 1;  # ???1-based??
                my $exon_number = $i + 1;

                # ???????
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
        warn "Unknown strand '$strand' in line: $line\n";
    }
}

# ??????
close($genpred_fh);
close($output_fh) if defined $output_file;

# ????
sub print_usage {
    print << "USAGE";
Usage: $0 [--output OUTPUT] [--genePredExt] INPUT

Description:
    This script processes a genePred or genePredExt file and outputs the annotation of each basepair
    in each transcript. The output includes gene id, transcript id, chromosome,
    strand, genomic position, transcript position, exon number, and total exons.
    
    INPUT is the path to the input genePred or genePredExt file.
    --output OUTPUT is the path to the output file. If not provided, the output will be
    written to standard output.
    --genePredExt    Specify that the input file is in genePredExt format.
    
Options:
    --help      Show this help message and exit.

Example:
    perl $0 --output output.txt input.gpd
    perl $0 --genePredExt --output output.txt input.gpdExt
    perl $0 input.gpd > output.txt

USAGE
}
