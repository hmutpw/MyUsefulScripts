#!/usr/bin/perl

use strict;
use Getopt::Long;

my $gencode_file = "-";
my $output_file = "";
my $bed_file = "";
my $help = 0;

GetOptions(
    "o|output:s" => \$output_file,
    "b|bed:s" => \$bed_file,
    "h|help" => \$help,
);

if ($help) {
    print "Usage: $0 [-b <output_bed_file>] [-o <output_file>|<STDOUT>] [<input_file>|<STDIN>]\n";
    print "  <input_file>   Input GTF file from GENCODE (optional, '-' for stdin)\n";
    print "  -b, --bed  Output bed file (optional)\n";
    print "  -o, --output  Output file (optional)\n";
    print "  -h, --help    Display this help message\n";
    exit;
}

$gencode_file = shift @ARGV if @ARGV;

if (!$gencode_file) {
    die "Please provide an input GTF file.\n";
}

if ($gencode_file eq "-") {
    # Read from standard input
    open(IN, "<&STDIN") or die "Can't read from STDIN.\n";
} else {
    # Open the specified file
    open(IN, "<$gencode_file") or die "Can't open $gencode_file.\n";
}

if ($output_file) {
    open(OUT, ">$output_file") or die "Can't open $output_file for writing.\n";
} else {
    open(OUT, ">&STDOUT") or die "Can't open STDOUT for writing.\n";
}

my @data;
my %all_column_names;

while (<IN>) {
    next if /^#/;

    chomp;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");

    my %row_data;

    my @attributes_data = split(";", $attributes);
    for my $attr (@attributes_data) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type = $1;
        my $c_value = $2;
        $c_value =~ s/\"//g;
        $row_data{$c_type} = $c_value;
        $all_column_names{$c_type} = 1;
    }

    my $id;
    if ($type eq "gene") {
        $id = $row_data{'gene_id'};
    } elsif ($type eq "transcript") {
        $id = $row_data{'transcript_id'};
    } else {
        $id = $row_data{'exon_id'};
    }

    if ($id eq ".") {
        if ($type eq "transcript") {
            $id = $row_data{'gene_id'};
        } else {
            $id = $row_data{'transcript_id'} || $row_data{'gene_id'};
        }
    }

    $row_data{'chr'} = $chr;
    $row_data{'source'} = $source;
    $row_data{'type'} = $type;
    $row_data{'start'} = $start;
    $row_data{'end'} = $end;
    $row_data{'score'} = $score;
    $row_data{'strand'} = $strand;
    $row_data{'phase'} = $phase;
    $row_data{'id'} = $id;

    push @data, \%row_data;
}

close(IN);

my @sorted_column_names = sort keys %all_column_names;
print OUT join("\t", 'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', @sorted_column_names) . "\n";

for my $row (@data) {
    print OUT join("\t", map { $row->{$_} // '.' } ('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', @sorted_column_names)) . "\n";
}

close(OUT);

if ($bed_file) {
    open(BED, ">$bed_file") or die "Can't open $bed_file for writing.\n";

    for my $row (@data) {
        my $id = $row->{'id'} || '.';
        print BED join("\t", $row->{'chr'}, $row->{'start'}, $row->{'end'}, $id, $row->{'type'}, $row->{'strand'}) . "\n";
    }
    close(BED);
}
