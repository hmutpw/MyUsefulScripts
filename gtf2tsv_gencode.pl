#!/usr/bin/perl
# usage: gtf2tsv_gencode.pl -i gencode.gtf (-o output.tsv)
use strict;
use Getopt::Long;

my $gencode_file = "";
my $output_file = "";

GetOptions(
    "i:s" => \$gencode_file,
    "o:s" => \$output_file,
);

if ($gencode_file eq "") {
    die "Please provide an input file using -i option.\n";
}

open(IN, "<$gencode_file") or die "Can't open $gencode_file.\n";

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
    # save values into column 1-8
    $row_data{'chr'} = $chr;
    $row_data{'source'} = $source;
    $row_data{'type'} = $type;
    $row_data{'start'} = $start;
    $row_data{'end'} = $end;
    $row_data{'score'} = $score;
    $row_data{'strand'} = $strand;
    $row_data{'phase'} = $phase;

    push @data, \%row_data;
}

close(IN);

# output column names
my @sorted_column_names = sort keys %all_column_names;
print OUT join("\t", 'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', @sorted_column_names) . "\n";

# output data
for my $row (@data) {
    print OUT join("\t", map { $row->{$_} // '.' } ('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', @sorted_column_names)) . "\n";
}

close(OUT);

if ($output_file) {
    print "Results written to $output_file\n";
}
