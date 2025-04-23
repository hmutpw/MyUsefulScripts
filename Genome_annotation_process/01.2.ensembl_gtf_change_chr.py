#!/usr/bin/env python3
import sys
import gzip
import re
import argparse
from collections import defaultdict

def convert_chromosome(chrom):
    """Convert only standard chromosomes to GENCODE format."""
    original = chrom
    # Regular expression pattern for chromosomes needing conversion
    if re.fullmatch(r'^(MT|X|Y|\d+)$', chrom, re.IGNORECASE):
        if chrom.upper() == 'MT':
            converted = 'chrM'
        else:
            converted = f'chr{chrom.upper()}'
        return converted, True
    return original, False

def smart_open(filename, mode):
    """Handle both compressed and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    return open(filename, mode)

def generate_report(conversion_map, conversion_counter):
    """Display conversion statistics in terminal."""
    print("\nChromosome Conversion Report:")
    print("{:<20} {:<20} {:<15}".format("Original", "Converted", "Entry Count"))
    print("-" * 55)
    
    for original, converted in sorted(conversion_map.items()):
        count = conversion_counter.get(converted, 0)
        print(f"{original:<20}  ->  {converted:<20} {count:>15,}")
    
    print("\nSummary Statistics:")
    print(f"Total converted chromosomes: {len(conversion_map)}")
    print(f"Total converted entries: {sum(conversion_counter.values()):,}")

def main():
    parser = argparse.ArgumentParser(
        description="Ensembl to GENCODE GTF converter with selective chromosome prefixing",
        epilog="Example: %(prog)s input.gtf.gz output.gtf.gz")
    
    parser.add_argument('input', help='Input GTF file (plain or gzipped)')
    parser.add_argument('output', help='Output GTF file (plain or gzipped)')
    args = parser.parse_args()

    conversion_map = {}
    conversion_counter = defaultdict(int)

    with smart_open(args.input, 'rt') as fin, \
         smart_open(args.output, 'wt') as fout:

        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                fout.write(line)
                continue

            original_chrom = fields[0]
            converted_chrom, changed = convert_chromosome(original_chrom)
            
            if changed:
                if original_chrom not in conversion_map:
                    conversion_map[original_chrom] = converted_chrom
                conversion_counter[converted_chrom] += 1
                fields[0] = converted_chrom
            
            fields[8] = fields[8].replace("_biotype", "_type")
            fout.write('\t'.join(fields) + '\n')

    if conversion_map:
        generate_report(conversion_map, conversion_counter)
    else:
        print("\nNo chromosome conversions were needed.")

if __name__ == "__main__":
    main()
    