#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import uuid
import gzip
import pysam
import random
import string
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Convert FASTQ (or FASTQ.GZ) file to BAM file")
    parser.add_argument("input_fastq", help="Path to the input FASTQ or FASTQ.GZ file")
    parser.add_argument("output_bam", help="Path to the output BAM file")
    parser.add_argument("--pb_version", default="5.0.0", help="pb version in the @HD header line, default is 5.0.0")
    parser.add_argument("--readtype", default="CCS", help="READTYPE in the @RG header line, default is CCS")
    parser.add_argument("--platform_mode", default="SEQUELII", choices=["SEQUELII", "SEQUEL", "PacBioRS", "REVIO"],
                        help="Platform mode in the @RG header line, default is SEQUELII")
    return parser.parse_args()

def fastq_reader(fastq_path):
    """
    Generator to read records from a FASTQ or FASTQ.GZ file one by one.
    Each iteration yields a tuple (header, sequence, plus, quality)
    """
    if fastq_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'  # Text mode
    else:
        open_func = open
        mode = 'r'

    with open_func(fastq_path, mode) as f:
        while True:
            header = f.readline().rstrip()
            if not header:
                break
            sequence = f.readline().rstrip()
            plus = f.readline().rstrip()
            quality = f.readline().rstrip()
            yield (header, sequence, plus, quality)

def parse_pu(header):
    """
    Parse the FASTQ header to extract PU.
    If it matches the pattern @m54082_181216_092744/4194411/ccs,
    extract 'm54082_181216_092744'. Otherwise, return None.
    """
    if header.startswith('@'):
        header = header[1:]
    
    match = re.match(r'^([^/]+)/\d+/ccs$', header)
    if match:
        return match.group(1)
    else:
       # match = re.match(r'^([^/]+)/\d+$', header)  # Adjusted to match non-CCS formats
        match = re.match(r'^([^/]+)(?:/[^ ]+)?(?:\s+\d+\s+length=\d+)?$', header) # Adjusted to match non-CCS formats

        if match:
            return match.group(1)
    return None

def generate_pu():
    """
    Generate a PU with three parts separated by underscores.
    The second part is the current date in YYYYMMDD format.
    Example: abcd1234_20240427_efgh5678
    """
    part1 = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    date_part = datetime.now().strftime('%Y%m%d')
    part3 = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    pu = f"{part1}_{date_part}_{part3}"
    return pu

def generate_rg_id():
    """
    Generate an 8-character alphanumeric RG ID.
    """
    return ''.join(random.choices(string.ascii_letters + string.digits, k=8))

def build_header(pb_version, readtype, platform_mode, rg_id, pu):
    """
    Construct the BAM file header including @HD and @RG lines.
    """
    header = {
        'HD': {
            'VN': '1.6',
            'SO': 'unknown',
            'pb': pb_version
        },
        'RG': [{
            'ID': rg_id,
            'PL': 'PACBIO',
            'DS': f'READTYPE={readtype}',
            'PU': pu,
            'PM': platform_mode
        }]
    }
    return header

def create_aligned_segment(header_str, seq, qual, rg_id, zm_i):
    """
    Create a pysam AlignedSegment from FASTQ read data.
    """
    a = pysam.AlignedSegment()
    a.query_name = header_str[1:] if header_str.startswith('@') else header_str
    a.flag = 4  # Unmapped
    a.reference_id = -1
    a.reference_start = 0  # POS = 0
    a.mapping_quality = 255
    a.cigar = ()
    a.next_reference_id = -1
    a.next_reference_start = 0  # PNEXT = 0
    a.template_length = 0  # TLEN = 0
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)

    # Add tags RG:Z: and zm:i:
    a.set_tag("RG", rg_id)
    a.set_tag("zm", zm_i)

    return a

def parse_zm_i(header):
    """
    Parse the FASTQ header to extract zm_i.
    If it matches the pattern @m54082_181216_092744/4194411/ccs,
    extract '4194411' as an integer. Otherwise, return None.
    """
    match = re.search(r'/(\d+)/ccs$', header)
    if match:
        try:
            zm_i = int(match.group(1))
            return zm_i
        except ValueError:
            return None
    else:
        match = re.search(r'/(\d+)$', header)  # Match for non-CCS patterns
        if match:
            try:
                zm_i = int(match.group(1))
                return zm_i
            except ValueError:
                return None
    return None

def main():
    args = parse_args()

    # Initialize FASTQ reader generator
    reader = fastq_reader(args.input_fastq)

    try:
        # Read the first read
        first_read = next(reader)
    except StopIteration:
        print("Error: FASTQ file is empty.", file=sys.stderr)
        sys.exit(1)

    first_header, first_seq, first_plus, first_qual = first_read

    # Extract PU from the first read's header
    pu = parse_pu(first_header)
    if pu is None:
        pu = generate_pu()

    # Generate RG:Z: ID
    rg_id = generate_rg_id()

    # Build BAM header
    header = build_header(args.pb_version, args.readtype, args.platform_mode, rg_id, pu)

    # Open BAM file for writing
    try:
        with pysam.AlignmentFile(args.output_bam, "wb", header=header) as bam_file:
            # Initialize zm_counter for reads that don't match the pattern
            zm_counter = 1

            # Process the first read
            zm_i = parse_zm_i(first_header)
            if zm_i is None:
                zm_i = zm_counter
                zm_counter += 1

            # Create and write the first alignment segment
            aligned_segment = create_aligned_segment(first_header, first_seq, first_qual, rg_id, zm_i)
            bam_file.write(aligned_segment)

            # Iterate through the rest of the reads
            for record in reader:
                header_str, seq, plus, qual = record
                zm_i = parse_zm_i(header_str)
                if zm_i is None:
                    zm_i = zm_counter
                    zm_counter += 1

                # Create and write the alignment segment
                aligned_segment = create_aligned_segment(header_str, seq, qual, rg_id, zm_i)
                bam_file.write(aligned_segment)

    except Exception as e:
        print(f"Error while writing BAM file: {e}", file=sys.stderr)
        sys.exit(1)

    print("Conversion completed successfully!", file=sys.stderr)

if __name__ == "__main__":
    main()
