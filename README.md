# MyUsefulScripts

### 1. 

### 2. gtf2tsv_gencode.pl
#### 2.1 The coordinates of bed files are 0-based, and the TSV files are 1-based. GTF/GFF files are 1-based (https://slowkow.com/notes/genomic-intervals/).
Usage: gtf2tsv_gencode.pl [-b <output_bed_file>] [-o <output_file>|<STDOUT>] [<input_file>|<STDIN>]
  <input_file>   Input GTF file from GENCODE (optional, '-' for stdin)
  -b, --bed  Output bed file (optional)
  -o, --output  Output file (optional)
  -h, --help    Display this help message

### 3. genpred2transcriptLoc.pl
#### 3.1. mapping transcript location to genomic location from genePred files

1. get genePred from gtf files:
```{shell}
gtfToGenePred -genePredExt in.gtf out.genePred
```
2. get localization files
Usage: genpred2transcriptLoc.pl [--output OUTPUT] [--genePredExt] INPUT

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
    perl genpred2transcriptLoc.pl --output output.txt input.gpd
    perl genpred2transcriptLoc.pl --genePredExt --output output.txt input.gpdExt
    perl genpred2transcriptLoc.pl input.gpd > output.txt


