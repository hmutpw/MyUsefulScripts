# MyUsefulScripts

### 1. 

### 2. gtf2tsv_gencode.pl
#### 2.1 The coordinates of bed files are 0-based, and the TSV files are 1-based. GTF/GFF files are 1-based (https://slowkow.com/notes/genomic-intervals/).
Usage: gtf2tsv_gencode.pl [-b <output_bed_file>] [-o <output_file>|<STDOUT>] [<input_file>|<STDIN>]
  <input_file>   Input GTF file from GENCODE (optional, '-' for stdin)
  -b, --bed  Output bed file (optional)
  -o, --output  Output file (optional)
  -h, --help    Display this help message
