# 1. replace genome fa
zcat Danio_rerio.GRCz11.with.SIRV.fa.gz | 
sed -E '
  /^>/ {
    s/^>([0-9]+)(\s|$)/>chr\1\2/;   # 1-25 → chr1-chr25
    s/^>MT(\s|$)/>chrM\1/;          # MT → chrM
    s/^>X(\s|$)/>chrX\1/;           # X → chrX
    s/^>Y(\s|$)/>chrY\1/;           # Y → chrY
  }' | 
gzip > Danio_rerio.GRCz11.with.SIRV.chr_renamed.fa.gz

# 2. replace gtf chr

zcat Danio_rerio.GRCz11.gtf.gz | 
awk -v OFS="\t" '{
  if ($0 ~ /^#/) {  
    print $0;
  } else {
    if ($1 == "MT") $1 = "chrM";
    else if ($1 == "X") $1 = "chrX";
    else if ($1 == "Y") $1 = "chrY";
    else if ($1 ~ /^[0-9]+$/) $1 = "chr"$1;
    print $0;
  }
}' | 
gzip > Danio_rerio.GRCz11.chr_renamed.gtf.gz

