#!/usr/bin/perl

#commond: perl 03.1.1.get_ERCC_transcript_fasta.pl ERCC92.gtf ERCC92.fa ERCC92.transcript.fa
open GTF_IN,"$ARGV[0]" or die $!;
open FASTA_IN,"$ARGV[1]" or die $!;
open OUTPUT, "> $ARGV[2]" or die $!;

# get transcript gene id relation.
my %Transcrpt;
while(<GTF_IN>){
        chomp;
        my @line=split(/\t/,$_);
        my $gene_id=$line[0];
        my @infor=split(/\"/,$line[8]);
        my $transcript_id=$infor[3];
        $Transcrpt{$gene_id}=$transcript_id;
}

#add transcript ids into fasta header.
while(<FASTA_IN>){
        chomp;
        if($_=~/^>/){
        my $line=$_;
        $line=~/^>(ERCC-\d+)/;
        my $gene_id=$1;
        my $transcript_id=$Transcrpt{$gene_id};
        my $newline=">".$transcript_id.'|'.$gene_id;
        print OUTPUT "$newline\n";
        }else{
        print OUTPUT "$_\n";
        }
}
