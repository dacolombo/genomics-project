#!/usr/bin/env bash

set -euxo pipefail

########################################
###########     Help     ###############
########################################
Help()
{
  echo "Description:"
  echo "  Pipeline to trim, align to a reference and assemble Oxford Nanopore reads."
  echo
  echo "Syntax:"
  echo "  AssemblyOxN.sh -i <input_fastq_file> -r <reference_fasta_file> -t <number_of_threads>"
  echo
  echo "Options:"
  echo "  -i <input_fastq_file>		input fastq file containing the OxN reads (can be in .fastq or .fastq.gz format)"
  echo "  -r <reference_fasta_file>       fasta file containing the reference genome"
  echo "  -t <number_of_threads>          number of threads to use"

}



########################################
###########     Main     ###############
########################################

#-------------#
# Get options #
#-------------#

while getopts ":i:r:t:h" opt; do
  case $opt in
    i) reads="$OPTARG"
    ;;
    r) reference="$OPTARG"
    ;;
    t) threads="$OPTARG"
    ;;
    h) Help
       exit
    ;;
    \?) echo "Invalid option -$OPTARG"
        exit
    ;;
  esac
done


if [ "$#" -eq 0 ]
then
  Help
  exit
fi


prefix=${reads%%.*}


#---------------------#
# Preprocessing reads #
#---------------------#

echo "Trimming reads..."
porechop --threads=$threads -i ${reads} -o ${prefix}_trimmed.fastq.gz


#---------#
# Mapping #
#---------#

echo "Mapping reads and sorting the bam file..."
minimap2 -t $threads -ax map-ont $reference ${prefix}_trimmed.fastq.gz | samtools view -@ $threads -b -o ${prefix}.bam -;
samtools sort -@ $threads -o ${prefix}.sort.bam ${prefix}.bam;
rm ${prefix}.bam;

echo "Extracting mapped reads..."
samtools view -@ $threads -b -F 4 ${prefix}.sort.bam -o ${prefix}_mapped.sort.bam;

echo "Converting mapped reads to fastq..."
bedtools bamtofastq -i ${prefix}_mapped.sort.bam -fq ${prefix}_mapped.fastq;
gzip ${prefix}_mapped.fastq;

echo "Computing coverage..."
samtools faidx ${reference}
bedtools genomecov -pc -d -ibam ${prefix}_mapped.sort.bam -g ${reference%%.*}.fai > ${prefix}_mapped.covbed.txt;

#----------#
# Assembly #
#----------#

echo "Performing the sampling...";
samples=(5000 10000 20000 30000);

mkdir samples;

for size in "${samples[@]}";
do
  seqtk sample -s1000 ${prefix}_mapped.fastq.gz $size > samples/${prefix}_$((size/1000))K.fastq;
  gzip samples/${prefix}_$((size/1000))K.fastq
done;

echo "Assembling with canu...";

genomeSize=$(tail -n +2 ${reference} | tr -d "\n" | wc -m);

mkdir canu;
cd canu;

for size in "${samples[@]}";
do
  canu -d "canu$((size/1000))K" -p "${prefix}_$((size/1000))K" genomeSize=${genomeSize} -nanopore ../samples/${prefix}_$((size/1000))K.fastq.gz;
  FastaSeqStats canu$((size/1000))K/${prefix}_$((size/1000))K.contigs.fasta > canu_$((size/1000))K_stats.txt;
done;


