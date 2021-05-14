#!/usr/bin/env bash

set -euxo pipefail;


########################################
###########     Help     ###############
########################################
Help()
{
  echo "Description:"
  echo "  Pipeline to trim, align to a reference and assemble Illumina reads."
  echo
  echo "Syntax:"
  echo "  AssemblyIllumina.sh -1 <input_fastq_file_1> -2 <input_fastq_file_2> -a <illumina_adapters_file> -r <reference_fasta_file> -t <number_of_threads>"
  echo
  echo "Options:"
  echo "  -1 <input_fastq_file_1>         input fastq file containing the forward Illumina reads (can be in .fastq or .fastq.gz format)"
  echo "  -2 <input_fastq_file_2>         input fastq file containing the reverse Illumina reads (can be in .fastq or .fastq.gz format)"
  echo "  -a <illumina_adapters_file>     input fasta file containing the illumina adapters to be trimmed from the reads"
  echo "  -r <reference_fasta_file>       fasta file containing the reference genome"
  echo "  -t <number_of_threads>          number of threads to use"

}



########################################
###########     Main     ###############
########################################



k_values=$(seq 17 8 99);
samples=(25000 55000 11000 192500 275000 410000 550000);



#-------------#
# Get options #
#-------------#

while getopts ":1:2:a:r:t:h" opt; do
  case $opt in
    1) reads1="$OPTARG"
    ;;
    2) reads2="$OPTARG"
    ;;
    a) adapters="$OPTARG"
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




#--------------#
# Inputs check #
#--------------#

# Print help and exit if no option is passed
if [ "$#" -eq 0 ];
then
  Help;
  exit;
fi;


# Print error message and exit if the reads are in a wrong format
if { [ "${reads1#*.}" != "fastq" ] && [ "${reads1#*.}" != "fastq.gz" ]; } || { [ "${reads2#*.}" != "fastq" ] && [ "${reads2#*.}" != "fastq.gz" ]; };
then
    echo "Wrong format of input reads";
    exit;
fi;


# Print error message and exit if the reference is in a wrong format
if [ "${reference##*.}" != "fasta" ] && [ "${reference##*.}" != "fa" ];
then
    echo "reference is not in fasta format";
    exit;
fi;



# Save prefix to be used to call files
prefix=${reads1%%_*};
prefix_reference=${reference%%.*};




#---------------------#
# Preprocessing reads #
#---------------------#

# Compress reads if not compressed
if [ "${reads1##*.}" != "gz" ];
then
    echo "Compressing reads...";
    gzip ${reads1};
    reads1=${reads1}.gz;

elif [ "${reads2##*.}" != "gz" ];
then
    echo "Compressing reads...";
    gzip ${reads2};
    reads2=${reads2}.gz;
fi;


# Trim adapters from reads
echo "Trimming reads...";
fastq-mcf -o ${prefix}_trimmed_1.fastq.gz -o ${prefix}_trimmed_2.fastq.gz ${adapters} ${reads1} ${reads2};




#---------#
# Mapping #
#---------#

# Build reference index
echo "Building reference index...";
bowtie2-build --threads ${threads} ${reference} ${prefix_reference};

# Align
echo "Mapping reads...";
bowtie2 -p ${threads} -x ${prefix_reference} -1 ${prefix}_trimmed_1.fastq.gz -2 ${prefix}_trimmed_2.fastq.gz | samtools view -@ ${threads} -b -o ${prefix}.bam -;

# Sort bam file by position
echo "Sorting reads...";
samtools sort -@ ${threads} -o ${prefix}.sort.bam ${prefix}.bam;
rm ${prefix}.bam;

# Reads extraction
echo "Extracting mapped reads...";
samtools view -@ ${threads} -b -F 4 -o ${prefix}_mapped.sort.bam ${prefix}.sort.bam;

# Convert mapped reads into fastq
echo "Converting mapped reads into fastq file...";
samtools fastq -@ ${threads} -1 ${prefix}_mapped_1.fastq.gz -2 ${prefix}_mapped_2.fastq.gz -c 6 -N ${prefix}_mapped.sort.bam;

# Compute the coverage
echo "Computing coverage...";
samtools faidx ${reference};
bedtools genomecov -pc -d -ibam ${prefix}_mapped.sort.bam -g ${prefix_reference}.fai > ${prefix}_mapped.covbed.txt;




#-------------------#
# Best k estimation #
#-------------------#

# Perform assembly with different k
for k in ${k_values};
do
    mkdir k$k;
    abyss-pe -C k$k name=${prefix} k=$k in="../${prefix}_mapped_1.fastq.gz ../${prefix}_mapped_2.fastq.gz";
done;

# Print table with results from all the k
echo "These are the assembly stats resulting from different k:";
/usr/lib/abyss/abyss-fac k*/*unitigs.fa;

printf "Input the best k value: ";
read k;




#----------#
# Sampling #
#----------#
echo "Preforming the sampling...";

mkdir samples;

for size in "${samples[@]}";
do
    size_name="$((size/1000))K";
    seqtk sample -s1000 ${prefix}_mapped_1.fastq.gz ${size} > samples/${prefix}_${size_name}_1.fastq;
    gzip samples/${prefix}_${size_name}_1.fastq;
    seqtk sample -s1000 ${prefix}_mapped_2.fastq.gz ${size} > samples/${prefix}_${size_name}_2.fastq;
    gzip samples/${prefix}_${size_name}_2.fastq;
done;




#----------#
# Assembly #
#----------#

mkdir Abyss;
cd Abyss;

echo "Assembling with ABySS...";
for size in "${samples[@]}";
do
    size_name="$((size/1000))K";
    mkdir ${size_name};
    cd ${size_name};
    abyss-pe name="${prefix}_${size_name}" k=$k in="../../samples/${prefix}_${size_name}_1.fastq.gz ../../samples/${prefix}_${size_name}_2.fastq.gz";
    cd ..;
done;




#---------#
# Results #
#---------#

/usr/lib/abyss/abyss-fac *K/*scaffolds.fa > assemblies_stats.txt;




#------#
# Exit #
#------#
cd ..;

