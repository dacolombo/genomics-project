# Genomics Project scripts
A collection of scripts used to perform a chloroplast genome assembly for a project of the Genomics course:

# AssemblyIllumina
This script performs trimming, mapping to a reference and assembly of the mapped reads. The reads used as input must be Illumina paired end reads.\
The script will use:
- **fastq-mcf** to perform adapters trimming and quality filtering;
- **bowtie2** to map the reads on the reference;
- **ABySS** to perform the assembly of the mapped reads.

## Arguments
The help of the script can be shown with `AssemblyIllumina.sh -h`.

The inputs needed to make the script function properly are:
- `-1`: file containing the forward Illumina reads. It has to be in either .fastq or .fastq.gz format;
- `-2`: file containing the reverse Illumina reads. It has to be in either .fastq or .fastq.gz format;
- `-a`: file containing the adapters possibly present in the reads. It has to be in .fasta format;
- `-r`: file containing the reference genome. It has to be in .fasta format;
- `-t`: integer indicating the number of threads to be used.

Right before the assembly, different k values are used to then ask the user what is the optimal one. Predefined values are inserted in the script in the variable `k_values` based on a specific project. To change them, modify *line 34* accordingly.
Sampling of the mapped reads is performed in order to obtain the optimal assembly. Predefined values are inserted in the script in the variable `samples` based on a specific project. To change them, modify *line 35* accordingly.


## Important Notes
More or less at half of its running time, the script will output a table showing the stats to decide the best k parameter to use. The user needs at this point to input the integer number corresponding to the best k and then press enter.


## Dependencies
This script uses these command line tools:
- gzip
- fastq-mcf
- bowtie2
- SAMTools
- BEDTools
- ABySS
- seqtk



# AssemblyOxN
This script performs trimming, mapping to a reference and assembly of the mapped reads. The reads used as input must be Oxford Nanopore reads.\
The script will use:
- **porechop** to perform adapters trimming and quality filtering;
- **minimap2** to map the reads on the reference;
- **canu** to perform the assembly of the mapped reads.

## Arguments
The help of the script can be shown with `AssemblyOxN.sh -h`.

The inputs needed to make the script function properly are:
- `-i`: file containing the Oxford Nanopore reads. It has to be in either .fastq or .fastq.gz format;
- `-r`: file containing the reference genome. It has to be in .fasta format;
- `-t`: integer indicating the number of threads to be used.

Sampling of the mapped reads is performed in order to obtain the optimal assembly. Predefined values are inserted in the script in the variable `samples` based on a specific project. To change them, modify *line 32* accordingly.

## Dependencies
This script uses these command line tools:
- gzip
- porechop
- minimap2
- SAMTools
- BEDTools
- seqtk
- canu
- FastaSeqStats


# SingleFastaTrim
A script to trim a sequence's end in a single sequence FASTA file based on a given position.\
This script was used to remove duplicated regions in the ends of the chloroplast assembly, which are quite common when assemblying circular genomes.
It uses standard bash and sed.

## Arguments
The help of the script can be shown with `Arguments.sh -h`.

The inputs needed to make the script function properly are:
- `-i`: file containing the assembly sequence. It has to be in .fasta format and have only one sequence in it.
- `-p`: integer indicating the position in the sequence to perform the trimming.
- `-s`: the side to trim away. Arguments accepted are "right", "left", "r", "l".

The count for the position in the sequence starts from 1 for both sides. The base at the position indicated is trimmed aswell.

Note: the script makes use of the /tmp folder in order to write temporary files.

