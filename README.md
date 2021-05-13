# Genomics Project scripts
A collection of scripts used to perform a chloroplast genome assembly for a project of the Genomics course:
- `SingleFastaTrim.sh`: a script to trim a sequence in a single sequence FASTA file based on a given position. This was used to remove duplicated regions in the ends of the assembly, which are quite common when assemblying circular genomes.
- `assembly_pipeline_oxn.sh`: a pipeline for assembling Oxford Nanopore reads. It's based on *porechop*, *minimap2* and *canu*.

The usage of the scripts can be seen by using the `-h` option for all the commands.
