#!/usr/bin/env bash

set -e

########################################
###########     Help     ###############
########################################
Help()
{
  echo "Description:"
  echo "Script to trim a sequence in a single-FASTA file given a position and the side to remove"
  echo
  echo "Syntax:"
  echo "SingleFastaTrim -i <input_fasta_file> -p <position> -s <side_to_trim>"
  echo
  echo "Options:"
  echo "-i <input_fasta_file>		input fasta file (must be in .fasta format)"
  echo "-p <position>                   position to trim the sequence. The count starts from 1 for both sides. The base at <position> is trimmed aswell"
  echo "-s <side_to_trim>               side to trim. Accepted values are 'right' 'left' 'r' 'l'"

}




########################################
###########     Main     ###############
########################################

#-------------#
# Get options #
#-------------#

while getopts ":i:p:s:h" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    p) position="$OPTARG"
    ;;
    s) side="$OPTARG"
    ;;
    h) Help
       exit
    ;;
    \?) echo "Invalid option -$OPTARG"
        exit
    ;;
  esac
done

prefix=${input%%.*}


#--------------#
# Main program #
#--------------#

header_count=$(grep -c ">" ${input})
file_format=${input##*.}

if [ $header_count -eq 1 ] && { [ $file_format = "fasta" ] || [ $file_format = "fa" ]; }
then

  if [ $side = "right" ] || [ $side = "r" ]
  then

    # trim to the right
    sequence=$(tail -n +2 ${input} | tr -d '\n')
    echo ${sequence::${position}} > /tmp/trimmed.txt

  elif [ $side = "left" ] || [ $side = "l" ]
  then

    # trim to the left
    sequence=$(tail -n +2 ${input} | tr -d '\n')
    echo ${sequence:$((position - 1))} > /tmp/trimmed.txt

  fi

else
  echo "Input file is not a single-FASTA file"
  exit

fi


# Add new line every 100 bases in the sequence
sed -e "s/.\{100\}/&\n/g" < /tmp/trimmed.txt > /tmp/trimmed_newlines.txt


#----------------------#
# Write result to file #
#----------------------#

# write header
head -n 1 ${input} > ${prefix}_trimmed.fasta
# write sequence
cat /tmp/trimmed_newlines.txt >> ${prefix}_trimmed.fasta


#-----------#
# Clean tmp #
#-----------#

rm /tmp/trimmed.txt
rm /tmp/trimmed_newlines.txt
