#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/"

#define default values
TESTING_DATA_FILE_DEFAULT=$scriptdir/../data/testing_data
RAW_DATA_FILE_DEFAULT=$scriptdir/../data/raw_data
OUTPUT_FILE_DEFAULT=$scriptdir/tmp/condel_input

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-t FILE    specify path of testing data file
-r FILE    specify path of raw data file
-o FILE    specify path of output file
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#parse param
while getopts "t:r:o:" OPTION; do
  case "$OPTION" in
    t)
      testing_data_file="$OPTARG"
      ;;
    r)
      raw_data_file="$OPTARG"
      ;;
    o)
      output_file="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

#setting default values:
: ${testing_data_file=$TESTING_DATA_FILE_DEFAULT}
: ${raw_data_file=$RAW_DATA_FILE_DEFAULT}
: ${output_file=$OUTPUT_FILE_DEFAULT}

#display parameter
cat <<EOF
configuration:
testing data file : $testing_data_file
raw data file     : $raw_data_file
output file       : $output_file
EOF


join -t $'\t' -1 1 -2 7 -o 2.1,2.2,2.3,2.5 $testing_data_file $raw_data_file > $output_file
