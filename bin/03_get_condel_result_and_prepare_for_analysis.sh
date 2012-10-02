#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/"

#define default values
RESULT_FILE_DEFAULT=$scriptdir/tmp/condel_score
TASK_LIST_FILE_DEFAULT=$scriptdir/tmp/task_list
TESTING_DATA_FILE_DEFAULT=$scriptdir/../data/testing_data

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-r FILE    specify path of result file
-t FILE    specify path of task list file
-T FILE    specify path of testing data file
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#parse param
while getopts "r:t:T:" OPTION; do
  case "$OPTION" in
    r)
      result_file="$OPTARG"
      ;;
    t)
      task_list_file="$OPTARG"
      ;;
    T)
      testing_data_file="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

#setting default values:
: ${result_file=$RESULT_FILE_DEFAULT}
: ${task_list_file=$TASK_LIST_FILE_DEFAULT}
: ${testing_data_file=$TESTING_DATA_FILE_DEFAULT}

#display parameter
cat <<EOF
configuration:
result file       : $result_file
task list file    : $task_list_file
testing data file : $testing_data_file
EOF

################################## start coding section ####################################
tmp_file=$scriptdir/tmp/tmp_file
tmp_condel_result=$scriptdir/tmp/tmp_condel_result
parsed_condel_result_file=$scriptdir/tmp/parsed_condel_result

#clear old data
if [ -e $tmp_condel_result ]
then
    rm $tmp_condel_result
fi
if [ -e $tmp_file ]
then
    rm $tmp_file
fi
if [ -e $parsed_condel_result_file ]
then
    rm $parsed_condel_result_file
fi

pad_length=2
pad=$(printf '%0.1s' "0"{1..2})

while read line
do
    curl -X get $line | grep -v "^#" >> $tmp_condel_result
done < $task_list_file
    
nawk -F $'\t' '{printf("%s:%s:%s\n", $2, $3, $15)}' < $tmp_condel_result >> $tmp_file

grep "^[0-9]" $tmp_file | awk -F ':' '{printf("%d\t%d\t%s\t%.3f\t%02d|%012d|%s|%.3f\t%02d|%012d|%s\n", $1, $2, $3, $4, $1, $2, $3, $4, $1, $2, $3)}' | sort -k5 -r | uniq -f5 | sort -k5 >> $parsed_condel_result_file
grep -v "^[0-9]" $tmp_file | awk -F ':' '{printf("%s\t%d\t%s\t%.3f\t%s|%012d|%s|%.3f\t%s|%012d|%s\n", $1, $2, $3, $4, $1, $2, $3, $4, $1, $2, $3)}' | sort -k5 -r | uniq -f5 | sort -k5 >> $parsed_condel_result_file

join -t $'\t' -1 6 -2 1 -o 0,2.2,2.3,2.4,2.5,2.6,2.7,1.4,2.8 $parsed_condel_result_file $testing_data_file > $result_file

