#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#scriptdir="/bubo/home/h10/jessada/private/combiner/"

#define default values
INPUT_FILE_DEFAULT=$scriptdir/tmp/condel_input
TASK_LIST_FILE_DEFAULT=$scriptdir/tmp/task_list

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i FILE    specify path of input file
-t FILE    specify path of task list file
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#parse param
while getopts "i:t:" OPTION; do
  case "$OPTION" in
    i)
      input_file="$OPTARG"
      ;;
    t)
      task_list_file="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

#setting default values:
: ${input_file=$INPUT_FILE_DEFAULT}
: ${task_list_file=$TASK_LIST_FILE_DEFAULT}

#display parameter
cat <<EOF
configuration:
input file       : $input_file
task list file   : $task_list_file
EOF

function split_and_submit_to_condel {
    tmp_condel_file_name=$1"_"`printf "%05d" "$2"`"_"`printf "%05d" "$3"`
    echo $tmp_condel_file_name
    sed -n $2","$3"p" $input_file > $tmp_condel_file_name
    curl -X PUT -T $tmp_condel_file_name http://bg.upf.edu/condel/taskService >> $4
}  

#clear old data
if [ -e $task_list_file ]
then
    rm $task_list_file
fi

#split condel input and then submit it to condel
file_size="$( wc -l $input_file | awk '{print $1}')"
begin_position=1
end_position=300
round=1
echo $file_size
while :
do
    echo "round" $round
    if [ $end_position -ge $file_size ]
    then
        split_and_submit_to_condel $input_file $begin_position $file_size $task_list_file
        break
    fi
    split_and_submit_to_condel $input_file $begin_position $end_position $task_list_file
    
    round=$(( round + 1 ))
    begin_position=$(( begin_position + 300 ))
    end_position=$(( end_position + 300 ))
done

