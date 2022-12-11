# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


for i in `seq 0  1 99`; do 
     Rscript run.R $i upstream  >/mnt/extStorage/constantin/windows/upstream/$i/alphabeta_output.txt &
     Rscript run.R $i gene  > /mnt/extStorage/constantin/windows/gene/$i/alphabeta_output.txt &
     Rscript run.R $i downstream > /mnt/extStorage/constantin/windows/downstream/$i/alphabeta_output.txt &
done 