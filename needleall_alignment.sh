#!/bin/bash
##SBATCH --account=congle
#SBATCH --partition=interactive
#SBATCH --time=48:00:00
#SBATCH --job-name="TATSalignmentR500_0727"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
##SBATCH --array=1
#SBATCH -o CRcount-%A-%a.out
#SBATCH -e CRcount-%A-%a.err
module purge
module load emboss
ref_path=/labs/congle/PRT/Nextseq_20191110/Data_Processing/ref_pre
echo '$1 = ' $1
base=/home/wjtang93
ls *.fa > $base"/"$1".txt"
mkdir $1
mv *.fa $1
mv $1".txt" $1
cd ${1}
for fq in $(cat $1".txt")
do
name2=$(echo "$fq" | cut -f 1 -d '.')
needleall $ref_path"/"$fq $base"/"$1"/"$fq -gapopen 13 -gapextend 0.5 -aformat3 sam -outfile $name2".needleall"
done
rm *.fa
