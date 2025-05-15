for i in chr1H chr2H chr3H chr4H chr5H chr6H chr7H;
do
echo $i

for m in {0..70};
do
echo $m

echo "#!/bin/bash

#SBATCH --job-name=${i}.${m}                              # Job name, will show up in squeue output
#SBATCH --auks=yes                                        # use Kerberos
#SBATCH --partition=cpu                                   # possible values: cpu, gpu
#SBATCH --ntasks=1                                        # Number of tasks (default=1)
#SBATCH --cpus-per-task=4                                 # number of cpus for this task (default=1)
#SBATCH --mem=80G                                         # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                         # Ensure that all cores are on one machine
#SBATCH --output=${i}.${m}.out                            # File to which standard out will be written
#SBATCH --error=${i}.${m}.err                             # File to which standard err will be written
#SBATCH --exclude=slurm-03,slurm-20
date

/opt/Bio/bcftools/1.15.1/bin/bcftools mpileup --threads 4 -a AD,DP,DV -b bamlist.txt -q 20 -Q 20 -d 1200000000 --ns 3332 -r ${i}:$(($m*10000000+1))-$(($m*10000000+10000000)) -Ou -f /filer-dg/agruppen/dg5/gaog/cultivars/MorexV3/200416_MorexV3_pseudomolecules.fasta 2> ${i}.${m}_mpileup.out | /opt/Bio/bcftools/1.15.1/bin/bcftools call --threads 4 -f GQ,GP -mv -Oz -o ${i}.${m}_raw.vcf.gz > ${i}.${m}_mpileup.out 2>>${i}.${m}_mpileup.out;echo finish >>${i}.${m}_mpileup.out

date"  > ${i}.${m}.zsh

sbatch ${i}.${m}.zsh

cd /filer-dg/agruppen/dg4/gaog/rawvcf4 # 工作目录

done 
done
