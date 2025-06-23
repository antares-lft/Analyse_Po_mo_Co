
#!/bin/bash

source /beegfs/data/gueguen/monpy/bin/activate
source /beegfs/home/alafitte/.bash_aliases

dir=/beegfs/home/alafitte/Pop_size/Internship/Script/Population_size/N100
outdir=/beegfs/data/alafitte/Pop_size/Internship/Results/Population_size/N100

mkdir -p ${outdir}

for num in {1..30}
do
mkdir -p ${outdir}/sim_${num}
echo "#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --partition=normal
#SBATCH --error=AA_Landscape_${num}.err
#SBATCH --output=AA_Landscape_${num}.out

set -e
set -u
set -o pipefail

hostname

echo "source $HOME/.bashrc\n" > ${outdir}/sim_${num}/slim_${num}.sh

python3  ${dir}/slim_parse.py -S ${dir}/Population_Splitting_Codon_Template.slim -N ${num} -n 100 -r 0.5 -m 1e-6 -l 1e3 -t ${dir}/bat_genes_complete_filtered-PhyML_tree.nhx -b 30000 -f ${dir}/AA_Fitness_Landscape" > ${outdir}/sim_${num}/slim_${num}.sh

cd ${outdir}/sim_${num}
sbatch slim_${num}.sh
cd ../..

done
