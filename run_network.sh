#!/usr/bin/env bash
# Executes three SLURM jobs, one for each step of network construction
# User should edit sbatch scripts for <> options
# author:   Kevin Chau
# date:     2018 02 06


# soft thresholding
RUN_SFT="runsft.sh"

# sbatch options
echo "#!/usr/bin/env bash\n" > $RUN_SFT
cat >> $RUN_SCRIPT <<- EOM
#SBATCH --job-name=net_construct
#SBATCH --output=net_con.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -t 10:00:00
#SBATCH -p shared

module load R
./1.sft.R <infile.txt> <sftimg.png (opt)> <outfile.RData (opt)>
EOM

SFT_JOB_ID=$(sbatch runsft.sh | awk '{ print $4 }')


# tom construction
RUN_TOM="runtom.sh"

# sbatch options
echo "#!/usr/bin/env bash\n" > $RUN_TOM
cat >> $RUN_TOM <<- EOM
#SBATCH --job-name=tom
#SBATCH --output=tom.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -t 48:00:00
#SBATCH -p large-shared
#SBATCH --mem=200G

module load R
./2.TOM.R <infile.txt> <sft.RData (opt)> <out.RData (opt)>
EOM

TOM_JOB_ID=$(sbatch --dependency=afterok:$SFT_JOB_ID | awk '{ print $4 }')


# network construction
RUN_WGCNA="runwgcna.sh"

# sbatch options
echo "#!/usr/bin/env bash\n" > $RUN_TOM
cat >> $RUN_WGCNA <<- EOM
#SBATCH --job-name=wgcna
#SBATCH --output=wgcna.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 48:00:00
#SBATCH -p compute
#SBATCH --mem=200G

module load R
./3.SingleBlockWGCNA.R <TOM.RData> <data.txt> <network.RData (opt)> <outdir (opt)>
EOM

WGCNA_JOB_ID=$(sbatch --dependency=afterok:$TOM_JOB_ID | awk '{ print $4 }')
