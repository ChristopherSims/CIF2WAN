#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --job-name=CS_W90
#SBATCH --mail-user=christopherssims95@gmail.com
#SBATCH --mail-type=ALL


#Select how logs get stored
mkdir $SLURM_JOB_ID
export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"
export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"


#Load Modules
ml load espresso


cd $SLURM_SUBMIT_DIR
# Create Log File
echo $SLURM_SUBMIT_DIRecho "JobID: $SLURM_JOB_ID" >> $debug_logs
echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs
echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs
echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs
echo  "Current working directory is `pwd`" >> $debug_logs


# Module debugging
module list >> $debug_logs
which mpirun >> $debug_logs


#Start Timestamp
date >> $benchmark_logs
echo "ulimit -l: " >> $benchmark_logs
ulimit -l >> $benchmark_logs


# run file
wannier90.x CS
sleep 3 
grep CONV CS.wout >> wconv.txt
echo "Program is finished with exit code $? at: `date`"


#End Timestamp
date >> $benchmark_logs
echo "ulimit -l" >> $benchmark_logs
ulimit -l >> $benchmark_logs


#Cleanup
mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/
mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/
rm -rf $SLURM_JOB_ID
