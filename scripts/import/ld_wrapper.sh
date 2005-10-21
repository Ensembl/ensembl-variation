#!/bin/sh

BASE_DIR=/ecs4/scratch5/dani/ld_data

# Get the first parameter number from the command line.
firstfile=$1
secondfile=$2

# Copy the input files locally
lsrcp ecs4a:${firstfile} /tmp/file1.$LSB_JOBID


# Run the job and save the exit status
   /nfs/acari/dr2/projects/src/branch-normalized-alleles/ensembl-variation/scripts/import/new_calc_genotypes \
	< /tmp/file1.$LSB_JOBID > \
	/tmp/output.$LSB_JOBID
    status=$?
    
# Copy the output files back, if it worked
    if [ $status -eq 0 ]; then
	lsrcp /tmp/output.$LSB_JOBID ecs4a:${secondfile}
    fi
    
# Delete the temporary files
    rm /tmp/file1.$LSB_JOBID \
	/tmp/output.$LSB_JOBID
    
# Finally return the job's real exit status to LSF
exit $status
