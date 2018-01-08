#!/usr/bin/env bash

module unload xalt
module load cmake3

mkdir build
cd build
cmake ..
make

# Save the current build directory and move to lustre
BUILD_DIR=$(pwd)
cd /lustre/atlas/scratch/atj/stf007

# Assign a unique name to the executable based upon the GitLab CI job ID
UNIQUE_EXE_NAME="DG_SparseGrid-${CI_JOB_ID}"

# Copy the executable to the lustre
cp ${BUILD_DIR}/DG_SparseGrid ${UNIQUE_EXE_NAME}

# Create the aprun command to run the tests
cat << EOF > aprun_command
aprun -n 1 ${UNIQUE_EXE_NAME}
EOF
chmod +x ./aprun_command

# Run through tests
qsub -ASTF007 -lwalltime=00:05:00,nodes=1 -I -x ./aprun_command