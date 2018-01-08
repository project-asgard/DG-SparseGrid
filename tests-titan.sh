#!/usr/bin/env bash

module unload xalt
module load cmake3

mkdir build
cd build
cmake ..
make

# Save the current build directory and move to lustre
BUILD_DIR=$(pwd)
RUN_DIR="/lustre/atlas/scratch/atj/stf007/DGSG-${CI_JOB_ID}"
mkdir -p ${RUN_DIR}
cd ${RUN_DIR}

# Copy the executable to the lustre
cp ${BUILD_DIR}/DG_SparseGrid ${RUN_DIR}

# Create the aprun command to run the tests
cat << EOF > aprun_command
#!/bin/bash -ex
cd ${RUN_DIR}
aprun -n 1 ./DG_SparseGrid
EOF
chmod +x ./aprun_command

# Run through tests
qsub -ASTF007 -lwalltime=00:05:00,nodes=1 -I -x ./aprun_command