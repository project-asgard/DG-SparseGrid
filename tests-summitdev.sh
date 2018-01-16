#!/usr/bin/env bash

module unload xalt
module load cmake

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
cat << EOF > jsrun_command
#!/bin/bash -ex
cd ${RUN_DIR}
jsrun -n1 ./DG_SparseGrid
EOF
chmod +x ./jsrun_command

# Run through tests
bsub -P STF007SUMMITDEV -W 30 -nnodes 1 -I $(pwd)/jsrun_command