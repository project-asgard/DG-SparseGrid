#!/usr/bin/env bash

module unload xalt
module load cmake3

mkdir build
cd build
cmake ..
make
