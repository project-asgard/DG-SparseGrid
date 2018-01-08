#!/usr/bin/env bash

module unload xalt
module load cmake

mkdir build
cd build
cmake ..
make
