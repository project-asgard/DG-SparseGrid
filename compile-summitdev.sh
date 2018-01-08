#!/usr/bin/env bash

module unload xalt

mkdir build
cd build
cmake ..
make