# Build and install AMGX
```
export BASE_DIR=`pwd`
module switch PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module unload xalt
module load cmake3
git clone https://github.com/NVIDIA/AMGX.git
cd AMGX
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$BASE_DIR/AMGX_INSTALL ..
```

# Build test application
```
cd $BASE_DIR
mkdir build
cd build
cmake -DAMGX_LIB_DIR=$BASE_DIR/AMGX_INSTALL/lib -DAMGX_INCLUDE_DIR=$BASE_DIR/AMGX_INSTALL/include ..
```
