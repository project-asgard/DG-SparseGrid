cmake_minimum_required(VERSION 3.9)
project(AMGX_test)

set(CMAKE_CXX_STANDARD 14)

add_executable(AMGX_test main.cpp)

# Add cusparse libraries
find_package(CUDA REQUIRED)
target_link_libraries(AMGX_test ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY})

# TODO - make AMGX findable
include_directories(${AMGX_INCLUDE_DIR})
target_link_libraries(AMGX_test ${AMGX_LIB_DIR}/libamgx.a)
