#!/bin/sh

# Copyright 2021, Oak Ridge National Laboratory.
# MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
# Author: Jieyang Chen (chenj3@ornl.gov)
# Date: April 2, 2021
# Script for building the example

set -x
set -e

module load rocm/5.3.0
module load PrgEnv-gnu/8.3.3
module load craype-accel-amd-gfx90a
ml cmake

export CC=cc
export CXX=CC

export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export GPU_TARGET=gfx908
export OMPI_CC=hipcc

# Setup MGARD installation dir


install_dir=/ccs/home/gongq/indir/frontier/
#install_dir=/lustre/orion/csc303/proj-shared/jieyang/MGARD/install-hip-frontier

rm -f build/CMakeCache.txt
mkdir -p build 
cmake -S .  -B ./build \
	    -Dmgard_ROOT=${install_dir} \
        -DIOCOMP_ENABLE_HIP=ON\
        -DCMAKE_PREFIX_PATH="${install_dir}"\
        -DCMAKE_C_COMPILER=hipcc \
        -DCMAKE_CXX_COMPILER=hipcc 
	  
cmake --build ./build
