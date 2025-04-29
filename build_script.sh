#!/bin/sh

# Script for building the example

set -x
set -e

module load rocm/6.3.1
module load PrgEnv-gnu
module load craype-accel-amd-gfx90a
ml cmake

export CC=cc
export CXX=CC

export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export GPU_TARGET=gfx908
export OMPI_CC=hipcc

# Setup MGARD installation dir

install_dir=/lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/MGARD/install-hip-frontier/
#install_dir=/lustre/orion/proj-shared/csc143/jieyang/MGARD/install-hip-frontier/
#zstd_dir=/lustre/orion/proj-shared/csc143/jieyang/MGARD/install-hip-frontier/lib/cmake/zstd/
#install_dir=/lustre/orion/cfd164/proj-shared/gongq/Software/MGARD/install-hip-frontier/
#install_dir=/ccs/home/gongq/indir/frontier/

rm -f build/CMakeCache.txt
mkdir -p build 
cmake -S .  -B ./build \
	    -Dmgard_ROOT=${install_dir} \
        -DCMAKE_HIP_ARCHITECTURES="gfx90a" \
        -DIOCOMP_ENABLE_HIP=ON\
        -DCMAKE_PREFIX_PATH="${install_dir}"\
        -DCMAKE_C_COMPILER=hipcc \
        -DCMAKE_CXX_COMPILER=hipcc 
	  
cmake --build ./build
