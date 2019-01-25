#!/bin/bash

# Run when first logged into BW

# BLUE WATERS SETUP ######

cd /scratch/sciteam/$USER
mkdir -p radical.pilot.sandbox
cd radical.pilot.sandbox
wget https://raw.githubusercontent.com/radical-cybertools/radical.pilot/devel/bin/radical-pilot-create-static-ve
sh ./radical-pilot-create-static-ve ve.ncsa.bw_aprun.0.50.7 bw

# Installing Gromacs on Blue Waters

# https://bluewaters.ncsa.illinois.edu/gromacs/

export GROMACS=$HOME/gromacs
export VERSION=gromacs-5.1.1

mkdir $GROMACS
cd $GROMACS
wget ftp://ftp.gromacs.org/pub/gromacs/${VERSION}.tar.gz
tar zxvf ${VERSION}.tar.gz
cd $VERSION

module swap PrgEnv-cray PrgEnv-gnu
module add fftw
module add cmake
module add boost

export CRAYPE_LINK_TYPE=dynamic
export CRAY_ADD_RPATH=yes
export CXX=CC
export CC=cc
export CMAKE_PREFIX_PATH=$FFTW_DIR/../
export FLAGS="-dynamic -O3 -march=bdver1 -ftree-vectorize -ffast-math -funroll-loops"

export INSTALL=$GROMACS/$VERSION/build-cpu
mkdir $INSTALL
cd $INSTALL

cmake ../ -DGMX_MPI=ON -DGMX_OPENMP=ON -DGMX_GPU=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON -DGMX_X11=OFF -DGMX_DOUBLE=OFF -DCMAKE_SKIP_RPATH=YES -DCMAKE_INSTALL_PREFIX=$INSTALL/.. -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_CXX_FLAGS="$FLAGS" -DGMX_CPU_ACCELERATION="AVX_128_FMA"

make -j4
make install

# Gromacs executable location $HOME/gromacs/gromacs-5.1.1/bin/gmx_mpi

# Installing CoCo-Md on BW



##########################
