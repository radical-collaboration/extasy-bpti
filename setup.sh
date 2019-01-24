#!bin/bash


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

# Installing CoCo-Md on BW



##########################



BW_LOGIN="ADD_BW_USERNAME_HERE"

# Creates a virtual env and sets up RP
# Assumes python is installed

sudo apt install virtualenv
virtualenv --system-site-packages $HOME/shared/ve
source $HOME/shared/ve/bin/activate
pip install radical.pilot
pip install radical.entk
 
export RADICAL_PILOT_DBURL="mongodb://franklinb:Dg6w*l08@ds133816.mlab.com:33816/frank01" #Using Franklins MongoDB
export RADICAL_VERBOSE=DEBUG
export RADICAL_PILOT_VERBOSE=DEBUG
export RADICAL_PILOT_AGENT_VERBOSE=DEBUG
export RMQ_HOSTNAME='two.radical-project.org'
export RMQ_PORT=33211
export RADICAL_ENMD_PROFILING=1
 
# viveks install gsissh script for Ubuntu (Works on 16.04!)

# INSTALL GLOBUS FIRST!!!!

## Installing "globus-proxy-utils" "globus-simple-ca"		

echo "deb http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" | sudo tee -a /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list
echo "deb-src http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" | sudo tee -a /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list		
sudo apt-get update		
sudo apt-get install globus-proxy-utils globus-simple-ca

# Installing "ca-policy-egi-core" - 1

cd
mkdir tmp
cd tmp
wget http://www.globus.org/ftppub/gt6/installers/repo/globus-toolkit-repo%5flatest%5fall.deb
dpkg -i globus-toolkit-repo_latest_all.deb
sudo apt-get update
sudo apt-get install gsi-openssh-clients globus-data-management-client myproxy globus-proxy-utils globus-simple-ca ca-certificates ca-digicertassuredidrootca-root ca-digicertgridca-1-classic ca-digicertgridca-1g2-classic-2015 ca-digicertgridrootca-root ca-digicertgridtrustca-classic ca-digicertgridtrustcag2-classic globus-gsi-cert-utils-progs libglobus-gsi-cert-utils0 ssl-cert 

## Installing "ca-policy-egi-core" - 2

curl https://dist.eugridpma.info/distribution/igtf/current/GPG-KEY-EUGridPMA-RPM-3 | sudo apt-key add -
cd /etc/apt/sources.list.d
sudo touch egi-igtf-core.list
sudo echo "deb http://repository.egi.eu/sw/production/cas/1/current egi-igtf core" >> egi-igtf-core.list
sudo apt-get update
sudo apt-get install ca-policy-egi-core

# Obtain Blue Waters Certificate

myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 -l $BW_LOGIN
