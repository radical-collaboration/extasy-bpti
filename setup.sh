#!bin/bash

BW_LOGIN="ADD_BW_USERNAME_HERE"

# Creates a virtual env and sets up RP
# Assumes python is installed

sudo apt install virtualenv
virtualenv --system-site-packages $HOME/ve
source $HOME/ve/bin/activate
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

## Installing "globus-proxy-utils" "globus-simple-ca"		

sudo echo "deb http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" >> /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list		
sudo echo "deb-src http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" >> /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list		
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
