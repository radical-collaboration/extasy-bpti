#!/bin/bash


BW_LOGIN="suess" # Replace with BW username

## Installing "globus-proxy-utils" "globus-simple-ca"           

echo "deb http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" | sudo tee -a /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list
echo "deb-src http://toolkit.globus.org/ftppub/gt6/stable/deb xenial contrib" | sudo tee -a /etc/apt/sources.list.d/globus-toolkit-6-stable-xenial.list
sudo apt-get update
sudo apt-get install globus-proxy-utils globus-simple-ca

## Installing "ca-policy-egi-core" - 2

sudo apt install curl
curl https://dist.eugridpma.info/distribution/igtf/current/GPG-KEY-EUGridPMA-RPM-3 | sudo apt-key add -
cd /etc/apt/sources.list.d
sudo touch egi-igtf-core.list
echo "deb http://repository.egi.eu/sw/production/cas/1/current egi-igtf core" | sudo tee -a egi-igtf-core.list
sudo apt-get update
sudo apt-get install ca-policy-egi-core

# Installing "ca-policy-egi-core" - 1

cd
mkdir tmp
cd tmp
wget http://www.globus.org/ftppub/gt6/installers/repo/globus-toolkit-repo%5flatest%5fall.deb
sudo dpkg -i globus-toolkit-repo_latest_all.deb
sudo apt-get update
sudo apt-get install gsi-openssh-clients globus-data-management-client myproxy globus-proxy-utils globus-simple-ca ca-certificates ca-digicertassuredidrootca-root ca-digicertgridca-1-classic ca-digicertgridca-1g2-classic-2015 ca-digicertgridrootca-root ca-digicertgridtrustca-classic ca-digicertgridtrustcag2-classic globus-gsi-cert-utils-progs libglobus-gsi-cert-utils0 ssl-cert

# Obtain Blue Waters Certificate

myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 -l $BW_LOGIN
