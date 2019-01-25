#!/bin/bash

# Run on Local (Xbow) Machine

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
