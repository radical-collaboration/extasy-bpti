# Instructions to run EnTK 0.7 scripts

## Installation and environment setup

* Virtualenv creation

```
virtualenv $HOME/myenv
source $HOME/myenv/bin/activate
```

* EnTK Installation: You can directly install from pypi

```
pip install radical.entk
```

You can check the version with ```radical-stack```, this
should print 0.7.* (where * is a number indicating versions
of patches).

NOTE: If your target machine is Blue Waters, note the radical
pilot version printed.

* RabbitMQ

A docker instance is already setup on radical.two with
hostname 'two.radical-project.org' and port '33211' that
you can specify to your EnTK script using two environment
variables:

```
export RMQ_HOSTNAME='two.radical-project.org'
export RMQ_PORT=33211
```

If you want to create your own, please see [EnTK docs](https://radicalentk.readthedocs.io/en/latest/install.html#installing-rabbitmq).

* MongoDB

You can create a new mongodb instance on [mlab](https://mlab.com/).
Steps broadly include:

1. Create a new db instance with approriate geographic location,
instance size, etc.
2. Once created, select the db and add an user to it. Specify
a user name and password.
3. Copy the db url and specify it as environment variable

```
export RADICAL_PILOT_DBURL='<mongo url>'
```

* Setup on Blue Waters

There is an additional step to perform on Blue Waters alone.
Log into Blue Waters and traverse to ```/scratch/sciteam/$USER```.
Create a folder ```radical.pilot.sandbox``` if it does not exist.
Next, you will be creating a virtualenv in this location. You
can follow these commands:

```
cd /scratch/sciteam/$USER
mkdir -p radical.pilot.sandbox
cd radical.pilot.sandbox
wget https://raw.githubusercontent.com/radical-cybertools/radical.pilot/devel/bin/radical-pilot-create-static-ve
sh ./radical-pilot-create-static-ve ve.ncsa.bw_aprun.0.50.7 bw
```

Note: '0.50.7' above is the version of RADICAL Pilot. This should
match the version that was printed above when you ran 
```radical-stack```.

## Executing your script

You can use the same cmd as in the previous versions:

```
python runme.py --Kconfig gmxcoco.wcfg --RPconfig bluewaters.rcfg
```

The default verbosity should give you updates about each stage of the
pipeline. If you want to increase the verbosity, you can specify
```RADICAL_ENTK_VERBOSE``` to ```INFO```.

## Issue reporting

If there are any issues/questions, please create a ticket in this 
Github repository. After investigating the issue, I can move it to
EnTK/RP/SAGA.