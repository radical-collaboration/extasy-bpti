
# tmux_auto_submit_job.sh will run any extasy script you specifiy  by making the approriate edits to it. 

# cat tmux_auto_submit_job.sh
#!/bin/bash 


for i in gmx*;
do
sess="$i"
dir="$i"
empty=" "
window=${sess}:0
tmux new-session -d -s  $sess
tmux send -t $sess.0 "cd $dir/" Enter
tmux send -t $sess.0 "export RADICAL_VERBOSE=DEBUG" Enter
tmux send -t $sess.0 "export RADICAL_PILOT_VERBOSE=DEBUG" Enter
tmux send -t $sess.0 "export RADICAL_PILOT_AGENT_VERBOSE=DEBUG" Enter
tmux send -t $sess.0 "source ~/Documents/feature_entk-0.7/myenv/bin/activate"  Enter
tmux send -t $sess.0 "export RMQ_HOSTNAME='two.radical-project.org'" Enter
tmux send -t $sess.0 "export RMQ_PORT=33211" Enter
# make sure to insert the link for your mongoDB information.
#you can get the link by logging into you mongoDB account and going to that DB instance
tmux send -t $sess.0 "export RADICAL_PILOT_DBURL='mongodb://<username>:<password>@ds151530.mlab.com:51530/frank02'" Enter
tmux send -t $sess.0 "export RADICAL_ENMD_PROFILING=1" Enter
tmux send -t $sess.0 "python runme.py --RPconfig bluewaters.rcfg --Kconfig gmxcoco.wcfg " Enter
done
# end script 

# this script needs to be run from the same dir as the dir containing the extasy scripts you wish to run. 
# for example given the following

hal9000@js-17-8:~/Documents/feature_entk-0.7/extasy-bpti$ ls
data          how_to_use_tmux_auto_submit_job.md  new_analysis  protein_pdb  remote_jupyter  Shaw_Data_Analysis
gmxcoco-bpti  meetings                            Papers        README.md    setup           tmux_auto_submit_job.sh
hal9000@js-17-8:~/Documents/feature_entk-0.7/extasy-bpti$ 

# note that tmux_auto_submit_job.sh is in the same dir as gmxcoco-bpti which is the target dir to run. 
# submit the job as follows

hal9000@js-17-8:~/Documents/feature_entk-0.7/extasy-bpti$ bash tmux_auto_submit_job.sh

# this will submit the job through tmux. It will title the tmux session the same as the dir you ar running. 

hal9000@js-17-8:~/Documents/feature_entk-0.7/extasy-bpti$ tmux ls
gmxcoco-bpti: 1 windows (created Sun Jan 13 19:41:08 2019) [130x31]


# listed here are some very useful Tmux links 
