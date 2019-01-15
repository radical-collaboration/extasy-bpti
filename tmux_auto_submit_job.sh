#!/bin/bash 


#p12b02_left_d3_itrs48_k12_100_k34_100

for i in gmx*; 
do 
	sess="$i" 
	dir="$i"
	empty=" "
	window=${sess}:0
	tmux new-session -d -s  $sess  
	tmux send -t $sess.0 "cd $dir/" Enter
	#tmux send -t $sess.0 " " Enter
	tmux send -t $sess.0 "export RADICAL_VERBOSE=DEBUG" Enter
	tmux send -t $sess.0 "export RADICAL_PILOT_VERBOSE=DEBUG" Enter
	tmux send -t $sess.0 "export RADICAL_PILOT_AGENT_VERBOSE=DEBUG" Enter
	tmux send -t $sess.0 "source ~/Documents/feature_entk-0.7/myenv/bin/activate"  Enter
	tmux send -t $sess.0 "export RMQ_HOSTNAME='two.radical-project.org'" Enter
	tmux send -t $sess.0 "export RMQ_PORT=33211" Enter
	tmux send -t $sess.0 "export RADICAL_PILOT_DBURL='mongodb://franklinb:Dg6w*l08@ds151530.mlab.com:51530/frank02'" Enter
	tmux send -t $sess.0 "export RADICAL_ENMD_PROFILING=1" Enter
	tmux send -t $sess.0 "python runme.py --RPconfig bluewaters.rcfg --Kconfig gmxcoco.wcfg " Enter 
done 
#tmux send -t $sess.0 "echo $\RADICAL_ENMD_PROFILING" Enter
#tmux send -t $sess.0 




