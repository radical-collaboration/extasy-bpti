from radical.entk import Pipeline, Stage, Task, AppManager
import os
import argparse
import glob
import sys
import imp
import json
import traceback
# ------------------------------------------------------------------------------
# Set default verbosity
if os.environ.get('RADICAL_ENTK_VERBOSE') == None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'


## Global variables

# Description of how the RabbitMQ process is accessible
# No need to change/set any variables if you installed RabbitMQ has a system
# process. If you are running RabbitMQ under a docker container or another
# VM, set "RMQ_HOSTNAME" and "RMQ_PORT" in the session where you are running
# this script.
hostname = os.environ.get('RMQ_HOSTNAME', 'localhost')
port = int(os.environ.get('RMQ_PORT', 5672))

# Share directory with ALL data
#shareDir="$SHARED"
#shareDir="staging://" # $SHARED is place holder and is replaced at runtime by "staging://" 
#https://github.com/radical-cybertools/radical.entk/blob/master/src/radical/entk/execution_plugin/staging/placeholders.py#L25	
#shareDir="/work/fbettenc/radical.pilot.sandbox/p13b01_left_d3_k12_1000_k34_1000" 
# # note tried without / before work and failed. diff err for /work/.. than work/..
#shareDir="/work/fbettenc/radical.pilot.sandbox/rp.session.js-17-187.jetstream-cloud.org.hal9000.017508.0005-pilot.0000/staging_area"

#shareDir="/work/fbettenc/p14b01_pool/staging_area"
shareDir='$SHARED'

# References to previous iterations, effects starting iteration of current run
# prev_sim_last_iter_to_use=48
# iterMod=iteration+prev_sim_last_iter_to_use
prev_sim_last_iter_to_use=0


# Retrieve basename and extension
outbase, ext = None, None

# NOTE: All stages, tasks, filenames are indexed from 0.

def generate_pipeline(index, iterations, ensemble_size):
#def generate_pipeline(index): #, iterations, ensemble_size):

    # Create a Pipeline object
    p = Pipeline()
    p.name = 'extasy-pipeline-%s' %index

    # Determine total iterations of all runs, current+historic
    total_iterations = prev_sim_last_iter_to_use + iterations

    for iter_cnt in range(prev_sim_last_iter_to_use, total_iterations):

        if iter_cnt != 0:

            
            # STAGE/KERNEL 1
            # Create a Stage object
            s1 = Stage()
            s1.name = 'grompp-stage-%s' %iter_cnt

            for t_cnt in range(ensemble_size):

                # Create a Task object
                t = Task()
                t.name = 'grompp-task-%s'%t_cnt

                # Same pre_exec as in the kernel def files for the specific target resource
                t.pre_exec = [  "export PATH=$PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                                "export GROMACS_LIB=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64",
                                "export GROMACS_INC=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/include",
                                "export GROMACS_BIN=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                                "export GROMACS_DIR=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu",
                                "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64"]

                # Same executable as in the kernel def files for the specific target resource                                
                t.executable = ["gmx grompp"]

                # Same flags arguments as in the kernel def file, but secondary arguments from the main file
                t.arguments = [ '-f','{0}'.format(os.path.basename(Kconfig.grompp_1_mdp)),
                                '-c','{0}'.format(os.path.basename(Kconfig.restr_file)),
                                '-p','{0}'.format(os.path.basename(Kconfig.top_file)),
                                '-o','min-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]

                # Copy link_input_data from main file, edit 'iterMod'->iter_cnt and 
                # 'instance'->t_cnt appropriately
                t.link_input_data = [   shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_1_mdp)),
                                        shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file)),
                                        shareDir+'/{0}'.format(os.path.basename(Kconfig.restr_file)),
                                        shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_1_itp_file)),
                                        shareDir+'/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iter_cnt-1,t_cnt,ext)
                                    ]

                # Copy copy_output_data from main file, edit 'iterMod'->iter_cnt and 
                # 'instance'->t_cnt appropriately
                t.copy_output_data = [  'min-{0}_{1}.tpr > '.format(iter_cnt,t_cnt) + 
                                        shareDir+'/min-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]

                # Add task to our Stage
                s1.add_tasks(t)

            # Add Stage to our pipeline
            p.add_stages(s1)

            # STAGE/ KERNEL 2
            s2 = Stage()
            s2.name = 'mdrun-stage-%s'%iter_cnt

            for t_cnt in range(ensemble_size):

                #Create a Task Object
                t = Task()
                t.name = 's2-task-%s'%t_cnt

                # Same pre_exec as in k def files for specific target resources 
                t.pre_exec = ["export PATH=$PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/bin",
                              "export GROMACS_LIB=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/lib64",
                              "export GROMACS_INC=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/include",
                              "export GROMACS_BIN=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/bin",
                              "export GROMACS_DIR=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu",
                              "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/lib64"]
                t.executable = ["gmx_mpi mdrun"] 

                # Number of cores for non-mpi Gromacs

                #from kernel mdrun.py - arguments = ['-deffnm','{0}'.format(self.get_arg("--deffnm="))]
                #from nwexgmx_v002.py - k2_min_kernel.arguments = ["--deffnm=min-{0}_{1}".format(iterMod-1,instance-1)]
                #NOTE - check for error here
                t.arguments = ["-deffnm=min-{0}_{1}".format(iter_cnt,t_cnt)]

                #k2_min_kernel.link_input_data = [shareDir+'/min-{0}_{1}.tpr > min-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
                t.link_input_data = [ shareDir+'/min-{0}_{1}.tpr > min-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]
                
                t.copy_output_data = ['min-{0}_{1}.gro >'.format(iter_cnt,t_cnt) +shareDir+'/min-{0}_{1}.gro'.format(iter_cnt,t_cnt)]
            
                s2.add_tasks(t)

            p.add_stages(s2)

            # STAGE/ KERNEL 3
            s3 = Stage()
            s3.name = 'grompp-stage-%s'%iter_cnt

            for t_cnt in range(ensemble_size):
                t = task()
                t.name = 'grompp-task-%s't_cnt

                #k3_prep_eq_kernel.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_mdp)),
                #                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file)),
                #                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.restr_file)),
                #                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_itp_file))]

                t.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_mdp)),
                                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file)),
                                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.restr_file)),
                                                     shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_itp_file))]

                #k3_prep_eq_kernel.link_input_data = k3_prep_eq_kernel.link_input_data + [shareDir+'/min-{0}_{1}.gro > min-{0}_{1}.gro'.format(iterMod-1,instance-1)]
                #k3_prep_eq_kernel.link_input_data = k3_prep_eq_kernel.link_input_data + [shareDir+'/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iterMod-2,instance-1,ext)]
                
                t.link_input_data = t.link_input_data + [shareDir+'/min-{0}_{1}.gro > min-{0}_{1}.gro'.format(iter_cnt,t_cnt)]
                                    # NOTE - ext - is this properly being assigned in this program? 
                t.link_input_data = t.link_input_data + [shareDir+'/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iter_cnt-1,t_cnt,ext)]
 
                #k3_prep_eq_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_2_mdp)),
                #                               "--ref={0}_{1}_{2}.{3}".format(outbase,iterMod-2,instance-1,ext),
                #                               "--top={0}".format(os.path.basename(Kconfig.top_file)),
                #                               "--gro=min-{0}_{1}.gro".format(iterMod-1,instance-1),
                #                               "--tpr=eq-{0}_{1}.tpr".format(iterMod-1,instance-1)]

                t.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_2_mdp)),
                                "--ref={0}_{1}_{2}.{3}".format(outbase,iter_cnt-1,t_cnt,ext),
                                "--top={0}".format(os.path.basename(Kconfig.top_file)),
                                "--gro=min-{0}_{1}.gro".format(iter_cnt,t_cnt),
                                "--tpr=eq-{0}_{1}.tpr".format(iter_cnt,t_cnt)]

                #k3_prep_eq_kernel.copy_output_data = ['eq-{0}_{1}.tpr > '.format(iterMod-1,instance-1) +shareDir+'/eq-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
                t.copy_output_data = ['eq-{0}_{1}.tpr > '.format(iter_cnt,t_cnt) +shareDir+'/eq-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]
                s3.add_tasks(t)
            p.add_stages(s3)


            # STAGE/ KERNEL 4
            s4 = Stage()
            s4.name = 'mdrun-stage-%s'%iter_cnt

            for t_cnt in range(ensemble_size):
                t = task()
                t.name = 'mdrun-task-%s'%t_cnt
                #k4_eq_kernel.link_input_data = [shareDir+'/eq-{0}_{1}.tpr > eq-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
                t.link_input_data = [shareDir+'/eq-{0}_{1}.tpr > eq-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]
                
                ## CHECK THIS 
                #k4_eq_kernel.cores = Kconfig.num_cores_per_sim_cu - VIVEK - does this change with that cpu config option 

                #k4_eq_kernel.arguments = ["--deffnm=eq-{0}_{1}".format(iterMod-1,instance-1)]
                t.arguments = ["--deffnm=eq-{0}_{1}".format(iter_cnt,t_cnt)]
                
                #k4_eq_kernel.copy_output_data = ['eq-{0}_{1}.gro > '.format(iter_cnt,t_cnt) +shareDir+'/eq-{0}_{1}.gro'.format(iter_cnt,t_cnt)]
                t.copy_output_data = ['eq-{0}_{1}.gro > '.format(iter_cnt,t_cnt) +shareDir+'/eq-{0}_{1}.gro'.format(iter_cnt,t_cnt)]

                s4.add_tasks(t)
            p.add_stages(s4)


        # TODO: Kernel 5 -> Stage 5
        # STAGE / KERNEL 5 
        s5 = Stage()
        s5.name = 'grompp-stage-%s'%iter_cnt

        for t_cnt in range(ensemble_size):

            t = Task()
            t.name = 's5-task-%s'%t_cnt

            k5_prep_sim_kernel.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_3_mdp)),
                                             shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file))]
            
            # Same pre_exec as in the kernel def files for the specific target resource
            t.pre_exec = [  "export PATH=$PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                            "export GROMACS_LIB=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64",
                            "export GROMACS_INC=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/include",
                            "export GROMACS_BIN=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                            "export GROMACS_DIR=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu",
                            "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64"]
            # Same executable as in the kernel def files for the specific target resource                                
            t.executable = ["gmx grompp"]
            
            if((iterMod-1)==0):
                #k5_prep_sim_kernel.link_input_data =  k5_prep_sim_kernel.link_input_data + [shareDir+'/{0}'.format(os.path.basename(Kconfig.initial_crd_file))]
                
                # Copy link_input_data from main file, edit 'iterMod'->iter_cnt and 
                # 'instance'->t_cnt appropriately
                t.link_input_data = t.link_input_data + [shareDir+'/{0}'.format(os.path.basename(Kconfig.initial_crd_file))]

                #k5_prep_sim_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_3_mdp)),
                #                               "--gro={0}".format(os.path.basename(Kconfig.initial_crd_file)),
                #                               "--top={0}".format(os.path.basename(Kconfig.top_file)),
                #                               "--tpr=md-{0}_{1}.tpr".format(iterMod-1,instance-1)]  

                                # Same flags arguments as in the kernel def file, but secondary arguments from the main file
                t.arguments = [ '-f','{0}'.format(os.path.basename(Kconfig.grompp_3_mdp)),
                                '-c','{0}'.format(os.path.basename(Kconfig.initial_crd_file)),
                                '-p','{0}'.format(os.path.basename(Kconfig.top_file)),
                                '-o','min-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]
            else:
                #k5_prep_sim_kernel.link_input_data =  k5_prep_sim_kernel.link_input_data + [shareDir+'/eq-{0}_{1}.gro > eq-{0}_{1}.gro'.format(iterMod-1,instance-1)]
                
                t.link_input_data = t.link_input_data + [shareDir+'/eq-{0}_{1}.gro > eq-{0}_{1}.gro'.format(iter_cnt,t_cnt)]
                
                #k5_prep_sim_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_3_mdp)),
                #                               "--gro=eq-{0}_{1}.gro".format(iterMod-1,instance-1),
                #                               "--top={0}".format(os.path.basename(Kconfig.top_file)),
                #                               "--tpr=md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
                
                # Same flags arguments as in the kernel def file, but secondary arguments from the main file
                t.arguments = [ '-f','{0}'.format(os.path.basename(Kconfig.grompp_3_mdp)),
                                '-c','eq-{0}_{1}.gro'.format(iter_cnt,t_cnt),
                                '-p','{0}'.format(os.path.basename(Kconfig.top_file)),
                                '-o','md-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]

            #k5_prep_sim_kernel.copy_output_data = ['md-{0}_{1}.tpr > $SHARED/md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]        
            #k5_prep_sim_kernel.copy_output_data = ['md-{0}_{1}.tpr > '.format(iterMod-1,instance-1) + shareDir +'/md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            s5.copy_output_data = ['md-{0}_{1}.tpr > '.format(iter_cnt,t_cnt) + shareDir +'/md-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]

            s5.add_tasks(t)

        p.add_stages(s5)



        # STAGE / KERNEL 6 

        s6 = Stage ()
        s6.name = 'mdrun-stage-%s'%iter_cnt

        for t_cnt in range(ensemble_size):

            #Create a Task Object
            t = Task()
            t.name = 's6-task-%s'%t_cnt

            # Same pre_exec as in k def files for specific target resources 
            t.pre_exec = ["export PATH=$PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/bin",
                          "export GROMACS_LIB=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/lib64",
                          "export GROMACS_INC=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/include",
                          "export GROMACS_BIN=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/bin",
                          "export GROMACS_DIR=/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu",
                          "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210_OMPI20151210-DYN/install-cpu/lib64"]
            t.executable = ["gmx_mpi mdrun"] 

            # Number of cores for non-mpi Gromacs
            #k6_sim_kernel.arguments = ["--deffnm=md-{0}_{1}".format(iterMod-1,instance-1)]    
            t.arguments = ["-deffnm=md-{0}_{1}".format(iter_cnt,t_cnt)]

            #k6_sim_kernel.link_input_data = [shareDir+'/md-{0}_{1}.tpr > md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            t.link_input_data = [ shareDir+'/md-{0}_{1}.tpr > md-{0}_{1}.tpr'.format(iter_cnt,t_cnt)]

            k6_sim_kernel.copy_output_data = ["md-{0}_{1}.gro > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.gro".format(iterMod-1,instance-1),
                                             "md-{0}_{1}.xtc > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]

            t.copy_output_data = ["md-{0}_{1}.gro > ".format(iter_cnt,t_cnt) +shareDir+"/md-{0}_{1}.gro".format(iter_cnt,t_cnt),
                                  "md-{0}_{1}.xtc > ".format(iter_cnt,t_cnt) +shareDir+"/md-{0}_{1}.xtc".format(iter_cnt,t_cnt)]
            
            s6.add_tasks(t)

        p.add_stages(s6)     



        # STAGE / KERNEL 7 
        s7 = Stage()
        s7.name = 'traj-stage-%s'%iter_cnt
        
        for t_cnt in range(ensemble_size):
            t = task()
            t.name = 'traj-task-%s'%t_cnt

            #k7_sim_kernel.link_input_data = [shareDir+"/md-{0}_{1}.gro > md-{0}_{1}.gro".format(iterMod-1,instance-1),
            #                                 shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
            
            t.link_input_data = [shareDir+"/md-{0}_{1}.gro > md-{0}_{1}.gro".format(iter_cnt,t_cnt),
                                             shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iter_cnt,t_cnt)]

            t.pre_exec = ["export PATH=$PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                                "export GROMACS_LIB=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64",
                                "export GROMACS_INC=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/include",
                                "export GROMACS_BIN=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/bin",
                                "export GROMACS_DIR=/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu",
                                "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/sciteam/gkd/gromacs/5.1.1/20151210-NO_MPI/install-cpu/lib64"]
            t.executable = ["/bin/bash"]

            #k7_sim_kernel.arguments = ["--echo=System",
            #                           "--f=md-{0}_{1}.gro".format(iterMod-1,instance-1),
            #                           "--s=md-{0}_{1}.tpr".format(iterMod-1,instance-1),
            #                           "--o=md-{0}_{1}_whole.gro".format(iterMod-1,instance-1),
            #                           "--pbc=whole"]
            
            # NOTE - check this one conversion was little tricky 
            t.arguments = ['-l','-c','echo System | gmx trjconv -f md-{0}_{1}.gro -o md-{0}_{1}_whole.gro -s md-{0}_{1}.tpr -pbc whole'.format(iter_cnt,t_cnt)]


            #k7_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.gro > $SHARED/md-{0}_{1}.gro".format(iterMod-1,instance-1)]        
            k7_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.gro > ".format(iter_cnt,t_cnt) +shareDir+"/md-{0}_{1}.gro".format(iter_cnt,t_cnt)]       
            t.copy_output_data = ["md-{0}_{1}_whole.gro > ".format(iter_cnt,t_cnt) +shareDir+"/md-{0}_{1}.gro".format(iter_cnt,t_cnt)]
            s7.add_tasks(t)

        p.add_stages(s7)


        #KERNEL / STAGE 8 
        s8 = Stage()
        s8.name = 'trajconv-stage-%s'iter_cnt

        for t_cnt in range(ensemble_size):
            t = Task()
            t.name = 'traj-task-%s'%t_cnt

            #k8_sim_kernel.link_input_data = [shareDir+"/md-{0}_{1}.xtc > md-{0}_{1}.xtc".format(iterMod-1,instance-1),
            #                              shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
            
            t.link_input_data = [shareDir+"/md-{0}_{1}.xtc > md-{0}_{1}.xtc".format(iterMod-1,instance-1),
                                          shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]

            k8_sim_kernel.arguments = ["--echo=System",
                                       "--f=md-{0}_{1}.xtc".format(iterMod-1,instance-1),
                                       "--s=md-{0}_{1}.tpr".format(iterMod-1,instance-1),
                                       "--o=md-{0}_{1}_whole.xtc".format(iterMod-1,instance-1),
                                       "--pbc=whole"]
            t.arguments = ['-l','-c','echo System | gmx trjconv -f md-{0}_{1}.xtc -o md-{0}_{1}_whole.xtc -s md-{0}_{1}.tpr -pbc whole'.format(iter_cnt,t_cnt)]
            
            #if(iterMod%Kconfig.nsave==0):
            if(iter_cnt%Kconfig.nsave==0):
                k8_sim_kernel.download_output_data = ["md-{0}_{1}_whole.xtc > output/iter{0}/md-{0}_{1}_whole.xtc".format(iterMod-1,instance-1)]
                t.download_output_data = ["md-{0}_{1}_whole.xtc > output/iter{0}/md-{0}_{1}_whole.xtc".format(iter_cnt,t_cnt)]

            #k8_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.xtc > $SHARED/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]        
            #k8_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.xtc > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]        
            t.copy_output_data = ["md-{0}_{1}_whole.xtc > ".format(iter_cnt,t_cnt) +shareDir+"/md-{0}_{1}.xtc".format(iter_cnt,t_cnt)]

            s8.add_tasks(t)

        p.add_stages(s8)


    return p        
    # # Add the Task to the Stage
    # s1.add_tasks(t1)

    # # Add Stage to the Pipeline
    # p.add_stages(s1)

    # # Create another Stage object to hold character count tasks
    # s2 = Stage()
    # s2.name = 's2'
    # s2_task_uids = []

    # for cnt in range(10):

    #     # Create a Task object
    #     t2 = Task()
    #     t2.name = 't%s' % (cnt + 1)
    #     t2.executable = ['/bin/bash']
    #     t2.arguments = ['-l', '-c', 'grep -o . output.txt | sort | uniq -c > ccount.txt']
    #     # Copy data from the task in the first stage to the current task's location
    #     t2.copy_input_data = ['$Pipline_%s_Stage_%s_Task_%s/output.txt' % (p.name, s1.name, t1.name)]

    #     # Add the Task to the Stage
    #     s2.add_tasks(t2)
    #     s2_task_uids.append(t2.name)

    # # Add Stage to the Pipeline
    # p.add_stages(s2)

    # # Create another Stage object to hold checksum tasks
    # s3 = Stage()
    # s3.name = 's3'

    # for cnt in range(10):

    #     # Create a Task object
    #     t3 = Task()
    #     t3.name = 't%s' % (cnt + 1)
    #     t3.executable = ['/bin/bash']
    #     t3.arguments = ['-l', '-c', 'sha1sum ccount.txt > chksum.txt']
    #     # Copy data from the task in the first stage to the current task's location
    #     t3.copy_input_data = ['$Pipline_%s_Stage_%s_Task_%s/ccount.txt' % (p.name, s2.name, s2_task_uids[cnt])]
    #     # Download the output of the current task to the current location
    #     t3.download_output_data = ['chksum.txt > chksum_%s.txt' % cnt]

    #     # Add the Task to the Stage
    #     s3.add_tasks(t3)

    # # Add Stage to the Pipeline
    # p.add_stages(s3)
    def analysis_stage(self, iteration, instance):
        '''
        function : Perform CoCo Analysis on the output of the simulation from the current iterMod. Using the .xtc
         files generated in all instances, generate .gro files (as many as the num_CUs) to be used in the next simulations. 
        

        coco :-

                Purpose : Runs CoCo analysis on a set of MD trajectory files in this case xtc files and generates several coordinates file to be

                Arguments : --grid           = Number of points along each dimension of the CoCo histogram
                            --dims           = The number of projections to consider from the input pcz file
                            --frontpoints    = Number of CUs
                            --topfile        = Topology filename
                            --mdfile         = MD Input filename
                            --output         = Output filename
                            --cycle          = Current iterMod number
                            --atom_selection = Selection of the biological part of the system we want to consider for analysis
        '''
        #shareDir="$SHARED"
        #shareDir="/work/fbettenc/radical.pilot.sandbox/rp.session.js-17-187.jetstream-cloud.org.hal9000.017508.0005-pilot.0000/staging_area"
    shareDir = "/work/fbettenc/p14b01_pool/staging_area"

        prev_sim_last_iter_to_use=48
        iterMod=iteration+prev_sim_last_iter_to_use

        k1_ana_kernel = Kernel(name="custom.coco")

        outbase, ext = os.path.basename(Kconfig.output).split('.')
        if ext == '':
            ext = '.pdb'        
        
        k1_ana_kernel.arguments = ["--grid={0}".format(Kconfig.grid),
                                   "--dims={0}".format(Kconfig.dims),
                                   "--frontpoints={0}".format(Kconfig.num_CUs),
                                   "--topfile=md-{0}_0.gro".format(iterMod-1),
                                   "--mdfile=*.xtc",
                                   "--output={0}_{1}_.gro".format(outbase,iterMod-1),
                                   "--atom_selection={0}".format(Kconfig.sel)]
       # k1_ana_kernel.cores = min(Kconfig.num_CUs,RPconfig.PILOTSIZE)
        k1_ana_kernel.cores = min(Kconfig.num_CUs*(iterMod+1),RPconfig.PILOTSIZE) # set to iterMod+1 bec at first iter coco analysis of k8 output so coco is iter ahead sort of

        print " "
    print "iter,iterMod,AnaCUcores = ",iteration,", ",iterMod,", ", k1_ana_kernel.cores
    print " "

    k1_ana_kernel.uses_mpi = True
        k1_ana_kernel.link_input_data = [shareDir+'/md-{1}_0.gro > md-{1}_0.gro'.format(iterMod,iterMod-1)]
        for iter in range(1,iterMod+1):
            for i in range(1,Kconfig.num_CUs+1):        
                k1_ana_kernel.link_input_data = k1_ana_kernel.link_input_data + [shareDir+'/md-{2}_{3}.xtc > md-{2}_{3}.xtc'.format(iter,i,iter-1,i-1)]
        
                
        k1_ana_kernel.copy_output_data = []
        for i in range(0,Kconfig.num_CUs):
            #k1_ana_kernel.copy_output_data += ["{0}_{1}_{2}.gro > $SHARED/{0}_{1}_{2}.gro".format(outbase,iterMod-1,i,ext)]
            k1_ana_kernel.copy_output_data += ["{0}_{1}_{2}.gro > ".format(outbase,iterMod-1,i,ext) +shareDir+"/{0}_{1}_{2}.gro".format(outbase,iterMod-1,i,ext)]


        k1_ana_kernel.download_output_data = ["coco.log > output/coco-iter{0}.log".format(iterMod-1)]   
        

        return [k1_ana_kernel]
        
    def post_loop(self):
        pass


if __name__ == '__main__':


    try: 
	parser = argparse.ArgumentParser()
        parser.add_argument('--RPconfig', help='link to Radical Pilot related configurations file')
        parser.add_argument('--Kconfig', help='link to Kernel configurations file')

        args = parser.parse_args()

        if args.RPconfig is None:
            parser.error('Please enter a RP configuration file')
            sys.exit(1)
        if args.Kconfig is None:
            parser.error('Please enter a Kernel configuration file')
            sys.exit(0)

        RPconfig = imp.load_source('RPconfig', args.RPconfig)
        Kconfig = imp.load_source('Kconfig', args.Kconfig)

        # Retrieve basename and extension
        outbase, ext = os.path.basename(Kconfig.output).split('.')
        if ext == '':
            ext = '.pdb'


        # Create Application Manager
        appman = AppManager(hostname=hostname, port=port)

        # Create a dictionary describe four mandatory keys:
        # resource, walltime, and cpus
        # resource is 'local.localhost' to execute locally
        res_dict = {

            'resource': RPconfig.REMOTE_HOST,
            'cpus': RPconfig.PILOTSIZE,
            'walltime':RPconfig.WALLTIME,
            'project': RPconfig.ALLOCATION, #project
            'queue': RPconfig.QUEUE,
            'access_schema': 'gsissh'
        }

        # Assign resource request description to the Application Manager
        appman.resource_desc = res_dict

        # Specify shared data
        appman.shared_data = [  Kconfig.initial_crd_file,
                                Kconfig.grompp_1_mdp,
                                Kconfig.grompp_2_mdp,
                                Kconfig.grompp_3_mdp,
                                Kconfig.grompp_1_itp_file,
                                Kconfig.grompp_2_itp_file,
                                Kconfig.top_file,
                                Kconfig.restr_file]

        # Assign the workflow as a set or list of Pipelines to the Application Manager
        # Note: The list order is not guaranteed to be preserved
        #def generate_pipeline(index, iterations, ensemble_size):
	extasy_pipeline = generate_pipeline(1,Kconfig.num_iterations,Kconfig.num_CUs)
        appman.workflow = [extasy_pipeline]

        # Run the Application Manager
        appman.run()

    except Exception as ex:
        print 'Current run failed with error: ', ex
        raise
