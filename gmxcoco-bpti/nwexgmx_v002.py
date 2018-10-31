#!/usr/bin/env python

__author__        = "The ExTASY project <ardita.shkurti@nottingham.ac.uk>"
__copyright__     = "Copyright 2015, http://www.extasy-project.org/"
__license__       = "MIT"
__use_case_name__ = "'Gromacs + CoCo' simulation-analysis (ExTASY)."


from radical.ensemblemd import Kernel
from radical.ensemblemd import EnsemblemdError
from radical.ensemblemd import SimulationAnalysisLoop
from radical.ensemblemd import ResourceHandle
from radical.ensemblemd.engine import get_engine

import imp
import argparse
import sys
import os
import json

# ------------------------------------------------------------------------------
# Set default verbosity

if os.environ.get('RADICAL_ENTK_VERBOSE') == None: #note ENTK was ENMD
    os.environ['RADICAL_ENTK_VERBOSE'] = 'INFO'  #note ENTK was ENMD

# ------------------------------------------------------------------------------
#Load all custom Kernels

from kernel_defs.grompp import grompp_Kernel
get_engine().add_kernel_plugin(grompp_Kernel)

from kernel_defs.mdrun import mdrun_Kernel
get_engine().add_kernel_plugin(mdrun_Kernel)

from kernel_defs.trjconv import trjconv_Kernel
get_engine().add_kernel_plugin(trjconv_Kernel)

from kernel_defs.coco import kernel_coco
get_engine().add_kernel_plugin(kernel_coco)

# ------------------------------------------------------------------------------
#

class Extasy_CocoGromacs_Static(SimulationAnalysisLoop):

    def __init__(self, maxiterMods, simulation_instances, analysis_instances):
        SimulationAnalysisLoop.__init__(self, maxiterMods, simulation_instances, analysis_instances)

    def pre_loop(self):
		pass

    def simulation_stage(self, iteration, instance):

	#shareDir="$SHARED"
	#shareDir="staging://" # $SHARED is place holder and is replaced at runtime by "staging://"
		#https://github.com/radical-cybertools/radical.entk/blob/master/src/radical/entk/execution_plugin/staging/placeholders.py#L25
        #shareDir="/work/fbettenc/radical.pilot.sandbox/p13b01_left_d3_k12_1000_k34_1000"
		# note tried without / before work and failed. diff err for /work/.. than work/..
	#shareDir="/work/fbettenc/radical.pilot.sandbox/rp.session.js-17-187.jetstream-cloud.org.hal9000.017508.0005-pilot.0000/staging_area"
        shareDir="/work/fbettenc/p14b01_pool/staging_area"
	prev_sim_last_iter_to_use=48
	iterMod=iteration+prev_sim_last_iter_to_use

	kernel_list = []

        outbase, ext = os.path.basename(Kconfig.output).split('.')
        if ext == '':
            ext = '.pdb'


        if ((iterMod-1)!=0):
            k1_prep_min_kernel = Kernel(name="custom.grompp")
            k1_prep_min_kernel.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_1_mdp)),
                                                  shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file)),
                                                  shareDir+'/{0}'.format(os.path.basename(Kconfig.restr_file)),
                                                  shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_1_itp_file))]
            #k1_prep_min_kernel.link_input_data = k1_prep_min_kernel.link_input_data + ['$PREV_ANALYSIS_INSTANCE_1/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iterMod-2,instance-1,ext)]
            k1_prep_min_kernel.link_input_data = k1_prep_min_kernel.link_input_data + [shareDir+'/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iterMod-2,instance-1,ext)]

            k1_prep_min_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_1_mdp)),
                                            "--ref={0}_{1}_{2}.{3}".format(outbase,iterMod-2,instance-1,ext),
                                            #"--ref={0}".format(os.path.basename(Kconfig.restr_file)),
                                            "--top={0}".format(os.path.basename(Kconfig.top_file)),
                                            "--gro={0}".format(os.path.basename(Kconfig.restr_file)),
                                            #"--gro={0}_{1}_{2}.{3}".format(outbase,iterMod-2,instance-1,ext),
                                            "--tpr=min-{0}_{1}.tpr".format(iterMod-1,instance-1)]
            #k1_prep_min_kernel.copy_output_data = ['min-{0}_{1}.tpr > $SHARED/min-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            k1_prep_min_kernel.copy_output_data = ['min-{0}_{1}.tpr > '.format(iterMod-1,instance-1) +shareDir+'/min-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            kernel_list.append(k1_prep_min_kernel)

            k2_min_kernel = Kernel(name="custom.mdrun")
            k2_min_kernel.link_input_data = [shareDir+'/min-{0}_{1}.tpr > min-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            k2_min_kernel.cores = Kconfig.num_cores_per_sim_cu
            k2_min_kernel.arguments = ["--deffnm=min-{0}_{1}".format(iterMod-1,instance-1)]
            #k2_min_kernel.copy_output_data = ['min-{0}_{1}.gro > $SHARED/min-{0}_{1}.gro'.format(iterMod-1,instance-1)]
            k2_min_kernel.copy_output_data = ['min-{0}_{1}.gro >'.format(iterMod-1,instance-1) +shareDir+'/min-{0}_{1}.gro'.format(iterMod-1,instance-1)]

            kernel_list.append(k2_min_kernel)

            k3_prep_eq_kernel = Kernel(name="custom.grompp")
            k3_prep_eq_kernel.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_mdp)),
                                                 shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file)),
                                                 shareDir+'/{0}'.format(os.path.basename(Kconfig.restr_file)),
                                                 shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_2_itp_file))]
            k3_prep_eq_kernel.link_input_data = k3_prep_eq_kernel.link_input_data + [shareDir+'/min-{0}_{1}.gro > min-{0}_{1}.gro'.format(iterMod-1,instance-1)]
            #k3_prep_eq_kernel.link_input_data = k3_prep_eq_kernel.link_input_data + ['$PREV_ANALYSIS_INSTANCE_1/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iterMod-2,instance-1,ext)]
            k3_prep_eq_kernel.link_input_data = k3_prep_eq_kernel.link_input_data + [shareDir+'/{0}_{1}_{2}.{3} > {0}_{1}_{2}.{3}'.format(outbase,iterMod-2,instance-1,ext)]
            k3_prep_eq_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_2_mdp)),
                                           "--ref={0}_{1}_{2}.{3}".format(outbase,iterMod-2,instance-1,ext),
                                           #"--ref=min-{0}_{1}.gro".format(iterMod-1,instance-1),
                                           "--top={0}".format(os.path.basename(Kconfig.top_file)),
                                           ##"--gro=min-{0}_{1}".format(iterMod-1,instance-1),
                                           "--gro=min-{0}_{1}.gro".format(iterMod-1,instance-1),
                                           #"--gro={0}_{1}_{2}.{3}".format(outbase,iterMod-2,instance-1,ext),
                                           "--tpr=eq-{0}_{1}.tpr".format(iterMod-1,instance-1)]
            #k3_prep_eq_kernel.copy_output_data = ['eq-{0}_{1}.tpr > $SHARED/eq-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            k3_prep_eq_kernel.copy_output_data = ['eq-{0}_{1}.tpr > '.format(iterMod-1,instance-1) +shareDir+'/eq-{0}_{1}.tpr'.format(iterMod-1,instance-1)]

            kernel_list.append(k3_prep_eq_kernel)

            k4_eq_kernel = Kernel(name="custom.mdrun")
            k4_eq_kernel.link_input_data = [shareDir+'/eq-{0}_{1}.tpr > eq-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
            k4_eq_kernel.cores = Kconfig.num_cores_per_sim_cu
            k4_eq_kernel.arguments = ["--deffnm=eq-{0}_{1}".format(iterMod-1,instance-1)]
            #k4_eq_kernel.copy_output_data = ['eq-{0}_{1}.gro > $SHARED/eq-{0}_{1}.gro'.format(iterMod-1,instance-1)]
            k4_eq_kernel.copy_output_data = ['eq-{0}_{1}.gro > '.format(iterMod-1,instance-1) +shareDir+'/eq-{0}_{1}.gro'.format(iterMod-1,instance-1)]

            kernel_list.append(k4_eq_kernel)

        k5_prep_sim_kernel = Kernel(name="custom.grompp")
        k5_prep_sim_kernel.link_input_data = [shareDir+'/{0}'.format(os.path.basename(Kconfig.grompp_3_mdp)),
                                             shareDir+'/{0}'.format(os.path.basename(Kconfig.top_file))]
        if((iterMod-1)==0):
            k5_prep_sim_kernel.link_input_data =  k5_prep_sim_kernel.link_input_data + [shareDir+'/{0}'.format(os.path.basename(Kconfig.initial_crd_file))]
            k5_prep_sim_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_3_mdp)),
                                           "--gro={0}".format(os.path.basename(Kconfig.initial_crd_file)),
                                           "--top={0}".format(os.path.basename(Kconfig.top_file)),
                                           "--tpr=md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
        else:
            k5_prep_sim_kernel.link_input_data =  k5_prep_sim_kernel.link_input_data + [shareDir+'/eq-{0}_{1}.gro > eq-{0}_{1}.gro'.format(iterMod-1,instance-1)]
            k5_prep_sim_kernel.arguments = ["--mdp={0}".format(os.path.basename(Kconfig.grompp_3_mdp)),
                                           "--gro=eq-{0}_{1}.gro".format(iterMod-1,instance-1),
                                           "--top={0}".format(os.path.basename(Kconfig.top_file)),
                                           "--tpr=md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
        #k5_prep_sim_kernel.copy_output_data = ['md-{0}_{1}.tpr > $SHARED/md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
        k5_prep_sim_kernel.copy_output_data = ['md-{0}_{1}.tpr > '.format(iterMod-1,instance-1) + shareDir +'/md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]

        kernel_list.append(k5_prep_sim_kernel)

        k6_sim_kernel = Kernel(name="custom.mdrun")
        k6_sim_kernel.link_input_data = [shareDir+'/md-{0}_{1}.tpr > md-{0}_{1}.tpr'.format(iterMod-1,instance-1)]
        k6_sim_kernel.cores = Kconfig.num_cores_per_sim_cu
        k6_sim_kernel.arguments = ["--deffnm=md-{0}_{1}".format(iterMod-1,instance-1)]
        #k6_sim_kernel.copy_output_data = ["md-{0}_{1}.gro > $SHARED/md-{0}_{1}.gro".format(iterMod-1,instance-1),
        #                                  "md-{0}_{1}.xtc > $SHARED/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]
        k6_sim_kernel.copy_output_data = ["md-{0}_{1}.gro > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.gro".format(iterMod-1,instance-1),
                                          "md-{0}_{1}.xtc > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]
        kernel_list.append(k6_sim_kernel)

        k7_sim_kernel = Kernel(name="custom.trjconv")
       # k7_sim_kernel.link_input_data = ["$SHARED/md-{0}_{1}.gro > md-{0}_{1}.gro".format(iterMod-1,instance-1),
       #                                  "$SHARED/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
        k7_sim_kernel.link_input_data = [shareDir+"/md-{0}_{1}.gro > md-{0}_{1}.gro".format(iterMod-1,instance-1),
                                         shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
        k7_sim_kernel.arguments = ["--echo=System",
                                   "--f=md-{0}_{1}.gro".format(iterMod-1,instance-1),
                                   "--s=md-{0}_{1}.tpr".format(iterMod-1,instance-1),
                                   "--o=md-{0}_{1}_whole.gro".format(iterMod-1,instance-1),
                                   "--pbc=whole"]
        #k7_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.gro > $SHARED/md-{0}_{1}.gro".format(iterMod-1,instance-1)]
        k7_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.gro > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.gro".format(iterMod-1,instance-1)]

        kernel_list.append(k7_sim_kernel)

        k8_sim_kernel = Kernel(name="custom.trjconv")
        #k8_sim_kernel.link_input_data = ["$SHARED/md-{0}_{1}.xtc > md-{0}_{1}.xtc".format(iterMod-1,instance-1),
        #                                 "$SHARED/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]
        k8_sim_kernel.link_input_data = [shareDir+"/md-{0}_{1}.xtc > md-{0}_{1}.xtc".format(iterMod-1,instance-1),
                                         shareDir+"/md-{0}_{1}.tpr > md-{0}_{1}.tpr".format(iterMod-1,instance-1)]

        k8_sim_kernel.arguments = ["--echo=System",
                                   "--f=md-{0}_{1}.xtc".format(iterMod-1,instance-1),
                                   "--s=md-{0}_{1}.tpr".format(iterMod-1,instance-1),
                                   "--o=md-{0}_{1}_whole.xtc".format(iterMod-1,instance-1),
                                   "--pbc=whole"]
        if(iterMod%Kconfig.nsave==0):
            k8_sim_kernel.download_output_data = ["md-{0}_{1}_whole.xtc > output/iter{0}/md-{0}_{1}_whole.xtc".format(iterMod-1,instance-1)]
        #k8_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.xtc > $SHARED/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]
        k8_sim_kernel.copy_output_data = ["md-{0}_{1}_whole.xtc > ".format(iterMod-1,instance-1) +shareDir+"/md-{0}_{1}.xtc".format(iterMod-1,instance-1)]

        kernel_list.append(k8_sim_kernel)


        return kernel_list


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

# ------------------------------------------------------------------------------
#
if __name__ == "__main__":

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

        # Create a new static execution context with one resource and a fixed
        # number of cores and runtime.

        cluster = ResourceHandle(
            resource=RPconfig.REMOTE_HOST,
            cores=RPconfig.PILOTSIZE,
            walltime=RPconfig.WALLTIME,
            username = RPconfig.UNAME, #username
            project = RPconfig.ALLOCATION, #project
            queue = RPconfig.QUEUE,
            database_url = RPconfig.DBURL,
	    access_schema = 'gsissh'
        )

        cluster.shared_data = [
                                Kconfig.initial_crd_file,
                                Kconfig.grompp_1_mdp,
                                Kconfig.grompp_2_mdp,
                                Kconfig.grompp_3_mdp,
                                Kconfig.grompp_1_itp_file,
                                Kconfig.grompp_2_itp_file,
                                Kconfig.top_file,
                                Kconfig.restr_file
                            ]

        cluster.allocate()

        coco_gromacs_static = Extasy_CocoGromacs_Static(maxiterMods=Kconfig.num_iterations, simulation_instances=Kconfig.num_CUs, analysis_instances=1)
        cluster.run(coco_gromacs_static)

        cluster.deallocate()

    except EnsemblemdError, er:

        print "The gromacs-coco ExTASY workflow not completed correctly due to an Ensemble MD Toolkit Error: {0}".format(str(er))
        raise # Just raise the execption again to get the backtrace
