from radical.entk import Pipeline, Stage, Task, AppManager
import os

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
port = os.environ.get('RMQ_PORT', 5672)

# Share directory with ALL data
#shareDir="$SHARED"
#shareDir="staging://" # $SHARED is place holder and is replaced at runtime by "staging://" 
#https://github.com/radical-cybertools/radical.entk/blob/master/src/radical/entk/execution_plugin/staging/placeholders.py#L25	
#shareDir="/work/fbettenc/radical.pilot.sandbox/p13b01_left_d3_k12_1000_k34_1000" 
# # note tried without / before work and failed. diff err for /work/.. than work/..
#shareDir="/work/fbettenc/radical.pilot.sandbox/rp.session.js-17-187.jetstream-cloud.org.hal9000.017508.0005-pilot.0000/staging_area"

shareDir="/work/fbettenc/p14b01_pool/staging_area"

# References to previous iterations, effects starting iteration of current run
# prev_sim_last_iter_to_use=48
# iterMod=iteration+prev_sim_last_iter_to_use
prev_sim_last_iter_to_use=0


# Retrieve basename and extension
outbase, ext = None, None

# NOTE: All stages, tasks, filenames are indexed from 0.

def generate_pipeline(index, iterations, ensemble_size):

    # Create a Pipeline object
    p = Pipeline()
    p.name = 'extasy-pipeline-%s' %index

    # Determine total iterations of all runs, current+historic
    total_iterations = prev_sim_last_iter_to_use + iterations

    for iter_cnt in range(prev_sim_last_iter_to_use, total_iterations):

        if iter_cnt != 0:

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

                # Number of cores for non-mpi Gromacs
                # t.cpu_reqs = {
                #           ‘processes’: 1,                # Default is 1
                #           ‘process_type’: None,          # Default is None, Other option is 'MPI'
                #           ‘threads_per_process’: Kconfig.num_cores_per_sim_cu,     # Default is 1
                #           ‘thread_type’: OpenMP}         # Default is None, Other option is 'OpenMP'


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

        # TODO: Kernel 5 -> Stage 5

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

    return p

if __name__ == '__main__':


    try: 

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
        extasy_pipeline = generate_pipeline('1')
        appman.workflow = [extasy_pipeline]

        # Run the Application Manager
        appman.run()

    except Exception as ex:
        print 'Current run failed with error: ', ex
        raise