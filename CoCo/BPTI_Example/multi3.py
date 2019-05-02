from xbowflow import xflowlib
from xbowflow.xflowlib import SubprocessKernel
from extasycoco import complement
import mdtraj as mdt
import os
import time
import subprocess

def restrained_md(tprfiles):
    """
    Run a Gromacs multi job via xflow, returning final coordinatest
    """
    nreps = len(tprfiles)
    num_pes = int(os.environ.get('PBS_NP'))
    multi = SubprocessKernel('aprun -n {} gmx_mpi mdrun -multi {} -s x.tpr'.format(num_pes, nreps))
    multi.set_inputs(['x{}.tpr'.format(i) for i in range(nreps)])
    multi.set_outputs(['confout{}.gro'.format(i) for i in range(nreps)])
    results = multi.run(*tprfiles)
    return results

def unrestrained_md(tprfiles):
    """
    Run a Gromacs multi job via xflow, returning trajectories
    """
    nreps = len(tprfiles)
    num_pes = int(os.environ.get('PBS_NP'))
    multi = SubprocessKernel('aprun -n {} gmx_mpi mdrun -multi {} -s x.tpr'.format(num_pes, nreps))
    multi.set_inputs(['x{}.tpr'.format(i) for i in range(nreps)])
    multi.set_outputs(['traj_comp{}.xtc'.format(i) for i in range(nreps)])
    results = multi.run(*tprfiles)
    return results

def unrestrained_grompp(crdfiles, itpfile, topfile, mdpfile):
    """
    Run multiple Grompp jobs
    """
    num_pes = int(os.environ.get('PBS_NP'))
    grompp = SubprocessKernel('aprun -n 1 gmx_mpi grompp -f x.mdp -c x.gro -p x.top -o x.tpr -maxwarn 1')
    grompp.set_inputs(['x.mdp', 'x.gro', 'x.top'])
    grompp.set_outputs(['x.tpr'])
    grompp.set_constant('x.mdp', mdpfile)
    grompp.set_constant('x.top', topfile)
    grompp.set_constant('posre.itp', itpfile)

    tprfiles = [grompp.run(c) for c in crdfiles]
    return tprfiles

def restrained_grompp(crdfiles, itpfile, topfile, mdpfile, targets):
    """
    Run multiple Grompp jobs, that include restraints
    """
    num_pes = int(os.environ.get('PBS_NP'))
    grompp = SubprocessKernel('aprun -n 1 gmx_mpi grompp -f x.mdp -c x.gro -p x.top -o x.tpr -r ref.gro -maxwarn 1')
    grompp.set_inputs(['x.mdp', 'x.gro', 'x.top', 'ref.gro'])
    grompp.set_outputs(['x.tpr'])
    grompp.set_constant('x.mdp', mdpfile)
    grompp.set_constant('x.top', topfile)
    grompp.set_constant('posre.itp', itpfile)

    tprfiles = [grompp.run(c, t) for c, t in zip(crdfiles, targets)]
    return tprfiles

def complement(trajfiles, topfile, gridsize=30, ndims=3, selection='all', npoints=1):
    """
    Run a coco job on a compute node.
    """
    trajlist = ' '.join(trajfiles)
    outlist = ' '.join(['coco_{}.pdb'.format(i) for i in range(npoints)])

    command = 'aprun -n 16 pyCoCo -i {} -s {} -t {} -d {} -g {} -o {} -l tmp.log'.format(trajlist, selection, topfile, ndims, gridsize, outlist)
    result = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    newtraj = mdt.load(outlist.split())
    return newtraj
    

if __name__ == '__main__':

    crdfile = 'inp_files/bpti.gro'
    topfile = 'inp_files/topol.top'
    mdpfile1 = 'inp_files/run1.mdp'
    mdpfile2 = 'inp_files/run2.mdp'
    mdpfile3 = 'inp_files/run3.mdp'

    selection = "'name CA and resid 1 to 55'"
    selection2 = 'name CA and resid 1 to 55'
    fc = 1000.0
    nreps = 20
    ncycles = 50

    crds = xflowlib.load(crdfile)
    top = xflowlib.load(topfile)
    mdp1 = xflowlib.load(mdpfile1)
    mdp2 = xflowlib.load(mdpfile2)
    mdp3 = xflowlib.load(mdpfile3)

    traj = mdt.load(crdfile)
    sel = traj.topology.select(selection2)
    with open('tmp.itp', 'w') as f:
        f.write(' [ position_restraints ]\n')
        for i in sel:
            f.write(' {:5d}  1  {} {} {}\n'.format(i + 1, fc, fc, fc))
    itp = xflowlib.load('tmp.itp')
    #os.remove('tmp.itp')

    all_trajfiles = [crdfile]
    crdlist = [crds] * nreps
    targets = crdlist
    cocolog = open('coco.log', 'w')
    for cycle in range(ncycles):
       md_start_time = time.time()
       tprs = restrained_grompp(crdlist, itp, top, mdp1, targets)
       mincrds = restrained_md(tprs)
       tprs = restrained_grompp(mincrds, itp, top, mdp2, targets)
       mdcrds = restrained_md(tprs)
       tprs = unrestrained_grompp(mdcrds, itp, top, mdp3)
       trajectories = unrestrained_md(tprs)
       md_end_time = time.time()
       print('Cycle {:3d}: Time for MD: {:8.0f} seconds.'.format(cycle, md_end_time - md_start_time))

       trajfile_names = ['bpti_{:04d}_{:04d}.xtc'.format(cycle, rep) for rep in range(nreps)]
       for r, t in zip(trajectories, trajfile_names):
           r.save(t)

       all_trajfiles += trajfile_names
       c = complement(all_trajfiles, crdfile, npoints=nreps, gridsize=30, selection=selection)
       targets = [c[i] for i in range(nreps)]
       with open('tmp.log') as f:
           cocolog.write(f.read())
       coco_end_time = time.time()
       print('           Time for CoCo analysis: {:8.0f} seconds.\n'.format(cycle, coco_end_time - md_end_time))
       cocolog.flush()

    #traj.save('bpti_all.xtc')
    cocolog.close()
