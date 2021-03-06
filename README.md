The aim of the project is to explore the performance of one of the ExTASY tools
– CoCo-MD – for advanced simulation when applied to a biologically-relevant test
case, the protein BPTI. Very long and computationally expensive conventional
simulations have shown how hard it is to properly explore the shape-space of
BPTI – can CoCo-MD do better? Using the ExTASY toolkit, you will run advanced
sampling simulations of BPTI on a range of supercomputer resources. Working with
biomoleciular simulation scientists in the UK, you will analyse the sampling
data produced, and optimise the simulation parameters to identify how best
performance can be achieved.

In the folder *Papers* is a copy of the 2010 Science paper by Shaw et al that 
describes their 1 millisecond simulation of BPTI.

In the folder *Shaw_Data_Analysis* is a copy of the trajectory file from that 
simulation, a copy of a PDB-format file of the sytructure of bpti, and an ipython 
notebook (and markdown version) showing how I have used PCA analysis to identify 
the 'rare' conformational state seen in that simulation. Also there is another
ipython notebook that shows how you can analyse a slug of trajectory data of the
type that this project should produce (in much greater volume!), and see if it
is getting anywhere near this rare state.

The tarball *gmxcoco.tgz* contains a first stab at the CoCo-MD workflow for this
project. It has been adapted from the examples used in the Edinburgh workshop
last year. It has been tested on Archer and runs OK, but will probably need 
tweaking for use on another HPC resource, and in any case the real thing will
need more replicates and cycles. A couple of README files are included which may 
help.
