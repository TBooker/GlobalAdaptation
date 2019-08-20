Global adaptation confounds the search for local adaptation
======

In this repo you'll find the scripts, config files and the plotting scripts that I used for our study **"Global adaptation confounds the search for local adaptation"** *link to preprint - add once it's submitted* .

There are three strands to this study: 
  1. [Simuations of a model of parapatry (or isolation with migration)](Parapatry/)
  2. [A population genetic model for the number of incomplete sweeps under panmixia](IncompleteSweeps/)
  3. [Simuations of a stepping-stone model](SteppingStone/)

Below I give an overview of each strand, but check out their individual repos for more detailed information.

*Because the simulation data files are large, I cannot add the raw data to the GitHub repo. I will include summary data, but the full simulation data I analysed will be deposited on Dryad. I'll add a link to Dryad soon.*

Simuating of a model of parapatry (or isolation with migration)
------
   I simulated adaptive evolution occurring in a populations structured according to a model of isolation with migration - or a two-deme island model. I performed the simulations in *SLiM* v3.2 making use of the tree-sequence recording feature. This feature records the coalescent history of your simulated population as a tree-sequence object. It is computationally efficient as neutral mutations (which do not influence the coalescent process *per se*) do not need to be tracked, reducing the computational burden. Adding neutral mutations to the stored tree-sequences once the simulation is complete is very similar to coalescent simulations (e.g. **ms** or **msprime**) so it takes very little time. 

The basic work-flow was as follows:
- Simulate adaptive evolution in structured populations in *SLiM*, storing the genealogy (and the beneficial mutations segregating in the population)
- Sprinkle neutral mutations onto the tree object using *PySLiM* and generate VCF files 
- Perform genome scans on the resulting VCF files using *VCFtools*

A population genetic model for the number of incomplete sweeps under panmixia
------
  
  The Mathematica file is deposited at REF_TO_DIR. This was the first project for which I needed to use Mathematica, so the code it is probably pretty horrible and inefficient for anyone who is used to that software. 

Also included in this directory is the plotting script I used to make the figure.

Simuations of a stepping-stone model
------
I simulated very large stepping-stone population; 500 demes of 5,000 individuals each (i.e. 2.5 million individuals). I chose to model a large population so that I could see advantageous mutations spread in the form of "Fisher waves". Modelling recurrent adaptive evolution in such a population was not feasible except with a two-locus genotype simulation.   

I simulated a simple 2-locus model, I applied the deterministic forces of selection, migration and recombination to genotype frequencies at two loci.  

We came up with a way of modelling recurrent adaptive evolution at many sites in the genome as follows: 
- Simulate genetic drift at a neutral locus for 100,000 generations, allele frequencies vary (100,000 replicates) 
- Examine the distribution of Fst from the neutral simulations to establish cut-off value
- Randomly draw a 

If there is anything incomplete, you find bugs or you want clarification on anything, please don't hesitate to get in touch (my_last_name [at] zoology [dot] ubc [dot] ca)

**A Disclaimer** If you want to use any of this code, keep an eye on path names. I am lazy and generally hard-code in paths and directories
