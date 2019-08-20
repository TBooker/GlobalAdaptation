
Simuating of a model of parapatry (or isolation with migration)
------
   I simulated adaptive evolution occurring in a populations structured according to a model of isolation with migration - or a two-deme island model. I performed the simulations in *SLiM* v3.2 making use of the tree-sequence recording feature. This feature records the coalescent history of your simulated population as a tree-sequence object. It is computationally efficient as neutral mutations (which do not influence the coalescent process *per se*) do not need to be tracked, reducing the computational burden. Adding neutral mutations to the stored tree-sequences once the simulation is complete is very similar to coalescent simulations (e.g. **ms** or **msprime**) so it takes very little time. 

The basic work-flow was as follows:
- Simulate adaptive evolution in structured populations in *SLiM*, storing the genealogy (and the beneficial mutations segregating in the population)
- Sprinkle neutral mutations onto the tree object using *PySLiM* and generate VCF files 
- Perform genome scans on the resulting VCF files using *VCFtools*

The *SLiM* config files are in the [configs/ repository](configs/)

When running SLiM at the command line you can specify parameters for your 
