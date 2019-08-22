
Simuating of a model of parapatry (or isolation with migration)
------

I simulated adaptive evolution occurring in a populations structured according to a model of isolation with migration - or a two-deme island model. I performed the simulations in *SLiM* v3.2 making use of the tree-sequence recording feature. This feature records the coalescent history of your simulated population as a tree-sequence object. It is computationally efficient as neutral mutations (which do not influence the coalescent process *per se*) do not need to be tracked, reducing the computational burden. Adding neutral mutations to the stored tree-sequences once the simulation is complete is very similar to coalescent simulations (e.g. **ms** or **msprime**) so it takes very little time. 

The *SLiM* config files are in the [configs/ repository](configs/)

When running SLiM at the command line you can specify parameters for your simulation as arguments. I make use of that to make parellisation easier. For example, you could use GNU parallel to perform 200 relicates each of simulations assuming migration rates of *Nm* = 1 or 10, with 1% or 0.1% of new selected mutations having a mean advantageous effect of *Nesa* = 100 or 10 as follows:

```
parallel "~/path/to/slim -d s_hat={1} -d M={2} -d pA={3} -d REP={4}" ::: 100 10 ::: 1 10 ::: 0.01 0.001 ::: $(seq 1 200)
```

Or, like I did, you could run your simulations on a nice big cluster (props to ComputeCanada) and submit an array job or the like.

Analysis workflow
------

The basic work-flow for analysing the simulation data was as follows:
- Sprinkle neutral mutations onto the tree object using *PySLiM* and generate VCF files
- Perform genome scans on the resulting VCF files using *VCFtools*
- Extract the segregating sites from the tree files and store as a new file

I also performed a big set of neutral simulations for both migration rates I tested. These were the same as the selection ones, except that I said that 


Making Figures 1, 2 and S1
------

Other than the "demo plots" (Figures 1, 2 and S1), population sizes were set to 5,000 dipoids per deme (i.e. *N* = 10,000). The difference in size was due to the fact that I was recording the population state every XX generations, so I ended up writing a large number of large files. I reduced the population size, but kept the population-scaled parameters the same, to reduce the  
