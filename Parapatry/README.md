
Simulating of a model of parapatry (or isolation with migration)
======

I simulated adaptive evolution occurring in a populations structured according to a model of isolation with migration - or a two-deme island model. I performed the simulations in *SLiM* v3.2 making use of the tree-sequence recording feature. This feature records the coalescent history of your simulated population as a tree-sequence object. It is computationally efficient as neutral mutations (which do not influence the coalescent process *per se*) do not need to be tracked, reducing the computational burden. Adding neutral mutations to the stored tree-sequences once the simulation is complete is very similar to coalescent simulations (e.g. **ms** or **msprime**) so it takes very little time. 

Additionally, I used the 'recapitation' feature of PySLiM, which uses the **msprime** to perform a coalescent simulation back in time from the start of the forward simulation to ensure coalescence. For this reason, I added 1,000 generations of 

*I've tried to make this repository obvious, but if it's not and you want me to clarify anything, don't hesitate to get in touch!*

Running the SLiMulations
------

Specific details of the simulations are described in the Methods section of the manuscript, but briefly, I simulated a two-deme model, with symmetrical migration but constant population size. An initial population of *N* = 10,000 individuals was split into two equally sized demes of 5,000 individuals. Simulated chromosomes had stretches of functional or "gene-like" sequence that were 5,000 bp long. These experienced mutations with fitness effects drawn from an exponential distribution of fitness effects (DFE). Recombination and mutation rates were uniform across the chromomsome.

For the purposes of comparison, I also ran a set of simulations modelling antagonistic pleiotropy, or spatially varying selection. In this case, the simulated chromosome was shortened and a single copy of a selected allele was added to the popualtion at a particular time. 

The *SLiM* config files are in the [configs/](configs/) repository.

When running SLiM at the command line you can specify parameters for your simulation as arguments. I make use of that to make parellisation easier. For example, you could use GNU parallel to perform 200 relicates each of simulations assuming migration rates of *Nm* = 1 or 10, with 1% or 0.1% of new selected mutations having a mean advantageous effect of *Nesa* = 100 or 10 as follows:

```
parallel "~/path/to/slim -d s_hat={1} -d M={2} -d pA={3} -d REP={4} /home/booker/configs/recap_template.txt" ::: 100 10 ::: 1 10 ::: 0.01 0.001 ::: $(seq 1 200)
```

Or, like I did, you could run your simulations on a nice big cluster (props to ComputeCanada) and submit an array job or the like. I used Cedar on Compute Canada, which runs the SLURM scheduler. An example submission script might look like:
```
#!/bin/bash
#SBATCH --array=1-200 # 200 Jobs
#SBATCH --job-name=TomSlim
#SBATCH --time=05:00:00
#SBATCH --mem=2500
#SBATCH --output=TomSlim.%A%a.out
#SBATCH --error=TomSlim.%A%a.err

/home/booker/bin/build/slim -d M=1 -d REP=$SLURM_ARRAY_TASK_ID -d s_hat=0.001 -d pA=0.001 /home/book
er/configs/recap_template.txt
```
I don't know if this would work under other scheduling software.


Analysis workflow
------

The basic work-flow for analysing the simulation data was as follows:
- Sprinkle neutral mutations onto the tree object using **PySLiM** and generate VCF files
- Perform genome scans on the resulting VCF files using **VCFtools**
- Extract the segregating sites from the tree files and store as a new file - note that only because I used treeSeq, there were only advantageous mutations in the SLiM portion of the simulations

The analysis I performed for each simulation replicate is in the file [bin/run.sh](bin/run.sh)

I also performed a big set of neutral simulations for both migration rates I tested. These were the same as the selection ones, except that no beneficial mutations occurred. From these simulations, I examined the distribution of Weir and Cockerham's estimator *Fst* calculated using **VCFtools**.


Making Figures 1, 2 and S1
------

Other than the "demo plots" (Figures 1, 2 and S1), population sizes were set to 5,000 dipoids per deme (i.e. *N* = 10,000). The difference in size was due to the fact that I was recording the population state every XX generations, so I ended up writing a large number of large files. I reduced the population size, but kept the population-scaled parameters the same, to reduce the file sizes I reduced population size.

I basically followed the analysis workflow outlined above and made Manhattan plots of *Fst*, colouring the points that passed the significance threshold defined from neutral simulations (Figures 1 and S1). The script to make those plots is [ref_to_plotting_script.R]().

Figure 2 was a little bit more involved. To make that figure, I divided the dataset into seven time slices and asked what the maximum individual *Fst* datapoint was during each of the slices. For each of the *Fst* hits, I obtained the allele frequency for linked beneficial mutations. The scripts for doing this are [ref_to_script.py]() 


Making Figures 3, S2 and S3
------

For a given parameter set, I had 2,000 simulation replicates. For each replicate, I had performed an *Fst* genome scan, so from each I extracted the analysis windows which overlapped "gene-like" sequences. Using the *Fst* cut-off values determined from the neutral simualtions, I classified windows as being either outliers or not. I then bootstrap the data by simulation replicate, generating samples and calculating the proportion of outlier analysis windows observed in that particular replicate. I did that 2,000 times for each dataset. The results of each of those were plotted as the violin plots shown in Figure 3. The R script [ADD PLOTTING SCRIPT](bin/PlottingScript.R), and the files that it links to, should generate Figure 3.

Figures S2 and S3 were a bit more complicated because I had to 

*Note that this approach may have been slightly conservative as I only considered 20 analysis windows per simulation replicate. There may be analysis windows that exhibit elevated Fst as a result of global sweeps*

Making Figure 4
------

To make the figure showing genetic diversity and *Fst*, I just extracted the top 100 analysis windows from the parameter set specified in the figure legend and calculated nucleotide diversity for the same windows. I repeated this with simulations modelling antagonistic pleiotropy.

The script that did that is [ref_to_script.py]().

Calculating ALPHA - Table S1
------
In order to calculate ALPHA, the proportion of substitutions fixed by positive selection, I did an extra set of simulations [link_to_configs](). I did these because ALPHA is **very** sensitive to what happens to nonsynonymous mutations that aren't beneficial. To calculate ALPHA, you need to know dN and dS, and you need to break dN into the number of substitutions attributable to positive selection and othewise. I set a DFE for harmful mutations at a proportion of sites (1 - pA) assuming the one for *Drosophila* estimated by Loewe and Charlesworth (2006) - a gamma distribution with mean -200 and shape parameter = 0.3.


The substitutions which occured in these simulations were then used to calculate ALPHA. The script [ref_to_script.py]() does that. If you are playing around with these, keep and eye on the mutation id variables from SLiM, my scripts are probably pretty fragile with respect to those. 
