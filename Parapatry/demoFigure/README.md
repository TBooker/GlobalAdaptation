Plotting Fst and allele frequency over time
======

Here's the workflow I used for making Figures 1, 2 and S1.

**NOTE** There is mix of Python2.7 and Python3.x used in this project. Remember that I'm using SLiM v3.x


I ran a set of 250 simulations, acting like each one was treated like a chromosome. 
The simulated chromosomes I was analysing are 2Mbp in size, so 250 * 2Mbp gives a genome size of around 500Mbp of similar to a lot of eukaryotes

Here's an example of how you might run the workflow for 5 simulations making use of GNU Parallel (again, you could run them in serial or on a cluster, whatever works for you):

Make a new directory, run the simulations using ```SLiM```,  then use the tree files to make VCF files, then run a VCFtools genome scan on the resulting data
```
mkdir M1_p0.0001
cd M1_p0.0001
parallel "slim -d M=1 -d REP={} -d pA=0.0001 ../../configs/timeCourse.slim" ::: $(seq 1 5)
cd ../
parallel "sh bin/run.sh {}" ::: $(ls -d M1_p0.0001/*_.txt)
```
*Keep an eye on all the paths to the individual files.*


Parse the allele frequencies for seg sites in all simulation replicates into a single file:
```
python bin/parseAlleleFreqs.py M1_p0.0001 M1_p0.0001.alleleFreqs.csv
```

Take all of the Fst scans and combine them into one file for each generation:
```
python bin/parseOutput.py M1_p0.0001 M1_p0.0001.alleleFreqs.csv
```

The previous two Python scripts are pretty straightforward. But the next one is a bit more complicated, it takes the files and:
1. Splits the data set into a number of chunks (specified using the ```num``` argument)
2. Finds the top *Fst* value in each chunk
3. Returns a file with the top 
4. Finds sites that are segregating at the loci that exhibit Fst peaks and returns their frequencies over time in a separate file

```
python python bin/getData.py -i M1_p0.0001_files/ --alleleFreq M1_p0.0001.alleleFreqs.csv  --output TESTTEST --num 7
```
I tried to comment that script a fair bit to make it clear what was going on.

With the resulting files, I make plots using the following R scripts:

[bin/demoPlot_S1.R](bin/demoPlot_S1.R) 
[bin/demoPlot_3part.R](bin/demoPlot_3part.R) 



