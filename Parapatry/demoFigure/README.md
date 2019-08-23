Plotting Fst and allele frequency over time
======

Here's the workflow I used for making Figures 1, 2 and S1.

**NOTE** There is mix of Python2.7 and Python3.x used in this project. Remember that I'm using SLiM v3.x


I ran a set of 250 simulations, acting like each one was treated like a chromosome. 
The simulated chromosomes I was analysing are 2Mbp in size, so 250 * 2Mbp gives a genome size of around 500Mbp of similar to a lot of eukaryotes

Here's an example of how you might run the workflow for 5 simulations making use of GNU Parallel (again, you could run them in serial or on a cluster, whatever works for you):

```
mkdir M1_p0.0001
cd M1_p0.0001
parallel "slim -d M=1 -d REP={} -d pA=0.0001 ../../configs/timeCourse.slim" ::: $(seq 1 5)
cd ../
parallel "sh bin/run.sh {}" ::: $(ls -d M1_p0.0001/*_.txt)
```


Parse the allele frequencies for seg sites in all sim. replicates into a single file
```
python bin/parseAlleleFreqs.py M1_p0.0001 M1_p0.0001.alleleFreqs.csv
```

Take all of the Fst scans and combine them into one file per generation
```
python bin/parseOutput.py M1_p0.0001 M1_p0.0001.alleleFreqs.csv
```


```
python python bin/getData.py -i M1_p0.0001_4_files/ --alleleFreq M1_p0.0001_4.alleleFreqs.csv  --output TESTTEST --num 7
```