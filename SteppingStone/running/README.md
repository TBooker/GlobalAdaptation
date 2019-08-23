Running the stepping-stone simulations
======

Here I lay out the steps with script names to try and make it obvious which scripts are used when. 

Initially, you make a bunch of neutral simulations. I made 100,000.

The script [SSdrifter.py](../bin/SSdrifter.py) simulates neutral burnin.

For example, you could run 100,000 reps each seeded with an initial allele frequency of 0.5 in one randomly chosen deme.

```
parallel "python bin/SSdrifter.py -k 500 -N 5000 -m 0.666 -p 0.5 -d 100000 -o neutralRuns/{}.txt" ::: $(seq 1 100000)
```

I then moved all the neutral runs into a new directory, giving each replicate a name that contains their allele frequency:

```
mkdir neutralRuns_formatted

python bin/renameSims.py --direc neutralRuns/ -N 5000 --output neutralRuns_formatted

```

Now you have a database of neutral runs, you can make a second database of selected ones, running 500 replicates for each value of *s* from 0.1 to 0.0001 in increments of 0.0001:
```
for i in $(seq 0.1000 -0.0001 0.0001)
do
    echo $i
    mkdir selectionDataBase/$i
    python bin/runSelectedRuns.py --direc neutralRuns_formatted/ -s $i -m 0.666 -c 0.0001 -o selectionDataBase/$i/ --batch --reps 500 --nproc 45
done
```
The ```nproc``` argument allows you to run the simulations in parallel using Python's multiprocessing module.

**Check paths in the individual files!!**

Now you have the database of selected runs, you can do the DFE sampling like so:

```
python bin/sampleFromDFE.py --direc neutralRuns_formatted/ --generation 100000 --output my_run.csv --Ua 1 --exp_mean 0.02 --sample_T 99999 --distance 400

In this example, I chose to simulate an advantageous mutation rate of ```Ua``` = 1, an exponential DFE with mean ```exp_mean```, reporting Fst for populations that were ```distance``` = 400 demes apart.

I did this 30 times for each parmeter combination I tested. I then plotted the resulting data using the Rscript [my_plot.R](../bin/myPlot.R)