Stepping Stone Simulations
======

I approximated the spread of beneficial mutations in continuous environments using a one-dimensional stepping-stone model with many demes. 

Doing a full simulation with neutral and selected sequence, randomly introduced mutations and recombination like I did for the paraptry case wasn't feasible, so I opted for a two-locus model instead. 

Check out the manuscript for full details of the simulation, but to understand the structure of what I did a bit of explanation is needed. Simulating the mutation and subsequent equilibration of neutral variants in a stepping-stone population of 2.5 million individuals would take a very long time, so I took a little bit of a short cut. The initial paper describing the stepping-stone model (Kimura and Weiss 1964) showed that the correlation in allele frequencies between adjacent demes can be very high. I wanted to simulate popualtions with a distribution of *Fst* at neutral sites which was similar to what is seen in real organisms. Often within species pair-wise *Fst* is less than 0.1, so I played around a bit and settled on a migration rate of *m* = 0.[666](https://www.youtube.com/watch?v=WxnN05vOuSM). Weiss and Kimura (1964) showed that the correlation in allele frequencies between immediately adjacent demes with such a migration rate is  

The population consisted of *k* = 500 demes each with *N* = 5,000 haploid individuals.






Here are the scripts for performing the stepping-stone simulations.
