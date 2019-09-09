## python3 script

## A script to get pi for analysis windowsimport glob, argparse
import pandas as pd
import numpy as np
import msprime, pyslim,  sys, subprocess, argparse, glob

def SFS_from_frequencies(frequencies, length,N):
	SFS = [0]*(N+1)
	for i in frequencies:
		if i > N:
			print( "SFS_from_frequencies: Error in your frequencies vector: One of the values is greater than the number of individuals\nThe offending value is: " + str(i) +" and the sample is "+str(N) )
			return
		SFS[i] += 1
	SFS[0] = length - len(frequencies)
	if sum(SFS) < length:
		print("SFS_from_frequencies: Error in your frequencies vector: Fewer items in the SFS than the length of the region")
		return
	if sum(SFS) > length:
		print("SFS_from_frequencies: Error in your frequencies vector: More items in the SFS than the length of the region")
		return
	return SFS
 
def pi(SFS,per_site = True):
	if sum(SFS) ==0: 
		return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	pi = sum([(1.0*i*(N-i)*(SFS[i]))/(binom) for i in range(N) if i != 0])
	if per_site == True:
		return pi/sum(SFS)
	else:
		return pi
 
def recapAndSprinkleTrees(Trees, output, sample = 50, just_vcf = False):

## I'll now graft a coalescent history for the period of shared
## history to the SLiM trees
	print('read trees')

	ts = pyslim.load(Trees)

	N = ts.get_sample_size()

## Recapitate, these parameters are the same as those I've used for the SLiM simulation
	print('recapitate trees')
#	recap = ts
	recap = ts.recapitate(recombination_rate=2.5e-7, Ne=1e4)
	print('recapitated')

## Extract the trees only for a random sample from each population
## If you use the whole population it slows everything down and chews up disc space
#	np.random.seed(10)

	SamplePop1 = np.random.choice(recap.samples()[:int(N/2)],sample, replace = False)
#	print('random sample 1:', SamplePop1)

	SamplePop2 = np.random.choice(recap.samples()[int(N/2):],sample, replace = False)
#	print('random sample 2:', SamplePop2)

	mySample = np.concatenate((SamplePop1, SamplePop2), axis=None)	
	
## Extract the trees only for the random sample
	ts2 = recap.simplify(samples = mySample)

##  Sprinkle mutations onto the trees for the sample
	sprinkled = msprime.mutate(ts2, rate= 2.5e-7,  keep=False)
	print('sprinkle')
	if just_vcf:
		print('writing output')
		with open(output, "w") as vcf_file:
			sprinkled.write_vcf(vcf_file, 2)
		return


def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--trees","-i", 
			required = True,
			dest = "trees",
			type =str,
			help = "The .trees file for the simulation, will detect whether it is zipped or not")
	parser.add_argument("--fst", 
			required = True,
			dest = "fst",
			help = "The fst file for simulation")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			help = "The name of the output file")
	args = parser.parse_args()

	if args.trees.endswith('z'):
		subprocess.Popen(['gunzip', args.trees])
		in_name = args.trees.strip('.gz')
		print ('!!!!!!!!!!!!!!!!', in_name)
	else:
		in_name = args.trees
	
	print (in_name)
	sampleN = 50

	recapAndSprinkleTrees(in_name, args.output+'.vcf',  sample = 50, just_vcf = True)
	#return

	flatten = lambda l: [item for sublist in l for item in sublist]
	
	positions = []
	positionData = {}
	for i in open(args.output+'.vcf', 'r'):
		if i.startswith('#'): continue
		
		x = i.strip().split()
		pos = int(x[1])
		genos = x[9:]
		p1 = flatten([j.split('|') for j in genos[:int(sampleN/2)] ])
		p2 = flatten([j.split('|') for j in genos[int(sampleN/2):] ])
		
		p1_count = p1.count('1')
		p2_count = p2.count('1')		

		variant = [p1_count, p2_count, p1_count + p2_count ]

		positions.append(pos)
		positionData[pos] = variant

	positions = np.array( positions )
	fst = pd.read_csv(args.fst, sep = '\t')

	count = 0
	for index, row in fst.iterrows():
		width = row.BIN_END - row.BIN_START + 1

		snps_in_window = positions[ ( positions >= row.BIN_START ) & ( positions <= row.BIN_END ) ]
		p1_alleles = [positionData[k][0] for  k in snps_in_window if positionData[k][0] != 0]
		p2_alleles = [positionData[k][1] for  k in snps_in_window if positionData[k][1] != 0]
		all_alleles = [positionData[k][2] for  k in snps_in_window]

		p1_sfs = SFS_from_frequencies( p1_alleles, width, sampleN)
		p2_sfs = SFS_from_frequencies( p2_alleles, width, sampleN)
		all_sfs = SFS_from_frequencies( all_alleles, width, sampleN*2)
		
		p1_pi = pi(p1_sfs)
		p2_pi = pi(p2_sfs)
		all_pi = pi(all_sfs)

## Add the calculated pi to the Fst file
		fst.loc[index,'p1_pi'] = pi(p1_sfs)
		fst.loc[index,'p2_pi'] = pi(p2_sfs)
		fst.loc[index,'all_pi'] = pi(all_sfs)
	
	
	fst.to_csv(args.output, index = False)

	subprocess.Popen(['rm', args.output+'.vcf'])

	subprocess.Popen(['gzip', in_name])
	

main()
