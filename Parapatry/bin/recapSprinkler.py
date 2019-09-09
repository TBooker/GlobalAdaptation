# Let's sprinkle mutations onto the subPop tree
import msprime, pyslim,  sys, subprocess, argparse, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

def SNPwindows(SNPs, width, jump):
	data = []
#	print(SNPs)
#	print([k for k in range(0, SNPs.shape[0], jump)])
	
	
	for i in range(0, SNPs.shape[0], jump):
#		print(i, i+width, SNPs.shape[0])
		if i + width > SNPs.shape[0]-1: continue
		temp = SNPs.iloc[i:i+width ,]
		estimate = temp.WCnum.sum()/ temp.WCden.sum()
		start = SNPs.iloc[int(i)].pos
		end = SNPs.iloc[int(i+width)].pos
		mid = SNPs.iloc[int(i+width/2)].pos
		temp = {'start':start, 'end':end, 'mid':mid, 'WC':estimate} 
		data.append(temp)
	return data

def PhysWindows(SNPs, width, jump):
	data = []

	windows = range(0, 2047499, jump)
	for start in windows:
		end = start + width -1
		temp = SNPs[(SNPs['pos']>= start) & (SNPs['pos'] <= end)]
		if temp.shape[0] == 0: 
			continue
		estimate = temp.WCnum.sum()/ temp.WCden.sum()
		mid = (start + end)/2
		temp = {'start':start, 'end':end, 'mid':mid, 'WC':estimate} 
		data.append(temp)

	return data


def wc2pops(p1,p2):
	p_bar = (p1 + p2)/2.
	s2 = (p1 - p_bar)**2 + (p2 - p_bar)**2 
	if p1 == p2 == 0.0 or p1 == p2 == 1.0:
		return 'Na'
	if s2 + p_bar*(1-p_bar) == -2:
		print(p1,p2) 
	return s2 , s2/2 + p_bar*(1-p_bar)


 
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
	print('random sample 1:', SamplePop1)

	SamplePop2 = np.random.choice(recap.samples()[int(N/2):],sample, replace = False)
	print('random sample 2:', SamplePop2)

	mySample = np.concatenate((SamplePop1, SamplePop2), axis=None)	
	
## Extract the trees only for the random sample
	ts2 = recap.simplify(samples = mySample)

##  Sprinkle mutations onto the trees for the sample
	sprinkled = msprime.mutate(ts2, rate= 2.5e-7, keep=False)
	print('sprinkle')
	if just_vcf:
		print('writing output')
		with open(output, "w") as vcf_file:
			sprinkled.write_vcf(vcf_file, 2)
		return
		
	else:
		with open(Trees+'.vcf', "w") as vcf_file:
			sprinkled.write_vcf(vcf_file, 1)

	WCnum_tot = 0
	WCdenom_tot = 0

	#SNPs = OrderedDict()
	SNPs = []
	for i in open(Trees+'.vcf', "r"):
		if i.startswith('#'):
			if i.startswith('#CHROM'):
				header = {}
				x = i.strip().split()
				for j in range(len(x)):
					header[x[j]] = j
			continue
		x = i.strip().split()
		pos = int(x[1])
	
		alleleFreqs = x[header['FORMAT']+1:]
		pop1 = list(map(int, alleleFreqs[:int(sample/2)]))
		pop2 = list(map(int, alleleFreqs[int(sample/2):]))
		p1 = sum(pop1)/len(list(pop1))
		p2 = sum(pop2)/len(list(pop2))

		if p1+p2 == 2: continue
		if sum(pop1 + pop2) > 2*sample:
			continue
		if 2 in pop1:
			continue
	
	## If the segregating site is at a frequency of 1 in the sample, o not include it in the analysis
		
		WCnum, WCden = wc2pops(p1,p2)

		SNPs.append({'pos':pos,'p1':p1,'p2':p2, 'WCnum':WCnum, 'WCden':WCden, 'WC':WCnum/WCden})

#		if WCnum/WCden > 0.5:
#			print(p1, p2,WCnum/WCden)
			
		WCnum_tot += WCnum
		WCdenom_tot += WCden
	if not just_vcf:
		subprocess.Popen(['rm', Trees+'.vcf'])
	return SNPs


def main():
	parser = argparse.ArgumentParser(description="This script takes Trees files from SLiMs and uses msprime/pyslim to add mutations and writes a VCF file")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the output files")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "Name the output file, give a prefix")
	parser.add_argument("--sample","-s", 
			required = False,
			dest = "sample",
			type = int,
			help = "How many individuals do you want to sample from each deme",
			default = 50)
	parser.add_argument("--just_vcf", 
			required = False,
			dest = "just_vcf",
			action = 'store_true',
			help = "Just make a VCF and then quit",
			default = False)			
	args = parser.parse_args()

	if args.just_vcf:
		recapAndSprinkleTrees(args.input, args.output, sample = args.sample, just_vcf = args.just_vcf)
	#	SNPs = recapAndSprinkleTrees(args.input, None, sample = args.sample)
	#	pd.DataFrame(SNPs).to_csv('snps'+args.output, index = False, sep = '\t')
		return
	else:
		SNPs = recapAndSprinkleTrees(args.input, None, sample = args.sample)
		
	#pd.DataFrame(SNPs).to_csv(args.output, index = False)
	byPhys = PhysWindows(pd.DataFrame(SNPs),1000,500)
	pd.DataFrame(byPhys).to_csv(args.output, index = False)
	
#	bySNPs = SNPwindows(pd.DataFrame(SNPs),500,5)
#	pd.DataFrame(bySNPs).to_csv(args.output, index = False)
	return



if __name__ == '__main__':
    main()
