# A script to ascertain the neutral distribution of Fst in the 2 deme simulations performed in SLiM
import glob, argparse
import pandas as pd
import numpy as np
from scipy.stats import binom

def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the input files. MAKE SURE IT'S A DIRECTORY OF NEUTRAL SIMULATIONS")
	parser.add_argument("--probs", 
			type = str,
			help = "specify the file of proportion of outliers",
			default = False,
			required = True)
	parser.add_argument("--quantiles", 
			type = str,
			help = "specify the file of quantiles for specfiying outliers",
			default = False,
			required = True)
	parser.add_argument("--window", 
			type = int,
			required = True,
			help = "Specify the window size that you want to analyse",
			default = False)

	args = parser.parse_args()
	count = 0
	
	
	geneStarts = np.array([300000])

	quantiles = {} 

	for k in open(args.quantiles):
		x = k.strip().split(' ')
		quantiles[str(int(float(x[0]))) + '_' + x[1]] = float(x[2])


	probs = {} 

## Read in the outlier proportions obtained from neutral simulations
	for k in open(args.probs):
		x = k.strip().split(' ')
		if int(x[1]) != args.window:
			print 'the window size you specified does not match the neutral data'
			return
		probs[x[2]] = float(x[0])

## Make dictionaries of the top-candidate thresholds for each number of SNPs from 1 to 2000
	thresh_999 = {i+1: binom.ppf(0.9999, i+1, probs['0.999']) for i in range(2000)}
	thresh_9999 = {i+1: binom.ppf(0.9999, i+1, probs['0.9999']) for i in range(2000)}
	thresh_99999 = {i+1: binom.ppf(0.9999, i+1, probs['0.99999']) for i in range(2000)}

	windowDat = {}
	files = 0
	fileList =  []


	for i in glob.glob(args.input + '/*.n50.weir.fst.gz'):
		print i
		nameStringRaw = i.split('/')[-1].split('_')

#		rep = nameStringRaw[ nameStringRaw.index('rep') + 1]
#		shat = nameStringRaw[ nameStringRaw.index('shat') + 1]
#		generation = nameStringRaw[ nameStringRaw.index('gen') + 1]
#		M = nameStringRaw[ nameStringRaw.index('M') + 1]
#		pA = nameStringRaw[ nameStringRaw.index('pA') + 1]
		
#		if generation == '4000':
#			continue
#		if pA == 0:
#			shat = 0

		files +=1

		fileList.append(i)
		
		p_999 = []
		p_9999 = []
		p_99999 = []

		starts = 1 + geneStarts - (args.window - 5000)/2
		repDF = pd.read_csv(i, sep = '\t').dropna()
		
		print repDF[repDF['WEIR_AND_COCKERHAM_FST'] == repDF['WEIR_AND_COCKERHAM_FST'].max()]

		print repDF[repDF['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.99999' ] ].copy()
		
#		print repDF
		continue
		for s in  starts:
			snps = repDF[(repDF['POS'] >= s)&(repDF['POS'] <= s + args.window)]
#			outliers = snps[snps['WEIR_AND_COCKERHAM_FST'] > 
#			print snps
			numSNPs = len(snps)

			print snps[snps['WEIR_AND_COCKERHAM_FST'] > 0.4 ].copy()

			outliers_999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.999' ] ].copy()
			outliers_9999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.9999' ] ].copy()
			outliers_99999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.99999' ] ].copy()

			p_999 = len(outliers_999)
			p_9999 = len(outliers_9999)
			p_99999 = len(outliers_99999)

			if len(outliers_99999) > 0:
#			if len(outliers_99999) > thresh_99999[numSNPs]:
				print '!', s 
				print outliers_999

				print len(outliers_999), thresh_999[numSNPs]
				print len(outliers_9999), thresh_9999[numSNPs]
				print len(outliers_99999), thresh_99999[numSNPs]
				print
			


		continue
		datPoint = {'shat':shat,
		'w_size':args.window,
		'generation':generation,
		'M':M,
		'pA': pA,
		'o0.99999':len(outliers_99999),
		'o0.9999':len(outliers_9999),
		'o0.999':len(outliers_999),
		'o0.99999_irs':len(outliers_99999_IRS),
		'o0.9999_irs':len(outliers_9999_IRS),
		'o0.999_irs':len(outliers_999_IRS)}
		windowDat[i] = datPoint

#		if files == 200: break
		
	print sum(p_999) / len(p_999), args.window, 0.999
	print sum(p_9999) / len(p_9999), args.window, 0.9999
	print sum(p_99999) / len(p_99999), args.window, 0.99999
	



main()
