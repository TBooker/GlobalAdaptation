# A script to ascertain the neutral distribution of Fst in the 2 deme simulations performed in SLiM
import glob, argparse
import pandas as pd
import numpy as np

def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the input files")
	parser.add_argument("--quantiles", 
			type = str,
			help = "specify the path to the file of quantiles",
			default = False)
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "A prefix to add to each of your output files")
	parser.add_argument("--antag", 
			dest = "antag",
			action = 'store_true',
			help = "Are you looking at the local adaptation sims?")

	args = parser.parse_args()
	count = 0
	
	geneStarts = np.array([47500, 152500, 257500, 362500, 467500, 572500, 677500, 782500, 887500, 992500, 1097500, 1202500, 1307500, 1412500, 1517500, 1622500, 1727500, 1832500, 1937500,2042500])
	if args.antag:
			geneStarts = np.array([300000]) 
	quantiles = {} 
	for k in open(args.quantiles):
		x = k.strip().split(' ')
		quantiles[str(int(float(x[0]))) + '_' + x[1]] = float(x[2])

	files = glob.glob(args.input + '/*.n50.weir.fst.gz')


	windowDat = {}
	count = 0
	SNPsPerGene = []
	
	ana = 0
	for i in files:
		ana +=1
	#	if ana == 20: break
	#	print i
		nameStringRaw = i.split('/')[-1].split('_')

		
		if 'windowed' not in i:
			w_size = 1
			
		else:
			w_size = float([k for k in nameStringRaw[-1].split('.') if k.startswith('w')][0][1:])
		if not args.antag:
			shat = nameStringRaw[ nameStringRaw.index('shat') + 1]
			generation = nameStringRaw[ nameStringRaw.index('gen') + 1]
			M = nameStringRaw[ nameStringRaw.index('M') + 1]
			repli = nameStringRaw[ nameStringRaw.index('rep') + 1]
			pA = nameStringRaw[ nameStringRaw.index('pA') + 1]
		elif args.antag:
			print nameStringRaw
			shat = nameStringRaw[ nameStringRaw.index('s') + 1]
			generation = nameStringRaw[ nameStringRaw.index('generation') + 1]
			M = nameStringRaw[ nameStringRaw.index('M') + 1]
			repli = nameStringRaw[ nameStringRaw.index('rep') + 1]
			pA = 1 

		if pA == 0:
			shat = 0

		if w_size == 1.:
			count +=1
			lociList = []

			starts = 1 + geneStarts - (10000 - 5000)/2
			ends   = 1 + geneStarts + 5000 + (10000 - 5000)/2
#			print starts
#			print ends
			
			rep = pd.read_csv(i, sep = '\t').dropna()
			#print w_size, starts
			for s, e in zip(starts, ends):
				sub = rep[(rep['POS'] >= s) &  (rep['POS'] < e)].copy()
				sub['shat'] = shat
				sub['generation'] = generation
				sub['locus'] = s
				sub['rep'] = repli
				sub['M'] = M
				sub['pA'] = pA
				lociList.append(sub)
			SNPsPerGene += lociList
#		if count == 100: break

	print 'analysed', ana,'files'
#	return
	allLoci = []
	for L in SNPsPerGene:
		nbSNPs = L.shape[0]
#		for q in [0.999, 0.9999, 0.99999]:
		for q in [0.99999]:
			locusDat = {}
			thresh = quantiles[ str(1) +'_' + str(q) ]
#			print str(1) +'_' + str(q), thresh
			outliers = L[L['WEIR_AND_COCKERHAM_FST'] > thresh].copy()
			nbOutliers = outliers.shape[0]
			locusDat['quantile'] = q
			locusDat['outliers'] = nbOutliers
			locusDat['SNPs'] = nbSNPs
			locusDat['generation'] = L['generation'].iloc[0]
			locusDat['M'] = L['M'].iloc[0]
			locusDat['rep'] = L['rep'].iloc[0]
			locusDat['locus'] = L['locus'].iloc[0]
			locusDat['shat'] = L['shat'].iloc[0]
			locusDat['pA'] = L['pA'].iloc[0]
			allLoci.append(locusDat)
	locusOutput = pd.DataFrame(allLoci)
	locusOutput.to_csv('SNPs_' + args.output, index = False)
#	return

	return

main() 
