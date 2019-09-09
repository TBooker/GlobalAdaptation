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
	parser.add_argument("--neutral", 
			action = "store_true",
			help = "Use this flag when analysing neutral sites",
			default = False)
	parser.add_argument("--quantiles", 
			type = str,
			help = "specify the path to the file of quantiles",
			default = False)
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "A prefix to add to each of your output files")

	args = parser.parse_args()
	count = 0
	
	geneStarts = np.array([300000])
	quantiles = {} 
	for k in open(args.quantiles):
		x = k.strip().split(' ')
		quantiles[str(int(float(x[0]))) + '_' + x[1]] = float(x[2])
	#print quantiles

	
	windowDat = {}
	count = 0
	SNPsPerGene = []
	files = 0
	for i in glob.glob(args.input + '/*weir.fst.gz'):
		print i
		files +=1
		nameStringRaw = i.split('/')[-1].split('_')
		if 'windowed' not in i:
			w_size = 1
			
		else:
			print nameStringRaw
			w_size = float([k for k in nameStringRaw[-1].split('.') if k.startswith('w')][0][1:])
		
		sRaw = [j for j in nameStringRaw if j.startswith('s')][0]
		print nameStringRaw.index(sRaw)+1
		
		s = float(nameStringRaw[nameStringRaw.index(sRaw)+1])

		generation = nameStringRaw[nameStringRaw.index([j for j in nameStringRaw if j.startswith('generation')][0])+1]
		M = nameStringRaw[nameStringRaw.index([j for j in nameStringRaw if j.startswith('M')][0])+1]
		print generation


		if w_size != 1.:
			count +=1
			starts = 1 + geneStarts - (w_size - 5000)/2

			
			rep = pd.read_csv(i, sep = '\t').dropna()
			#print w_size, starts
			sub = rep[rep['BIN_START'].isin(starts)].copy()
			sub['s'] = s	
#			print(shat)
			sub['w_size'] = w_size
			sub['generation'] = generation
			sub['M'] = M
			if w_size not in windowDat.keys():
				windowDat[w_size] = [sub]
			else:
				windowDat[w_size].append(sub)
		elif w_size == 1.:
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
				sub['s'] = s
				sub['generation'] = generation
				sub['locus'] = s
				sub['M'] = M
				lociList.append(sub)
			SNPsPerGene += lociList
#		if count == 100: break

#	print 'analysed', files,'files'
#	return
	allLoci = []
	for L in SNPsPerGene:
		nbSNPs = L.shape[0]
		for q in [0.999, 0.9999, 0.99999]:
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
			locusDat['s'] = L['s'].iloc[0]
			allLoci.append(locusDat)
#			print locusDat
	locusOutput = pd.DataFrame(allLoci)
	locusOutput.to_csv('SNPs_' + args.output, index = False)
#	return
	windowedAnalyses = []
	for i in  windowDat.keys():
		stuff = pd.concat(windowDat[i])
		dataEntry = {}
		nbWindows = stuff.shape[0]
		
		for q in [0.999, 0.9999, 0.99999]:

			thresh = quantiles[ str(int(i)) +'_' + str(q) ]
			temp = stuff[stuff['WEIGHTED_FST'] > thresh].copy()
			temp['q'] = q
			windowedAnalyses.append(temp)
	windowOutput = pd.concat(windowedAnalyses)

	windowOutput.to_csv('windows_' + args.output, index = False)

	return

main() 
