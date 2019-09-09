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
			help = "The directory containing the input files. MAKE SURE IT'S A DIRECTORY OF NEUTRAL SIMULATIONS")
	parser.add_argument("--quantiles", 
			type = str,
			help = "specify the path to the file of quantiles",
			default = False,
			required = True)
	parser.add_argument("--window", 
			type = int,
			required = True,
			help = "Specify the window size that you want to analyse",
			default = False)

	args = parser.parse_args()
	count = 0
	
	geneStarts = np.array([47500, 152500, 257500, 362500, 467500, 572500, 677500, 782500, 887500, 992500, 1097500, 1202500, 1307500, 1412500, 1517500, 1622500, 1727500, 1832500, 1937500,2042500])

	quantiles = {} 

	for k in open(args.quantiles):
		x = k.strip().split(' ')
		quantiles[str(int(float(x[0]))) + '_' + x[1]] = float(x[2])

#	print quantiles
#	return
	windowDat = {}
	files = 0
	fileList =  []

	p_999_num = 0
	p_999_den = 0
	p_9999_num = 0
	p_9999_den = 0
	p_99999_num = 0
	p_99999_den = 0
	
	for i in glob.glob(args.input + '/*.n50.weir.fst.gz'):
#		print i
		nameStringRaw = i.split('/')[-1].split('_')

		rep = nameStringRaw[ nameStringRaw.index('rep') + 1]
		shat = nameStringRaw[ nameStringRaw.index('shat') + 1]
		generation = nameStringRaw[ nameStringRaw.index('gen') + 1]
		M = nameStringRaw[ nameStringRaw.index('M') + 1]
		pA = nameStringRaw[ nameStringRaw.index('pA') + 1]
		
		if generation == '4000':
			continue
		if pA == 0:
			shat = 0

		files +=1

		fileList.append(i)
		

		starts = 1 + geneStarts - (args.window - 5000)/2

		repDF = pd.read_csv(i, sep = '\t').dropna()
		
		for s in  starts:
			snps = repDF[(repDF['POS'] >= s)&(repDF['POS'] <= s +args.window)]
#			outliers = snps[snps['WEIR_AND_COCKERHAM_FST'] > 
#			print snps
			numSNPs = len(snps)
#			print snps[snps['WEIR_AND_COCKERHAM_FST'] > 0.389138 ].copy()

			outliers_999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.999' ] ].copy()
			outliers_9999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.9999' ] ].copy()
			outliers_99999 = snps[snps['WEIR_AND_COCKERHAM_FST'] > quantiles[ '1_0.99999' ] ].copy()
			
			if len(outliers_999) > 0:
				p_999_num += len(outliers_999)
				p_999_den += numSNPs
			if len(outliers_9999) > 0:
				p_9999_num += len(outliers_9999)
				p_9999_den += numSNPs
			if len(outliers_99999) > 0:
				p_99999_num += len(outliers_99999)
				p_99999_den += numSNPs
			
		if files == 200: break
		
	print float(p_999_num) / p_999_den, args.window, 0.999
	print float(p_9999_num) / p_9999_den, args.window, 0.9999
	print float(p_99999_num) / p_99999_den, args.window, 0.99999
	



main()
