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
			default = False,
			required = True)
	parser.add_argument("--window", 
			type = int,
			required = True,
			help = "Specify the window size that you want to analyse",
			default = False)
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "A prefix to add to each of your output files")

	args = parser.parse_args()
	count = 0
	
	geneStarts = np.array([47500, 152500, 257500, 362500, 467500, 572500, 677500, 782500, 887500, 992500, 1097500, 1202500, 1307500, 1412500, 1517500, 1622500, 1727500, 1832500, 1937500,2042500])

	quantiles = {} 

	for k in open(args.quantiles):
		x = k.strip().split(' ')
		quantiles[str(int(float(x[0]))) + '_' + x[1]] = float(x[2])

	windowDat = {}
	count = 0
	files = 0
	fileList =  []

	for i in glob.glob(args.input + '/*.windowed.weir.fst.gz'):
#		print i
		nameStringRaw = i.split('/')[-1].split('_')

		if 'windowed' not in i:
			w_size = 1
		else:
			w_size = float([k for k in nameStringRaw[-1].split('.') if k.startswith('w')][0][1:])
		if w_size != args.window: continue

		
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
		

		starts = 1 + geneStarts - (w_size - 5000)/2

#		print i.split('x')
		alleleFreqs = pd.read_csv(i.split('x')[0] +'xt.segSites.txt')
		alleleFreqs['BIN_START'] = [starts[p] for p in np.array(alleleFreqs['locus'])-1 ] 
#		print alleleFreqs
		
			
		rep = pd.read_csv(i, sep = '\t').dropna()

		sub = rep[rep['BIN_START'].isin(starts)].copy()
#		print sub
		
		result = pd.merge(sub, alleleFreqs, how='left', on='BIN_START')
		result.to_csv('YOYO.csv')
		result['Ns'] = result['s'] * 10000
		result = result[result['Ns'] > 1].copy()
#		print result
		result.drop_duplicates('BIN_START', inplace=True)
#		print result
		
		result = result[(result['q'] >= 0.2) & (result['q'] <= 0.8)].copy()

		
		outliers_999 = sub[sub['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.999) ] ].copy()
		outliers_9999 = sub[sub['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.9999) ] ].copy()
		outliers_99999 = sub[sub['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.99999) ] ].copy()

		outliers_999_IRS = result[result['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.999) ] ].copy()
		outliers_9999_IRS = result[result['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.9999) ] ].copy()
		outliers_99999_IRS = result[result['WEIGHTED_FST'] > quantiles[ str(int(w_size)) +'_' + str(0.99999) ] ].copy()


		if files == 1:
			header = [shat, w_size, M, pA]
		datPoint = {'shat':shat,
		'w_size':w_size,
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
	
#		break
#	print files

	fileListArray = np.array(fileList)
	listLen = len(fileList)
	
	
	pointEst = pd.DataFrame([windowDat[k] for k in fileList])
	pointEst99999 = float(pointEst['o0.99999'].sum()) / (listLen*20)
	pointEst9999 = float(pointEst['o0.9999'].sum()) / (listLen*20)
	pointEst999 = float(pointEst['o0.999'].sum()) / (listLen*20)
	
	pointEst99999_irs = float(pointEst['o0.99999_irs'].sum()) / (listLen*20)
	pointEst9999_irs = float(pointEst['o0.9999_irs'].sum()) / (listLen*20)
	pointEst999_irs = float(pointEst['o0.999_irs'].sum()) / (listLen*20)
	output =  [header + ['est', pointEst999, pointEst9999, pointEst99999, pointEst999_irs, pointEst9999_irs, pointEst99999_irs] ]

	for b in range(2000):
		bootSample = np.random.choice( fileListArray, size= listLen, replace=True)
		bootDF = pd.DataFrame([ windowDat[k] for k in bootSample])
		prop99999 = float(bootDF['o0.99999'].sum()) / (listLen*20)
		prop9999 = float(bootDF['o0.9999'].sum()) / (listLen*20)
		prop999 = float(bootDF['o0.999'].sum()) / (listLen*20)
		prop99999_irs = float(bootDF['o0.99999_irs'].sum()) / (listLen*20)
		prop9999_irs = float(bootDF['o0.9999_irs'].sum()) / (listLen*20)
		prop999_irs = float(bootDF['o0.999_irs'].sum()) / (listLen*20)
		output.append( header + [b, prop999, prop9999, prop99999, prop999_irs, prop9999_irs, prop99999_irs] )
	pd.DataFrame(output, columns = ['shat','w_size','M','pA','rep','o999','o9999','o99999', 'o999_irs','o9999_irs','o99999_irs']).to_csv(args.output, index = False)
	
	return

main() 
