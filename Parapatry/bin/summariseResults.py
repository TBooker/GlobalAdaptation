import recapSprinkler, glob, re
import pandas as pd
from multiprocessing import Pool

def analyseReplicate(trees_file):
	rep = trees_file.split('/')[-1].split('.')[0].split('_')[1]
	gen = trees_file.split('/')[-1].split('.')[1].split('_')[1]
	SNPs = recapSprinkler.recapAndSprinkleTrees(trees_file, None)
	return [pd.DataFrame(SNPs), rep, gen]

def main():
	
	windowSize = 100
	Step =  50
	
	for i in glob.glob('*/'):
		if i.endswith('.csv/'):
			continue
		else:
			print(i)

			numbers = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", i)
			pA = float(numbers[0])
			M = float(numbers[1])
			s_hat = float(numbers[2])
		if s_hat ==0: continue
		if M !=5.: continue
		p = Pool(3)
		files_to_analyse = glob.glob(i+'/*txt')

		data = p.map(analyseReplicate,files_to_analyse)
		reps = []
		for datum in data:
			rep = datum[1]
			gen = datum[2]
			analysis = pd.DataFrame( recapSprinkler.PhysWindows(datum[0],100000,10000) )
			
			#SNPwindows(datum[0], windowSize, Step) )
			analysis['gen'] = gen
			analysis['rep'] = rep
			reps.append(analysis)
		final = pd.concat(reps)
		final['s_hat'] = s_hat
		final['pA'] = pA
		outstring = 'shat'+str(s_hat)+'_pA'+str(pA)+'_M'+str(M)+'_w'+str(windowSize) +'_j'+str(Step)+'.csv'
		final.to_csv(outstring, index = False)
		print(outstring)

if __name__ == '__main__':
    main()
