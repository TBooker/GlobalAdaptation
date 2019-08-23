import numpy as np
import pandas as pd
import glob, sys

allelesFreqs = []
for i in glob.glob(sys.argv[1]+'/*segSites.txt'):
	x = i.split('/')[-1].split('_')

#	print allF[(allF['q'] > 0.3) & (allF['q'] < 0.7)]

		
	temp = pd.read_csv(i)
	temp['REP'] = x[1]
	temp['gen'] = x[3]
	allelesFreqs.append(temp)

pd.concat(allelesFreqs).to_csv(sys.argv[2]+'.alleleFreqs.csv', index = False)