import numpy as np
import pandas as pd
import glob, sys

starts = np.array([45000,150000, 360000,255000,360000,465000,570000,675000,780000,885000,990000,1095000,1200000,1305000,1410000,1515000,1620000,1725000,1830000,1935000,2040000])


allFreqs = []

count = 0 
for i in glob.glob(sys.argv[1]+'/*w10000.*'):

	x = i.split('/')[-1].split('_')
		
	temp = pd.read_csv(i, sep = '\t')
	temp['REP'] = x[1]
	temp['gen'] = x[3]

	if int(x[3]) == int(sys.argv[2]):
		count += 1
                if count ==1:
                        temp.to_csv('gen'+str(sys.argv[2])+'.' + sys.argv[3], index = False)

                else:
                        temp.to_csv('gen'+str(sys.argv[2])+'.' + sys.argv[3], mode = 'a',header = False, index = False)
                       
#pd.concat(allFreqs).to_csv('alleleFreqs.csv', index = False)
