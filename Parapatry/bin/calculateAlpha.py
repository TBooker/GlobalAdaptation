import pandas as pd
import glob, sys
from collections import Counter
SubsData = {}

for i in glob.glob(sys.argv[1] + '/*.txt'):
	subCounter = Counter()
	for j in open(i):	
		x = j.strip().split(' ')
		if len( x ) < 8:
			continue
		fixationGen = int(x[8])
		s = float(x[4])
		mutType = x[2]
		if fixationGen <= 200000: continue

		if mutType == 'm1':
			subCounter['m1'] += 1 
		elif mutType == 'm2':
			if s * 20000 <= 1:
				subCounter['m1'] += 1 
			elif s * 20000 >= 1:
				subCounter['m2'] += 1 
	SubsData[i] = subCounter

final = []
for k in SubsData.keys():
	if SubsData[k]['m1'] == 0:
		continue
	FileData = {}
	label = k.split('/')[-1].split('.t')[0].split('_')
	print 
	print label
	for b,d in zip(label[::2], label[1::2]):
		FileData[b] = d
	FileData['m1'] = SubsData[k]['m1']
	FileData['m2'] = SubsData[k]['m2']
	FileData['alpha'] = float(SubsData[k]['m2']) / (SubsData[k]['m1'] + SubsData[k]['m2'] )
	final.append(FileData)

pd.DataFrame(final).to_csv('alphaEstimates.csv', index = False)
		