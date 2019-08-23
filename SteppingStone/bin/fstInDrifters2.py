import pandas as pd
import argparse, glob
import numpy as np

def wc2pops(p1,p2):
	p_bar = (p1 + p2)/2.
	s2 = (p1 - p_bar)**2 + (p2 - p_bar)**2 
	if p1 == p2 == 0.0 or p1 == p2 == 1.0:
		return 'Na'
	return s2 /( s2/2 + p_bar*(1-p_bar))

def wc(alleleFreqs):
	p_bar = sum(alleleFreqs)/len(alleleFreqs)
	if p_bar == 1.0:
		return None, None
	s2 = sum([(p - p_bar)**2 for p in alleleFreqs])/len(alleleFreqs)
	return s2 /( s2/2 + p_bar*(1-p_bar))


def parseData(data):#, pop1, pop2, index = True):
	csv =  pd.read_csv(data,header = None,sep = '\t')
#	csv =  pd.read_csv(data,header = None,compression = 'gzip',sep = '\t')
#	csv['pA'] = (csv[2] + csv[3])/(csv[2] + csv[3] + csv[4] + csv[5])
	csv['pB'] = (csv[3])/(csv[0] + csv[1] + csv[2] + csv[3])
	wcPop = wc( list(csv['pB']) )
	meanP = csv['pB'].mean()
	d50 = csv.iloc[[49]]
	d100 = csv.iloc[[99]]
	d150 = csv.iloc[[149]]
	d200 = csv.iloc[[199]]
	d300 = csv.iloc[[299]]
	d350 = csv.iloc[[349]]
	d400 = csv.iloc[[399]]
	d450 = csv.iloc[[449]]

	wc100 = wc2pops(list(d200['pB'])[0],list(d300['pB'])[0])
	wc200 = wc2pops(list(d150['pB'])[0],list(d350['pB'])[0])
	wc300 = wc2pops(list(d100['pB'])[0],list(d400['pB'])[0])
	wc400 = wc2pops(list(d50['pB'])[0],list(d450['pB'])[0])

	return [wcPop, wc100, wc200,wc300,wc400,meanP]

def main():
	parser = argparse.ArgumentParser(description="Stepping stone model")
 	parser.add_argument("-i	", 
 			required = True,
 			dest = "inputDir",
 			type =str, 
 			help = "the directory with the drift runs")
 	parser.add_argument("-o	", 
 			required = True,
 			dest = "output",
 			type =str, 
 			help = "name the output file")
	args = parser.parse_args()

	count = 0
	WCnum = 0
	WCdenom = 0
	outputList = []
	for f in glob.glob(args.inputDir+'/*txt'):
		count +=1
		output = parseData(f)
		outputList.append(output)
#		if type(num) == str or type(denom) == str:continue
#		WCnum += num
#		WCdenom +=denom
#		print WCnum/WCdenom, num/denom
		
	#	if count == 2000:break
	pd.DataFrame(outputList, columns = ['pop','d100','d200','d300','d400','meanP']).to_csv(args.output,index = False)

main()
