import glob, argparse
import numpy as np
import pandas as pd

def yieldPops(directory):
	for g in glob.glob(directory+'/*'):
		yield np.array([map(float, i.strip().split() ) for i in open(g)])

def wc2pops(p1,p2):
	p_bar = (p1 + p2)/2.
	s2 = (p1 - p_bar)**2 + (p2 - p_bar)**2 
	if p1 == p2 == 0.0:
		return 'Na'
	return s2 / (s2/2 + p_bar*(1-p_bar))

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
	output = []
	for i in yieldPops(args.inputDir):
		p0 = i.T[3]/sum(i[0])
		fst_l250 = wc2pops(p0[124],p0[374])
		fst_l200 = wc2pops(p0[149],p0[349])
		fst_l150 = wc2pops(p0[174],p0[324])
		fst_l100 = wc2pops(p0[199],p0[299])
		fst_l50 = wc2pops(p0[224],p0[274])
		output.append([p0.mean(), fst_l250, fst_l200, fst_l150, fst_l100, fst_l50])
		
	pd.DataFrame(output, columns = ['p0','l250','l200','l150','l100','l50']).to_csv(args.output, index = False)
		
if '__name__':
	main()
	
