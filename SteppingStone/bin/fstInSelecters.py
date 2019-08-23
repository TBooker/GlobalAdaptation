### This script will calculate the Fst between different points in the landscape for all timepoints for a run. Using this, we can then get something like this:

## generation	Fst400	Fst300	Fst200	Fst100

import argparse

def gimmeGens(inputFile):
	current = []
	ident = '0'
	for i in open(inputFile):
		x = i.strip().split(',')
		z = x[:2] + map(float, x[2:])
		if x[0] != ident:
			yield current
			current = []
		current.append(z)
		ident = x[0]  

def wc2pops(pb1,pb2):
	p1 = (pb1[2] + pb1[4]) / sum(pb1[2:])
	p2 = (pb2[2] + pb2[4]) / sum(pb2[2:]) 
#	print p1, p2
	p_bar = (p1 + p2)/2.
	s2 = (p1 - p_bar)**2 + (p2 - p_bar)**2 
	if p1 == p2 == 0.0 or p1 == p2 == 1.0:
		return 'Na'
	return s2 / (s2/2 + p_bar*(1-p_bar))

def analyseGen(gen):
	fst_100 = wc2pops(gen[3], gen[4])
	fst_200 = wc2pops(gen[2], gen[5])
	fst_300 = wc2pops(gen[1], gen[6])
	fst_400 = wc2pops(gen[0], gen[7])
	return map(str, [fst_100,fst_200,fst_300,fst_400] )

def main():
	parser = argparse.ArgumentParser(description="this script will calulate Fst between pairs of populations from simulation output. It will write the output to  a file with the name INPUT.FST.csv")
	parser.add_argument("-i", 
			required = True,
			dest = "input",
			type = str, 
			help = "The input file that you want to get Fst for")
	args = parser.parse_args()
	output = open(args.input+'.fst.csv','w')
	for i in gimmeGens(args.input):
		name = args.input.split('.')[1] ## BREAKABLE!!
		line = analyseGen( i )
	#	print ','.join([args.input, i[0][0]] + line)
		output.write( ','.join([args.input, i[0][0]] + line) + '\n')
	output.close() 
#		genotypes = map(float, i[2:])

main()
