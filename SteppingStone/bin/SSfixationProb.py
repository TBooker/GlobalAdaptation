import argparse
from SSfunctions import * 
def main():
	parser = argparse.ArgumentParser(description="Stepping stone model")
	parser.add_argument("-i	", 
			required = True,
			dest = "inputFile",
			type =str, 
			help = "the drift seed fle that you want to use")
	parser.add_argument("-s", 
			required = True,
			dest = "s",
			type = float, 
			help = "selection in favour of the derived allele")
	parser.add_argument("-m", 
			required = True,
			dest = "m",
			type =float, 
			help = "the probability of migration")
	parser.add_argument("-c", 
			required = False,
			dest = "c",
			type =float, 
			help = "the probability of recombination between the selected and neutral loci [0.0]",
			default = 0)
	parser.add_argument("-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "the name you want to give to the output file")

	from tom import brace

	args = parser.parse_args()
	startingPop  = loadPop(args.inputFile, gzip = True)
	N = startingPop[0].sum()
	k = len(startingPop)
	
	migArray = np.array([ args.m/2, 1-args.m, args.m/2])
	sArray = np.array([1+args.s, 1+args.s, 1, 1])

	pop = mutate(startingPop.copy())

	generation = 0
	lost = False
	data = {}

	while not lost:
		selectedAlleleCount = pop.T[0].sum() + pop.T[1].sum()
		if selectedAlleleCount ==0:
			lost = 'lost'
		elif selectedAlleleCount/float(N*k) >= 0.99999:
			fixed = True
			lost = 'fixed'
		data[generation] = pop
		pop = churnPopVector(pop, k, N, migArray, sArray, args.c, burnin= False)
		generation += 1
		print generation, selectedAlleleCount
	output = open(args.output, 'w')
	if lost == 'fixed':
		output.write('1\n')
	elif lost == 'lost':
		output.write('0\n')
	output.close()
if '__name__':
	main()
