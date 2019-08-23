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
	parser.add_argument("--store", 
			default = False,
			dest = "store",
			action = 'store_true', 
			help = "Use this flag and the program will save the population state a t")
	parser.add_argument("--gzip", 
			default = False,
			dest = "gzip",
			action = 'store_true', 
			help = "Use this flag if the input file is gzipped")
#	from tom import brace

	args = parser.parse_args()
	if args.gzip:
		startingPop  = loadPop(args.inputFile, gzip = True)
	else:
		startingPop  = loadPop(args.inputFile)
		
	N = startingPop[0].sum()
	k = len(startingPop)
	
	migArray = np.array([ args.m/2, 1-args.m, args.m/2])
	sArray = np.array([1+args.s, 1+args.s, 1, 1])

	pop = mutate(startingPop.copy())

	fixed = False
	generation = 0
	lost = 0
	#if args.store:
	data = {}
#	np.random.seed(0)

	while not fixed:

		selectedAlleleCount = pop.T[0].sum() + pop.T[1].sum()
#		print generation, selectedAlleleCount /(N*k)
		if selectedAlleleCount ==0:
			lost +=1
			pop = mutate(startingPop.copy())
			generation = 0
			#if args.store:
			data = {}

		elif selectedAlleleCount/float(N*k) >= 0.99999:
#			print generation
			fixed = True
	## for testing:
	#	if generation == 50:
	#		break
	#	print generation
	#	print pop
		data[generation] = pop
#		popC = churnPop(pop, k, N, migArray, sArray, args.c, burnin= False)
	#	print '############################\n\n'
		pop = churnPopVector(pop, k, N, migArray, sArray, args.c, burnin= False)
	#	brace()
		generation += 1
	if args.store:
		writeLotsOfSelectedPops(data, args.output)
	else:
		writeLotsOfSelectedPopsGenotypes(data, args.output + '_' + str(round(1./(1+lost),6)) +'_.txt')

#	print 'lost',lost,'advantageous mutations',round(1./lost,6)
# Every ith generation, store the population state, and re-initialise from that point if the allele is lost

if '__name__':
	main()
