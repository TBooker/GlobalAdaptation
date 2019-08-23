import argparse
from SSfunctions import * 
from tom import brace
def main():
	parser = argparse.ArgumentParser(description="Stepping stone model")
	parser.add_argument("-k", 
			required = True,
			dest = "k",
			type = int, 
			help = "number of demes")
	parser.add_argument("-N", 
			required = True,
			dest = "N",
			type = int, 
			help = "population size")
#	parser.add_argument("-s", 
#			required = True,
#			dest = "s",
#			type = float, 
#			help = "selection in favour of the derived allele")
# 	parser.add_argument("-t", 
# 			required = False,
# 			dest = "t",
# 			type =int, 
# 			default = None,
# 			help = "The number of generations to run the simulation for, the default is to run the simulation until there is complete fixation")
	parser.add_argument("-m", 
			required = True,
			dest = "m",
			type =float, 
			help = "the probability of migration")
	parser.add_argument("-p", 
			required = True,
			dest = "p0",
			type =float, 
			help = "the initial frequnecy of the neutral allele")
 	parser.add_argument("-d", 
 			required = True,
 			dest = "drift",
 			type = int, 
 			help = "Here you can specify the number of generations that you want the simulation to drift for")
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
 	parser.add_argument("--seed", 
 			required = False,
 			dest = "seed",
 			type = int, 
 			help = "Give a random number seed to the simulator, FOR TESTING",
 			default = -123)

	args = parser.parse_args()
	## For testing
	if args.seed != -123:
		np.random.seed(args.seed)

	if args.p0 == 99:
		args.p0 = np.random.random()*1.0
		args.output = str(args.p0)[:5]+'.'+args.output
		pop = initialisePop(args.N,args.k,args.p0)
	else:
		pop = initialisePop(args.N,args.k,args.p0)
#	print pop.T[3]
#	print gimmeDemes(pop, 1, args.k)
#	display = []
	worker = True
	attempt = 0
	migArray = np.array([ args.m/2., 1-args.m, args.m/2.])
	sArray = np.array([1, 1, 1, 1])
	if args.store:
		data = {}
	while worker:
		for i in range(args.drift):
			if args.store:
				if i % 100 == 0:
					data[i] = pop
			if i % 1000 == 0: print i
			
			pop = churnPopVector(pop,  args.k, args.N, migArray, sArray, 0)
			if np.sum(pop, axis=0)[3] == 0:
				pop = initialisePop(args.N,args.k,args.p0)
				print 'restarting'
				break
			if i == args.drift -1:
				worker = False
	if args.store:
		writeLotsOfPops( data , args.output)
	else:
		writePop( pop , args.output)

if '__name__':
	main()
