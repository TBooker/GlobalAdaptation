import glob, os, argparse, random
from multiprocessing import Pool

def main():
	parser = argparse.ArgumentParser(description="Run positive selection simulations")
	parser.add_argument("--direc", 
			required = True,
			dest = "direc",
			type = str, 
			help = "The path to a directory containing the neutral burn-in simulations")
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
			required = True,
			dest = "c",
			type =float, 
			help = "the recombination rate")
	parser.add_argument("--reps", 
			required = True,
			dest = "reps",
			type = int, 
			help = "the number of replicates you want to run")
        parser.add_argument("--nproc",
                        required = False,
                        dest = "nproc",
                        type = int,
			default = 4,
                        help = "the number of threads you want to use [4]")
	parser.add_argument("-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "the prefix you want to give to each of the output files")
 	parser.add_argument("--batch", 
 			default = False,
 			dest = "batch",
 			action = 'store_true', 
 			help = "Use this flag when saving files to directory")

			
	args = parser.parse_args()
	
	myfiles = [f for f in glob.glob(args.direc+'/*')]
	
	commands = []
	
	for i in range(args.reps):
		selectedNeturalRun = random.choice(myfiles)
		selectedNeturalRunName = selectedNeturalRun.split('/')[-1]
		if selectedNeturalRunName.endswith('.gz'):
			outputName = args.output + '.'+str(i) +'.' + selectedNeturalRunName.split('.gz')[0]

		else:
			outputName = args.output + '.'+str(i) +'.' + selectedNeturalRunName
		if args.batch:
			outputName = args.output + str(i) +'_' + str(args.s) 
		pythonString = 'python /home/booker/FishersWaveFst/simulations/bin/SSselecter.py -i '+selectedNeturalRun + ' --gzip -s ' + str(args.s) + ' -m ' + str(args.m) + ' -c ' + str(args.c) + ' -o ' + outputName
		commands.append(pythonString)

#	return
	print commands
	
	p = Pool(args.nproc)
	p.map(os.system, commands)


if '__name__':
	main()
