import pandas as pd
import sys, os, glob, argparse

def main():
	parser = argparse.ArgumentParser(description="A script that renames a bunch of neutral simulations so that they have a searchable name")
	parser.add_argument("--direc", 
			required = True,
			dest = "direc",
			type = str, 
			help = "The path to a directory containing the neutral burn-in simulations")

        parser.add_argument("-N",
                        required = True,
                        dest = "N",
                        type = int,
                        help = "Nnumber of individuals per deme")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "the directory where you want to deposit the renamed files")
	args = parser.parse_args()

	numDigits = 6 # The number of digits that you want to retain for the file name
	
	files = glob.glob(args.direc + '/*')
	counter = 0
	for i in files:
		counter +=1
		temp = pd.read_csv(i, header = None, sep = '\t')
		mean = (temp[3]/args.N).mean()
#		print str(round(mean, 5))
		name = str(round(mean, 5))+'.'+str(counter)+'.txt'
#		print name
		print 'cp ' + i + ' ' + args.output + '/'+ name
		os.system('cp ' + i + ' ' + args.output + '/'+ name)

main()
