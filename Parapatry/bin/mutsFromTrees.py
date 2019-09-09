# Let's sprinkle mutations onto the subPop tree
import msprime, pyslim,  sys, subprocess, argparse, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from collections import Counter

 
def MutsFromTrees(Trees):

	print('read trees')

	ts = pyslim.load(Trees).simplify()

	variants =[]
#	for  variant in pyslim.extract_mutation_metadata(ts.tables):
#		print(variant)
#	return

	for variant, mut, node in zip(ts.variants(), pyslim.extract_mutation_metadata(ts.tables), pyslim.extract_node_metadata(ts.tables) ):
		
#		print ()
#		print (variant)
#		print (mut)
#		print (node)

		freq = sum(variant.genotypes)/len(variant.genotypes)
#		print (freq)

		if freq >= 1:continue
		if len(mut) > 1: pass
		variants.append([variant.site.position,freq, mut[0].selection_coeff])

	return variants



def main():
	parser = argparse.ArgumentParser(description="This script takes Trees files from SLiMs and uses msprime/pyslim to add mutations and writes a VCF file")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the output files")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "Name the output file, give a prefix")
#	parser.add_argument("--sample","-s", 
#			required = False,
#			dest = "sample",
#			type = int,
#			help = "How many individuals do you want to sample from each deme",
#			default = 50)
#	parser.add_argument("--just_vcf", 
#			required = False,
#			dest = "just_vcf",
#			action = 'store_true',
#			help = "Just make a VCF and then quit",
#			default = False)			
	args = parser.parse_args()
	
	variantsRaw = MutsFromTrees(args.input)

	variants = pd.DataFrame(variantsRaw, columns =  ['POS','q', 's'])
#	print(args.input)
#	print(variants)

	geneStarts = np.array([47500, 152500, 257500, 362500, 467500, 572500, 677500, 782500, 887500, 992500, 1097500, 1202500, 1307500, 1412500, 1517500, 1622500, 1727500, 1832500, 1937500])

	geneEnds =  geneStarts + 4999
	counter = 0 
	newVars = []
	for i, j in zip(geneStarts, geneEnds):
		counter += 1 
		temp= variants[(variants['POS'] >= i) & (variants['POS'] <= j)].copy()
		temp['locus'] = counter
		newVars.append(temp)
	pd.concat(newVars).to_csv(args.output, index =False)

if __name__ == '__main__':
    main()
