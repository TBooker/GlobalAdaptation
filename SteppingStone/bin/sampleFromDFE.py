### This script will simulate the occurence of advantageous mutations arising in a species.
import numpy as np
import pandas as pd
import random, glob, argparse, subprocess
from collections import OrderedDict
#from tom import brace

def wc2pops(p1,p2):
	p_bar = (p1 + p2)/2.
	s2 = (p1 - p_bar)**2 + (p2 - p_bar)**2 
	if p1 == p2 == 0.0 or p1 == p2 == 1.0:
		return 'Na'
	return s2 / (s2/2 + p_bar*(1-p_bar))

def expDFEdraw(mean):
	
### Generate a random number from an exponential distribution with 
### mean given as an argument of the function
#	return np.random.choice([0.1,0.01,0.001],  p=[0.0, 0.01, 0.99])
	s = None
	while s == None:
		s = np.random.exponential(scale=mean)
#		print '!', s
		if s < 0.0001:
			return 0.0
		else:
			return s
#			s =  round(s,4)
#			if len(str(s)) == 6:
#				return s
#			elif len(str(s)) == 5:
	#			print s
#				return "{:5.4f}".format(s)
		
def DFEdraw():
## For testing
	s = np.random.choice([0.1,0.01,0.001],  p=[0.0, 0.01, 0.99])
	return s
	
def fixation(s):
#	print s	
#	Fst = 0.132
	FixProb = 2*s#*(1-Fst) ## Just for testing, the actual fixation probability has Fst built into it

	test = np.random.random()
	if test <= FixProb:
#		print test, 2*s
		return True
	else:
		return False

def fileGrabber(s, simulationDirectory):
#	print simulationDirectory + '/'+str(s)+'/*'
#	print simulationDirectory + '/'+str(s)+'/*'
	files = glob.glob(simulationDirectory + '/'+str(s)+'/*')
#	print s
#	print files [:2]
	return random.choice(files)

def parseData(data, gen, pop1, pop2, s, mut, index = True):
	csv =  pd.read_csv(data,header = None)
	csv['pA'] = (csv[2] + csv[3])/(csv[2] + csv[3] + csv[4] + csv[5])
	csv['pB'] = (csv[2] + csv[4])/(csv[2] + csv[3] + csv[4] + csv[5])
	d1 = csv[csv[1] == pop1].copy()
	d2 = csv[csv[1] == pop2].copy()
	temp = np.array([ [wc2pops(i,j),k+gen, s, mut] for i,j,k in zip(d1['pB'],d2['pB'],d1[0]) ])
	if index:
		return pd.DataFrame(temp, columns = ['fst','generation', 's','mut']).set_index('generation')
	else:
		return pd.DataFrame(temp, columns = ['fst','generation','s','mut'])

def parseData2(data, gen, pop1, pop2, s, mut, index = True):
	csv =  pd.read_csv(data,header = None)
	csv['pA'] = (csv[2] + csv[3])/(csv[2] + csv[3] + csv[4] + csv[5])
	csv['pB'] = (csv[2] + csv[4])/(csv[2] + csv[3] + csv[4] + csv[5])
	d1 = csv[csv[1] == pop1].copy()
	d2 = csv[csv[1] == pop2].copy()
	wc = np.array([ wc2pops(i,j)for i,j in zip(d1['pB'],d2['pB']) ])
	dat = pd.DataFrame({'WC' : wc,
		's' : np.full_like(wc, s),
		'mut' :np.full_like(wc, mut),
		'generation' : np.array(d1[0] + gen)})

#	print 'generation', dat['generation']
#	dat['WC'] = np.array([ wc2pops(i,j)for i,j in zip(d1['pB'],d2['pB']) ])
#	dat['s']  = np.full_like(dat['WC'], s)
#	dat['mut'] = np.full_like(dat['WC'], mut)
#	dat['generation'] = np.array(d1[0] + gen)

	return dat

def main():
#	Run the simulation, for an arbitary number of generations
	parser = argparse.ArgumentParser(description="A script that renames a bunch of neutral simulations so that they have a searchable name")
	parser.add_argument("--direc", "-d", 
			required = True,
			dest = "direc",
			type = str, 
			help = "The directory that contains sub-directories with ")
	parser.add_argument("--generations","-g", 
			required = True,
			dest = "generations",
			type = int, 
			help = "The number of generations to run for")
	parser.add_argument("--output", "-o", 
			required = True,
			dest = "output",
			type = str, 
			help = "Give a name to the output file, it will be in CSV format")
	parser.add_argument("--Ua", 
			required = True,
			dest = "Ua",
			type = float, 
			help = "What is the mean number of new advantageous mutations per generatio?")
	parser.add_argument("--exp_mean", 
			required = False,
			dest = "exp_mean",
			type = float, 
			help = "Give the mean of an exponential DFE",
			default = 0)
	parser.add_argument("--sampleT", 
			required = False,
			dest = "sampleT",
			type = int, 
			help = "Give the time at which you want to sample the population",
			default = 0)	
	parser.add_argument("--store", 
			default = False,
			dest = "store",
			action = 'store_true', 
			help = "Use this flag and the program will save the trajectories of all mutations, otherwise it will only save the mutations that are segregating when ")
	parser.add_argument("--distance", 
			required = True,
			dest = "distance",
			type = int, 
			help = "The distance between sampled demes")
	args = parser.parse_args()

	generations = args.generations
	fixed = 0
	arose = 0
#	simulationDirectory = '/home/booker/work/FishersWaveFst/simulations/selectionTestFast/'
	simulationDirectory = args.direc
	muts_to_add = OrderedDict()
	
	if args.distance == 100:
		d1 = 200 
		d2 = 300
	elif args.distance == 200:
		d1 = 150 
		d2 = 350
	elif args.distance == 300:
		d1 = 100 
		d2 = 400
	elif args.distance == 400:
		d1 = 50 
		d2 = 450
	elif args.distance not in [100,200,300,400]:
		print 'Need to give a distance from [100,200,300,400]'
		return
	
	for g in range(generations):
		if g% 10000 == 0:
			print g
		## Start with new mutations per generation
		mutations = np.random.poisson(args.Ua) # The value here is the expected number of advantageous mutations that arise per generation
		for mutation in range(mutations):
			arose +=1

			if args.exp_mean == 0:
				s = DFEdraw()
			else:
				s = expDFEdraw(args.exp_mean)
			if s >0.1:
				s = 0.1

			if len(str(round(s,4))) == 6:
					t = round(s,4)
			elif len(str(round(s,4))) == 5:
				t = "{:5.4f}".format(round(s,4))
			elif s == 0:
				continue
		#	t = 0.01 ## For testing
#			print t
			x = fileGrabber(t,simulationDirectory)
#			print x
			FixProb = x.split('/')[-1].split('_')[2]
			if fixation(float(FixProb)):
				fixed += 1
	
				muts_to_add[arose] = [x,g,s]
					
				## Here we draw an advantageous mutation from the appropriate directory
	print fixed, arose, float(fixed)/arose

	count = 0
	if not  args.store:
		output = []
	for i in muts_to_add.keys():
#		print i
		print muts_to_add[i]
		fileName = muts_to_add[i][0]
		selectionCoefficient = muts_to_add[i][2]
		generation = muts_to_add[i][1]
		line = subprocess.check_output(['tail', '-1', fileName])
		if int(line.split(',')[0]) + generation < args.sampleT:
			continue
		else:
			count +=1

		tempCSV = parseData2(fileName, generation, d1, d2, selectionCoefficient, i, index = False)
		tempCSV['lab'] = count
		if args.store:
			if count == 1:
				DATAOUT = tempCSV.copy()
			else:
				DATAOUT = DATAOUT.append(tempCSV)
		else:
			output.append( tempCSV[tempCSV['generation'] == args.sampleT ].copy() )
	if args.store:
		DATAOUT.to_csv(args.output, index = False)
	else:
		pd.concat(output).to_csv(args.output, index = False)
	return

main()
