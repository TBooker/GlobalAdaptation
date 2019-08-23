# -*- coding: utf-8 -*-
import sys, glob, pandas as pd, argparse
from tom import brace

def main():
	parser = argparse.ArgumentParser(description="This script takes a directory full of Fst-per-generation files and generates files which can be used to make a plot like Figure 1")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the Fst files")
	parser.add_argument("--alleleFreq","-a", 
			required = True,
			dest = "alleleFreq",
			type =str,
			help = "Name the file of allele frequencies, output from parseAlleleFreqs.py")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "Name the output file, several files will be generated with different prefixes, so don't specify path!")
	parser.add_argument("--num", 
			required = False,
			dest = "num",
			type = int,
			help = "How many traces do you want to extract? [3]",
			default = 3)

	args = parser.parse_args()

### Here's a script that will do the following:

	#	1. For each generation, get the top 10 Fst hits
	#	2. For each hit, get the trace of Fst over time
	#	3. For each hit get the allele frequencies at the overlapping locus
	#	4. Combine and result in a dataframe with Fst over time
	#	5. Add zeros and ones to the dataframe before and after mutation and fixation to make the plot make more sense
	#	6. Make a separate file with X randomly chosen trajectories

	## These are the location of the gene-like regions in the simulations
	starts = [47500, 152500, 257500, 362500, 467500, 572500, 677500, 782500, 887500, 992500, 1097500, 1202500, 1307500, 1412500, 1517500, 1622500, 1727500, 1832500, 1937500, 2042500]

	# Make a place to store the data

	data = []

	# read in the allele frequency data
	alleleFreqs = pd.read_csv(args.alleleFreq)

	# make a id variable for easy searching later
	alleleFreqs['MutId'] = alleleFreqs['POS'].map(str) + '_' + alleleFreqs['REP'].map(str)

	#brace()

	count = 0
	zz = []
	# Loop over each of the processed files in the chosen directory
	for i in glob.glob(args.input + '/gen*'):
		print i

#		count += 1 
#		if count == 10:
#				break
	## read in the per-generation files
		temp = pd.read_csv(i)
		temp = temp[ temp['BIN_END'] < 2130001].copy()

	## get the max value of FST for each replicate this generation
		top5index = temp.groupby(['REP'])['WEIGHTED_FST'].transform(max) == temp['WEIGHTED_FST']
	
	## get the top 5 values of FST this generation
		largest = temp[ top5index ].nlargest(5, 'WEIGHTED_FST').copy()
	## Extract the allele frequency data for this generation
		genAlleleFreqs = alleleFreqs[alleleFreqs['gen'] == temp['gen'].max()].copy()

	## Make some dummy variables
		largest['diff'] = 999999999
		largest['MutId'] = '999999999_'
		largest['firstGen'] = 999999999

	## iterating over the top 5 Fst hits

		for index, row in largest.iterrows():
			pos = (row.BIN_END + row.BIN_START) / 2 

	## get all segregating sites corresponding to the replicate the Fst hit came from
			genRepAlleleFreqs = genAlleleFreqs[ genAlleleFreqs['REP'] == row.REP].copy()
			if len(genRepAlleleFreqs) == 0:
				continue
	## get the distance between the Fst hit and the analysis window
			genRepAlleleFreqs['diff'] =  genRepAlleleFreqs['POS'] - pos
	## get the closest segregating site to the centre of the analysis window
			genRepAlleleFreqs = genRepAlleleFreqs[abs(genRepAlleleFreqs['diff']) ==  abs(genRepAlleleFreqs['diff']).min() ].copy()

	## get the distance 

			diff = list(genRepAlleleFreqs['diff'])[0]
	## if the distance is greater than 12500 bp away, let's not consider it linked 
			if abs(diff) >12500: continue

			MutId = list(genRepAlleleFreqs['MutId'])[0]
			firstGen = list(genRepAlleleFreqs['gen'])[0]
	## add the relavant stats to the datafame, so that we can track the mutation freqs through time
		
			largest.at[index, 'diff'] = diff
			largest.at[index, 'MutId'] = MutId
			largest.at[index, 'firstGen'] = firstGen


		largestSkim = largest[largest['diff'] != 999999999].copy()
		data.append( largestSkim )
	
	

	top5s = pd.concat(data)

	top5s['id'] = top5s['REP'].map(str) + '_' + top5s['BIN_START'].map(str)

	top5ids = list(top5s['id'])

	count = 0

	tracer = []


	# Loop over each of the processed files in the chosen directory
	for i in glob.glob(args.input + '/gen*'):
		print i
#		count += 1 
#		if count >= 10:
#			break

	# Read in the file
		temp = pd.read_csv(i)
		temp['id'] = temp['REP'].map(str) + '_' + temp['BIN_START'].map(str)
	# extract the rows corresponding to the top5s from above	
		matches = temp[temp["id"].isin(top5ids)].copy()

		print 'matches'
		matches['MutId'] = 'XX_'
		matches['MutGen'] = 'XX_'

		for index, row in matches.iterrows():
			mutID = top5s[top5s.id == row.id].iloc[0]
		
			matches.at[index, 'MutId'] = mutID.MutId
			matches.at[index, 'firstGen'] = mutID.firstGen

		tracer.append( matches )

	# Make one big DF from all the little ones
	tracerOut = pd.concat(tracer)

	# Write the Fst trace to file
	tracerOut.to_csv('tracer.'+args.output)

	# merge the Fst trace and the alleleFreq trace dataframes.

	df_merge_col = pd.merge(tracerOut,alleleFreqs, how='left', on=['MutId', 'gen'])

	print df_merge_col['gen'].min(), df_merge_col['gen'].max()

	timeSpan = df_merge_col['gen'].max() - df_merge_col['gen'].min()
	startGen = df_merge_col['gen'].min()
	timeSlice = int( timeSpan/args.num )

	print timeSpan, timeSlice
	end = startGen
	Maxxers = []
	
	for slice in range(args.num):
		start = end 
		end = end + timeSlice 
		print start, end

		FstSlice = df_merge_col[(df_merge_col['gen'] >= start) & (df_merge_col['gen'] <= end)]

		Maxxers.append( FstSlice.ix[FstSlice['WEIGHTED_FST'].idxmax()].id )

		print FstSlice.ix[FstSlice['WEIGHTED_FST'].idxmax()].id
		
	MaxMatch = df_merge_col[df_merge_col["id"].isin(Maxxers)].copy()


		
	## Psuedo-code
	#	Using the first and final generations get the total time of the simulation 
	#	Divide the total time up into N chunks
	#	For each chunk:
	#		grab the max Fst, and retain the ID
	#
	#	From the list of IDs, make a new DF (will be much smaller than the other)
	#	write to file
	#
	#

	## Add zeroes and ones before and after the adv.mut, respectively
	
	
	

	# Write the combined trace to file
	df_merge_col.to_csv('freqs.' + args.output, index = False)


	# Write the combined trace to file - from the sample
	MaxMatch.to_csv('freqSlice.' + args.output, index = False)

	sys.exit()


main()








