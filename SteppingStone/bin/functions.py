import numpy as np
import itertools as it
#np.random.seed(0)

## The array of genotypes corresponds to the types:
## np.array([n1,n2,n3,n3])
##
## where n1 = the number of AB
## where n2 = the number of Ab
## where n3 = the number of aB
## where n4 = the number of ab
##
## A/a is the selected locus and B/b is the neutral marker

def selectionAltVec(arr, sArray):
	N = arr[1].sum()
	selUnscaled = (sArray*arr)/N
	output = (selUnscaled.T/ np.sum(selUnscaled, axis=1)).T
	output[0] = 0
	output[-1] = 0
	return output


## Gets the deterministic change in genotype frequnecies due to 
## selection
 
def selectionAlt(arr, sArray):
	N = arr[1].sum()
	selUnscaled = (sArray*arr)/N
	output = (selUnscaled.T/ np.sum(selUnscaled, axis=1)).T
	return output
	
## This function takes a Python iterable (list, array etc.) and return
## groups of three items. Eg. threes([1,2,3,4,5,6,7,8,9]) gives:
## [1,2,3],[2,3,4],[3,4,5]...
## This is used to get the nearest neighbours for the migration step
	
def threes(iterator):
#    "s -> (s0,s1,s2), (s1,s2,s3), (s2, s3,4), ..."
    a, b, c = it.tee(iterator, 3)
    next(b, None)
    next(c, None)
    next(c, None)
    return zip(a, b, c)

## This function calucaltes that deterministic change in 
## genotype frequencies under nearest neghtbour dispersal at rate m 
## The migration rate information is captured in the migArray
## This is generated once at the start of a simultation to save time.

## This function takes as input, three arrays: the l.h.s pop, the focal 
## pop and the r.h.s. pop. The migration into the pop from either side 
## is calculated for each of the genotypes, as is the number of 
## genotypes that stay put. This is all done using array multiplication

def mig(mm,migArray):
	if np.isnan(mm[0][0]):
		mm[0] = np.zeros(4)
	if np.isnan(mm[2][0]):
		mm[2] = np.zeros(4)
#	print '**'
#	print mm.T
	X = mm.T*migArray
#	print X
	Z = np.array([X[0][0]+X[0][1]+X[0][2],
					X[1][0]+X[1][1]+X[1][2],
					X[2][0]+X[2][1]+X[2][2],
					X[3][0]+X[3][1]+X[3][2]])
#	print Z
	return Z/(Z[0]+Z[1]+Z[2]+Z[3])

def migVec(mm,migArray):
#	from tom import brace
#	print '!!\n'	
#	staying = mm.T[:, 1:-1]
#	leftSide = mm.T[:, :-2]
#	rightSide = mm.T[:, 2:]
#	print 'start',mm.T
#	print 'staying',staying
#	print 'lhs',leftSide
#	print 'rhs',rightSide
#	print 'migration vector', mm
#	print 'new',staying*migArray[1] + migArray[0]*leftSide + migArray[2]*rightSide
#	brace()
	new = (mm.T[:, 1:-1] * migArray[1] + migArray[0]*mm.T[:, :-2] + migArray[2]*mm.T[:, 2:])

#	print new
#	print new.T.sum(axis=1)
	return (new/ new.T.sum(axis=1)).T

## This function recombines the genotypes present in the population. 
## Note that if the advantageous mutation is absent, there will be no
## effect of recombination so we don't recombine in the drift burn-in
def recombine(pop,c):
	for i in range(len(pop)):
		if np.isnan(pop[i][0]):continue
		pA = pop[i][0] + pop[i][1]
		pB = pop[i][0] + pop[i][2]
		if pA == 0: continue
		pop[i] = np.array([(1.-c)*pop[i][0] + c*pA*pB, (1.-c)*pop[i][1] + c*pA*(1-pB), (1.-c)*pop[i][2] + c*(1-pA)*pB, (1.-c)*pop[i][3] + c*(1-pA)*(1-pB)] )
	return pop
	
## This is the main function of the simulation. Each generation is a 
## churn. The steps of the churn are:
##		Apply selection (deterministic) -- (if necessary)
## 		Recombine (deterministic) -- (if necessary)
##		Migrate (deterministic)
##		Multinomial sampling (stochastic)

def churnPopVector(pop, k, N, migArray, sArray, c, burnin = True):
#	from tom import brace
#	print 'start', pop
	if burnin == False:
		pop = selectionAltVec(pop, sArray)
		if c > 0:
			pop = recombine(pop, c)
	else:
		pop = pop/N
	

	pop = migVec(pop, migArray)
	
	output =  np.array([[0,0,0,0]]+[np.random.multinomial(N, t) for t in pop]+[[0,0,0,0]], dtype=float)

	return  output


## This is the main function of the simulation. Each generation is a 
## churn. The steps of the churn are:
##		Apply selection (deterministic) -- (if necessary)
## 		Recombine (deterministic) -- (if necessary)
##		Migrate (deterministic)
##		Multinomial sampling (stochastic)

## For a given allele frequency, this function returns a correlated 
## value, incorporating some random error.

def correlatedValueA4(p0,r):
	sel =(1+ np.random.normal(0,np.sqrt( 1 - (r*r) )))
	p1 = p0*r
	return p1*sel

## This function simply rounds values up or down
def roundValue(x):
	if x > 1:
		return 1
	elif x < 0:
		return 0
	else:
		return x

## This function seeds the allele freuqnecy across the population
	
def initialiseAlleleFreqs(s1,s2,p0,r):
	p1 = p0
	p2 = p0 
	for k in range(1,len(s1)):
		p1 = roundValue( correlatedValueA4(p1,r) )
		s1[k] = p1
	s1[0] = p0
	for k in range(1,len(s2)):
		p2 = roundValue( correlatedValueA4(p2,r) )
		s2[k] = p2
	s2[0] = p0
	return s1[::-1] +s2[1:]

## This function creates a population array, seeding a neutral mutation 
## at a random position and then spreading it across the
## rest of the population using correlated calues

def initialisePop(N,k, p0):
	r = 0.99
	d1 = int(k*np.random.random())
	popStart = [0]*k
	s1 = popStart[:d1+1]
	s2 =popStart[d1:]

	alleleFreqList = [float(int(i*N)) for i in initialiseAlleleFreqs(s1,s2,p0,r)]
	pop =[]
	pop.append([0,0,0,0])
	for q in alleleFreqList:
		pop.append([0,0,N-q, q])
#		pop = np.array([[0,0,0,N]]*(k-1) + [[1,0,0,N-1]])
	pop.append([0,0,0,0])
	return np.array(pop)

## This function dumps the current state of the population to a 
## simple file, with a column for each genotype

def writePop(pop, name):
	out = open(name, 'w')
	for i in pop:
		if sum(i) == 0:continue # Note that the first and last positions in the array are zeroes, there's no point saving those files
		out.write( '\t'.join(map(str,i)) + '\n' )
	out.close()

## This function takes a large number of population states and writes 
## them to a file with the columns: generation, allele frequency, deme

## NEED TO ADD A FUNCTION LIKE THIS THAT WILL TAKE BOTH THE NEUTRAL 
## AND THE SELECTED ALLELE FREQUNECIES

def writeLotsOfPops(pops, name, N):
	out = open(name, 'w')
	
	for k in pops.keys():
		ind = 0
		for i in pops[k].T[2]:
			if ind == 0:
				ind+=1
				continue
			elif ind == len(pops[k])-1:
				pass
			else:
				out.write(','.join( [ str(k), str(i/N), str(ind) ] ) + '\n')
				ind +=1

#			out.write( '\t'.join(map(str,i)) + '\n' )
	out.close()

## This function takes a large number of population states and writes 
## them to a file with the columns: generation, allele frequency, deme

## NEED TO ADD A FUNCTION LIKE THIS THAT WILL TAKE BOTH THE NEUTRAL 
## AND THE SELECTED ALLELE FREQUNECIES

def writeLotsOfSelectedPops(pops, name):
	out = open(name, 'w')	
	for k in pops.keys():
		pA = (pops[k].T[0] + pops[k].T[1] )/ sum(pops[k][1]) 
		pB = (pops[k].T[0] + pops[k].T[2] )/ sum(pops[k][1]) 
	#	print pA
		
		for i in range(len(pA)):  
			if i == 0:
				continue
			elif i == len(pops[k])-1:
				continue
			else:
				out.write(','.join( [ str(k), str(pA[i]),str(pB[i]), str(i) ] ) + '\n')

	out.close()
## This function takes a large number of population states and writes 
## them to a file with the columns: generation, allele frequency, deme

## NEED TO ADD A FUNCTION LIKE THIS THAT WILL TAKE BOTH THE NEUTRAL 
## AND THE SELECTED ALLELE FREQUNECIES

def writeLotsOfSelectedPopsGenotypes(pops, name):
	out = open(name, 'w')	
	for k in pops.keys():
#		print k, pops[k]
		for i in [50, 100,150,200,250,300,350,400,450]:
			line = map(str, [k, i] + list(pops[k][i])	)
			#print line
			out.write(','.join(line) + '\n')

	out.close()

## This function loads the population from a file that contains
## the population state after a period of burn-in

def loadPop(inputFile, gzip = False):
	if gzip:
		import gzip
		pop = np.array([map(float, i.strip().split() ) for i in gzip.open(inputFile)])
	else:
		pop = np.array([map(float, i.strip().split() ) for i in open(inputFile)])
	return pop

## This function introduces an advantageous mutation into 
## the population at a randomly selected location.
def mutate(pop):
	d1 = int(pop.shape[0]*np.random.random())
	alleleFreq = pop[d1][2]/pop[0].sum()	
	test = np.random.random()
	if test <= alleleFreq:
		pop[d1][0]+=1
		pop[d1][2]-=1
	else:
		pop[d1][1]+=1
		pop[d1][3]-=1
	zero = np.array([[0.,0.,0.,0.]])	
	return np.concatenate((zero,pop, zero))

## DEPRECATED

def churnPop(pop, k, N, migArray, sArray, c, burnin = True):
#	from tom import brace
#	np.random.seed(0)
#	print 'churnStart',pop

#	print 'From the churn',pop.T[2]
	print 'churn',pop
	if burnin == False:
		pop = selectionAlt(pop, sArray)
		if c > 0:
			pop = recombine(pop, c)
	else:
		pop = pop/N
#	print np.random.random()
#	print pop
#	brace()
## I know, I know, this is a really horrible list comprehension
## It takes the populaiton, applies the migration step, performs the 
## multinomial sampling and adds the caps back on to the 
## end of the popuatlio
	print 'chunking', np.array([[0,0,0,0]]+[mig(np.array(t), migArray) for t in threes(pop)]+[[0,0,0,0]], dtype=float).T

	output =  np.array([[0,0,0,0]]+[np.random.multinomial(N, mig(np.array(t), migArray)) for t in threes(pop)]+[[0,0,0,0]], dtype=float)
#	print 'churn2',output
	return output
