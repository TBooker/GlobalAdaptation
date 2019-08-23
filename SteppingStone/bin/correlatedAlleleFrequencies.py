import numpy as np

def correlatedValueA1(p0,r):
	p1 = p0*r + numpy.random.normal(0,numpy.sqrt( 1 - (r*r) ))
	return p1	

def correlatedValueA2(p0,r):
	p1 = p0*r + numpy.random.normal(0,numpy.sqrt( 1 - (r*r) ))
	if p1 <= 0:
		return 0.0
	elif p1 >= 1:
		return 1.0
	else:
		return p1
		
def correlatedValueA3(p0,r):
	p1 = p0*r + numpy.random.normal(0,numpy.sqrt( 1 - (r*r) ))
	return 1./(1.+ numpy.exp(-p1) )	
		
def correlatedValueA4(p0,r):
	sel =(1+ np.random.normal(0,np.sqrt( 1 - (r*r) )))
	p1 = p0*r
	return p1*sel
def roundValue(x):
	if x > 1:
		return 1
	elif x < 0:
		return 0
	else:
		return x
	
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

def main():
	
	k = 300 
	p0 = 0.01
	r = 0.98
	d1 = int(k*np.random.random())
	pop = [0]*k
	s1 = pop[:d1+1]
	s2 =pop[d1:]
	alleleFreqList = initialiseAlleleFreqs(s1,s2,p0,r)
	print [int(i*1000) for i in alleleFreqList]


main()
