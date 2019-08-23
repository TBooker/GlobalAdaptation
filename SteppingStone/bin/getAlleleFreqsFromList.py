import sys 
for i in open(sys.argv[1]):
	if 'e' in i.split('.')[0]:
		print i.split('.')[0]
	else:
		print '.'.join(i.split('.')[:2])
	
	
	
