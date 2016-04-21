import os

E = [5*10**3,16*10**3, 62.5*10**3]

for e in E:
	for ref in [1,2,3]:
		string = 'mpirun -np 4 python BM.py %s %s'%(ref,e)
		print string
		os.system(string)
