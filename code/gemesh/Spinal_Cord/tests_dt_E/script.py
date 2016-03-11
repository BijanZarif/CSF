from FSI_Solver import *
import time
set_log_active(False)
t0 = time.time()

rho_f = (1./1000)		# g/mm
nu_f = (0.658)			# mm**2/s
mu_f = (nu_f*rho_f)		# g/(mm*s)


Youngs = [62.5*10**3, 0.7*10**6, 10**6]

Pr = 0.479

delta_t = [0.0002, 0.0001]


for E in Youngs:
	A = FSI_Solver(rho_f, nu_f, E, Pr)
	A.make_functions('longer_top_spinal.xml')
	for dt in delta_t:
		print '\n Running with dt = %g and E = %g \n' %(dt,E)
		folder = 'RESULTS/dt_%d_E_%d'%(int(dt*1E4),E/10**2)
		sc = int(20*(delta_t[0]/dt))								# img every 10 millisecond. 
		A.run(dt,5,folder,save_count=sc,solver='default')  
	
		print 'Ended run for dt = %g and E = %g \n' %(dt,E)
