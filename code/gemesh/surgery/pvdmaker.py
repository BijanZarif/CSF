u = 'velocity'


outfile = open('%s.pvd'%u,'w')

n = 2111
outfile.write('<?xml version="1.0"?>\n')
outfile.write('<VTKFile type="Collection" version="0.1">\n')
outfile.write('  <Collection>\n')
for i in range(n):
	outfile.write('    <DataSet timestep="%d" part="0" file="%s%.6d.vtu" />\n'%(i,u,i))
outfile.write('  </Collection>\n')
outfile.write('</VTKFile>')

outfile.close()
