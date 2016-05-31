#!/usr/bin/python
 
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
 
import sys, time, os
import pymol
#import build_seq
 
pymol.finish_launching()
 
##
# Read User Input
#print(sys.argv[0])
archivo=sys.argv[1]
f = open(archivo)

#Crea un diccionario con el id de la secuancia como llave y la secuencia como valor 
seqs={}
for line in f:
	line=line.rstrip()
	if  line.startswith('>') :#"line[0]=='>': # or line.startswith('>')			
		words=line.split()
		name=words[0][1:]
		seqs[name]=''
	else : # sequence, not header
		seqs[name] = seqs[name] + line	

#Cierra el archivo
f.close()
'''
staticStructurePath = os.path.abspath(sys.argv[1])
staticStructureName = staticStructurePath.split('/')[-1].split('.')[0]
mobileStructurePath = os.path.abspath(sys.argv[2])
mobileStructureName = mobileStructurePath.split('/')[-1].split('.')[0]
''' 

# Load Structures
 
#pymol.cmd.load(mobileStructurePath, mobileStructureName)
#pymol.cmd.load(staticStructurePath, staticStructureName)
 
# CEAlign STATIC, MOBILE
# CEAlign produces same alignment as the complex SUPER below
pymol.cmd.do("run build_seq.py ") # Import Module
i=sys.argv[2]
for numId,secuencia in seqs.items():
	#i=i+1
	#pymol.cmd.do("create object"+str(i)+" , (all)")
	pymol.cmd.do("build_seq "+secuencia+", ss=helix")
	#pymol.cmd.do("build_seq QGAADLESLGQYFEEMKTKLIQDMTE, ss=helix")
	time.sleep(2) # Dunno why, but if I don't wait, structures do not align properly..

	# Save Superimposition
	# save(file, selection, state (0 default), format)
	#pymol.cmd.save("/home/proyecto/Documentos/resultados/"+numId+".pdb", "ALL", -1, 'pdb')
	n=str(i)
	pymol.cmd.save("/home/proyecto/Documentos/resultados/temp"+n+".pdb", "ALL", -1, 'pdb')

	comando="perl -T /home/proyecto/Documentos/PDB/pymol/rename.pl -infile /home/proyecto/Documentos/resultados/temp"+n+".pdb -tochain A"
	#perl -T /home/proyecto/Documentos/PDB/pymol/rename.pl -infile /home/proyecto/Documentos/resultados/temp1.pdb -tochain A
	resultado = os.system(comando)

	time.sleep(1)
	
	comando="python /home/proyecto/Documentos/discotope-1.1/discotope.py -f '/home/proyecto/Documentos/resultados/temp"+n+".pdb' -chain A  >> /home/proyecto/Documentos/resultados/resultadoDiscotope"+n+".txt"
	resultado = os.system(comando)

	time.sleep(1)

	#pymol.cmd.do("delete object"+str(i))
 
# Get out!
pymol.cmd.quit()